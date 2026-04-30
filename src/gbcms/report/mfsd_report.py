"""
mFSD Per-Variant HTML Report Generator.

Generates a standalone, interactive HTML report with per-variant fragment size
distributions, CH-vs-ctDNA fragment origin signals, and summary statistics.
Uses Plotly.js for interactive histograms with normalized KDE density overlays
and a STRiDE-inspired design system.
"""

import json
import logging
import math
from pathlib import Path
from typing import Any

from gbcms.io.output import CH_GENES

logger = logging.getLogger(__name__)


# ── Fragment Origin Signal Classification ────────────────────────────────────
def _classify_origin(
    hugo: str,
    sub_nuc_enrichment: float,
    ks_pval: float,
    alt_count: int,
    min_alt: int,
) -> tuple[str, str]:
    """Classify variant as TUMOR-LIKE, CH-LIKE, AMBIGUOUS, or INSUFFICIENT.

    Returns (signal_label, explanation).
    """
    if alt_count < min_alt:
        return "INSUFFICIENT", f"ALT fragment count ({alt_count}) below threshold ({min_alt})."

    is_ch_gene = hugo.upper() in CH_GENES if hugo else False
    enriched = not math.isnan(sub_nuc_enrichment) and sub_nuc_enrichment > 1.3
    sig_ks = not math.isnan(ks_pval) and ks_pval < 0.05
    low_enrich = not math.isnan(sub_nuc_enrichment) and sub_nuc_enrichment < 1.2
    nonsig_ks = math.isnan(ks_pval) or ks_pval > 0.05

    if enriched and sig_ks and not is_ch_gene:
        return "TUMOR-LIKE", (
            f"Sub-nucleosomal enrichment={sub_nuc_enrichment:.2f} (>1.3), "
            f"KS p={ks_pval:.2e} (<0.05), non-CH gene."
        )
    if is_ch_gene and low_enrich and nonsig_ks:
        return "CH-LIKE", (
            f"Known CH gene ({hugo}), enrichment={sub_nuc_enrichment:.2f} (<1.2), "
            f"KS p={'NA' if math.isnan(ks_pval) else f'{ks_pval:.2e}'} (>0.05)."
        )
    return "AMBIGUOUS", "Mixed signals — does not clearly fit TUMOR-LIKE or CH-LIKE criteria."


def _safe_float(v: Any) -> float:
    """Convert value to float, returning NaN for non-numeric."""
    try:
        f = float(v)
        return f
    except (ValueError, TypeError):
        return float("nan")


# ── KDE Helper (pure Python, no scipy) ───────────────────────────────────────
def _compute_kde(
    sizes: list[int],
    grid_start: int = 50,
    grid_end: int = 500,
    grid_points: int = 200,
    bandwidth: float | None = None,
) -> tuple[list[float], list[float]]:
    """Compute a Gaussian KDE density estimate over a fixed grid.

    Uses Silverman's rule-of-thumb bandwidth when not specified.
    Returns (x_values, density_values) suitable for a Plotly scatter trace.
    Returns empty lists if fewer than 2 data points.

    Args:
        sizes: Raw fragment sizes (integers).
        grid_start: Left edge of the evaluation grid (bp).
        grid_end: Right edge of the evaluation grid (bp).
        grid_points: Number of evenly spaced evaluation points.
        bandwidth: KDE bandwidth (bp). If None, uses Silverman's rule.

    Returns:
        Tuple of (x_grid, density) as lists of floats.
    """
    n = len(sizes)
    if n < 2:
        logger.debug("KDE skipped: fewer than 2 data points (n=%d).", n)
        return [], []

    # Silverman's rule-of-thumb: h = 0.9 * min(σ, IQR/1.34) * n^(-1/5)
    if bandwidth is None:
        mean = sum(sizes) / n
        var = sum((x - mean) ** 2 for x in sizes) / n
        std = var**0.5
        sorted_s = sorted(sizes)
        q1 = sorted_s[n // 4]
        q3 = sorted_s[(3 * n) // 4]
        iqr = q3 - q1
        spread = min(std, iqr / 1.34) if iqr > 0 else std
        bandwidth = 0.9 * spread * (n**-0.2)
        # Guard against degenerate bandwidth
        if bandwidth < 1.0:
            bandwidth = 5.0
            logger.debug("KDE bandwidth clamped to 5.0 (degenerate spread).")

    # Build evaluation grid
    step = (grid_end - grid_start) / (grid_points - 1)
    x_grid = [grid_start + i * step for i in range(grid_points)]

    # Evaluate Gaussian kernel at each grid point
    inv_bw = 1.0 / bandwidth
    norm_factor = 1.0 / (n * bandwidth * (2.0 * math.pi) ** 0.5)
    density = []
    for x in x_grid:
        total = 0.0
        for s in sizes:
            z = (x - s) * inv_bw
            total += math.exp(-0.5 * z * z)
        density.append(total * norm_factor)

    return x_grid, density


def generate_mfsd_report(
    parquet_path: Path,
    maf_path: Path,
    output_path: Path,
    min_alt: int = 3,
    max_variants: int = 20,
    sample_name: str = "",
) -> Path:
    """Generate an interactive HTML mFSD report.

    Args:
        parquet_path: Path to the .fsd.parquet file with raw fragment sizes.
        maf_path: Path to the MAF output with mFSD columns.
        output_path: Path to write the HTML report.
        min_alt: Minimum ALT fragments for inclusion.
        max_variants: Maximum variants (-1 = unlimited).
        sample_name: Sample identifier for the report header.

    Returns:
        Path to the generated HTML report.

    Raises:
        FileNotFoundError: If parquet or MAF files don't exist.
        RuntimeError: If report generation fails.
    """
    if not parquet_path.exists():
        raise FileNotFoundError(f"Parquet file not found: {parquet_path}")
    if not maf_path.exists():
        raise FileNotFoundError(f"MAF file not found: {maf_path}")

    logger.info("Generating mFSD report: parquet=%s, maf=%s", parquet_path, maf_path)

    # ── Load data ────────────────────────────────────────────────────────────
    try:
        import pyarrow.parquet as pq
    except ImportError:
        raise RuntimeError(
            "pyarrow is required for mFSD report generation. " "Install with: pip install pyarrow"
        ) from None

    table = pq.read_table(parquet_path)
    parquet_df = table.to_pydict()
    n_variants_parquet = len(parquet_df.get("chrom", []))
    logger.info("Loaded %d variants from parquet", n_variants_parquet)

    # Load MAF for Hugo_Symbol and mFSD stats
    import csv

    maf_data: dict[str, dict[str, str]] = {}
    with open(maf_path) as f:
        # Skip comment lines
        lines = [line for line in f if not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        key = f"{row.get('Chromosome', '')}:{row.get('Start_Position', '')}:{row.get('Reference_Allele', '')}:{row.get('Tumor_Seq_Allele2', '')}"
        maf_data[key] = row
    logger.info("Loaded %d variants from MAF", len(maf_data))

    # ── Build variant records ────────────────────────────────────────────────
    variants = []
    for i in range(n_variants_parquet):
        chrom = parquet_df["chrom"][i]
        pos = parquet_df["pos"][i]
        ref = parquet_df["ref"][i]
        alt = parquet_df["alt"][i]
        ref_sizes = (
            [int(x) for x in parquet_df["ref_sizes"][i]] if parquet_df["ref_sizes"][i] else []
        )
        alt_sizes = (
            [int(x) for x in parquet_df["alt_sizes"][i]] if parquet_df["alt_sizes"][i] else []
        )

        # Match to MAF row
        maf_key = f"{chrom}:{pos}:{ref}:{alt}"
        maf_row = maf_data.get(maf_key, {})
        hugo = maf_row.get("Hugo_Symbol", "")

        # Get mFSD stats from MAF
        sub_nuc_enrichment = _safe_float(maf_row.get("mfsd_sub_nuc_enrichment", "nan"))
        ks_pval = _safe_float(maf_row.get("mfsd_pval_alt_ref", "nan"))
        alt_mean = _safe_float(maf_row.get("mfsd_alt_mean", "nan"))
        ref_mean = _safe_float(maf_row.get("mfsd_ref_mean", "nan"))
        delta = _safe_float(maf_row.get("mfsd_delta_alt_ref", "nan"))
        alt_llr = _safe_float(maf_row.get("mfsd_alt_llr", "nan"))
        sub_nuc_ref = _safe_float(maf_row.get("mfsd_sub_nuc_ref_frac", "nan"))
        sub_nuc_alt = _safe_float(maf_row.get("mfsd_sub_nuc_alt_frac", "nan"))

        alt_count = len(alt_sizes)
        ref_count = len(ref_sizes)

        signal, explanation = _classify_origin(
            hugo, sub_nuc_enrichment, ks_pval, alt_count, min_alt
        )

        variants.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "hugo": hugo,
                "ref_sizes": ref_sizes,
                "alt_sizes": alt_sizes,
                "alt_count": alt_count,
                "ref_count": ref_count,
                "alt_mean": alt_mean,
                "ref_mean": ref_mean,
                "delta": delta,
                "alt_llr": alt_llr,
                "ks_pval": ks_pval,
                "sub_nuc_enrichment": sub_nuc_enrichment,
                "sub_nuc_ref": sub_nuc_ref,
                "sub_nuc_alt": sub_nuc_alt,
                "signal": signal,
                "explanation": explanation,
                "is_ch_gene": hugo.upper() in CH_GENES if hugo else False,
            }
        )

    # Filter by min_alt and sort by alt_count descending
    variants = [v for v in variants if v["alt_count"] >= min_alt]
    variants.sort(key=lambda v: v["alt_count"], reverse=True)
    if max_variants > 0:
        variants = variants[:max_variants]

    logger.info(
        "Report will include %d variants (min_alt=%d, max_variants=%d)",
        len(variants),
        min_alt,
        max_variants,
    )

    if not variants:
        logger.warning("No variants pass min_alt=%d filter. Report will be empty.", min_alt)

    # ── Generate HTML ────────────────────────────────────────────────────────
    html = _build_html(variants, sample_name, parquet_path.name, min_alt)
    output_path.write_text(html, encoding="utf-8")
    logger.info("mFSD report written: %s (%d variants)", output_path, len(variants))
    return output_path


def _fmt_val(v: float, precision: int = 2) -> str:
    """Format float for display, handling NaN."""
    if math.isnan(v):
        return "N/A"
    return f"{v:.{precision}f}"


def _fmt_pval(v: float) -> str:
    """Format p-value for display."""
    if math.isnan(v):
        return "N/A"
    if v < 0.001:
        return f"{v:.2e}"
    return f"{v:.4f}"


def _signal_badge(signal: str) -> str:
    """HTML badge for fragment origin signal."""
    colors = {
        "TUMOR-LIKE": ("#e74c3c", "#fff"),
        "CH-LIKE": ("#3498db", "#fff"),
        "AMBIGUOUS": ("#f39c12", "#000"),
        "INSUFFICIENT": ("#95a5a6", "#fff"),
    }
    bg, fg = colors.get(signal, ("#95a5a6", "#fff"))
    return f'<span class="badge" style="background:{bg};color:{fg}">{signal}</span>'


def _build_variant_nav_html(n_variants: int) -> str:
    """Build the sticky variant navigator bar HTML.

    Only rendered when there are 2 or more variants. For single-variant
    or empty reports, returns an empty string to avoid unnecessary UI chrome.

    Args:
        n_variants: Total number of variants in the report.

    Returns:
        HTML string for the navigator bar, or empty string.
    """
    if n_variants < 2:
        logger.debug("Variant navigator skipped: only %d variant(s).", n_variants)
        return ""

    logger.debug("Rendering variant navigator for %d variants.", n_variants)
    return """
  <div class="variant-nav" id="variant-nav">
    <div class="nav-controls">
      <button class="nav-btn" id="nav-prev" title="Previous variant (←)">← Prev</button>
      <select class="nav-dropdown" id="nav-dropdown" aria-label="Select variant"></select>
      <button class="nav-btn" id="nav-next" title="Next variant (→)">Next →</button>
      <span class="nav-counter" id="nav-counter"></span>
    </div>
    <button class="nav-mode-toggle" id="nav-mode-toggle" title="Toggle between Show All and Focus mode">Focus</button>
  </div>"""


def _build_html(variants: list[dict], sample_name: str, parquet_name: str, min_alt: int) -> str:
    """Build the complete HTML report string."""

    # Build per-variant cards with embedded Plotly data
    cards_html = []
    plotly_calls = []
    for idx, v in enumerate(variants):
        div_id = f"hist_{idx}"
        label = f"{v['hugo']} " if v["hugo"] else ""
        label += f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}"

        ch_tag = ' <span class="ch-tag">CH Gene</span>' if v["is_ch_gene"] else ""

        # Include data attributes for navigator dropdown labeling
        nav_label = f"{label} [{v['signal']}]"
        card = f"""
    <div class="card" id="variant-{idx}" data-variant-index="{idx}" data-variant-label="{nav_label}">
      <div class="card-header">
        <div class="card-title">{label}{ch_tag}</div>
        <div>{_signal_badge(v['signal'])}</div>
      </div>
      <div class="stats-row">
        <div class="stat"><span class="stat-label">REF n</span><span class="stat-value">{v['ref_count']}</span></div>
        <div class="stat"><span class="stat-label">ALT n</span><span class="stat-value">{v['alt_count']}</span></div>
        <div class="stat"><span class="stat-label">REF mean</span><span class="stat-value">{_fmt_val(v['ref_mean'], 1)} bp</span></div>
        <div class="stat"><span class="stat-label">ALT mean</span><span class="stat-value">{_fmt_val(v['alt_mean'], 1)} bp</span></div>
        <div class="stat"><span class="stat-label">Δ(ALT−REF)</span><span class="stat-value">{_fmt_val(v['delta'], 1)} bp</span></div>
        <div class="stat"><span class="stat-label">KS p</span><span class="stat-value">{_fmt_pval(v['ks_pval'])}</span></div>
        <div class="stat"><span class="stat-label">LLR</span><span class="stat-value">{_fmt_val(v['alt_llr'])}</span></div>
        <div class="stat"><span class="stat-label">Sub-nuc enrich.</span><span class="stat-value">{_fmt_val(v['sub_nuc_enrichment'])}</span></div>
      </div>
      <div class="plot-container" id="{div_id}"></div>
      <div class="interpretation">
        <strong>Fragment Origin Signal:</strong> {v['signal']} — {v['explanation']}
      </div>
    </div>"""
        cards_html.append(card)

        # ── Build Plotly traces: histogram bars + KDE density overlays ────
        ref_data = json.dumps(v["ref_sizes"])
        alt_data = json.dumps(v["alt_sizes"])

        # Compute KDE density curves (Python-side; embedded as JSON arrays)
        ref_kde_x, ref_kde_y = _compute_kde(v["ref_sizes"])
        alt_kde_x, alt_kde_y = _compute_kde(v["alt_sizes"])
        ref_kde_x_json = json.dumps(ref_kde_x)
        ref_kde_y_json = json.dumps(ref_kde_y)
        alt_kde_x_json = json.dumps(alt_kde_x)
        alt_kde_y_json = json.dumps(alt_kde_y)

        # KDE traces are on secondary y-axis (yaxis='y2') so they don't
        # fight with histogram count scale. Only added if ≥2 data points.
        kde_traces = ""
        if ref_kde_x:
            kde_traces += f""",
      {{x: {ref_kde_x_json}, y: {ref_kde_y_json}, type: 'scatter', mode: 'lines',
       name: 'REF density', yaxis: 'y2',
       line: {{color: currentTheme === 'dark' ? '#90CAF9' : '#0D47A1', width: 2}},
       showlegend: true}}"""
        if alt_kde_x:
            kde_traces += f""",
      {{x: {alt_kde_x_json}, y: {alt_kde_y_json}, type: 'scatter', mode: 'lines',
       name: 'ALT density', yaxis: 'y2',
       line: {{color: currentTheme === 'dark' ? '#FFCC80' : '#BF360C', width: 2, dash: 'dash'}},
       showlegend: true}}"""

        plotly_calls.append(f"""
    Plotly.newPlot('{div_id}', [
      {{x: {ref_data}, type: 'histogram', name: 'REF', opacity: 0.55,
       marker: {{color: currentTheme === 'dark' ? '#64B5F6' : '#1565C0'}},
       xbins: {{start: 50, end: 500, size: 5}} }},
      {{x: {alt_data}, type: 'histogram', name: 'ALT', opacity: 0.55,
       marker: {{color: currentTheme === 'dark' ? '#FFB74D' : '#E65100'}},
       xbins: {{start: 50, end: 500, size: 5}} }}{kde_traces}
    ], {{
      barmode: 'overlay',
      title: {{text: '{label}', font: {{size: 14}}}},
      xaxis: {{title: 'Fragment Size (bp)', range: [50, 450],
               zeroline: false, mirror: false}},
      yaxis: {{title: 'Count', side: 'left',
               zeroline: false, rangemode: 'tozero'}},
      yaxis2: {{title: 'Density', side: 'right', overlaying: 'y',
               showgrid: false, zeroline: false, showline: false,
               rangemode: 'tozero',
               tickfont: {{size: 10}}, titlefont: {{size: 11}} }},
      legend: {{x: 0.75, y: 0.95}},
      margin: {{t: 40, b: 50, l: 55, r: 55}},
      paper_bgcolor: 'rgba(0,0,0,0)',
      plot_bgcolor: 'rgba(0,0,0,0)',
      font: {{color: currentTheme === 'dark' ? '#e0e0e0' : '#333', family: 'Inter, sans-serif'}}
    }}, {{responsive: true}});""")

    # Summary stats
    n_tumor = sum(1 for v in variants if v["signal"] == "TUMOR-LIKE")
    n_ch = sum(1 for v in variants if v["signal"] == "CH-LIKE")
    n_ambig = sum(1 for v in variants if v["signal"] == "AMBIGUOUS")
    n_insuff = sum(1 for v in variants if v["signal"] == "INSUFFICIENT")

    ch_gene_list = ", ".join(sorted(CH_GENES))

    return f"""<!DOCTYPE html>
<html lang="en" data-theme="dark">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>mFSD Report — {sample_name or parquet_name}</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
:root {{
  --bg-page: #0f1117; --bg-card: #1a1d27; --bg-card-hover: #22263a;
  --text-primary: #e8eaed; --text-secondary: #9aa0a6;
  --accent: #64B5F6; --accent-warm: #FFB74D;
  --border: #2d3140; --radius: 12px;
}}
[data-theme="light"] {{
  --bg-page: #f5f6fa; --bg-card: #ffffff; --bg-card-hover: #f0f2f8;
  --text-primary: #1a1d27; --text-secondary: #5f6368;
  --accent: #1565C0; --accent-warm: #E65100;
  --border: #dadce0;
}}
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: 'Inter', sans-serif; background: var(--bg-page); color: var(--text-primary); line-height: 1.6; }}
.hero {{
  background: linear-gradient(135deg, #1a237e 0%, #0d47a1 50%, #01579b 100%);
  padding: 32px 40px; color: #fff;
}}
[data-theme="light"] .hero {{ background: linear-gradient(135deg, #42a5f5 0%, #1e88e5 50%, #1565c0 100%); }}
.hero h1 {{ font-size: 1.8rem; font-weight: 700; margin-bottom: 4px; }}
.hero .subtitle {{ font-size: 0.95rem; opacity: 0.85; font-weight: 300; }}
.hero .badges {{ display: flex; gap: 12px; margin-top: 12px; flex-wrap: wrap; }}
.hero .metric-badge {{
  background: rgba(255,255,255,0.15); border-radius: 8px; padding: 6px 14px;
  font-size: 0.85rem; backdrop-filter: blur(4px);
}}
.metric-badge strong {{ font-weight: 600; }}
.toolbar {{
  display: flex; justify-content: space-between; align-items: center;
  padding: 12px 40px; background: var(--bg-card); border-bottom: 1px solid var(--border);
}}
.theme-toggle {{
  background: var(--border); border: none; color: var(--text-primary);
  padding: 6px 14px; border-radius: 6px; cursor: pointer; font-family: inherit; font-size: 0.85rem;
}}
.theme-toggle:hover {{ opacity: 0.8; }}
.container {{ max-width: 1200px; margin: 0 auto; padding: 24px 20px; }}
.card {{
  background: var(--bg-card); border: 1px solid var(--border); border-radius: var(--radius);
  padding: 20px 24px; margin-bottom: 20px; transition: background 0.2s;
}}
.card:hover {{ background: var(--bg-card-hover); }}
.card-header {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px; }}
.card-title {{ font-size: 1.1rem; font-weight: 600; font-family: 'JetBrains Mono', monospace; }}
.badge {{
  display: inline-block; padding: 3px 10px; border-radius: 4px;
  font-size: 0.75rem; font-weight: 600; letter-spacing: 0.5px;
}}
.ch-tag {{
  background: #3498db; color: #fff; padding: 2px 8px; border-radius: 4px;
  font-size: 0.7rem; font-weight: 600; margin-left: 8px; vertical-align: middle;
}}
.stats-row {{ display: flex; gap: 16px; flex-wrap: wrap; margin-bottom: 16px; }}
.stat {{ display: flex; flex-direction: column; min-width: 90px; }}
.stat-label {{ font-size: 0.72rem; color: var(--text-secondary); text-transform: uppercase; letter-spacing: 0.5px; }}
.stat-value {{ font-size: 0.95rem; font-weight: 500; font-family: 'JetBrains Mono', monospace; }}
.plot-container {{ width: 100%; height: 340px; margin-bottom: 12px; }}
.interpretation {{
  font-size: 0.85rem; color: var(--text-secondary); padding: 10px 14px;
  background: rgba(100,181,246,0.05); border-radius: 8px; border-left: 3px solid var(--accent);
}}
.section-title {{ font-size: 1.3rem; font-weight: 600; margin: 28px 0 16px; }}
.caveat {{
  background: rgba(243,156,18,0.1); border: 1px solid rgba(243,156,18,0.3);
  border-radius: 8px; padding: 16px 20px; margin: 20px 0; font-size: 0.88rem;
}}
.caveat strong {{ color: #f39c12; }}
.methodology {{
  background: var(--bg-card); border: 1px solid var(--border); border-radius: var(--radius);
  padding: 20px 24px; margin: 20px 0; font-size: 0.88rem;
}}
.methodology h3 {{ margin-bottom: 8px; font-size: 1rem; }}
.methodology code {{ font-family: 'JetBrains Mono', monospace; font-size: 0.82rem; background: rgba(100,181,246,0.1); padding: 1px 5px; border-radius: 3px; }}
.gene-list {{ font-family: 'JetBrains Mono', monospace; font-size: 0.78rem; word-break: break-word; color: var(--text-secondary); }}
.summary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin: 16px 0; }}
.summary-card {{
  background: var(--bg-card); border: 1px solid var(--border); border-radius: 8px;
  padding: 16px; text-align: center;
}}
.summary-card .big {{ font-size: 2rem; font-weight: 700; font-family: 'JetBrains Mono', monospace; }}
.summary-card .label {{ font-size: 0.78rem; color: var(--text-secondary); text-transform: uppercase; }}
/* ── Variant Navigator ───────────────────────────────────────────────────── */
.variant-nav {{
  display: flex; justify-content: space-between; align-items: center;
  background: var(--bg-card); border: 1px solid var(--border); border-radius: var(--radius);
  padding: 10px 16px; margin: 16px 0; position: sticky; top: 0; z-index: 100;
  backdrop-filter: blur(8px); -webkit-backdrop-filter: blur(8px);
}}
.nav-controls {{ display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }}
.nav-btn {{
  background: var(--border); border: none; color: var(--text-primary);
  padding: 6px 14px; border-radius: 6px; cursor: pointer; font-family: inherit;
  font-size: 0.82rem; font-weight: 500; transition: background 0.15s, opacity 0.15s;
}}
.nav-btn:hover {{ opacity: 0.8; }}
.nav-btn:disabled {{ opacity: 0.35; cursor: not-allowed; }}
.nav-dropdown {{
  background: var(--bg-page); color: var(--text-primary); border: 1px solid var(--border);
  border-radius: 6px; padding: 6px 10px; font-family: 'JetBrains Mono', monospace;
  font-size: 0.78rem; max-width: 360px; cursor: pointer;
}}
.nav-counter {{
  font-size: 0.82rem; color: var(--text-secondary); font-family: 'JetBrains Mono', monospace;
  min-width: 60px; text-align: center;
}}
.nav-mode-toggle {{
  background: var(--accent); border: none; color: #fff;
  padding: 6px 16px; border-radius: 6px; cursor: pointer; font-family: inherit;
  font-size: 0.82rem; font-weight: 600; transition: background 0.15s;
}}
.nav-mode-toggle:hover {{ opacity: 0.85; }}
.nav-mode-toggle.active {{ background: #e74c3c; }}
/* Card states for Focus mode */
.card.dimmed {{
  opacity: 0.15; transform: scale(0.98); pointer-events: none;
  transition: opacity 0.3s, transform 0.3s;
}}
.card.active-card {{
  border-color: var(--accent); box-shadow: 0 0 0 2px rgba(100,181,246,0.3);
  transition: border-color 0.3s, box-shadow 0.3s;
}}
.footer {{
  text-align: center; padding: 24px 20px; margin-top: 32px;
  border-top: 1px solid var(--border); font-size: 0.82rem;
  color: var(--text-secondary);
}}
.footer a {{
  color: var(--accent); text-decoration: none; font-weight: 500;
}}
.footer a:hover {{ text-decoration: underline; }}
@media print {{
  .hero {{ print-color-adjust: exact; -webkit-print-color-adjust: exact; }}
  .theme-toggle, .toolbar, .variant-nav {{ display: none; }}
  .card {{ page-break-inside: avoid; opacity: 1 !important; transform: none !important; pointer-events: auto !important; }}
  .card.dimmed {{ opacity: 1 !important; transform: none !important; pointer-events: auto !important; }}
}}
</style>
</head>
<body>

<div class="hero">
  <h1>mFSD Per-Variant Report</h1>
  <div class="subtitle">{sample_name or parquet_name} — Fragment Size Distribution Analysis</div>
  <div class="badges">
    <span class="metric-badge"><strong>{len(variants)}</strong> variants</span>
    <span class="metric-badge"><strong>{n_tumor}</strong> tumor-like</span>
    <span class="metric-badge"><strong>{n_ch}</strong> CH-like</span>
    <span class="metric-badge"><strong>{n_ambig}</strong> ambiguous</span>
    <span class="metric-badge"><strong>{n_insuff}</strong> insufficient</span>
    <span class="metric-badge">min ALT ≥ <strong>{min_alt}</strong></span>
  </div>
</div>

<div class="toolbar">
  <span style="font-size:0.85rem;color:var(--text-secondary)">Generated by gbcms mFSD</span>
  <button class="theme-toggle" onclick="toggleTheme()">🌓 Toggle Theme</button>
</div>

<div class="container">

  <div class="summary-grid">
    <div class="summary-card"><div class="big" style="color:#e74c3c">{n_tumor}</div><div class="label">Tumor-like</div></div>
    <div class="summary-card"><div class="big" style="color:#3498db">{n_ch}</div><div class="label">CH-like</div></div>
    <div class="summary-card"><div class="big" style="color:#f39c12">{n_ambig}</div><div class="label">Ambiguous</div></div>
    <div class="summary-card"><div class="big" style="color:#95a5a6">{n_insuff}</div><div class="label">Insufficient</div></div>
  </div>

  <div class="caveat">
    <strong>⚠️ Important:</strong> Fragment size alone cannot definitively distinguish CH from ctDNA.
    Some genes (e.g., TP53) can be both CH-driven and tumor-driven. Paired WBC sequencing remains
    the gold standard for CH exclusion. The Fragment Origin Signal is <strong>suggestive, not diagnostic</strong>,
    and should be interpreted in clinical context.
  </div>

  {_build_variant_nav_html(len(variants))}

  <h2 class="section-title">Per-Variant Analysis</h2>
  {''.join(cards_html) if cards_html else '<p style="color:var(--text-secondary)">No variants passed the minimum ALT fragment filter.</p>'}

  <h2 class="section-title">Methodology</h2>
  <div class="methodology">
    <h3>Fragment Origin Signal Classification</h3>
    <p>Each variant is classified based on three signals:</p>
    <ul style="margin:8px 0 8px 20px">
      <li><strong>TUMOR-LIKE:</strong> Sub-nucleosomal enrichment &gt;1.3, KS p&lt;0.05, and <em>not</em> in the CH gene set.</li>
      <li><strong>CH-LIKE:</strong> Known CH gene, enrichment &lt;1.2, and KS p&gt;0.05 (ALT distribution mirrors REF).</li>
      <li><strong>AMBIGUOUS:</strong> Mixed criteria — does not clearly fit either category.</li>
      <li><strong>INSUFFICIENT:</strong> ALT fragment count below threshold (<code>min_alt={min_alt}</code>).</li>
    </ul>
    <h3>Sub-nucleosomal Enrichment</h3>
    <p>Ratio of ALT fragments &lt;150bp to REF fragments &lt;150bp. ctDNA tends to show enrichment
    (ratio &gt;1.0) due to tumor-derived fragments being shorter. CH mirrors background cfDNA.</p>
    <h3>CH Gene Set ({len(CH_GENES)} genes)</h3>
    <p class="gene-list">{ch_gene_list}</p>
    <p style="margin-top:8px;font-size:0.82rem;color:var(--text-secondary)">
    Sources: Steensma et al. (2015), Jaiswal et al. (2014), Bolton et al. (2020).</p>
  </div>

  <div class="methodology">
    <h3>Glossary</h3>
    <ul style="margin:8px 0 8px 20px;font-size:0.85rem">
      <li><strong>REF n / ALT n:</strong> Fragment count classified as reference or alternate allele.</li>
      <li><strong>Δ(ALT−REF):</strong> Difference in mean fragment size (bp). Negative = ALT shorter.</li>
      <li><strong>KS p:</strong> Kolmogorov-Smirnov test p-value comparing ALT vs REF size distributions.</li>
      <li><strong>LLR:</strong> Log-likelihood ratio (positive = tumor-like fragment profile).</li>
      <li><strong>Sub-nuc enrich.:</strong> ALT(&lt;150bp fraction) / REF(&lt;150bp fraction).</li>
      <li><strong>Density overlay:</strong> Gaussian KDE (Silverman bandwidth) normalized to unit area, plotted on the right y-axis. Enables direct shape comparison between REF and ALT regardless of fragment count differences.</li>
    </ul>
  </div>
</div>

<footer class="footer">
  Made by <a href="https://github.com/rhshah" target="_blank" rel="noopener noreferrer">Ronak Shah (@rhshah)</a> using Antigravity
</footer>

<script>
let currentTheme = localStorage.getItem('mfsd-theme') ||
  (window.matchMedia('(prefers-color-scheme: light)').matches ? 'light' : 'dark');
document.documentElement.setAttribute('data-theme', currentTheme);

function toggleTheme() {{
  currentTheme = currentTheme === 'dark' ? 'light' : 'dark';
  document.documentElement.setAttribute('data-theme', currentTheme);
  localStorage.setItem('mfsd-theme', currentTheme);
  renderPlots();
}}

function renderPlots() {{
{''.join(plotly_calls)}
}}
renderPlots();

/* ── Variant Navigator Controller ────────────────────────────────────────── */
(function() {{
  'use strict';
  var nav = document.getElementById('variant-nav');
  if (!nav) return;  // No navigator for 0-1 variants

  var cards = document.querySelectorAll('.card[data-variant-index]');
  var total = cards.length;
  if (total < 2) return;  // Safety: should not happen since nav is only rendered for ≥2

  var dropdown = document.getElementById('nav-dropdown');
  var counter  = document.getElementById('nav-counter');
  var prevBtn  = document.getElementById('nav-prev');
  var nextBtn  = document.getElementById('nav-next');
  var modeBtn  = document.getElementById('nav-mode-toggle');

  var currentIndex = 0;
  var focusMode = false;  // 'Show All' is the default

  // Populate dropdown from card data attributes
  for (var i = 0; i < total; i++) {{
    var opt = document.createElement('option');
    opt.value = i;
    opt.textContent = (i + 1) + '. ' + cards[i].getAttribute('data-variant-label');
    dropdown.appendChild(opt);
  }}

  function updateUI() {{
    // Update counter and dropdown selection
    counter.textContent = (currentIndex + 1) + ' of ' + total;
    dropdown.value = currentIndex;

    // Update prev/next button states (wrap-around, so never truly disabled)
    prevBtn.disabled = false;
    nextBtn.disabled = false;

    // Apply card visibility based on mode
    for (var i = 0; i < total; i++) {{
      cards[i].classList.remove('active-card', 'dimmed');
      if (focusMode) {{
        if (i === currentIndex) {{
          cards[i].classList.add('active-card');
        }} else {{
          cards[i].classList.add('dimmed');
        }}
      }} else {{
        // Show All mode — only highlight active
        if (i === currentIndex) {{
          cards[i].classList.add('active-card');
        }}
      }}
    }}

    // Update mode button text and style
    if (focusMode) {{
      modeBtn.textContent = 'Show All';
      modeBtn.classList.add('active');
    }} else {{
      modeBtn.textContent = 'Focus';
      modeBtn.classList.remove('active');
    }}
  }}

  function goTo(index) {{
    currentIndex = ((index % total) + total) % total;  // Wrap-around
    updateUI();
    // Smooth scroll to the active card
    cards[currentIndex].scrollIntoView({{ behavior: 'smooth', block: 'center' }});
    // Trigger Plotly resize for the newly visible plot (deferred to after scroll)
    setTimeout(function() {{
      var plotDiv = cards[currentIndex].querySelector('.plot-container');
      if (plotDiv && typeof Plotly !== 'undefined') {{
        Plotly.Plots.resize(plotDiv);
      }}
    }}, 400);
  }}

  // Event listeners
  prevBtn.addEventListener('click', function() {{ goTo(currentIndex - 1); }});
  nextBtn.addEventListener('click', function() {{ goTo(currentIndex + 1); }});
  dropdown.addEventListener('change', function() {{ goTo(parseInt(this.value, 10)); }});
  modeBtn.addEventListener('click', function() {{
    focusMode = !focusMode;
    updateUI();
    // Scroll to current card when entering Focus mode
    if (focusMode) {{
      cards[currentIndex].scrollIntoView({{ behavior: 'smooth', block: 'center' }});
    }}
  }});

  // Keyboard navigation: ← and → arrow keys
  document.addEventListener('keydown', function(e) {{
    // Don't intercept when typing in an input/select
    if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') return;
    if (e.key === 'ArrowLeft')  {{ e.preventDefault(); goTo(currentIndex - 1); }}
    if (e.key === 'ArrowRight') {{ e.preventDefault(); goTo(currentIndex + 1); }}
  }});

  // Initial state: highlight first card, no dimming (Show All default)
  updateUI();
}})();
</script>
</body>
</html>"""
