import subprocess
import time
from pathlib import Path

# Configuration
VCF_INPUT = Path("~/Downloads/C-U1MJT8-L001-d.DONOR22-TP.combined-variants_anno.vcf").expanduser()
REF_FASTA = Path(
    "/Users/shahr2/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Downloads/03Jan2023/nucleo_qc_generation/test_data/Homo_sapiens_assembly19.fasta"
)
OUT_DIR = Path("~/Downloads/gbcms_benchmarks").expanduser()
BAM_LIST = OUT_DIR / "bams.txt"

BAMS = [
    (
        "C-6WTKCL-L001-d",
        "/Users/shahr2/Downloads/C-6WTKCL-L001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam",
    ),
    (
        "C-6WTKCL-N001-d",
        "/Users/shahr2/Downloads/C-6WTKCL-N001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam",
    ),
]


def setup():
    if not OUT_DIR.exists():
        OUT_DIR.mkdir(parents=True)

    with open(BAM_LIST, "w") as f:
        for name, path in BAMS:
            f.write(f"{name}\t{path}\n")


def run_benchmark(threads, output_format):
    cmd = [
        "gbcms",
        "dna",
        "--variants",
        str(VCF_INPUT),
        "--fasta",
        str(REF_FASTA),
        "--bam-list",
        str(BAM_LIST),
        "--output-dir",
        str(OUT_DIR),
        "--format",
        output_format,
        "--suffix",
        f".bench_t{threads}_{output_format}",
    ]

    if threads != "ALL":
        cmd.extend(["--threads", str(threads)])

    print(f"Running: Threads={threads}, Format={output_format}...", end="", flush=True)

    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end_time = time.time()

    duration = end_time - start_time

    if result.returncode != 0:
        print(f" FAILED! ({duration:.2f}s)")
        print(result.stderr)
    else:
        print(f" DONE ({duration:.2f}s)")

    return duration


def main():
    setup()

    print("Starting Benchmark")
    print(f"Input: {VCF_INPUT}")
    print(f"BAMs: {len(BAMS)} samples")
    print("-" * 40)

    results = []

    # Test cases
    thread_counts = [1, 2, 3, 4, "ALL"]
    formats = ["vcf", "maf"]

    for fmt in formats:
        for t in thread_counts:
            duration = run_benchmark(t, fmt)
            results.append({"threads": t, "format": fmt, "time": duration})

    print("-" * 40)
    print("Summary:")
    print(f"{'Format':<10} {'Threads':<10} {'Time (s)':<10}")
    print("-" * 30)
    for r in results:
        print(f"{r['format']:<10} {str(r['threads']):<10} {r['time']:.2f}")


if __name__ == "__main__":
    main()
