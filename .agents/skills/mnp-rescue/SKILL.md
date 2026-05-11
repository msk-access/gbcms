# MNP Rescue & Diagnostic Flags

## MNP Rescue Threshold

`--rescue-mnp-threshold` (0.0–1.0) controls MNP candidacy for the rescue pass.

### Boundary Semantics
- **1.0** (default): Permissive — rescue all MNPs. C++ compatible.
- **0.5**: Conservative — rescue only sparse MNPs (≤50% discriminating positions).
- **0.0**: Disabled — no MNP rescue eligibility. Still emits `MNP_DISC_RATIO` diagnostics.

### Logic (`pipeline.py`)
```python
ratio = disc_positions / mnp_length  # MNP_DISC_RATIO
eligible = ratio <= threshold         # MNP_RESCUE_ELIGIBLE
```
For `threshold=0.0`: ratio is always >0 for valid MNPs, so no MNPs are eligible.

## Diagnostic Flags

`gbcms_diagnostic` and `gbcms_rescue` are **strongly-typed Rust fields** on `PreparedVariant`, not dynamic Python attributes.

Format: semicolon-separated key=value strings:
```
MNP_DISC_RATIO=0.333;MNP_RESCUE_ELIGIBLE=true
```

## Key Files
- `pipeline.py`: `_compute_diagnostics()`, threshold gating logic
- `rust/src/types.rs`: `gbcms_diagnostic`, `gbcms_rescue` field definitions
- `src/gbcms/_rs.pyi`: type stubs (must include these fields)
- `tests/test_diagnostic_flags.py`: threshold-based gating tests
- `tests/test_rescue_mnp.py`: MNP rescue pass validation
