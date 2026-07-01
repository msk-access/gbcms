# Contributing to gbcms

Thank you for your interest in contributing to GetBaseCounts! This document provides guidelines and instructions for contributing.

## Development Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/msk-access/gbcms.git
   cd gbcms
   ```

2. **Install uv** (if not already installed)
   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

3. **Install development dependencies**
   ```bash
   # Install Rust (if not installed)
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   
   # Build and install in development mode
   pip install maturin
   maturin develop --release
   ```



## Development Workflow

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=gbcms --cov-report=html

# Run specific test file
pytest tests/test_accuracy.py

# Run with verbose output
pytest -v
```

### Code Quality

```bash
# Format code
black src/ tests/

# Lint code
ruff check src/ tests/

# Fix linting issues automatically
ruff check --fix src/ tests/

# Type checking
mypy src/
```

### Building and Testing Docker

```bash
# Build Docker image
docker build -t gbcms:latest .
```

## Code Style

- Follow PEP 8 style guidelines
- Use type hints for all function signatures
- Write docstrings for all public functions and classes
- Keep functions focused and single-purpose
- Maximum line length: 100 characters

## Testing Guidelines

- Write tests for all new features
- Maintain or improve code coverage
- Use descriptive test names
- Use fixtures for common test setup
- Test edge cases and error conditions

### Counting Engine Changes

Changes to the counting engine (`rust/src/counting/` — especially `engine.rs` and
`variant_checks.rs`) require additional verification:

- **Unit tests**: `cargo test` must pass (all existing tests + new tests for the change)
- **BAM-slice regression**: Run the BAM-slice regression suite and verify no unintended
  count shifts (see the developer guide for the current slice set)
- **Binned↔legacy parity**: any change to read classification / filtering / fragment
  consensus / fetch-window logic must be mirrored in the legacy `count_single_variant`,
  or the parity tests fail (see `.agents/rules/architecture.md`)
- **Variant-type coverage**: If modifying a specific variant type (SNP/MNP/Indel/Complex),
  ensure the regression MAF includes representative variants of that type
- **MnpResult paths**: Any changes to MNP handling must test all 5 `MnpResult` variants:
  `Ref`, `Alt`, `LowQuality`, `ThirdAllele`, `Structural`

## Pull Request Process

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**
   - Write code following the style guidelines
   - Add tests for new functionality
   - Update documentation as needed

3. **Run tests and linters** — the Python *and* Rust gates (`make test`/`make lint` cover
   Python only; the Rust lint/tests must be run separately, or you'll ship un-linted Rust):
   ```bash
   make test && make lint                                   # Python: pytest + ruff/black/mypy
   cd rust && cargo clippy --all-targets -- -D warnings && cargo test && cd ..
   ```
   The load-bearing invariants (binned↔legacy parity, one base-quality gate across
   alignment backends, the 4 counting-test invariants) are in
   [`.agents/rules/architecture.md`](.agents/rules/architecture.md) — read them before
   touching the engine.

4. **Commit your changes**
   ```bash
   git add .
   git commit -m "Description of your changes"
   ```

5. **Push to your fork**
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Create a Pull Request**
   - Provide a clear description of the changes
   - Reference any related issues
   - Ensure all CI checks pass

## Reporting Issues

When reporting issues, please include:

- Python version
- Operating system
- Steps to reproduce
- Expected behavior
- Actual behavior
- Error messages or logs

## Feature Requests

We welcome feature requests! Please:

- Check if the feature already exists
- Provide a clear use case
- Describe the expected behavior
- Consider submitting a PR if you can implement it

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Focus on constructive feedback
- Help others learn and grow

## Questions?

Feel free to open an issue for questions or reach out to the maintainers.

Thank you for contributing!
