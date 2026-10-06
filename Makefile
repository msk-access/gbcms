.PHONY: help install install-dev install-all test test-cov test-workflows lint format clean verify setup docker-build docker-run

help:
	@echo "Available commands:"
	@echo ""
	@echo "Installation:"
	@echo "  install       - Install package"
	@echo "  install-dev   - Install with dev dependencies"
	@echo "  install-all   - Install with all dependencies (including Ray)"
	@echo "  setup         - Complete setup and test"
	@echo "  verify        - Verify installation"
	@echo ""
	@echo "Testing:"
	@echo "  test          - Run unit tests"
	@echo "  test-cov      - Run tests with coverage"
	@echo "  test-workflows - Run end-to-end workflow tests"
	@echo ""
	@echo "Code Quality:"
	@echo "  lint          - Run linters"
	@echo "  format        - Format code"
	@echo "  clean         - Clean build artifacts"
	@echo ""
	@echo "Docker:"
	@echo "  docker-build  - Build Docker image"
	@echo "  docker-run    - Run example in Docker"

install:
	uv pip install .

install-dev:
	uv pip install -e . --group dev

install-all:
	uv pip install -e ".[all]" --group dev

setup:
	@chmod +x scripts/setup_and_test.sh
	@bash scripts/setup_and_test.sh

verify:
	@python3 scripts/verify_installation.py

test:
	pytest -v

test-cov:
	pytest --cov=gbcms --cov-report=html --cov-report=term-missing

test-workflows:
	@chmod +x scripts/test_vcf_workflow.sh scripts/test_maf_workflow.sh
	@echo "Running VCF workflow test..."
	@bash scripts/test_vcf_workflow.sh
	@echo ""
	@echo "Running MAF workflow test..."
	@bash scripts/test_maf_workflow.sh

lint:
	ruff check src/ tests/
	mypy src/ || true

format:
	black src/ tests/
	ruff check --fix src/ tests/

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache/
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete

docker-build:
	docker build -t gbcms:latest .
	@echo ""
	@echo "✅ Docker image built successfully!"
	@echo "Test it with: docker run --rm gbcms:latest version"

docker-run:
	@echo "Example Docker run (mount your data directory):"
	@echo "docker run -v /path/to/data:/data gbcms:latest count run \\"
	@echo "  --fasta /data/reference.fa \\"
	@echo "  --bam sample1:/data/sample1.bam \\"
	@echo "  --vcf /data/variants.vcf \\"
	@echo "  --output /data/output.txt"

docker-clean:
	docker rmi gbcms:latest gbcms:test 2>/dev/null || true
	@echo "✅ Docker images removed"
