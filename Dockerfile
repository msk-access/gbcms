# syntax=docker/dockerfile:1

# Dockerfile for gbcms
# 
# Build: docker build -t gbcms .
# Run:   docker run gbcms --help

# Stage 1: Builder (with Rust compilation)
FROM python:3.11-bookworm AS builder

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

# Install system dependencies for Rust compilation
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    git \
    cmake \
    clang \
    libclang-dev \
    llvm-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

ENV LIBCLANG_PATH="/usr/lib/llvm-14/lib"
ENV BINDGEN_EXTRA_CLANG_ARGS="-I/usr/lib/llvm-14/lib/clang/14.0.6/include"

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Cargo network resilience for CI: crates.io downloads intermittently fail with
# "[16] Error in the HTTP2 framing layer" (cargo's HTTP/2 multiplexing on a flaky
# connection). Force HTTP/1.1, retry, and fetch the git dependency (wfa2lib-rs)
# via the system git CLI rather than cargo's builtin fetcher.
ENV CARGO_NET_RETRY=10 \
    CARGO_HTTP_MULTIPLEXING=false \
    CARGO_NET_GIT_FETCH_WITH_CLI=true

# Install maturin with patchelf support
RUN pip install --no-cache-dir "maturin[patchelf]"

# Copy project files (selective copy for better layer caching)
COPY pyproject.toml README.md LICENSE ./
COPY rust/ rust/
COPY src/ src/

# Build unified wheel with maturin (includes both Python and Rust)
# Don't use --manifest-path; it's in pyproject.toml and ensures correct wheel name
# --no-default-features drops the `legacy-parity` feature so the shipped wheel does not
# export the per-variant `count_bam` parity oracle (test-only). See rust/Cargo.toml.
RUN maturin build --release --no-default-features --out /app/dist

# Stage 2: Runtime (slim image)
FROM python:3.11-slim-bookworm

# OCI Labels for GitHub Container Registry
LABEL org.opencontainers.image.title="gbcms"
LABEL org.opencontainers.image.description="Python Get Base Count Multi Sample - high-performance variant counting"
LABEL org.opencontainers.image.url="https://github.com/msk-access/gbcms"
LABEL org.opencontainers.image.source="https://github.com/msk-access/gbcms"
LABEL org.opencontainers.image.vendor="MSK-ACCESS"
LABEL org.opencontainers.image.licenses="AGPL-3.0"
LABEL org.opencontainers.image.authors="Ronak Shah <shahr2@mskcc.org>"

WORKDIR /app

# Install runtime dependencies for htslib
# procps and bash are required for Nextflow task metrics and shell execution
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4 \
    libssl3 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    procps \
    bash \
    && rm -rf /var/lib/apt/lists/*

# Copy and install the unified wheel
COPY --from=builder /app/dist/*.whl /app/dist/
RUN pip install --no-cache-dir /app/dist/*.whl

# Verify installation
RUN python -c "from gbcms import _rs; import gbcms; print(f'gbcms {gbcms.__version__} ready')"

ENTRYPOINT ["gbcms"]
CMD ["--help"]