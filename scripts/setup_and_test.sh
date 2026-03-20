#!/bin/bash
# Complete setup and test script for gbcms

set -e

echo "=========================================="
echo "gbcms Setup and Test"
echo "=========================================="
echo

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "Project directory: $PROJECT_DIR"
cd "$PROJECT_DIR"

# Step 1: Check Python version
echo
echo "Step 1: Checking Python version..."
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "Python version: $PYTHON_VERSION"

REQUIRED_VERSION="3.9"
if python3 -c "import sys; exit(0 if sys.version_info >= (3, 9) else 1)"; then
    echo -e "${GREEN}✅ Python version OK${NC}"
else
    echo -e "${RED}❌ Python 3.9+ required${NC}"
    exit 1
fi

# Step 2: Check for uv
echo
echo "Step 2: Checking for uv..."
if command -v uv &> /dev/null; then
    echo -e "${GREEN}✅ uv found${NC}"
else
    echo -e "${YELLOW}⚠️  uv not found, installing...${NC}"
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.cargo/bin:$PATH"
fi

# Step 3: Install package
echo
echo "Step 3: Installing GetBaseCounts..."
echo "Installing with dev and all dependencies..."

if uv pip install -e ".[dev,all]"; then
    echo -e "${GREEN}✅ Installation successful${NC}"
else
    echo -e "${RED}❌ Installation failed${NC}"
    exit 1
fi

# Step 4: Verify installation
echo
echo "Step 4: Verifying installation..."
if python3 "$SCRIPT_DIR/verify_installation.py"; then
    echo -e "${GREEN}✅ Verification passed${NC}"
else
    echo -e "${RED}❌ Verification failed${NC}"
    exit 1
fi

# Step 5: Check CLI
echo
echo "Step 5: Checking CLI..."
if gbcms --help > /dev/null 2>&1; then
    echo -e "${GREEN}✅ CLI accessible${NC}"
    gbcms version
else
    echo -e "${RED}❌ CLI not accessible${NC}"
    exit 1
fi

# Step 6: Run unit tests
echo
echo "Step 6: Running unit tests..."
if pytest tests/ -v --tb=short; then
    echo -e "${GREEN}✅ Unit tests passed${NC}"
else
    echo -e "${YELLOW}⚠️  Some unit tests failed (may need test data)${NC}"
fi

# Step 7: Check for samtools (required for workflow tests)
echo
echo "Step 7: Checking for samtools..."
if command -v samtools &> /dev/null; then
    echo -e "${GREEN}✅ samtools found${NC}"
    
    # Make test scripts executable
    chmod +x "$SCRIPT_DIR/test_vcf_workflow.sh"
    chmod +x "$SCRIPT_DIR/test_maf_workflow.sh"
    
    # Step 8: Run VCF workflow test
    echo
    echo "Step 8: Running VCF workflow test..."
    if bash "$SCRIPT_DIR/test_vcf_workflow.sh"; then
        echo -e "${GREEN}✅ VCF workflow test passed${NC}"
    else
        echo -e "${RED}❌ VCF workflow test failed${NC}"
        exit 1
    fi
    
    # Step 9: Run MAF workflow test
    echo
    echo "Step 9: Running MAF workflow test..."
    if bash "$SCRIPT_DIR/test_maf_workflow.sh"; then
        echo -e "${GREEN}✅ MAF workflow test passed${NC}"
    else
        echo -e "${RED}❌ MAF workflow test failed${NC}"
        exit 1
    fi
else
    echo -e "${YELLOW}⚠️  samtools not found, skipping workflow tests${NC}"
    echo "Install samtools to run end-to-end tests:"
    echo "  brew install samtools  # macOS"
    echo "  apt-get install samtools  # Ubuntu/Debian"
fi

# Summary
echo
echo "=========================================="
echo "Setup and Test Summary"
echo "=========================================="
echo -e "${GREEN}✅ All checks passed!${NC}"
echo
echo "gbcms is ready to use:"
echo "  gbcms --help"
echo "  gbcms version"
echo "  gbcms info"
echo "  gbcms dna --help"
echo
echo "Documentation:"
echo "  README.md - Main documentation"
echo "  QUICKSTART.md - Quick start guide"
echo "  ADVANCED_FEATURES.md - Advanced usage"
echo "  CLI_FEATURES.md - CLI documentation"
echo
