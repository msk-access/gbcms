#!/bin/bash
# Test Docker build and functionality

set -e

echo "=========================================="
echo "Docker Build and Test Script"
echo "=========================================="
echo

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_DIR"

# Test 1: Build production image
echo "Test 1: Building production Docker image..."
if docker build -t gbcms:test-build . ; then
    echo -e "${GREEN}✅ Production image built successfully${NC}"
else
    echo -e "${RED}❌ Failed to build production image${NC}"
    exit 1
fi

echo

# Test 2: Verify installation
echo "Test 2: Verifying installation..."
if docker run --rm gbcms:test-build version ; then
    echo -e "${GREEN}✅ Installation verified${NC}"
else
    echo -e "${RED}❌ Installation verification failed${NC}"
    exit 1
fi

echo

# Test 3: Check help command
echo "Test 3: Testing help command..."
if docker run --rm gbcms:test-build --help > /dev/null 2>&1 ; then
    echo -e "${GREEN}✅ Help command works${NC}"
else
    echo -e "${RED}❌ Help command failed${NC}"
    exit 1
fi

echo



echo

# Test 5: Check Ray availability
echo "Test 5: Checking Ray availability..."
if docker run --rm gbcms:test-build python -c "import ray; print('Ray version:', ray.__version__)" ; then
    echo -e "${GREEN}✅ Ray is available${NC}"
else
    echo -e "${YELLOW}⚠️  Ray not available (optional)${NC}"
fi

echo

# Test 6: Check Numba availability
echo "Test 6: Checking Numba availability..."
if docker run --rm gbcms:test-build python -c "import numba; print('Numba version:', numba.__version__)" ; then
    echo -e "${GREEN}✅ Numba is available${NC}"
else
    echo -e "${RED}❌ Numba not available (required)${NC}"
    exit 1
fi

echo

# Test 7: Check samtools
echo "Test 7: Checking samtools..."
if docker run --rm gbcms:test-build samtools --version | head -1 ; then
    echo -e "${GREEN}✅ samtools is available${NC}"
else
    echo -e "${RED}❌ samtools not available${NC}"
    exit 1
fi

echo

# Test 8: Build test image
echo "Test 8: Building test Docker image..."
if docker build -f Dockerfile.test -t gbcms:test-image . ; then
    echo -e "${GREEN}✅ Test image built successfully${NC}"
else
    echo -e "${RED}❌ Failed to build test image${NC}"
    exit 1
fi

echo

# Test 9: Check image sizes
echo "Test 9: Checking image sizes..."
echo "Production image:"
docker images gbcms:test-build --format "Size: {{.Size}}"
echo "Test image:"
docker images gbcms:test-image --format "Size: {{.Size}}"
echo -e "${GREEN}✅ Image sizes checked${NC}"

echo

# Test 10: Test with docker-compose
echo "Test 10: Testing docker-compose..."
if command -v docker-compose &> /dev/null; then
    if docker-compose build > /dev/null 2>&1 ; then
        echo -e "${GREEN}✅ docker-compose build works${NC}"
    else
        echo -e "${YELLOW}⚠️  docker-compose build failed${NC}"
    fi
else
    echo -e "${YELLOW}⚠️  docker-compose not installed${NC}"
fi

echo

# Summary
echo "=========================================="
echo "Summary"
echo "=========================================="
echo -e "${GREEN}✅ All critical tests passed!${NC}"
echo
echo "Docker images are ready to use:"
echo "  docker run --rm gbcms:test-build version"
echo "  docker run --rm gbcms:test-build --help"
echo
echo "To tag and use:"
echo "  docker tag gbcms:test-build gbcms:latest"
echo "  docker run --rm -v \$(pwd)/data:/data gbcms:latest dna ..."
echo

# Cleanup option
read -p "Clean up test images? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    docker rmi gbcms:test-build gbcms:test-image
    echo -e "${GREEN}✅ Test images removed${NC}"
fi

echo
echo "=========================================="
echo "Docker test complete!"
echo "=========================================="
