#!/bin/bash
# End-to-end test for VCF workflow

set -e  # Exit on error

echo "=========================================="
echo "gbcms VCF Workflow Test"
echo "=========================================="
echo

# Create test directory
TEST_DIR=$(mktemp -d)
echo "Test directory: $TEST_DIR"
cd "$TEST_DIR"

# Create test reference
echo "Creating test reference..."
cat > reference.fa << 'EOF'
>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
EOF

# Index reference
echo "Indexing reference..."
samtools faidx reference.fa

# Create test VCF
echo "Creating test VCF..."
cat > variants.vcf << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=120>
##contig=<ID=chr2,length=120>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10	.	A	T	.	PASS	.
chr1	20	.	C	G	.	PASS	.
chr1	30	.	GA	G	.	PASS	.
chr2	15	.	G	GC	.	PASS	.
EOF

# Create test BAM
echo "Creating test BAM..."
samtools view -bS > sample1.bam << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:120
@SQ	SN:chr2	LN:120
read1	99	chr1	5	60	20M	=	50	65	ATCGATCGATCGATCGATCG	IIIIIIIIIIIIIIIIIIII
read2	99	chr1	15	60	20M	=	60	65	CGATCGATCGATCGATCGAT	IIIIIIIIIIIIIIIIIIII
read3	99	chr1	25	60	20M	=	70	65	GATCGATCGATCGATCGATC	IIIIIIIIIIIIIIIIIIII
read4	99	chr2	10	60	20M	=	55	65	GCTAGCTAGCTAGCTAGCTA	IIIIIIIIIIIIIIIIIIII
EOF

# Index BAM
echo "Indexing BAM..."
samtools index sample1.bam

# Run gbcms
echo
echo "Running gbcms..."
export PYTHONPATH=$PYTHONPATH:$(pwd)/../../src
python3 -m gbcms.cli dna \
    --fasta reference.fa \
    --bam sample1.bam \
    --variants variants.vcf \
    --output-dir output \
    --format vcf

# Check output
echo
echo "Checking output..."
if [ -f output/sample1.vcf ]; then
    echo "✅ Output file created"
    echo
    echo "First few lines of output:"
    head -n 15 output/sample1.vcf
    
    # Count variants (excluding header)
    VARIANT_COUNT=$(grep -v "^#" output/sample1.vcf | wc -l)
    echo
    echo "Variants processed: $VARIANT_COUNT"
    
    if [ "$VARIANT_COUNT" -gt 0 ]; then
        echo "✅ VCF workflow test PASSED"
        EXIT_CODE=0
    else
        echo "❌ No variants in output"
        EXIT_CODE=1
    fi
else
    echo "❌ Output file not created"
    EXIT_CODE=1
fi

# Cleanup
echo
echo "Cleaning up..."
cd /
rm -rf "$TEST_DIR"

echo
echo "=========================================="
echo "Test complete"
echo "=========================================="

exit $EXIT_CODE
