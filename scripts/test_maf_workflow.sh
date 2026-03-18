#!/bin/bash
# End-to-end test for MAF workflow

set -e  # Exit on error

echo "=========================================="
echo "gbcms MAF Workflow Test"
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

# Create test MAF
echo "Creating test MAF..."
cat > variants.maf << 'EOF'
Hugo_Symbol	Chromosome	Start_Position	End_Position	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2	Tumor_Sample_Barcode	Matched_Norm_Sample_Barcode	t_ref_count	t_alt_count	n_ref_count	n_alt_count	Variant_Classification
GENE1	chr1	10	10	A	T		Tumor1	Normal1	10	5	15	0	Missense_Mutation
GENE2	chr1	20	20	C	G		Tumor1	Normal1	8	7	12	1	Missense_Mutation
GENE3	chr1	30	31	GA	G		Tumor1	Normal1	12	3	18	0	Frame_Shift_Del
GENE4	chr2	15	15	G	GC		Tumor1	Normal1	9	6	14	1	Frame_Shift_Ins
EOF

# Create test BAMs
echo "Creating test BAMs..."

# Tumor BAM
samtools view -bS > Tumor1.bam << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:120
@SQ	SN:chr2	LN:120
read1	99	chr1	5	60	20M	=	50	65	ATCGATCGATCGATCGATCG	IIIIIIIIIIIIIIIIIIII
read2	99	chr1	15	60	20M	=	60	65	CGATCGATCGATCGATCGAT	IIIIIIIIIIIIIIIIIIII
read3	99	chr1	25	60	20M	=	70	65	GATCGATCGATCGATCGATC	IIIIIIIIIIIIIIIIIIII
read4	99	chr2	10	60	20M	=	55	65	GCTAGCTAGCTAGCTAGCTA	IIIIIIIIIIIIIIIIIIII
EOF

# Normal BAM
samtools view -bS > Normal1.bam << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:120
@SQ	SN:chr2	LN:120
read5	99	chr1	5	60	20M	=	50	65	ATCGATCGATCGATCGATCG	IIIIIIIIIIIIIIIIIIII
read6	99	chr1	15	60	20M	=	60	65	CGATCGATCGATCGATCGAT	IIIIIIIIIIIIIIIIIIII
read7	99	chr1	25	60	20M	=	70	65	GATCGATCGATCGATCGATC	IIIIIIIIIIIIIIIIIIII
read8	99	chr2	10	60	20M	=	55	65	GCTAGCTAGCTAGCTAGCTA	IIIIIIIIIIIIIIIIIIII
EOF

# Index BAMs
echo "Indexing BAMs..."
samtools index Tumor1.bam
samtools index Normal1.bam

# Create BAM file-of-files
echo "Creating BAM file-of-files..."
cat > bam_files.txt << 'EOF'
Tumor1	Tumor1.bam
Normal1	Normal1.bam
EOF

# Run gbcms (MAF format)
echo
echo "Running gbcms (MAF format)..."
export PYTHONPATH=$PYTHONPATH:$(pwd)/../../src
python3 -m gbcms.cli dna \
    --fasta reference.fa \
    --bam-list bam_files.txt \
    --variants variants.maf \
    --output-dir output \
    --format maf

# Check outputs
echo
echo "Checking outputs..."
EXIT_CODE=0

if [ -f output/Tumor1.maf ]; then
    echo "✅ Tumor1 MAF output created"
    echo
    echo "First few lines of Tumor1 output:"
    head -n 3 output/Tumor1.maf
    
    VARIANT_COUNT=$(tail -n +2 output/Tumor1.maf | wc -l)
    echo
    echo "Variants in Tumor1 MAF: $VARIANT_COUNT"
else
    echo "❌ Tumor1 MAF output not created"
    EXIT_CODE=1
fi

echo

if [ -f output/Normal1.maf ]; then
    echo "✅ Normal1 MAF output created"
else
    echo "❌ Normal1 MAF output not created"
    EXIT_CODE=1
fi

if [ $EXIT_CODE -eq 0 ]; then
    echo
    echo "✅ MAF workflow test PASSED"
else
    echo
    echo "❌ MAF workflow test FAILED"
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
