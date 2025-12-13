#!/bin/bash
#SBATCH --job-name=align_8reps
#SBATCH --output=align_8reps_%j.out
#SBATCH --error=align_8reps_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --partition=courses

set -euo pipefail

echo "================================"
echo "Alignment job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "================================"

module purge
module load bwa/0.7.18
module load samtools/1.21

REF="/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/reference/reference.fasta"
echo "Reference: $REF"

# Build BWA index if missing
if [ ! -f "${REF}.bwt" ]; then
    echo "[INFO] BWA index not found. Building..."
    bwa index "$REF"
fi

# Configuration
THREADS_PER_SAMPLE=4
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-1}
MAX_JOBS=$(( TOTAL_CPUS / THREADS_PER_SAMPLE ))

# Ensure at least 1 job can run
[ "$MAX_JOBS" -lt 1 ] && MAX_JOBS=1

echo "[INFO] Using $THREADS_PER_SAMPLE threads per sample"
echo "[INFO] Max parallel samples: $MAX_JOBS"

# Function to align one sample
align_one() {
    local R1=$1
    local REF=$2
    local THREADS=$3

    local SAMPLE=${R1%_R1.fq}
    local R2=${SAMPLE}_R2.fq

    local BAM_UNSORTED=${SAMPLE}.unsorted.bam
    local BAM_SORTED=${SAMPLE}.sorted.bam

    echo "----------------------------------------"
    echo "[ALIGN] Sample: $SAMPLE"

    if [ ! -f "$R2" ]; then
        echo "[ERROR] R2 not found: $R2"
        return 1
    fi

    echo "[BWA-MEM] Aligning..."
    bwa mem -t "$THREADS" "$REF" "$R1" "$R2" \
        | samtools view -b -o "$BAM_UNSORTED" -

    echo "[SAMTOOLS] Sorting..."
    samtools sort -@ "$THREADS" -o "$BAM_SORTED" "$BAM_UNSORTED"
    rm -f "$BAM_UNSORTED"

    echo "[SAMTOOLS] Indexing..."
    samtools index "$BAM_SORTED"

    # Quick verification
    if samtools quickcheck "$BAM_SORTED" 2>/dev/null; then
        local read_count=$(samtools view -c "$BAM_SORTED" 2>/dev/null)
        echo "[SUCCESS] $SAMPLE: ${read_count} reads"
    else
        echo "[WARNING] $SAMPLE: BAM file may be corrupted"
    fi

    echo ""
}

export -f align_one

# Get the list of R1 files
R1_FILES=( synthetic_data*_R1.fq )

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "[ERROR] No synthetic_data*_R1.fq files found."
    exit 1
fi

echo "[INFO] Found ${#R1_FILES[@]} R1 files"
echo ""

# Start all alignments in parallel (8 jobs × 4 threads = 32 cores)
for R1 in "${R1_FILES[@]}"; do
    align_one "$R1" "$REF" "$THREADS_PER_SAMPLE" &
done

# Wait for all background jobs to finish
wait

echo "================================"
echo "Alignment finished at: $(date)"
echo ""

# List and verify all BAM files
echo "Created BAM files:"
if ls -1 *.sorted.bam 2>/dev/null; then
    echo ""
    echo "Verifying BAM files..."
    for bam in *.sorted.bam; do
        if samtools quickcheck "$bam" 2>/dev/null; then
            echo "✅ $bam: OK"
        else
            echo "❌ $bam: CORRUPTED"
        fi
    done
else
    echo "No BAM files found!"
fi

echo "================================"
