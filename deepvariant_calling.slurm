#!/bin/bash
#SBATCH --job-name=dv_batch
#SBATCH --partition=courses
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/logs/dv_%j.out
#SBATCH --error=/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/logs/dv_%j.err

SAMPLE='normal_chr20_rep1'
REFERENCE='/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/reference/reference.fasta'
BAM="/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/synthetic/bam/${SAMPLE}.bam"
OUTDIR='/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/results/deepvariant'

mkdir -p ${OUTDIR}

singularity exec docker://google/deepvariant:1.2.0 /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=${REFERENCE} \
  --reads=${BAM} \
  --output_vcf=${OUTDIR}/${SAMPLE}_DV.vcf \
  --num_shards=8
