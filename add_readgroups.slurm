#!/bin/bash
#SBATCH --job-name=add_rg
#SBATCH --partition=courses
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/logs/add_rg_%j.out
#SBATCH --error=/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/logs/add_rg_%j.err

SAMPLE=$1

INBAM="/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/bam/${SAMPLE}.sorted.bam"
OUTBAM="/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/bam/${SAMPLE}.bam"

echo "Adding read groups to: ${SAMPLE}"

singularity exec --bind /scratch docker://broadinstitute/gatk:4.2.3.0 gatk AddOrReplaceReadGroups \
  -I ${INBAM} \
  -O ${OUTBAM} \
  -RGID 1 \
  -RGLB lib1 \
  -RGPL ILLUMINA \
  -RGPU unit1 \
  -RGSM ${SAMPLE}

echo "Creating index..."
singularity exec --bind /scratch docker://broadinstitute/gatk:4.2.3.0 gatk BuildBamIndex -I ${OUTBAM}

echo "Done!"
