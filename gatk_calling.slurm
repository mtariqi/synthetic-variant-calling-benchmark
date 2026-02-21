#!/bin/bash
#SBATCH --job-name=gatk_new
#SBATCH --partition=courses
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/logs/gatk_%j.out
#SBATCH --error=/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/logs/gatk_%j.err

SAMPLE='tumor2'

REFERENCE='/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/reference/reference.fasta'
BAM="/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/bam/${SAMPLE}.bam"
OUTDIR='/scratch/islam.mdtar/vbench_group/variant-calling-benchmark/results/gatk'

mkdir -p ${OUTDIR}

# Call variants
singularity exec --bind /scratch docker://broadinstitute/gatk:4.2.3.0 gatk --java-options "-Xmx12g" HaplotypeCaller \
  -R ${REFERENCE} \
  -I ${BAM} \
  -O ${OUTDIR}/${SAMPLE}_simple.vcf

# Select SNPs
singularity exec --bind /scratch docker://broadinstitute/gatk:4.2.3.0 gatk SelectVariants \
  -R ${REFERENCE} \
  -V ${OUTDIR}/${SAMPLE}_simple.vcf \
  -select-type SNP \
  -O ${OUTDIR}/${SAMPLE}_SNPs.vcf

# Filter SNPs
singularity exec --bind /scratch docker://broadinstitute/gatk:4.2.3.0 gatk VariantFiltration \
  -R ${REFERENCE} \
  -V ${OUTDIR}/${SAMPLE}_SNPs.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.0" --filter-name "SOR3" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -O ${OUTDIR}/${SAMPLE}_SNPs_filtered.vcf

# Select INDELs
singularity exec --bind /scratch docker://broadinstitute/gatk:4.2.3.0 gatk SelectVariants \
  -R ${REFERENCE} \
  -V ${OUTDIR}/${SAMPLE}_simple.vcf \
  -select-type INDEL \
  -O ${OUTDIR}/${SAMPLE}_INDELs.vcf

# Filter INDELs
singularity exec --bind /scratch docker://broadinstitute/gatk:4.2.3.0 gatk VariantFiltration \
  -R ${REFERENCE} \
  -V ${OUTDIR}/${SAMPLE}_INDELs.vcf \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -O ${OUTDIR}/${SAMPLE}_INDELs_filtered.vcf

# Merge
singularity exec --bind /scratch docker://broadinstitute/gatk:4.2.3.0 gatk MergeVcfs \
  -I ${OUTDIR}/${SAMPLE}_SNPs_filtered.vcf \
  -I ${OUTDIR}/${SAMPLE}_INDELs_filtered.vcf \
  -O ${OUTDIR}/${SAMPLE}_HARDFILTER.vcf
