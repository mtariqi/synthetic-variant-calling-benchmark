# 🧬 **Synthetic Variant Calling Benchmark**

### *Controlled Somatic Mutation Framework for Evaluating Variant Callers*

**Md Tariqul Islam (mtariqi)**, **Atra Alimoradian**, **Raghad Al-Ampudi** -- Bioinformatics Department, Northeastern University

---

## 🔖 **Badges**

<p align="left">

<a href="#"><img src="https://img.shields.io/badge/Python-3.10+-blue.svg"></a> <a href="#"><img src="https://img.shields.io/badge/DeepVariant-v1.2.0-orange.svg"></a> <a href="#"><img src="https://img.shields.io/badge/GATK-4.4-blue.svg"></a> <a href="#"><img src="https://img.shields.io/badge/BWA--MEM-0.7.18-green.svg"></a> <a href="#"><img src="https://img.shields.io/badge/HPC-SLURM-critical.svg"></a> <a href="#"><img src="https://img.shields.io/badge/License-MIT-green.svg"></a> <a href="#"><img src="https://img.shields.io/badge/Zenodo-DOI%20Ready-yellow.svg"></a> <a href="#"><img src="https://img.shields.io/badge/Reproducible%20Research-Yes-success.svg"></a>

</p>

---
[![Build Status](https://img.shields.io/badge/CI-Pending-lightgrey)]()
[![Coverage](https://img.shields.io/badge/Coverage-Not%20Measured-lightgrey)]()
[![Docker Hub](https://img.shields.io/badge/Docker%20Hub-ready-blue)]()
[![Pre-commit](https://img.shields.io/badge/pre--commit-enabled-success)]()
[![Snakemake](https://img.shields.io/badge/Snakemake-Workflow-green)]()
[![Nextflow](https://img.shields.io/badge/Nextflow-SLURM%20Profile-purple)]()


<img width="1131" height="767" alt="image" src="https://github.com/user-attachments/assets/b38911ee-31e4-40bc-b37a-c3708468eeff" />



# 📘 **Abstract**

Variant calling accuracy depends heavily on genome coverage, read depth, sequencing noise, and the mutational landscape. Public benchmark datasets (e.g., GIAB) are extremely large (250–800 GB) and were blocked on our HPC.

To address this, we developed a **synthetic tumor–normal somatic mutation framework**, enabling controlled benchmarks of:

* **GATK HaplotypeCaller**
* **Google DeepVariant**

Our approach compares:

### ✔ Dataset-1: *Concentrated reads*

5,000 read pairs from a **600 bp chr1 window** → **>8,000× depth** → **successful variant calling**

### ✔ Dataset-2 & Dataset-3: *Diluted reads*

Same 5,000 read pairs from **entire hg38 genome** → **<0.001× depth** → **no variants detected**

This repository contains the complete pipeline, scripts, figures, and results, demonstrating why *coverage depth is the single most important factor in mutation detection*.

---
# 📚 Related Work & Comparison to Published Benchmarking Studies

Reliable benchmarking of variant callers typically depends on very large, high-quality truth sets such as the Genome in a Bottle (GIAB) consortium. One of the most comprehensive variant calling benchmark studies is:

Barbitoff et al., “Systematic benchmarking of variant calling pipelines for clinical diagnostics” (2022)
doi: 10.1186/s13073-022-01057-7

🔬 What the Barbitoff et al. Study Did

```
| Aspect                     | Barbitoff et al. (2022)                                                   | Our Work                                           |
| -------------------------- | ------------------------------------------------------------------------- | -------------------------------------------------- |
| **Data scale**             | Full GIAB sequencing datasets (40–80 GB BAM per sample)                   | Synthetic reads (400–800 KB BAM per sample)        |
| **Samples**                | 14 human genomes                                                          | 1 synthetic tumor–normal pair (expandable)         |
| **Variants per sample**    | ~20,000–25,000                                                            | ~180–200                                           |
| **Aligners tested**        | BWA-MEM, Bowtie2, Isaac, NovoAlign                                        | BWA-MEM (Future: Bowtie2)                          |
| **Variant callers tested** | 9 callers (GATK, DeepVariant, Strelka2, Octopus, FreeBayes, Clair3, etc.) | 2 callers: GATK HaplotypeCaller & DeepVariant      |
| **Filtering strategies**   | CNN filtering, VQSR, hard filters                                         | Hard filters only (no VQSR due to synthetic truth) |
| **Truth set**              | GIAB high-confidence regions                                              | Synthetic somatic truth injected programmatically  |
```

<img width="586" height="761" alt="image" src="https://github.com/user-attachments/assets/63e6ff6b-6e40-419d-b7d6-a75daca5d7c5" />



# 🧬 Why Our Project Is Novel and Scientifically Valuable

Unlike GIAB-based studies, which rely on naturally occurring germline variation, our project introduces:

✔ Synthetic Somatic Mutation Injection

A controlled variant landscape can be generated deterministically, allowing:

Known true positives / true negatives

Precise mutation rates

Exact tumor vs. normal comparisons

This produces a perfect ground truth, enabling unbiased evaluation of variant callers.

✔ Coverage-Controlled Experiment Design

Our study demonstrates a key biological principle:

Variant calling accuracy collapses not because callers fail, but because coverage depth is insufficient.

Dataset-1 and Dataset-2/3 highlight this phenomenon vividly.

✔ Lightweight, HPC-Friendly Pipeline

Unlike GIAB workflows (250+ GB), our synthetic datasets:

Download instantly

Run on any HPC

Enable rapid benchmarking without license restrictions

This makes the project especially appealing to:

Research labs

Teaching HPC courses

Diagnostic pipeline developers

Employers evaluating your bioinformatics engineering skills




# 🧱 **Architectural Design**

## 🔧 **Pipeline Architecture (High-Level)**

```
 ┌────────────────────┐
 │  Reference Genome   │  (hg38 or chr1 slice)
 └─────────┬──────────┘
           │
           ▼
 ┌────────────────────┐
 │ Synthetic Generator │  Inject somatic mutations
 │ (Normal + Tumor)    │  Create FASTQ (5,000 pairs)
 └─────────┬──────────┘
           │
           ▼
 ┌────────────────────┐
 │   Alignment Layer   │  BWA-MEM
 │  (SLURM parallel)   │  BAM + Index
 └─────────┬──────────┘
           │
           ▼
 ┌────────────────────┐
 │ Variant Calling     │  GATK HC + DeepVariant
 │ (Pipeline A & B)    │
 └─────────┬──────────┘
           │
           ▼
 ┌────────────────────┐
 │ Benchmark Layer     │  Compare tumor vs. truth
 │  Precision/Recall   │  Concordance metrics
 └────────────────────┘
```

---

## 🧬 **Synthetic Strategy Visualization**

Place this file in: `figures/strategy_comparison.png`

```
```
<img width="658" height="633" alt="image" src="https://github.com/user-attachments/assets/f2da39a2-1201-44ee-ae07-59f5da3a1e98" />

```

Add reference to image (GitHub will render it):

```markdown
```
<img width="610" height="470" alt="image" src="https://github.com/user-attachments/assets/e37a9c9b-314c-401c-ab12-15be5eec3316" />

```

```

# 🎯 **Key Insight: Why Dataset-1 Works & Dataset-2 Fails**

### Dataset-1 (Success)

* Reads concentrated in **600 bp region**
* **91%** mapping to chr1
* Depth **> 8,000×**
* Variant callers detect true somatic mutations

### Datasets-2 & 3 (Failure)

* 5,000 reads scattered across **3.1 billion bases**
* Coverage drops to:

[
\text{Depth} < 0.001\times
]

* No position has enough reads → callers cannot detect mutations

---

# 🧪 **Methods**

## 1. Synthetic Data Generation


<img width="1681" height="910" alt="image" src="https://github.com/user-attachments/assets/ad04096e-5857-4e20-a020-a0b4099a9b07" />


* Extract 600 bp window: `chr1:100000–100600`
* Inject somatic SNPs (2%)
* Generate normal & tumor paired-end FASTQ

Script:
`scripts/synthetic_generator.py`

---

## 2. HPC Job (SLURM)

`slurm/generate_reads.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=syn_reads
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem=16G

module load python/3.10
python scripts/synthetic_generator.py --num_reads 5000
```

---

## 3. Alignment (BWA-MEM)

```bash
bwa mem -t 16 reference.fasta sample_R1.fastq.gz sample_R2.fastq.gz |
    samtools sort -o sample.sorted.bam
samtools index sample.sorted.bam
```

---

## 4. Variant Calling Pipelines

### **Pipeline A — GATK HaplotypeCaller**

```bash
gatk HaplotypeCaller \
   -R reference.fasta \
   -I tumor.sorted.bam \
   -O tumor.gatk.raw.vcf
```

### **Pipeline B — DeepVariant**

```bash
run_deepvariant \
   --model_type=WGS \
   --ref=reference.fasta \
   --reads=tumor.sorted.bam \
   --output_vcf=tumor.dv.vcf
```

---

# 📊 **Results**

## Variant Counts

| Sample           | GATK | DeepVariant | Δ   |
| ---------------- | ---- | ----------- | --- |
| synthetic_NORMAL | 179  | 197         | +18 |
| synthetic_TUMOR  | 383  | 408         | +25 |

DeepVariant consistently detects more true variants.

---

# 📁 **Repository Structure**

```bash
synthetic-variant-calling-benchmark/
│
├── data/
│   ├── reference/
│   └── fastq/
│
├── scripts/
│   ├── synthetic_generator.py
│   ├── evaluate_metrics.py
│   └── pipeline/
│       ├── run_bwa.sh
│       ├── run_gatk.sh
│       └── run_deepvariant.sh
│
├── slurm/
│   ├── generate_reads.slurm
│   ├── align.slurm
│   └── call_variants.slurm
│
├── results/
│   ├── vcfs/
│   ├── metrics/
│   └── figures/
│
├── environment.yml
├── config.json
├── .gitignore
├── LICENSE
└── README.md
```

---

# ⚙️ **Reproducibility**

Install environment:

```bash
conda env create -f environment.yml
conda activate synthetic-vc
```

Run full pipeline:

```bash
bash run_pipeline.sh
```

---

# 🔖 **License**

MIT License is included for Zenodo archiving.

---

# 📚 **Citable DOI (Zenodo Ready)**

Once you publish:

```text
@dataset{tariq2025syntheticvc,
  author = {Md Tariqul Islam}, {Atra Alimoradian}, {Raghad Al-Ampudi}
  title  = {Synthetic Variant Calling Benchmark},
  year   = 2025,
  doi    = {10.xxxx/zenodo.xxxxxx},
  url    = {},
}
```

---
### References

1. Barbitoff, Y.A., Khmelkova, D.N., Shcherbakova, I.V., et al.  
   *Systematic benchmarking of variant calling pipelines for clinical diagnostics.*  
   Genome Medicine 14, 10 (2022).  
   https://doi.org/10.1186/s13073-022-01057-7

2. Zook, J.M., et al.  
   *Open-access next-generation sequencing resources for variant analysis.*  
   GIAB Consortium, NIST.

3. Poplin, R., et al.  
   *A universal SNP and small-indel variant caller using deep neural networks.*  
   Nature Biotechnology (2018). *(DeepVariant)*



# 🙏 **Acknowledgements**

This work was completed as part of the Northeastern University Variant Benchmarks Group. Special thanks to Atra and Raghad for insights into somatic variant modeling and HPC genomics workflows.

---



