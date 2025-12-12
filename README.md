```
███████╗██╗   ██╗███╗   ██╗████████╗██╗   ██╗███████╗
██╔════╝██║   ██║████╗  ██║╚══██╔══╝██║   ██║██╔════╝
███████╗██║   ██║██╔██╗ ██║   ██║   ██║   ██║█████╗  
╚════██║██║   ██║██║╚██╗██║   ██║   ██║   ██║██╔══╝  
███████║╚██████╔╝██║ ╚████║   ██║   ╚██████╔╝███████╗
╚══════╝ ╚═════╝ ╚═╝  ╚═══╝   ╚═╝    ╚═════╝ ╚══════╝
```
<object type="image/svg+xml" data="assets/workflow.svg" width="100%"></object>

----
``

# 🧬 **Synthetic Variant Calling Benchmark**

### *Controlled Somatic Mutation Framework for Evaluating Variant Callers*

**Md Tariqul Islam (mtariqi)**, **Atra Alimoradian**, **Raghad Al-Ampudi** -- Bioinformatics Department, Northeastern University

---

![Python](https://img.shields.io/badge/Python-3.10+-blue)
![Singularity](https://img.shields.io/badge/Container-Singularity-purple)
![Nextflow](https://img.shields.io/badge/Nextflow-SLURM-darkgreen)
![Snakemake](https://img.shields.io/badge/Snakemake-Workflow-orange)
![GATK](https://img.shields.io/badge/GATK-4.2.3.0-yellow)
![DeepVariant](https://img.shields.io/badge/DeepVariant-1.2.0-red)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow)
![Zenodo](https://img.shields.io/badge/DOI-Coming%20Soon-blue)

---



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

```

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

```
```
<svg width="1080" height="620" viewBox="0 0 1080 620" xmlns="http://www.w3.org/2000/svg">

  <!-- Background -->
  <rect width="100%" height="100%" fill="#ffffff"/>

  <!-- Styles -->
  <style>
    .box {
      fill: #f5f5f5;
      stroke: #333333;
      stroke-width: 2;
      rx: 8;
      ry: 8;
    }
    .title {
      font-family: Helvetica, Arial, sans-serif;
      font-size: 20px;
      font-weight: 600;
      fill: #111111;
    }
    .label {
      font-family: Helvetica, Arial, sans-serif;
      font-size: 15px;
      fill: #222222;
    }
    .arrow {
      stroke: #000000;
      stroke-width: 2.5;
      marker-end: url(#arrowhead);
    }
  </style>

  <!-- Arrowhead -->
  <defs>
    <marker id="arrowhead" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto">
      <polygon points="0 0, 12 6, 0 12" fill="#000000"/>
    </marker>
  </defs>

  <!-- FASTQ Input -->
  <rect x="60" y="80" width="220" height="90" class="box"/>
  <text x="170" y="125" text-anchor="middle" class="title">FASTQ Inputs</text>
  <text x="170" y="147" text-anchor="middle" class="label">Tumor + Normal</text>

  <!-- Arrow 1 -->
  <line x1="280" y1="125" x2="360" y2="125" class="arrow"/>

  <!-- Nextflow -->
  <rect x="360" y="60" width="260" height="130" class="box"/>
  <text x="490" y="105" text-anchor="middle" class="title">Nextflow Pipeline</text>
  <text x="490" y="130" text-anchor="middle" class="label">DSL2 Modules</text>
  <text x="490" y="150" text-anchor="middle" class="label">Process Orchestration</text>

  <!-- Arrow 2 -->
  <line x1="620" y1="125" x2="700" y2="125" class="arrow"/>

  <!-- Singularity -->
  <rect x="700" y="60" width="260" height="130" class="box"/>
  <text x="830" y="105" text-anchor="middle" class="title">Singularity Containers</text>
  <text x="830" y="130" text-anchor="middle" class="label">DeepVariant / GATK</text>
  <text x="830" y="150" text-anchor="middle" class="label">BWA / Minimap2</text>

  <!-- Arrow 3 (down from Singularity) -->
  <line x1="830" y1="190" x2="830" y2="260" class="arrow"/>

  <!-- HPC Cluster -->
  <rect x="650" y="260" width="360" height="130" class="box"/>
  <text x="830" y="305" text-anchor="middle" class="title">HPC Cluster (SLURM)</text>
  <text x="830" y="330" text-anchor="middle" class="label">Multi-core Execution</text>
  <text x="830" y="350" text-anchor="middle" class="label">Resource-managed Jobs</text>

  <!-- Arrow 4 (from HPC to results) -->
  <line x1="830" y1="390" x2="830" y2="460" class="arrow"/>

  <!-- Results -->
  <rect x="650" y="460" width="360" height="130" class="box"/>
  <text x="830" y="505" text-anchor="middle" class="title">Benchmark Outputs</text>
  <text x="830" y="530" text-anchor="middle" class="label">Sorted BAM / VCF</text>
  <text x="830" y="550" text-anchor="middle" class="label">Precision / Recall</text>
  <text x="830" y="570" text-anchor="middle" class="label">Synthetic Truth Match</text>

  <!-- Left Branch: Preprocessing -->
  <line x1="490" y1="190" x2="490" y2="260" class="arrow"/>

  <rect x="360" y="260" width="260" height="130" class="box"/>
  <text x="490" y="305" text-anchor="middle" class="title">Preprocessing</text>
  <text x="490" y="330" text-anchor="middle" class="label">Alignment → Sorting</text>
  <text x="490" y="350" text-anchor="middle" class="label">Indexing → QC</text>

  <!-- Arrow Back to HPC -->
  <line x1="620" y1="325" x2="650" y2="325" class="arrow"/>

</svg><svg width="1080" height="620" viewBox="0 0 1080 620" xmlns="http://www.w3.org/2000/svg">

  <!-- Background -->
  <rect width="100%" height="100%" fill="#ffffff"/>

  <!-- Styles -->
  <style>
    .box {
      fill: #f5f5f5;
      stroke: #333333;
      stroke-width: 2;
      rx: 8;
      ry: 8;
    }
    .title {
      font-family: Helvetica, Arial, sans-serif;
      font-size: 20px;
      font-weight: 600;
      fill: #111111;
    }
    .label {
      font-family: Helvetica, Arial, sans-serif;
      font-size: 15px;
      fill: #222222;
    }
    .arrow {
      stroke: #000000;
      stroke-width: 2.5;
      marker-end: url(#arrowhead);
    }
  </style>

  <!-- Arrowhead -->
  <defs>
    <marker id="arrowhead" markerWidth="12" markerHeight="12" refX="10" refY="6" orient="auto">
      <polygon points="0 0, 12 6, 0 12" fill="#000000"/>
    </marker>
  </defs>

  <!-- FASTQ Input -->
  <rect x="60" y="80" width="220" height="90" class="box"/>
  <text x="170" y="125" text-anchor="middle" class="title">FASTQ Inputs</text>
  <text x="170" y="147" text-anchor="middle" class="label">Tumor + Normal</text>

  <!-- Arrow 1 -->
  <line x1="280" y1="125" x2="360" y2="125" class="arrow"/>

  <!-- Nextflow -->
  <rect x="360" y="60" width="260" height="130" class="box"/>
  <text x="490" y="105" text-anchor="middle" class="title">Nextflow Pipeline</text>
  <text x="490" y="130" text-anchor="middle" class="label">DSL2 Modules</text>
  <text x="490" y="150" text-anchor="middle" class="label">Process Orchestration</text>

  <!-- Arrow 2 -->
  <line x1="620" y1="125" x2="700" y2="125" class="arrow"/>

  <!-- Singularity -->
  <rect x="700" y="60" width="260" height="130" class="box"/>
  <text x="830" y="105" text-anchor="middle" class="title">Singularity Containers</text>
  <text x="830" y="130" text-anchor="middle" class="label">DeepVariant / GATK</text>
  <text x="830" y="150" text-anchor="middle" class="label">BWA / Minimap2</text>

  <!-- Arrow 3 (down from Singularity) -->
  <line x1="830" y1="190" x2="830" y2="260" class="arrow"/>

  <!-- HPC Cluster -->
  <rect x="650" y="260" width="360" height="130" class="box"/>
  <text x="830" y="305" text-anchor="middle" class="title">HPC Cluster (SLURM)</text>
  <text x="830" y="330" text-anchor="middle" class="label">Multi-core Execution</text>
  <text x="830" y="350" text-anchor="middle" class="label">Resource-managed Jobs</text>

  <!-- Arrow 4 (from HPC to results) -->
  <line x1="830" y1="390" x2="830" y2="460" class="arrow"/>

  <!-- Results -->
  <rect x="650" y="460" width="360" height="130" class="box"/>
  <text x="830" y="505" text-anchor="middle" class="title">Benchmark Outputs</text>
  <text x="830" y="530" text-anchor="middle" class="label">Sorted BAM / VCF</text>
  <text x="830" y="550" text-anchor="middle" class="label">Precision / Recall</text>
  <text x="830" y="570" text-anchor="middle" class="label">Synthetic Truth Match</text>

  <!-- Left Branch: Preprocessing -->
  <line x1="490" y1="190" x2="490" y2="260" class="arrow"/>

  <rect x="360" y="260" width="260" height="130" class="box"/>
  <text x="490" y="305" text-anchor="middle" class="title">Preprocessing</text>
  <text x="490" y="330" text-anchor="middle" class="label">Alignment → Sorting</text>
  <text x="490" y="350" text-anchor="middle" class="label">Indexing → QC</text>

  <!-- Arrow Back to HPC -->
  <line x1="620" y1="325" x2="650" y2="325" class="arrow"/>

</svg>

```

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



