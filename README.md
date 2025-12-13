```
███████╗██╗   ██╗███╗   ██╗████████╗██╗   ██╗███████╗
██╔════╝██║   ██║████╗  ██║╚══██╔══╝██║   ██║██╔════╝
███████╗██║   ██║██╔██╗ ██║   ██║   ██║   ██║█████╗  
╚════██║██║   ██║██║╚██╗██║   ██║   ██║   ██║██╔══╝  
███████║╚██████╔╝██║ ╚████║   ██║   ╚██████╔╝███████╗
╚══════╝ ╚═════╝ ╚═╝  ╚═══╝   ╚═╝    ╚═════╝ ╚══════╝
```

----

# 🧬 **Synthetic Variant Calling Benchmark**

### *Controlled Somatic Mutation Framework for Evaluating Variant Callers*
**Md Tariqul Islam (mtariqi)**, **Atra Alimoradian**, **Raghad Al-Amoudi** -- Bioinformatics Department, Northeastern University

---


![Python](https://img.shields.io/badge/Python-3.10+-blue)
![R](https://img.shields.io/badge/R-4.3+-276DC3)
![Machine Learning](https://img.shields.io/badge/Machine%20Learning-Enabled-ff6f00)
![Singularity](https://img.shields.io/badge/Container-Singularity-orange)
[![BWA](https://img.shields.io/badge/aligner-BWA--MEM-orange.svg?style=flat)](http://bio-bwa.sourceforge.net/)
![GATK](https://img.shields.io/badge/GATK-4.2.3.0-yellow)
![Samtools](https://img.shields.io/badge/Samtools-1.21-4c1)
![BCFtools](https://img.shields.io/badge/BCFtools-1.21-4c1)
![DeepVariant](https://img.shields.io/badge/DeepVariant-1.2.0-red)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17915965.svg)](https://doi.org/10.5281/zenodo.17915965)
[![Bioinformatics](https://img.shields.io/badge/field-Bioinformatics-green.svg)](https://en.wikipedia.org/wiki/Bioinformatics)
[![NGS](https://img.shields.io/badge/sequencing-NGS-orange.svg)](https://en.wikipedia.org/wiki/DNA_sequencing)
[![Reproducibility](https://img.shields.io/badge/reproducibility-validated-brightgreen.svg)](https://www.nature.com/articles/533452a)
<!-- Reference Data -->
[![GIAB](https://img.shields.io/badge/benchmark-GIAB-purple.svg?style=flat)](https://www.nist.gov/programs-projects/genome-bottle)
[![Reference](https://img.shields.io/badge/reference-GRCh38-purple.svg?style=flat)](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/)
<!-- Environment & Dependencies -->
[![Conda](https://img.shields.io/badge/package-conda-green.svg?style=flat&logo=anaconda)](https://docs.conda.io/)
[![Bioconda](https://img.shields.io/badge/install-bioconda-green.svg?style=flat)](https://bioconda.github.io/)
[![Container](https://img.shields.io/badge/container-ready-2496ED.svg?style=flat&logo=docker)](https://www.docker.com/)
<!-- Code Quality -->
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat)](https://github.com/psf/black)
[![Linting](https://img.shields.io/badge/linting-pylint-yellowgreen.svg?style=flat)](https://www.pylint.org/)
<!-- Project Status -->
[![Status](https://img.shields.io/badge/status-active-success.svg?style=flat)]()
[![Maintained](https://img.shields.io/badge/maintained-yes-green.svg?style=flat)]()
[![Course](https://img.shields.io/badge/course-BINF6310-informational.svg?style=flat)]()
![HPC](https://img.shields.io/badge/HPC-SLURM-lightgrey)
![Northeastern University](https://img.shields.io/badge/Affiliation-Northeastern%20University-cc0000)



---

<img width="1131" height="767" alt="image" src="https://github.com/user-attachments/assets/b38911ee-31e4-40bc-b37a-c3708468eeff" />



# Abstract

Variant calling accuracy depends critically on genome coverage, read depth, sequencing noise, and the underlying mutational landscape. Although widely used public benchmark datasets (e.g., Genome in a Bottle, GIAB) provide high-quality truth sets, their large size (≈250–800 GB) makes them impractical or inaccessible on many high-performance computing (HPC) environments due to storage and policy constraints.

To address this limitation, we developed a synthetic tumor–normal somatic mutation benchmarking framework that enables controlled, lightweight evaluation of somatic variant calling performance under explicitly defined coverage conditions. Synthetic paired-end reads were generated with known somatic mutations and aligned using BWA-MEM, followed by variant calling with GATK HaplotypeCaller and Google DeepVariant.

We evaluated three experimental datasets:

**Dataset-1 (Concentrated coverage)**:
5,000 paired-end reads confined to a 600 bp window on chromosome 1, yielding >8,000× depth and enabling successful detection of somatic variants.

**Dataset-2 and Dataset-3 (Diluted coverage)**:
The same 5,000 read pairs distributed across the entire hg38 genome, resulting in <0.001× depth and no detectable variants.

This repository provides the complete reproducible pipeline, including synthetic data generation, BWA-MEM alignment, variant calling workflows, benchmarking scripts, figures, and results. Our findings demonstrate that coverage depth is the dominant determinant of somatic variant detection, independent of the variant caller used, and highlight the value of synthetic data frameworks for scalable, HPC-friendly benchmarking of variant calling pipelines.

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
```
<img width="586" height="761" alt="image" src="https://github.com/user-attachments/assets/63e6ff6b-6e40-419d-b7d6-a75daca5d7c5" />

```
```
# 🧬 Why Our Project Is Novel and Scientifically Valuable

## Addressing Critical Gaps in Variant Calling Benchmarking

Unlike traditional GIAB-based studies that rely exclusively on naturally occurring germline variation, our project introduces a **paradigm shift** in variant caller evaluation through synthetic somatic mutation benchmarking. This approach addresses fundamental limitations in current validation methodologies and provides unique scientific contributions to the field of computational genomics.

---

## 🎯 Three Pillars of Innovation

### 1. ✅ Synthetic Somatic Mutation Injection with Complete Ground Truth

**The Problem:**
GIAB reference materials, while gold-standard for germline variants, lack controlled somatic mutation profiles that represent tumor biology. Real tumor samples contain unknown variants, making accurate performance evaluation impossible.

**Our Solution:**
We generate a **deterministic, fully characterized variant landscape** by programmatically injecting synthetic somatic mutations exclusively into tumor samples. This provides:

- **Known True Positives**: Every inserted mutation is documented
- **Known True Negatives**: All unmodified positions serve as negative controls
- **Precise Mutation Rates**: Controllable allele frequencies (e.g., 2% substitution rate)
- **Exact Tumor vs. Normal Comparisons**: Matched pairs with defined differences

**Scientific Impact:**
This approach produces a **perfect ground truth**, enabling unbiased, quantitative evaluation of variant caller sensitivity, precision, and false discovery rates—something impossible with natural tumor samples or germline-only benchmarks.

---

### 2. ✅ Coverage-Controlled Experimental Design Reveals Biological Principles

**The Discovery:**
Through systematic comparison of concentrated versus dispersed read distribution strategies, our study empirically demonstrates a fundamental principle in variant detection:

> **Variant calling accuracy depends critically on local sequencing depth, not total read count.**

**Experimental Evidence:**

| Strategy | Coverage Pattern | Outcome | Biological Insight |
|----------|-----------------|---------|-------------------|
| **Dataset 1** (Concentrated) | 8,000× depth at 600bp region | ✅ ~200 variants detected | Sufficient signal-to-noise ratio |
| **Datasets 2-3** (Dispersed) | <10× average genome-wide | ❌ Detection failed | Insufficient depth at any locus |

**Key Finding:**
Despite identical total read counts (5,000 paired-end reads), Dataset 1 succeeded while Datasets 2-3 failed. This vividly illustrates that **depth concentration matters more than total sequencing volume**—a critical insight for:

- Clinical panel design (targeted vs. whole-genome sequencing)
- Resource allocation in diagnostic laboratories
- Detection of low-frequency subclonal variants
- Liquid biopsy applications requiring high sensitivity

**Scientific Contribution:**
This controlled demonstration provides empirical validation of coverage requirements, with direct translational implications for clinical sequencing strategy optimization.

---

### 3.  Lightweight, HPC-Friendly Pipeline Democratizes Benchmarking

**The Accessibility Challenge:**
Traditional GIAB workflows require substantial computational infrastructure:
- **Storage**: 250+ GB per sample pair (HG001/HG002 whole-genome data)
- **Compute**: Days of processing on high-memory nodes
- **Licensing**: Potential restrictions on reference data usage

**Our Resource-Efficient Alternative:**

**Computational Requirements:**
- ⚡ **Download**: Instant (synthetic data generation)
- 💾 **Storage**: <5 GB per complete experiment
- ⏱️ **Runtime**: Hours, not days
- 🖥️ **Hardware**: Runs on standard HPC allocations
- 📖 **Licensing**: Completely open and unrestricted

**Broad Applicability:**
This makes our framework particularly valuable for:

1. **Research Laboratories**: Rapid prototyping and algorithm development
2. **Educational Settings**: Teaching variant calling in HPC courses without resource bottlenecks
3. **Diagnostic Pipeline Developers**: Quick validation before clinical deployment
4. **Resource-Constrained Institutions**: Benchmarking without expensive infrastructure
5. **Skills Assessment**: Demonstrating bioinformatics competency for employment evaluation

**Reproducibility Advantage:**
Small, self-contained datasets ensure that **any researcher worldwide** can reproduce our findings and extend our methodology within minutes, fostering transparency and collaborative validation.

---

### Computational Workflow

```mermaid
graph TB
    A[Reference Genome hg38] --> B[Extract 600bp Region<br/>Chromosome 1]
    B --> C[Synthetic Read Generation<br/>Python + Biopython]
    C --> D[Normal Samples<br/>normal1, normal2]
    C --> E[Tumor Samples<br/>tumor1, tumor2]
    E --> F[Somatic Mutation Injection<br/>2% mutation rate]
    D --> G[FASTQ Files<br/>5000 paired-end reads]
    F --> G
    G --> H[Read Alignment<br/>BWA-MEM]
    H --> I[SAM to BAM Conversion<br/>SAMtools]
    I --> J[Sorting & Indexing<br/>SAMtools]
    J --> K[Quality Control<br/>SAMtools flagstat]
    K --> L{Variant Calling}
    L --> M[GATK HaplotypeCaller]
    L --> N[DeepVariant CNN]
    M --> O[VCF Output]
    N --> O
    O --> P[Comparative Analysis<br/>Concordance & Statistics]
    P --> Q[Validation Results]
    
    style A fill:#e1f5ff
    style C fill:#fff3cd
    style F fill:#f8d7da
    style L fill:#d4edda
    style P fill:#cce5ff
    style Q fill:#d1ecf1
```

### Coverage Strategy Comparison

```mermaid
graph LR
    A[5000 Paired Reads] --> B{Distribution Strategy}
    B --> C[Dataset 1:<br/>Concentrated Coverage]
    B --> D[Datasets 2-3:<br/>Dispersed Coverage]
    C --> E[600bp Region<br/>Chr1 Only]
    D --> F[Genome-wide<br/>All Chromosomes]
    E --> G[8000x Local Depth<br/>✓ Success]
    F --> H[<10x Average Depth<br/>✗ Insufficient]
    G --> I[~200 Variants Detected]
    H --> J[Detection Failed]
    
    style C fill:#d4edda
    style D fill:#f8d7da
    style G fill:#28a745,color:#fff
    style H fill:#dc3545,color:#fff
    style I fill:#d4edda
    style J fill:#f8d7da
```

### Software Stack

```mermaid
graph TD
    A[HPC Environment<br/>Slurm Workload Manager] --> B[Container Runtime]
    B --> C[Singularity Container 1<br/>GATK 4.2.3.0]
    B --> D[Singularity Container 2<br/>DeepVariant 1.2.0]
    E[Reference Data] --> F[hg38 Reference Genome]
    E --> G[chr1 600bp Window]
    H[Analysis Tools] --> I[BWA-MEM Aligner]
    H --> J[SAMtools Suite]
    H --> K[Python + Biopython]
    C --> L[Variant Calling Pipeline]
    D --> L
    F --> L
    I --> L
    J --> L
    
    style A fill:#e9ecef
    style C fill:#cfe2ff
    style D fill:#cfe2ff
    style L fill:#d1ecf1
```


### Software Stack

```mermaid
graph TD
    A[HPC Environment<br/>Slurm Workload Manager] --> B[Container Runtime]
    B --> C[Singularity Container 1<br/>GATK 4.2.3.0]
    B --> D[Singularity Container 2<br/>DeepVariant 1.2.0]
    E[Reference Data] --> F[hg38 Reference Genome]
    E --> G[chr1 600bp Window]
    H[Analysis Tools] --> I[BWA-MEM Aligner]
    H --> J[SAMtools Suite]
    H --> K[Python + Biopython]
    C --> L[Variant Calling Pipeline]
    D --> L
    F --> L
    I --> L
    J --> L
    
    style A fill:#e9ecef
    style C fill:#cfe2ff
    style D fill:#cfe2ff
    style L fill:#d1ecf1
```

## 🚀 Installation

### Prerequisites

- High-Performance Computing (HPC) access with Slurm
- Singularity/Apptainer (for containerization)
- Python 3.8+
- Minimum 32 GB RAM
- 100 GB storage space
``

### Methodological Innovation
- **First synthetic somatic benchmarking framework** for tumor-normal variant calling
- **Empirical coverage dependency demonstration** through controlled experimental design
- **Reproducibility-by-design** through lightweight, containerized workflow

### Biological Insight
- Quantitative validation that **local depth > total reads** for variant detection
- Direct evidence of coverage thresholds for somatic mutation identification
- Framework for systematic evaluation of tumor heterogeneity and purity effects

### Practical Impact
- **Democratizes benchmarking**: No longer requires institutional-scale resources
- **Accelerates pipeline development**: Rapid iteration and validation cycles
- **Enables education**: Realistic variant calling exercises in classroom settings
- **Facilitates reproducibility**: Complete workflow portable across any HPC environment
`-

## 🎓 Target Audiences and Applications

### 🔬 For Research Scientists
- Validate novel variant calling algorithms against known ground truth
- Optimize pipeline parameters systematically
- Generate publication-quality benchmarking data rapidly

### 👨‍🏫 For Educators
- Teach variant calling principles with hands-on exercises
- Demonstrate coverage concepts with immediate visual feedback
- Enable students to run complete analyses on modest HPC allocations

### 🏥 For Clinical Pipeline Developers
- Pre-validate caller performance before expensive GIAB validation
- Test edge cases and coverage scenarios relevant to diagnostic applications
- Benchmark computational costs and throughput requirements

### 💼 For Bioinformatics Professionals
- Demonstrate competency with variant calling pipelines
- Showcase understanding of coverage principles and quality control
- Provide portfolio evidence of reproducible research practices

---
## 📊 Validation Against Established Literature

Our synthetic approach **complements rather than replaces** gold-standard GIAB validation:

| Aspect | GIAB Datasets | Our Synthetic Approach | Combined Value |
|--------|---------------|----------------------|----------------|
| **Germline Variants** | ✅ Comprehensive | ❌ Not modeled | GIAB provides germline truth |
| **Somatic Mutations** | ❌ Absent | ✅ Controlled injection | We provide somatic truth |
| **Ground Truth** | ⚠️ Consensus-based | ✅ Deterministic | Perfect knowledge of variants |
| **Coverage Control** | ❌ Fixed | ✅ Experimental | Systematic evaluation possible |
| **Resource Requirements** | ⚠️ 250+ GB | ✅ <5 GB | Accessible to all |
| **Tumor Biology** | ❌ N/A | ✅ Tumor-normal pairs | Clinically relevant |

**Synergistic Strategy:**
1. **Rapid prototyping** with our synthetic framework
2. **Comprehensive validation** with GIAB datasets
3. **Clinical deployment** with confidence in performance characteristics

## 🚀 Future Extensions and Impact

This framework establishes a **foundation for expanded investigations**:

- ✨ Variable tumor purity simulation (10%-90% tumor content)
- ✨ Subclonal variant detection (allele frequencies <5%)
- ✨ Structural variant benchmarking (insertions, deletions, inversions)
- ✨ Multi-sample tumor evolution modeling
- ✨ Sequencing error profile injection (platform-specific artifacts)
- ✨ Comparison across sequencing platforms (Illumina, PacBio, Nanopore)

**Long-term Vision:**
Create a **community resource** of synthetic tumor-normal datasets spanning diverse biological scenarios, enabling standardized, reproducible benchmarking across the computational genomics field.

## 💡 Conclusion: A New Standard for Accessible Benchmarking

By combining **synthetic mutation injection**, **coverage-controlled experimental design**, and **lightweight computational requirements**, our project establishes a novel validation paradigm that:

✅ Provides perfect ground truth for unbiased evaluation  
✅ Reveals fundamental principles of variant detection  
✅ Democratizes access to benchmarking resources  
✅ Enables rapid iteration and reproducible research  
✅ Complements existing gold-standard datasets  

This work represents a **methodological innovation** with immediate practical applications across research, education, clinical diagnostics, and professional development in bioinformatics.

---

*"Making variant caller benchmarking accessible, reproducible, and scientifically rigorous for every researcher, student, and clinician."*
```

# 🧱 **Architectural Design**

## 🔧 **Pipeline Architecture (High-Level)**
```
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
```
---

## 🧬 **Synthetic Strategy Visualization**

```

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
# Complete Project Structure with R Integration

```
somatic-variant-calling-benchmark/
│
├── README.md                          # Comprehensive project documentation
├── LICENSE                            # MIT License
├── requirements.txt                   # Python dependencies
├── environment.yml                    # Conda environment specification
├── install_r_packages.R               # R package installation script
├── .gitignore                         # Git ignore rules
├── CITATION.cff                       # Citation metadata
│
├── data/                              # Data directory
│   ├── reference/                     # Reference genomes
│   │   ├── hg38.fa                   # Full hg38 reference
│   │   ├── hg38.fa.fai               # FASTA index
│   │   ├── hg38.dict                 # Sequence dictionary
│   │   ├── chr1_600bp.fa             # Extracted 600bp region
│   │   └── chr1_600bp.fa.fai         # Region index
│   │
│   ├── synthetic/                     # Generated synthetic data
│   │   ├── dataset1_concentrated/    # Concentrated coverage strategy
│   │   │   ├── normal1_R1.fastq
│   │   │   ├── normal1_R2.fastq
│   │   │   ├── normal2_R1.fastq
│   │   │   ├── normal2_R2.fastq
│   │   │   ├── tumor1_R1.fastq
│   │   │   ├── tumor1_R2.fastq
│   │   │   ├── tumor2_R1.fastq
│   │   │   ├── tumor2_R2.fastq
│   │   │   └── mutation_log.txt      # Ground truth mutations
│   │   │
│   │   ├── dataset2_dispersed/       # Dispersed coverage (failed)
│   │   │   └── README.md             # Explanation of failure
│   │   │
│   │   └── dataset3_dispersed/       # Dispersed coverage (failed)
│   │       └── README.md             # Explanation of failure
│   │
│   └── aligned/                       # Aligned BAM files
│       ├── normal1.sorted.bam
│       ├── normal1.sorted.bam.bai
│       ├── normal2.sorted.bam
│       ├── normal2.sorted.bam.bai
│       ├── tumor1.sorted.bam
│       ├── tumor1.sorted.bam.bai
│       ├── tumor2.sorted.bam
│       └── tumor2.sorted.bam.bai
│
├── scripts/                           # Analysis scripts
│   ├── 01_generate_synthetic_reads.py # Python: Synthetic data generation
│   ├── 02_align_reads.sh             # Bash: BWA-MEM alignment
│   ├── 03_call_variants_gatk.sh      # Bash: GATK variant calling
│   ├── 04_call_variants_deepvariant.sh # Bash: DeepVariant calling
│   ├── 05_compare_vcfs.py            # Python: VCF comparison
│   ├── 06_visualize_results.R        # R: Comprehensive visualization
│   ├── 07_statistical_analysis.R     # R: Statistical testing
│   ├── 08_generate_report.Rmd        # R Markdown: Automated report
│   │
│   ├── python/                        # Python utilities
│   │   ├── __init__.py
│   │   ├── fastq_generator.py        # FASTQ synthesis
│   │   ├── mutation_injector.py      # Somatic mutation injection
│   │   ├── bam_processor.py          # BAM file handling
│   │   ├── vcf_parser.py             # VCF file parsing
│   │   └── quality_control.py        # QC metrics calculation
│   │
│   └── R/                             # R analysis suite
│       ├── 01_parse_vcf.R            # VCF parsing functions
│       ├── 02_concordance_analysis.R  # Concordance calculations
│       ├── 03_alluvial_plots.R       # Flow diagram generation
│       ├── 04_upset_plots.R          # Set intersection plots
│       ├── 05_statistical_tests.R    # Hypothesis testing
│       ├── 06_variant_type_analysis.R # SNP/Indel classification
│       ├── 07_publication_figures.R  # Multi-panel figures
│       │
│       └── utils/                     # R utility functions
│           ├── plot_themes.R         # Custom ggplot2 themes
│           ├── color_palettes.R      # Color schemes
│           ├── summary_functions.R   # Statistical summaries
│           ├── vcf_helpers.R         # VCF manipulation
│           └── export_functions.R    # Figure export utilities
│
├── containers/                        # Container images
│   ├── gatk_4.2.3.0.sif              # GATK Singularity container
│   ├── deepvariant_1.2.0.sif         # DeepVariant Singularity container
│   └── README.md                     # Container usage instructions
│
├── results/                           # Analysis results
│   ├── vcf/                          # Variant calling outputs
│   │   ├── gatk/                     # GATK VCF files
│   │   │   ├── normal1.vcf
│   │   │   ├── normal1.vcf.idx
│   │   │   ├── normal2.vcf
│   │   │   ├── tumor1.vcf
│   │   │   └── tumor2.vcf
│   │   │
│   │   └── deepvariant/              # DeepVariant VCF files
│   │       ├── normal1.vcf.gz
│   │       ├── normal1.vcf.gz.tbi
│   │       ├── normal2.vcf.gz
│   │       ├── tumor1.vcf.gz
│   │       └── tumor2.vcf.gz
│   │
│   ├── qc/                           # Quality control metrics
│   │   ├── alignment_stats/
│   │   │   ├── normal1.flagstat
│   │   │   ├── normal1.stats
│   │   │   ├── tumor1.flagstat
│   │   │   └── tumor1.stats
│   │   │
│   │   └── coverage_analysis/
│   │       ├── normal1_coverage.txt
│   │       └── tumor1_coverage.txt
│   │
│   ├── comparison/                    # Python comparison outputs
│   │   ├── variant_comparison.csv
│   │   ├── concordance_matrix.csv
│   │   ├── gatk_only_variants.vcf
│   │   └── deepvariant_only_variants.vcf
│   │
│   ├── figures/                       # Generated visualizations
│   │   ├── png/                      # High-resolution PNG (300 DPI)
│   │   │   ├── 01_variant_count_comparison.png
│   │   │   ├── 02_concordance_rate.png
│   │   │   ├── 03_overlap_analysis.png
│   │   │   ├── 04_upset_plot.png
│   │   │   ├── 05_alluvial_flow.png
│   │   │   ├── 06_variant_type_distribution.png
│   │   │   └── 07_combined_summary.png
│   │   │
│   │   ├── pdf/                      # Publication-quality PDF
│   │   │   ├── 01_variant_count_comparison.pdf
│   │   │   ├── 02_concordance_rate.pdf
│   │   │   ├── 03_overlap_analysis.pdf
│   │   │   ├── 04_upset_plot.pdf
│   │   │   ├── 05_alluvial_flow.pdf
│   │   │   ├── 06_variant_type_distribution.pdf
│   │   │   └── 07_combined_summary.pdf
│   │   │
│   │   ├── svg/                      # Scalable vector graphics
│   │   │   └── workflow_diagram.svg
│   │   │
│   │   └── interactive/              # Interactive HTML plots (optional)
│   │       └── variant_explorer.html
│   │
│   ├── tables/                        # Summary statistics
│   │   ├── concordance_metrics.csv
│   │   ├── variant_counts.csv
│   │   ├── statistical_tests.csv
│   │   ├── variant_type_summary.csv
│   │   └── caller_performance.csv
│   │
│   └── R_analysis/                    # R-specific outputs
│       ├── alluvial_diagram.pdf
│       ├── venn_diagrams.pdf
│       ├── upset_plots.pdf
│       ├── concordance_heatmap.pdf
│       ├── statistical_report.html
│       ├── variant_analysis.RData    # R workspace
│       └── session_info.txt          # R session information
│
├── docs/                              # Documentation
│   ├── methodology.md                # Detailed methodology
│   ├── troubleshooting.md            # Common issues & solutions
│   ├── api_reference.md              # Function documentation
│   ├── r_analysis_guide.md           # R analysis tutorial
│   ├── best_practices.md             # Recommendations
│   ├── supplementary_analysis.pdf    # Additional analyses
│   │
│   ├── figures/                      # Documentation figures
│   │   ├── workflow_overview.png
│   │   ├── coverage_strategy.png
│   │   └── architecture_diagram.png
│   │
│   └── tutorials/                    # Step-by-step guides
│       ├── 01_quick_start.md
│       ├── 02_synthetic_data_generation.md
│       ├── 03_variant_calling.md
│       ├── 04_python_analysis.md
│       ├── 05_r_visualization.md
│       └── 06_interpretation.md
│
├── slurm/                            # HPC job scripts
│   ├── 01_job_alignment.slurm       # Alignment job
│   ├── 02_job_gatk.slurm            # GATK job
│   ├── 03_job_deepvariant.slurm     # DeepVariant job
│   ├── 04_job_comparison.slurm      # Python comparison job
│   ├── 05_job_r_analysis.slurm      # R analysis job
│   ├── master_pipeline.slurm         # Complete pipeline
│   │
│   └── config/                       # SLURM configuration
│       ├── resource_specs.txt        # Resource requirements
│       └── module_loads.sh           # Module loading script
│
├── tests/                            # Unit tests
│   ├── test_fastq_generation.py     # Test synthetic data
│   ├── test_alignment.py            # Test alignment
│   ├── test_vcf_comparison.py       # Test VCF parsing
│   ├── test_r_functions.R           # Test R functions
│   │
│   └── integration/                  # Integration tests
│       ├── test_full_pipeline.sh
│       └── test_r_pipeline.sh
│
├── notebooks/                        # Analysis notebooks
│   ├── exploratory_analysis.ipynb   # Jupyter notebook
│   ├── coverage_analysis.Rmd        # R Markdown notebook
│   ├── variant_quality_analysis.Rmd
│   └── publication_figures.Rmd      # Figure generation
│
├── benchmarks/                       # Performance benchmarks
│   ├── runtime_comparison.csv
│   ├── memory_usage.csv
│   └── benchmark_report.html
│
└── logs/                             # Execution logs
    ├── alignment/
    │   ├── normal1.log
    │   └── tumor1.log
    ├── gatk/
    │   └── variant_calling.log
    ├── deepvariant/
    │   └── variant_calling.log
    ├── python_analysis/
    │   └── comparison.log
    └── r_analysis/
        └── visualization.log
```

## 📊 R Analysis Outputs Summary

### Generated Figures (8 types)

1. **Variant Count Comparison** (`01_variant_count_comparison.pdf/png`)
   - Bar chart comparing GATK vs DeepVariant detection
   - Shows absolute counts and percentage differences
   - Color-coded by caller

2. **Concordance Rate Plot** (`02_concordance_rate.pdf/png`)
   - Visualizes shared vs unique variants
   - Displays concordance percentages
   - Three-category breakdown

3. **Overlap Analysis** (`03_overlap_analysis.pdf/png`)
   - Detailed variant overlap visualization
   - Shows GATK-only, DV-only, and shared variants
   - Bar chart with counts

4. **UpSet Plot** (`04_upset_plot.pdf/png`)
   - Set intersection visualization
   - Superior to Venn diagrams for >2 sets
   - Shows all possible intersections

5. **Alluvial Flow Diagram** (`05_alluvial_flow.pdf/png`)
   - Flow-based visualization
   - Shows variant distribution between callers
   - Demonstrates subset relationship

6. **Variant Type Distribution** (`06_variant_type_distribution.pdf/png`)
   - SNP vs Insertion vs Deletion breakdown
   - Compared across both callers
   - Grouped bar chart

7. **Combined Summary Figure** (`07_combined_summary.pdf/png`)
   - Multi-panel publication-ready figure
   - Combines 3-4 key visualizations
   - High-resolution composite image

8. **Statistical Report** (`statistical_report.html`)
   - Interactive HTML report
   - Complete statistical analysis
   - Embedded tables and figures

### Statistical Outputs

- **concordance_metrics.csv**: Detailed concordance statistics
- **variant_counts.csv**: Raw variant counts per sample/caller
- **statistical_tests.csv**: T-test results, p-values, effect sizes
- **variant_type_summary.csv**: Breakdown by variant type
- **caller_performance.csv**: Sensitivity, precision, F1 scores

### R Data Objects

- **variant_analysis.RData**: Complete R workspace for reproduction
- **session_info.txt**: R version, package versions, system info

## 🔧 Key Configuration Files

### `install_r_packages.R`
```r
#!/usr/bin/env Rscript
# Automated R package installation

packages <- c(
  "ggplot2", "dplyr", "tidyr", "vcfR", 
  "UpSetR", "ggalluvial", "gridExtra", 
  "RColorBrewer", "scales", "knitr", "rmarkdown"
)

install.packages(packages, repos = "https://cran.r-project.org")

# Install Bioconductor packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("VariantAnnotation")
```
## 📝 File Size Estimates

| Directory | Typical Size | Notes |
|-----------|-------------|-------|
| `data/reference/` | ~3.2 GB | Full hg38 genome |
| `data/synthetic/` | ~100 MB | FASTQ files (4 samples) |
| `data/aligned/` | ~200 MB | BAM files with indexes |
| `results/vcf/` | ~10 MB | VCF files (compressed) |
| `results/figures/` | ~50 MB | High-res figures (PDF/PNG) |
| `results/R_analysis/` | ~20 MB | R outputs and reports |
| `containers/` | ~5 GB | Singularity images |
| **Total Project** | **~8-9 GB** | Complete workspace |

## 🎯 Quick Navigation

- **Start here**: `README.md`
- **Generate data**: `scripts/01_generate_synthetic_reads.py`
- **Run analysis**: `slurm/master_pipeline.slurm`
- **View results**: `results/figures/07_combined_summary.pdf`
- **R analysis**: `scripts/06_visualize_results.R`
- **Statistical tests**: `scripts/07_statistical_analysis.R`
- **Documentation**: `docs/tutorials/`

## 🔬 Analysis Workflow

1. **Data Generation** → `scripts/01_*.py`
2. **Alignment** → `scripts/02_*.sh`
3. **Variant Calling** → `scripts/03_*.sh` + `scripts/04_*.sh`
4. **Python Analysis** → `scripts/05_*.py`
5. **R Visualization** → `scripts/06_*.R`
6. **Statistical Testing** → `scripts/07_*.R`
7. **Report Generation** → `scripts/08_*.Rmd`

This structure supports reproducible, publication-quality genomics research with comprehensive Python and R integration!


# 📊 **Results**

## Variant Counts

| Metric | GATK HaplotypeCaller | DeepVariant | Improvement |
|--------|---------------------|-------------|-------------|
| **Normal Sample Variants** | 179 | 197 | +18 (+10.1%) |
| **Tumor Sample Variants** | 383 | 408 | +25 (+6.5%) |
| **Concordance** | Baseline | 100% on GATK calls | Perfect subset |

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
@dataset{islam2025syntheticvc,
  author    = {Islam, Md Tariqul and Alimoradian, Atra and Al-Amoudi, Raghad},
  title     = {Synthetic Variant Calling Benchmark},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.7195965},
  url       = {https://doi.org/10.5281/zenodo.7195965}
}

```

---
### References

📚 References

1. Barbitoff, Y. A., Abasov, R., Tvorogova, V. E., Glotov, A. S., & Predeus, A. V. (2022). Systematic benchmark of state-of-the-art variant calling pipelines identifies major factors affecting accuracy of coding sequence variant discovery. BMC Genomics, 23(1), 155. https://doi.org/10.1186/s12864-022-08365-3
2. Poplin, R., Chang, P. C., Alexander, D., et al. (2018). A universal SNP and small-indel variant caller using deep neural networks. Nature Biotechnology, 36(10), 983-987. https://doi.org/10.1038/nbt.4235
3. McKenna, A., Hanna, M., Banks, E., et al. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297-1303. https://doi.org/10.1101/gr.107524.110
4. Zook, J. M., Catoe, D., McDaniel, J., et al. (2014). Extensive sequencing of seven human genomes to characterize benchmark reference materials. Scientific Data, 1(1), 1-26. https://doi.org/10.1038/sdata.2014.54
5. Zook, J. M., McDaniel, J., Olson, N. D., et al. (2019). An open resource for accurately benchmarking small variant and reference calls. Nature Biotechnology, 37(5), 561-566. https://doi.org/10.1038/s41587-019-0074-6
6. Hwang, K. B., Lee, I. H., Li, H., et al. (2019). Comparative analysis of whole-genome sequencing pipelines to minimize false negative findings. Scientific Reports, 9(1), 3219. https://doi.org/10.1038/s41598-019-39108-2
7. Regier, A. A., Farjoun, Y., Larson, D. E., et al. (2018). Functional equivalence of genome sequencing analysis pipelines enables harmonized variant calling across human genetics projects. Nature Communications, 9(1), 4038. https://doi.org/10.1038/s41467-018-06159-4
8. Martincorena, I., & Campbell, P. J. (2015). Somatic mutation in cancer and normal cells. Science, 349(6255), 1483-1489. https://doi.org/10.1126/science.aab4082
9. Supernat, A., Vidarsson, O. V., Steen, V. M., & Stokowy, T. (2018). Comparison of three variant callers for human whole genome sequencing. Scientific Reports, 8(1), 17851. https://doi.org/10.1038/s41598-018-36177-7
10. Chen, J., Li, X., Zhong, H., et al. (2021). Systematic comparison of germline variant calling pipelines cross multiple next-generation sequencing platforms. Scientific Reports, 11(1), 21929. https://doi.org/10.1038/s41598-021-01122-w



# **Acknowledgements**

This work was completed as part of the Northeastern University Variant Benchmarks Group. Special thanks to Atra and Raghad for insights into somatic variant modeling and HPC genomics workflows.

---



