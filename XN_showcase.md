# Benchmarking AI-Based and Traditional Variant Calling Pipelines Using Synthetic Data


**Raghad Al-Amoudi**, MSc Bioinformatics · **Atra Alimoradian**, MSc Bioinformatics · **Md Tariqul Islam**, MSc Bioinformatics
Faculty Advisor: Prof. Oyeronke Ayansola. BINF 6310: Introduction to Computational Methods in Bioinformatics

---

## TL;DR

We reproduced the central finding of [Barbitoff et al. (2022)](https://doi.org/10.1186/s12864-022-08365-3) using lightweight synthetic data on a resource-constrained HPC cluster. **DeepVariant detects 8–14% more variants than GATK HaplotypeCaller** across matched synthetic NORMAL and TUMOR samples, and confirms 100% of GATK's filtered calls. The directional finding; AI-based variant calling outperforms traditional methods, held even at small scale.

---

## Key Results

| Sample | GATK (PASS) | DeepVariant | Shared | GATK-only | DV-only | Concordance |
|--------|------------:|------------:|-------:|----------:|--------:|------------:|
| Synthetic NORMAL | 172 | 197 | 172 | **0** | +25 (+14.5%) | 87.3% |
| Synthetic TUMOR  | 377 | 408 | 377 | **0** | +31 (+8.2%)  | 92.4% |

GATK raw counts were 179 (NORMAL) and 383 (TUMOR) before hard filtering. Every variant GATK kept was also called by DeepVariant.

---

## Methods

- **Synthetic data:** 5,000 paired-end reads confined to a 600 bp window on chr1 → ~8,000× local coverage. GIAB samples (HG001–HG007) were blocked by HPC IT restrictions, so synthetic reads with injected indel mutations were generated for NORMAL and TUMOR samples.
- **Alignment:** BWA-MEM → sorted, indexed BAM with read groups.
- **GATK HaplotypeCaller v4.2.3.0:** `AddOrReplaceReadGroups → HaplotypeCaller → SelectVariants → VariantFiltration → MergeVcfs`. Hard filters from Table 3 of Barbitoff et al. (BQSR skipped because no compatible dbSNP/Mills resources for synthetic reference).
- **DeepVariant v1.2.0:** single command with `--model_type=WGS --num_shards=8`.
- **Containers:** both callers via Singularity (Docker requires root and is unavailable on shared HPC).
- **Concordance:** `bcftools isec` → GATK-only, DeepVariant-only, shared.

Full SLURM scripts: [`add_readgroups.slurm`](add_readgroups.slurm), [`gatk_variant_calling.slurm`](gatk_variant_calling.slurm), [`deepvariant_calling.slurm`](deepvariant_calling.slurm).

---

## What We Learned Beyond the Paper

1. **Read concentration matters more than total read count.** Two follow-up datasets used the same 5,000 reads but spread them across the full GRCh38 genome, coverage dropped below 0.001× and *zero variants* were detected. 
2. **BQSR cannot be run on synthetic references** without matched dbSNP/Mills resource files. Skipping it preserves a fair head-to-head comparison since both callers receive the same uncalibrated input.
3. **GATK applies strictly conservative filtering.** Across both samples, the GATK-only set was empty, DeepVariant captured every variant GATK kept, plus additional calls. The difference is driven by recall, not by either tool missing each other's calls.

---

## What We Cannot Claim

- Exact F1, precision, or recall scores; synthetic data has no GIAB-style truth set, so absolute accuracy metrics are not available.
- Performance across diverse genomic contexts (GC content, repetitive regions, exon boundaries).
- Generalization across ancestry or WES vs. WGS; these would require the full 14 GIAB samples from the original study.

---

## Contact & Code

- **Raghad Al-Amoudi** - `alamoudi.r@northeastern.edu`
- **Atra Alimoradian** - `alimoradian.a@northeastern.edu`
- **Md Tariqul Islam** - `islam.mdtar@northeastern.edu`
- **Full repository:** [github.com/mtariqi/synthetic-variant-calling-benchmark](https://github.com/mtariqi/synthetic-variant-calling-benchmark)
- **Full project README** (extended methods, synthetic data framework, coverage analysis): [README.md](README.md)

---

## References

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
