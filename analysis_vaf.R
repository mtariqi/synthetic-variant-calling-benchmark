#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(dplyr)
  library(ggplot2)
})

vcf_dir <- "results_multi"
vcf_files <- list.files(vcf_dir, pattern = "\\.(vcf|vcf.gz)$", full.names = TRUE)
if (length(vcf_files) == 0) stop("No VCF files found in results_multi/")

read_vcf_to_df <- function(path) {
  vcf <- readVcf(path)
  sample_name <- sub("\\.(vcf|vcf.gz)$", "", basename(path))
  rng <- rowRanges(vcf)
  ad  <- geno(vcf)$AD
  ref_counts <- sapply(ad, function(x) x[1])
  alt_counts <- sapply(ad, function(x) x[2])
  df <- data.frame(
    sample = sample_name,
    chrom  = as.character(seqnames(rng)),
    pos    = start(rng),
    ref    = as.character(ref(vcf)),
    alt    = as.character(unlist(alt(vcf))),
    refCount = ref_counts,
    altCount = alt_counts
  )
  df$depth <- df$refCount + df$altCount
  df$VAF   <- ifelse(df$depth > 0, df$altCount / df$depth, NA_real_)
  df
}

vcf_df <- bind_rows(lapply(vcf_files, read_vcf_to_df))
vcf_df$class <- ifelse(grepl("^NORMAL", vcf_df$sample), "NORMAL", "TUMOR")
vcf_df_filt <- filter(vcf_df, depth >= 5)

print(vcf_df_filt %>% group_by(class) %>%
        summarise(n_sites = n(),
                  mean_VAF = mean(VAF, na.rm = TRUE),
                  median_VAF = median(VAF, na.rm = TRUE)))

p <- ggplot(vcf_df_filt, aes(x = VAF, fill = class)) +
  geom_density(alpha = 0.4) +
  xlim(0, 1) +
  theme_minimal() +
  labs(title = "Tumor vs Normal VAF distributions",
       x = "Variant Allele Frequency (VAF)",
       y = "Density")
ggsave("VAF_density_tumor_normal.pdf", p, width = 6, height = 4)
