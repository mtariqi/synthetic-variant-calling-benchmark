##############################################
# R_Analysis_Figures
##############################################

library(data.table)
library(ggplot2)
library(gridExtra)
library(ggVennDiagram)

setwd("~/Desktop/vaf_python")

# Color palette
colors <- list(
  gatk = "#E91E63",        # Hot Pink
  deepvariant = "#2196F3", # Blue
  normal = "#9C27B0",      # Purple
  tumor = "#FF5722",       # Orange
  shared = "#4CAF50",      # Green
  additional = "#FFC107"   # Yellow
)

##############################################
# CORRECTED DATA (FILTERED COUNTS)
##############################################

# Key data - FILTERED GATK counts
normal_gatk <- 172      # NOT 179
normal_dv <- 197
normal_shared <- 172    # All GATK confirmed by DV
normal_dv_only <- 197 - 172  # = 25, NOT 18

tumor_gatk <- 377       # NOT 383
tumor_dv <- 408
tumor_shared <- 377     # All GATK confirmed by DV
tumor_dv_only <- 408 - 377   # = 31, NOT 25

##############################################
# FIGURE 1 — VARIANT COUNTS
##############################################

# Create data
mapping_df <- data.frame(
  Sample = c("Normal", "Tumor"),
  PercentMapped = c(100, 99)
)

# Create the plot with sky blue for Normal
p1 <- ggplot(mapping_df, aes(x = Sample, y = PercentMapped, fill = Sample)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("Normal" = "#87CEEB", "Tumor" = "#C73E3E")) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20), 
                     labels = paste0(seq(0, 100, 20), "%")) +
  labs(title = "Alignment QC: Mapping Rate per Sample\n(Consistent High Rates)") +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "none"
  )

plot(p1)

# Also save it
png("~/Desktop/FIGURE1_TwoBars.png", width = 10, height = 7, units = "in", res = 300)
print(p1)
dev.off()

cat("File saved with sky blue for Normal!\n")

##############################################
# FIGURE 2 — CONCORDANCE
##############################################

concordance_df <- data.frame(
  Sample = c("Synthetic NORMAL", "Synthetic TUMOR"),
  Concordance = c(
    (normal_shared / normal_dv) * 100,    # 172/197 = 87.3%
    (tumor_shared / tumor_dv) * 100       # 377/408 = 92.4%
  ),
  Label = c(
    paste0(round((normal_shared / normal_dv) * 100, 1), "%\n(", normal_shared, "/", normal_dv, ")"),
    paste0(round((tumor_shared / tumor_dv) * 100, 1), "%\n(", tumor_shared, "/", tumor_dv, ")")
  )
)

p2 <- ggplot(concordance_df, aes(x = Sample, y = Concordance, fill = Sample)) +
  geom_col(color = "black", width = 0.6, linewidth = 1.2) +
  geom_text(aes(label = Label), vjust = -0.3, 
            size = 6, fontface = "bold") +
  scale_fill_manual(values = c("Synthetic NORMAL" = colors$normal, 
                               "Synthetic TUMOR" = colors$tumor)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 16) +
  labs(
    title = "Pipeline Concordance Rate",
    subtitle = "% of DeepVariant calls confirmed by GATK",
    y = "Concordance (%)",
    x = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "gray40"),
    legend.position = "none",
    axis.text = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_line(color = "gray90")
  )

ggsave("Figure2_Concordance_CORRECTED.pdf", p2, width = 8, height = 6, dpi = 300)
cat("✓ Figure 2 (corrected) saved\n")

##############################################
# FIGURE 3 — OVERLAP
##############################################

overlap_df <- data.frame(
  Sample = rep(c("Synthetic NORMAL", "Synthetic TUMOR"), each = 2),
  Category = rep(c("Shared (Both Callers)", "DeepVariant Only"), 2),
  Count = c(normal_shared, normal_dv_only, tumor_shared, tumor_dv_only)
)

overlap_df$Category <- factor(overlap_df$Category, 
                              levels = c("Shared (Both Callers)", "DeepVariant Only"))

p3 <- ggplot(overlap_df, aes(x = Sample, y = Count, fill = Category)) +
  geom_col(color = "black", width = 0.6, linewidth = 1.2) +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            size = 7, fontface = "bold") +
  scale_fill_manual(
    values = c("Shared (Both Callers)" = colors$shared,
               "DeepVariant Only" = colors$additional)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic(base_size = 16) +
  labs(
    title = "Variant Detection Breakdown",
    subtitle = "100% of GATK variants confirmed by DeepVariant",
    y = "Number of Variants",
    x = NULL,
    fill = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "gray40"),
    legend.position = "top",
    legend.text = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_line(color = "gray90")
  )

ggsave("Figure3_Overlap_CORRECTED.pdf", p3, width = 9, height = 6, dpi = 300)
cat("✓ Figure 3 (corrected) saved\n")

##############################################
# FIGURE 4 — VENN DIAGRAMS
##############################################

# Install/load packages
if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  install.packages("ggalluvial")
}
library(ggplot2)
library(ggalluvial)

# Data with FILTERED counts
alluvial_final <- data.frame(
  Sample    = c("NORMAL", "NORMAL", "TUMOR", "TUMOR"),
  Detection = c("Both Callers", "DeepVariant Only",
                "Both Callers", "DeepVariant Only"),
  Count     = c(172, 25, 377, 31),
  Flow_ID   = c("Normal_Both", "Normal_DV", "Tumor_Both", "Tumor_DV")
)

# Verify format
stopifnot(is_alluvia_form(alluvial_final, axes = 1:2, weight = Count))

# Manual label positions to avoid overlap
label_positions <- data.frame(
  x = c(1.5, 1.5, 1.5, 1.5),
  y = c(533, 417, 288, 16),  # Manually positioned
  label = c("172", "25", "377", "31")
)

# Create plot
p_alluvial <- ggplot(
  alluvial_final,
  aes(axis1 = Sample, axis2 = Detection, y = Count)
) +
  geom_alluvium(
    aes(fill = Flow_ID), 
    width = 1/25,
    alpha = 0.8,
    curve_type = "sigmoid"
  ) +
  geom_stratum(
    width = 1/5,
    fill = "#FFF9E6",
    color = "#333333",
    linewidth = 1.2
  ) +
  geom_text(
    stat = "stratum",
    aes(label = after_stat(stratum)),
    size = 4,
    fontface = "bold",
    color = "#1a1a1a"
  ) +
  # MANUAL LABELS - positioned to avoid overlap
  geom_label(
    data = label_positions,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold",
    fill = "#FFFEF7",
    label.padding = unit(0.2, "lines"),
    linewidth = 0.3,
    color = "#1a1a1a",
    alpha = 0.95
  ) +
  
  # YOUR ORIGINAL COLORS
  scale_fill_manual(
    name = NULL,
    values = c(
      "Normal_Both"  = "#2E7D32",
      "Normal_DV"    = "#FFA726",
      "Tumor_Both"   = "#1976D2",
      "Tumor_DV"     = "#D32F2F"
    ),
    labels = c(
      "Normal_Both"  = "NORMAL to Both Callers",
      "Normal_DV"    = "NORMAL to DeepVariant Only",
      "Tumor_Both"   = "TUMOR to Both Callers",
      "Tumor_DV"     = "TUMOR to DeepVariant Only"
    )
  ) +
  
  scale_x_discrete(
    limits = c("Sample", "Detection"),
    expand = c(0.15, 0.15)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 600, 200),
    limits = c(0, 620),
    expand = c(0, 0)
  ) +
  
  theme_minimal(base_size = 16) +
  labs(
    title    = "Variant Detection Flow: DeepVariant Captures All GATK Variants",
    subtitle = "All GATK variants confirmed by DeepVariant | DeepVariant shows higher sensitivity",
    y        = "Number of Variants",
    x        = NULL,
    caption  = "Flow counts represent filtered variants. 'Both Callers' = variants detected by both GATK and DeepVariant.\n'DeepVariant Only' = additional variants unique to DeepVariant, demonstrating its increased sensitivity."
  ) +
  theme(
    plot.background = element_rect(fill = "#FDFCF9", color = NA),
    panel.background = element_rect(fill = "#FDFCF9", color = NA),
    panel.spacing = unit(0, "lines"),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 17, 
                                   color = "#1a1a1a",
                                   margin = margin(b = 5)),
    plot.subtitle   = element_text(hjust = 0.5, size = 11.5, color = "#3a3a3a",
                                   margin = margin(b = 20)),
    plot.caption    = element_text(hjust = 0, size = 10, color = "#4a4a4a",
                                   margin = margin(t = 15), lineheight = 1.3),
    legend.position = "bottom",
    legend.background = element_rect(fill = "#FDFCF9", color = NA),
    legend.text     = element_text(size = 10.5, lineheight = 1.2, color = "#1a1a1a"),
    legend.key.size = unit(0.7, "cm"),
    axis.text.y     = element_text(size = 12, face = "bold", color = "#1a1a1a"),
    axis.title.y    = element_text(face = "bold", size = 14, color = "#1a1a1a"),
    panel.grid      = element_blank(),
    axis.text.x     = element_text(size = 13, face = "bold", color = "#1a1a1a",
                                   margin = margin(t = 10)),
    plot.margin     = margin(20, 20, 20, 10),
    axis.ticks.length.y = unit(0, "pt")
  ) +
  coord_cartesian(clip = "off")

print(p_alluvial)

ggsave("Figure4_Alluvial_Final.pdf", p_alluvial, 
       width = 12, height = 8, dpi = 300)

ggsave("Figure4_Alluvial_Final.png", p_alluvial, 
       width = 12, height = 8, dpi = 300, bg = "#FDFCF9")

##############################################
# SUMMARY TABLE 
##############################################

summary_df <- data.frame(
  Sample = c("Synthetic NORMAL", "Synthetic TUMOR"),
  `GATK (filtered)` = c(normal_gatk, tumor_gatk),
  `DeepVariant` = c(normal_dv, tumor_dv),
  `Shared` = c(normal_shared, tumor_shared),
  `DV Additional` = c(
    paste0("+", normal_dv_only, " (+", round((normal_dv_only/normal_gatk)*100, 1), "%)"),
    paste0("+", tumor_dv_only, " (+", round((tumor_dv_only/tumor_gatk)*100, 1), "%)")
  ),
  `Concordance` = c(
    paste0(round((normal_shared/normal_dv)*100, 1), "%"),
    paste0(round((tumor_shared/tumor_dv)*100, 1), "%")
  ),
  check.names = FALSE
)

print(summary_df)

cat("\n==========================================\n")
cat("✓ ALL CORRECTED FIGURES GENERATED!\n")
cat("==========================================\n")
cat("\nCORRECTED Key Findings:\n")
cat("• NORMAL: DeepVariant detected", normal_dv_only, "more variants (+", 
    round((normal_dv_only/normal_gatk)*100, 1), "%)\n")
cat("• TUMOR: DeepVariant detected", tumor_dv_only, "more variants (+", 
    round((tumor_dv_only/tumor_gatk)*100, 1), "%)\n")
cat("• Concordance: 87.3% (NORMAL), 92.4% (TUMOR)\n")
cat("• All GATK variants confirmed by DeepVariant\n")
cat("==========================================\n")