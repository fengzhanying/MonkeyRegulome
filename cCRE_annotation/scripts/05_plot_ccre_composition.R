#!/usr/bin/env Rscript

library(ggplot2)
library(scales)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1 && args[1] == "--help") {
  cat("Usage: Rscript scripts/05_plot_ccre_composition.R <class_stats.txt> <output_dir>\n")
  quit(status = 0)
}
if (length(args) == 2) {
  stats_file <- args[1]
  out_dir <- args[2]
} else if (length(args) == 0 && Sys.getenv("CCRE_WORKDIR") != "") {
  stats_file <- file.path(Sys.getenv("CCRE_WORKDIR"), "04.classification", "cCRE_class_stats.txt")
  out_dir <- dirname(stats_file)
} else {
  stop("Supply class_stats.txt and output_dir, or set CCRE_WORKDIR. Use --help.")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

plot_df <- read.table(stats_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (ncol(plot_df) != 2) stop("Expected a two-column class-count table.")
colnames(plot_df) <- c("Class", "Count")
if (anyDuplicated(plot_df$Class) || any(!is.finite(plot_df$Count)) ||
    any(plot_df$Count < 0) || any(plot_df$Count %% 1 != 0)) {
  stop("Invalid class counts.")
}

class_name_map <- c(
  "PLS" = "Promoter",
  "pELS" = "Proximal enhancer",
  "dELS" = "Distal enhancer",
  "CA-H3K4me3" = "CA-H3K4me3",
  "CA-CTCF" = "CA-CTCF",
  "CA" = "CA"
)

if (!all(plot_df$Class %in% names(class_name_map))) stop("Unknown cCRE class.")
plot_df$Class <- class_name_map[plot_df$Class]
total_ccre <- sum(plot_df$Count)
if (total_ccre == 0) stop("The total class count must be positive.")
plot_df$Percentage <- plot_df$Count / total_ccre

write.csv(
  plot_df,
  file.path(out_dir, "Monkey_cCRE_Composition_Stats.csv"),
  row.names = FALSE
)

ccre_colors <- c(
  "Promoter" = "#FF0000",
  "Proximal enhancer" = "#FFA500",
  "Distal enhancer" = "#FFD700",
  "CA-H3K4me3" = "#FFB6C1",
  "CA-CTCF" = "#00BFFF",
  "CA" = "#00C957"
)

class_order <- c(
  "Promoter",
  "Proximal enhancer",
  "Distal enhancer",
  "CA-H3K4me3",
  "CA-CTCF",
  "CA"
)

plot_df$Class <- factor(plot_df$Class, levels = rev(class_order))
plot_df$Label <- paste0(
  format(round(plot_df$Count / 1000, 1), nsmall = 1),
  "K (",
  round(plot_df$Percentage * 100, 1),
  "%)"
)

p <- ggplot(plot_df, aes(x = Count / 1000, y = Class, fill = Class)) +
  geom_bar(stat = "identity", width = 0.65, color = NA) +
  geom_text(aes(label = Label), hjust = -0.08, size = 3.5, color = "black") +
  scale_fill_manual(values = ccre_colors) +
  scale_x_continuous(
    labels = comma_format(),
    expand = expansion(mult = c(0, 0.3)),
    name = "No. of macaque cCREs (x1000)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4)
  ) +
  labs(
    title = "Composition of macaque cCREs",
    subtitle = paste0("Total: ", format(total_ccre, big.mark = ","), " cCREs")
  )

ggsave(file.path(out_dir, "Monkey_cCRE_Composition.pdf"),
       p, width = 6.5, height = 2.5, useDingbats = FALSE)
ggsave(file.path(out_dir, "Monkey_cCRE_Composition.png"),
       p, width = 6.5, height = 2.5, dpi = 300)
