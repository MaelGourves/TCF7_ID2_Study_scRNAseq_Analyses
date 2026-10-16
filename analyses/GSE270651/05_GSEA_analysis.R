# ==============================================================================
# COPYRIGHT & LICENSE
# ==============================================================================
# Copyright (C) 2026 Gourvès Maël
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# ==============================================================================
# PURPOSE
# ==============================================================================
# This script performs differential-expression ranking and GSEA for two
# comparisons among GMM-defined TCF7/ID2 populations in GSE270651:
#
#   1. ID2+TCF7+ versus ID2-TCF7+
#   2. ID2+TCF7+ versus ID2+TCF7-
#
# Positive log2 fold changes and positive normalized enrichment scores (NES)
# indicate enrichment in ID2+TCF7+ cells. Negative values indicate enrichment
# in the comparison population.
#
# Required local files:
#   - Datasets_GMM_corrected_ID2TCF7_for_Volcano_and_GSEA.rds
#   - 04_GSEA_functions.R
#   - data/My_SelectedPathwaysClonotypes_GSEA.gmt
#
# The corrected RDS is the exact analysis object used in the original
# successful workflow. Differential-expression tables are calculated in
# memory; pre-generated Qlucore-style .txt files are not required.
#
# Cell-level tests are exploratory. Cells from the same patient/time point are
# not independent biological replicates.
# ==============================================================================


# ==============================================================================
# 0. PACKAGES, INPUTS, AND OUTPUT DIRECTORIES
# ==============================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(fgsea)

source("04_GSEA_functions.R")

input_file <- "Datasets_GMM_corrected_ID2TCF7_for_Volcano_and_GSEA.rds"

stopifnot(file.exists(input_file))

DE_OUTDIR <- file.path("output", "DE")
VOLCANO_OUTDIR <- file.path("output", "volcano")
GSEA_OUTDIR <- file.path("output", "GSEA")

dir.create(DE_OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(VOLCANO_OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GSEA_OUTDIR, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. LOAD THE EXACT CORRECTED ANALYSIS OBJECT
# ==============================================================================

df_clean <- readRDS(input_file)

stopifnot(
  is.list(df_clean),
  length(df_clean) >= 3
)

expr_genes_x_cells <- df_clean[[1]]
gene_info <- df_clean[[2]]
cell_anno <- df_clean[[3]]

stopifnot(
  nrow(gene_info) == nrow(expr_genes_x_cells),
  nrow(cell_anno) == ncol(expr_genes_x_cells),
  all(c("Gene", "GeneSymbol") %in% colnames(gene_info)),
  "ID2_TCF7_expression" %in% colnames(cell_anno)
)


# ==============================================================================
# 2. SHARED ANALYSIS FUNCTIONS
# ==============================================================================

run_mannwhitney_de <- function(
    df_list,
    reference_group,
    target_group,
    group_col = "ID2_TCF7_expression",
    expr_scale = c("log2", "linear"),
    pseudocount = 1e-6,
    fdr_method = "BH"
) {
  expr_scale <- match.arg(expr_scale)

  expression <- df_list[[1]]
  genes <- df_list[[2]]
  cells <- df_list[[3]]

  stopifnot(
    ncol(expression) == nrow(cells),
    group_col %in% colnames(cells)
  )

  reference_index <- which(cells[[group_col]] == reference_group)
  target_index <- which(cells[[group_col]] == target_group)

  if (length(reference_index) < 3 || length(target_index) < 3) {
    stop(
      "Too few cells: ", reference_group, " = ", length(reference_index),
      "; ", target_group, " = ", length(target_index)
    )
  }

  reference_mean <- rowMeans(
    expression[, reference_index, drop = FALSE],
    na.rm = TRUE
  )
  target_mean <- rowMeans(
    expression[, target_index, drop = FALSE],
    na.rm = TRUE
  )

  if (expr_scale == "log2") {
    log2_fold_change <- target_mean - reference_mean
    fold_change <- 2^log2_fold_change
  } else {
    fold_change <- (target_mean + pseudocount) /
      (reference_mean + pseudocount)
    log2_fold_change <- log2(fold_change)
  }

  p_value <- apply(expression, 1, function(values) {
    reference_values <- values[reference_index]
    target_values <- values[target_index]

    reference_values <- reference_values[is.finite(reference_values)]
    target_values <- target_values[is.finite(target_values)]

    if (length(reference_values) < 3 || length(target_values) < 3) {
      return(NA_real_)
    }

    suppressWarnings(
      wilcox.test(
        target_values,
        reference_values,
        paired = FALSE,
        exact = FALSE
      )$p.value
    )
  })

  data.frame(
    Gene = genes$Gene,
    GeneSymbol = genes$GeneSymbol,
    FoldChange = fold_change,
    log2FC = log2_fold_change,
    p_value = p_value,
    q_value = p.adjust(p_value, method = fdr_method),
    stringsAsFactors = FALSE
  )
}


make_volcano <- function(
    de_table,
    reference_label,
    target_label,
    q_cutoff = 0.05,
    log2fc_cutoff = 1,
    label_n = 25
) {
  plot_data <- de_table %>%
    mutate(
      q_plot = pmax(ifelse(is.na(q_value), 1, q_value), 1e-300),
      negLog10Q = -log10(q_plot),
      Significant = abs(log2FC) >= log2fc_cutoff & q_value <= q_cutoff,
      Direction = case_when(
        Significant & log2FC > 0 ~ target_label,
        Significant & log2FC < 0 ~ reference_label,
        TRUE ~ "Not significant"
      )
    )

  label_data <- plot_data %>%
    filter(Significant, !is.na(GeneSymbol), GeneSymbol != "") %>%
    arrange(q_value) %>%
    slice_head(n = label_n)

  colors <- c(
    "Not significant" = "grey75",
    setNames("#8080FF", reference_label),
    setNames("#41F0AE", target_label)
  )

  ggplot(plot_data, aes(log2FC, negLog10Q)) +
    geom_point(aes(color = Direction), alpha = 0.6, size = 1.5) +
    geom_vline(
      xintercept = c(-log2fc_cutoff, log2fc_cutoff),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(q_cutoff),
      linetype = "dashed"
    ) +
    geom_text_repel(
      data = label_data,
      aes(label = GeneSymbol),
      size = 3,
      max.overlaps = Inf,
      segment.color = "grey50"
    ) +
    scale_color_manual(values = colors) +
    labs(
      title = paste(target_label, "versus", reference_label),
      x = paste0("log2 fold change (", target_label, " / ", reference_label, ")"),
      y = "-log10(BH-adjusted P value)",
      color = NULL
    ) +
    theme_classic(base_size = 12)
}


build_gsea_ranks <- function(de_table, comparison_label) {
  de_table %>%
    # Preserve the filtering used by the original volcano/GSEA workflow.
    filter(grepl("^ENSG", Gene)) %>%
    transmute(
      gene = GeneSymbol,
      stat = log2FC,
      stat_type = "log2FC",
      source = "DE",
      label = comparison_label,
      p_value = p_value,
      q_value = q_value
    ) %>%
    filter(!is.na(gene), gene != "", is.finite(stat)) %>%
    group_by(gene) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(stat))
}


recode_pathway_names <- function(fgsea_table) {
  fgsea_table %>%
    mutate(
      pathway = recode(
        pathway,
        "GSE24081_CONTROLLER_VS_PROGRESSOR_HIV_SPECIFIC_CD8_TCELL_UP" =
          "HIV-specific CD8+ T cells from controllers",
        "GSE24081_CONTROLLER_VS_PROGRESSOR_HIV_SPECIFIC_CD8_TCELL_DN" =
          "HIV-specific CD8+ T cells from progressors",
        "Central memory (from Pace et al. 2018)" = "Central memory",
        "Effectors (from Pace et al. 2018)" = "Effectors",
        "Naive (from Szabo et al. 2019)" = "Naive",
        "Co-inhibitory molecules (from Chihara et al. 2018)" =
          "Co-inhibitory molecules",
        "Terminal exhausted (from Miller et al. 2019)" =
          "Terminal exhausted"
      )
    )
}


darken_color <- function(color, factor = 1.3) {
  rgb_values <- grDevices::col2rgb(color) / 255
  rgb_values <- rgb_values / factor
  grDevices::rgb(rgb_values[1], rgb_values[2], rgb_values[3])
}


plot_gsea_with_direction <- function(
    fgsea_table,
    reference_label,
    target_label,
    reference_color,
    target_color,
    title
) {
  if (nrow(fgsea_table) == 0) {
    stop("No pathways passed the selected GSEA thresholds.")
  }

  plot_data <- fgsea_table %>%
    mutate(
      NES = as.numeric(NES),
      padj = as.numeric(padj),
      leadingEdge_n = as.integer(leadingEdge_n),
      negLog10FDR = -log10(pmax(padj, 1e-300)),
      Direction = ifelse(NES >= 0, target_label, reference_label)
    ) %>%
    filter(is.finite(NES), is.finite(padj), is.finite(leadingEdge_n)) %>%
    arrange(NES) %>%
    mutate(pathway = factor(pathway, levels = unique(pathway)))

  direction_colors <- c(
    setNames(reference_color, reference_label),
    setNames(target_color, target_label)
  )

  x_limit <- ceiling(max(abs(plot_data$NES), na.rm = TRUE) * 1.1)
  pathway_count <- nrow(plot_data)

  ggplot(plot_data, aes(NES, pathway)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(
      aes(
        size = leadingEdge_n,
        alpha = negLog10FDR,
        color = Direction
      )
    ) +
    scale_color_manual(values = direction_colors, guide = "none") +
    scale_size_continuous(
      name = "Leading-edge genes",
      range = c(2.5, 5.5)
    ) +
    scale_alpha_continuous(
      name = "-log10(FDR)",
      range = c(0.35, 1)
    ) +
    # Direction arrows and labels above the plotting panel.
    annotate(
      "segment",
      x = 0, xend = -x_limit,
      y = pathway_count + 0.75, yend = pathway_count + 0.75,
      arrow = grid::arrow(length = grid::unit(0.12, "in"), type = "closed"),
      linewidth = 0.7,
      color = reference_color
    ) +
    annotate(
      "segment",
      x = 0, xend = x_limit,
      y = pathway_count + 0.75, yend = pathway_count + 0.75,
      arrow = grid::arrow(length = grid::unit(0.12, "in"), type = "closed"),
      linewidth = 0.7,
      color = target_color
    ) +
    annotate(
      "text",
      x = -x_limit / 2,
      y = pathway_count + 1.15,
      label = reference_label,
      fontface = "bold",
      color = reference_color,
      size = 3.2
    ) +
    annotate(
      "text",
      x = x_limit / 2,
      y = pathway_count + 1.15,
      label = target_label,
      fontface = "bold",
      color = target_color,
      size = 3.2
    ) +
    coord_cartesian(
      xlim = c(-x_limit, x_limit),
      ylim = c(0.5, pathway_count + 1.35),
      clip = "off"
    ) +
    labs(
      x = "Normalized enrichment score (NES)",
      y = NULL
    ) +
    theme_classic(base_size = 12) +
    guides(
      alpha = guide_legend(
        order = 1,
        title.position = "left",
        title.hjust = 0,
        nrow = 1,
        byrow = TRUE
      ),
      size = guide_legend(
        order = 2,
        title.position = "left",
        title.hjust = 0,
        nrow = 1,
        byrow = TRUE
      )
    ) +
    theme(
      axis.text.y = element_text(size = 9),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = "left",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 22, r = 8, b = 5, l = 8)
    )
}


run_comparison <- function(
    reference_group,
    target_group = "ID2+TCF7+",
    file_prefix,
    reference_color,
    pathways,
    gsea_fdr = 0.05,
    gsea_abs_nes = 1.5
) {
  comparison_label <- paste(target_group, "vs", reference_group)
  message("Running comparison: ", comparison_label)

  de_table <- run_mannwhitney_de(
    df_list = df_clean,
    reference_group = reference_group,
    target_group = target_group,
    expr_scale = "log2"
  )

  de_file <- file.path(DE_OUTDIR, paste0(file_prefix, "_DE.tsv"))
  write.table(
    de_table,
    de_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  volcano_plot <- make_volcano(
    de_table = de_table,
    reference_label = reference_group,
    target_label = target_group
  )

  ggsave(
    file.path(VOLCANO_OUTDIR, paste0(file_prefix, "_volcano.svg")),
    volcano_plot,
    width = 7,
    height = 5
  )

  rank_table <- build_gsea_ranks(de_table, comparison_label)
  write.csv(
    rank_table,
    file.path(GSEA_OUTDIR, paste0(file_prefix, "_ranked_genes.csv")),
    row.names = FALSE
  )

  fgsea_result <- run_mixed_fgsea(rank_table, pathways)
  raw_fgsea_file <- export_fgsea_csv_raw(
    fgsea_result,
    GSEA_OUTDIR,
    paste0(file_prefix, "_GSEA")
  )

  fgsea_filtered <- coerce_fgsea_types(fgsea_result) %>%
    filter(
      !is.na(padj),
      padj <= gsea_fdr,
      abs(NES) >= gsea_abs_nes
    ) %>%
    recode_pathway_names()

  # CSV files cannot store list-columns. Flatten the leading-edge genes only
  # in the export copy, while retaining the list-column in fgsea_filtered.
  fgsea_filtered_export <- fgsea_filtered
  if (
    "leadingEdge" %in% colnames(fgsea_filtered_export) &&
      is.list(fgsea_filtered_export$leadingEdge)
  ) {
    fgsea_filtered_export$leadingEdge <- vapply(
      fgsea_filtered_export$leadingEdge,
      function(genes) paste(genes, collapse = ";"),
      FUN.VALUE = character(1)
    )
  }

  write.csv(
    fgsea_filtered_export,
    file.path(GSEA_OUTDIR, paste0(file_prefix, "_GSEA_filtered.csv")),
    row.names = FALSE
  )

  gsea_plot <- plot_gsea_with_direction(
    fgsea_table = fgsea_filtered,
    reference_label = reference_group,
    target_label = target_group,
    reference_color = reference_color,
    target_color = darken_color("#41F0AE", 1.3),
    title = paste("GSEA:", comparison_label)
  )

  ggsave(
    file.path(GSEA_OUTDIR, paste0(file_prefix, "_GSEA.svg")),
    gsea_plot,
    width = 8,
    height = 4.5
  )

  invisible(list(
    comparison = comparison_label,
    de = de_table,
    volcano = volcano_plot,
    ranks = rank_table,
    fgsea = fgsea_result,
    fgsea_filtered = fgsea_filtered,
    gsea_plot = gsea_plot,
    raw_fgsea_file = raw_fgsea_file
  ))
}


# ==============================================================================
# 3. LOAD THE EXACT CUSTOM GENE SETS USED IN THE ORIGINAL ANALYSIS
# ==============================================================================

gmt_file <- file.path(
  "data",
  "My_SelectedPathwaysClonotypes_GSEA.gmt"
)

if (!file.exists(gmt_file)) {
  stop(
    "Missing custom gene-set file: ", gmt_file, "\n",
    "Place the exact GMT used for the original analysis in the data folder."
  )
}

pathways_custom <- fgsea::gmtPathways(gmt_file)

required_pathways <- c(
  "GSE24081_CONTROLLER_VS_PROGRESSOR_HIV_SPECIFIC_CD8_TCELL_UP",
  "GSE24081_CONTROLLER_VS_PROGRESSOR_HIV_SPECIFIC_CD8_TCELL_DN",
  "Central memory (from Pace et al. 2018)",
  "Effectors (from Pace et al. 2018)",
  "Terminal exhausted (from Miller et al. 2019)",
  "Co-inhibitory molecules (from Chihara et al. 2018)",
  "IL6-JAK/STAT3 signalling"
)

missing_pathways <- setdiff(required_pathways, names(pathways_custom))
if (length(missing_pathways)) {
  stop(
    "The GMT is not the gene-set file used in the original analysis. ",
    "Missing pathways: ", paste(missing_pathways, collapse = ", ")
  )
}


# ==============================================================================
# 4. COMPARISON 1 — ID2+TCF7+ VERSUS ID2-TCF7+
# ==============================================================================

comparison_1 <- run_comparison(
  reference_group = "ID2-TCF7+",
  target_group = "ID2+TCF7+",
  file_prefix = "ID2high_TCF7high_vs_ID2low_TCF7high",
  reference_color = darken_color("#8080FF", 1.05),
  pathways = pathways_custom
)

p_de_comparison_1 <- comparison_1$gsea_plot
p_de_comparison_1

ggsave(
  filename = "ID2+TCF7+_vs_ID2-TCF7+.svg",
  plot = p_de_comparison_1,
  device = "svg",
  width = 8,
  height = 4.5,
  units = "in"
)

# ==============================================================================
# 5. COMPARISON 2 — ID2+TCF7+ VERSUS ID2+TCF7-
# ==============================================================================

comparison_2 <- run_comparison(
  reference_group = "ID2+TCF7-",
  target_group = "ID2+TCF7+",
  file_prefix = "ID2high_TCF7high_vs_ID2high_TCF7low",
  reference_color = darken_color("#FFA040", 1.15),
  pathways = pathways_custom
)

p_de_comparison_2 <- comparison_2$gsea_plot
p_de_comparison_2

ggsave(
  filename = "ID2+TCF7+_vs_ID2+TCF7-.svg",
  plot = p_de_comparison_2,
  device = "svg",
  width = 8,
  height = 4.5,
  units = "in"
)

# ==============================================================================
# 6. REPRODUCIBILITY INFORMATION
# ==============================================================================

writeLines(
  capture.output(sessionInfo()),
  file.path("output", "sessionInfo_GSEA.txt")
)

message("Finished both comparisons. Outputs are available under output/.")
