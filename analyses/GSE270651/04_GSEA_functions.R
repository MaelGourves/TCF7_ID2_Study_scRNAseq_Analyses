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
#
# This script defines helper functions which may be used by `05_GSEA_analysis.R` to:
#
# - Retrieve Reactome, Gene Ontology Biological Process, and Hallmark
#   gene sets through msigdbr
# - Construct ranked gene vectors
# - Run fgsea
# - Plot significant enrichment results
# - Export and import fgsea result tables
#
# This script defines functions only. It does not load the expression
# dataset or run a GSEA comparison by itself.
# ==============================================================================

library(dplyr)
library(fgsea)
library(msigdbr)
library(stringr)
library(ggplot2)


# Gene set loaders (Reactome + GO BP + Hallmark) and mixer

get_msig_sets <- function(type = c("reactome", "go_bp", "hallmark"),
                          species = "Homo sapiens") {
  type <- match.arg(type)
  
  if (type == "reactome") {
    m <- msigdbr(species = species, collection = "C2", subcollection = "CP:REACTOME")
    sets <- split(m$gene_symbol, paste0("Reactome ", str_replace_all(m$gs_name, "^REACTOME_", "") %>% str_replace_all("_"," ")))
  }
  
  if (type == "go_bp") {
    m <- msigdbr(species = species, collection = "C5", subcollection = "GO:BP")
    sets <- split(m$gene_symbol, paste0("GO ", m$gs_name %>% str_replace_all("^GOBP_", "") %>% str_replace_all("_"," ")))
  }
  
  if (type == "hallmark") {
    m <- msigdbr(species = species, collection = "H")
    sets <- split(m$gene_symbol, paste0("Hallmark ", m$gs_name %>% str_replace_all("^HALLMARK_", "") %>% str_replace_all("_"," ")))
  }
  
  # ensure unique genes per set
  sets <- lapply(sets, unique)
  return(sets)
}

# Combine multiple gene-set lists into one named list
combine_gene_sets <- function(...) {
  sets_list <- list(...)
  sets <- do.call(c, sets_list)
  # drop any empty sets
  sets <- sets[lengths(sets) > 0]
  return(sets)
}

# Build ranks from the standardized rank table

make_ranks_from_ranktable <- function(rank_tbl) {
  stopifnot(all(c("gene","stat") %in% colnames(rank_tbl)))
  
  x <- rank_tbl %>%
    filter(!is.na(gene), gene != "", is.finite(stat)) %>%
    group_by(gene) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(stat))
  
  ranks <- x$stat
  names(ranks) <- x$gene
  ranks
}

# Run fgsea on mixed gene sets (+ optional name filtering)

# drop_regex <- "(?i)HIV|SARS|INFECTION|VIRAL|ADME|DRUG|PARACETAMOL" IF FILTERING IS REQUIRED

run_mixed_fgsea <- function(rank_tbl,
                            pathways,
                            nperm = 10000,
                            minSize = 10,
                            maxSize = 500,
                            drop_name_regex = NULL) { # FILTERING IF REQUIRED
  ranks <- make_ranks_from_ranktable(rank_tbl)
  
  fg <- fgsea(
    pathways = pathways,
    stats = ranks,
    nperm = nperm,
    minSize = minSize,
    maxSize = maxSize
  ) %>%
    arrange(padj) %>%
    mutate(
      leadingEdge_n = lengths(leadingEdge),
      negLog10FDR = -log10(padj),
      Direction = ifelse(NES >= 0, "Positive", "Negative")
    )
  
  if (!is.null(drop_name_regex)) {
    fg <- fg %>% filter(!str_detect(pathway, drop_name_regex))
  }
  
  fg
}

# Dotplot (top pos + top neg)

plot_fgsea_dot_significant <- function(
    fg_res,
    fdr_max = 0.05,
    max_terms_per_source = 50,   # set to Inf to have everything
    title = "Mixed GSEA",
    subtitle = NULL
) {
  # keep only significant
  df <- fg_res %>%
    dplyr::filter(!is.na(padj), padj <= fdr_max) %>%
    dplyr::mutate(source = stringr::str_extract(pathway, "^(Reactome|GO|Hallmark)")) %>%
    dplyr::filter(!is.na(source)) %>%
    dplyr::arrange(padj)
  
  if (nrow(df) == 0) stop("No pathways pass the FDR cutoff.")
  
  # optional: cap number shown per source (keeps GO from dominating)
  if (is.finite(max_terms_per_source)) {
    df <- df %>%
      dplyr::group_by(source, Direction) %>%   # Direction is already in your fg_res
      dplyr::slice_head(n = max_terms_per_source) %>%
      dplyr::ungroup()
  }
  
  # nicer label: remove prefix for printing, but keep source in a separate column
  df <- df %>%
    dplyr::mutate(
      pathway_clean = stringr::str_replace(pathway, "^(Reactome|GO|Hallmark)\\s+", "")
    )
  
  ggplot(df, aes(x = NES, y = reorder(pathway_clean, NES))) +
    geom_point(aes(size = leadingEdge_n, alpha = negLog10FDR)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    coord_cartesian(clip = "off") +
    scale_y_discrete(expand = expansion(add = 0.8)) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      title = title,
      subtitle = if (!is.null(subtitle)) subtitle else paste0("All significant pathways (FDR ≤ ", fdr_max, ")"),
      size = "Leading-edge genes",
      alpha = "-log10(FDR)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.margin = margin(t = 10, r = 25, b = 10, l = 55),  # big left margin
      axis.text.y = element_text(size = 7, lineheight = 0.95, margin = margin(r = 10))
    )
}

# =========================
# Export / import helpers for fgsea tables
# =========================

export_fgsea_csv_raw <- function(fg_df, outdir, prefix) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  fg_out <- fg_df
  
  # Flatten list-columns for CSV export
  if ("leadingEdge" %in% colnames(fg_out)) {
    fg_out$leadingEdge <- vapply(
      fg_out$leadingEdge,
      function(x) paste(x, collapse = ";"),
      FUN.VALUE = character(1)
    )
  }
  
  out_file <- file.path(outdir, paste0(prefix, "_raw.csv"))
  utils::write.csv(fg_out, out_file, row.names = FALSE, quote = TRUE)
  
  message("Raw fgsea results written to: ", out_file)
  invisible(out_file)
}


import_fgsea_csv <- function(file) {
  df <- utils::read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Rebuild leadingEdge as list-column if present
  if ("leadingEdge" %in% colnames(df)) {
    df$leadingEdge <- strsplit(ifelse(is.na(df$leadingEdge), "", df$leadingEdge), ";", fixed = TRUE)
  }
  
  df
}


coerce_fgsea_types <- function(fg_df) {
  out <- fg_df
  
  num_cols <- intersect(c("pval", "padj", "ES", "NES", "size", "leadingEdge_n", "negLog10FDR"), colnames(out))
  int_cols <- intersect(c("size", "leadingEdge_n"), colnames(out))
  
  for (nm in num_cols) {
    out[[nm]] <- suppressWarnings(as.numeric(out[[nm]]))
  }
  
  for (nm in int_cols) {
    out[[nm]] <- suppressWarnings(as.integer(out[[nm]]))
  }
  
  if ("leadingEdge" %in% colnames(out) && !is.list(out$leadingEdge)) {
    out$leadingEdge <- strsplit(ifelse(is.na(out$leadingEdge), "", out$leadingEdge), ";", fixed = TRUE)
  }
  
  if (!"negLog10FDR" %in% colnames(out) && "padj" %in% colnames(out)) {
    out$negLog10FDR <- -log10(pmax(out$padj, 1e-300))
  }
  
  if (!"Direction" %in% colnames(out) && "NES" %in% colnames(out)) {
    out$Direction <- ifelse(out$NES >= 0, "Positive", "Negative")
  }
  
  out
}


plot_fgsea_from_csv <- function(
    fg_df,
    title = "Mixed GSEA",
    subtitle = NULL,
    fdr_max = NULL,
    max_terms = Inf
) {
  df <- coerce_fgsea_types(fg_df)
  
  needed <- c("pathway", "NES", "padj", "leadingEdge_n")
  stopifnot(all(needed %in% colnames(df)))
  
  df <- df %>%
    dplyr::filter(
      !is.na(pathway),
      is.finite(NES),
      is.finite(padj),
      is.finite(leadingEdge_n)
    )
  
  if (!is.null(fdr_max)) {
    df <- df %>% dplyr::filter(padj <= fdr_max)
  }
  
  if (nrow(df) == 0) stop("No pathways available to plot after filtering.")
  
  if (!"negLog10FDR" %in% colnames(df)) {
    df$negLog10FDR <- -log10(pmax(df$padj, 1e-300))
  }
  
  if (!"Direction" %in% colnames(df)) {
    df$Direction <- ifelse(df$NES >= 0, "Positive", "Negative")
  }
  
  if (is.finite(max_terms) && nrow(df) > max_terms) {
    df <- df %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = max_terms)
  }
  
  df <- df %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(pathway = factor(pathway, levels = unique(pathway)))
  
  ggplot(df, aes(x = NES, y = pathway)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point(aes(size = leadingEdge_n, alpha = negLog10FDR, color = Direction)) +
    scale_color_manual(values = c("Positive" = "firebrick3", "Negative" = "dodgerblue3")) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = "Leading-edge genes",
      alpha = "-log10(FDR)",
      color = "Direction"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )
}