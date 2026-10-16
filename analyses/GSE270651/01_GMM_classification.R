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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ==============================================================================
# DATA SOURCES & INPUT REQUIREMENTS
# ==============================================================================
# Original Data Source:
# Gene Expression Omnibus (GEO) Accession: GSE270651
# Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270651
# Original File: GSE270651_FilteredCells_Normalized_sceobject_Processed.rds.gz
#
# Original Authors: White, E., Papagno, L., Samri, A. et al. Clonal succession after prolonged antiretroviral therapy rejuvenates CD8+ T cell responses against HIV-1. 
# Nat Immunol 25, 1555–1564 (2024). https://doi.org/10.1038/s41590-024-01931-9
#
# Processed Input File for this Script:
# This script requires 'Datasets.rds', a Qlucore-exported file with specific structure.
# The original GEO file cannot be used directly without conversion.
#
# For reproducibility, the processed file is archived on Zenodo:
# DOI: 10.5281/zenodo.20854108
# Download: 
#
# Instructions:
# 1. Download 'Datasets.rds' from Zenodo
# 2. Place it in the same directory as this script
# 3. Run the script
# ==============================================================================


# Packages needed
library(mclust)
library(dplyr)
library(ggplot2)

# Importing the dataset
df <- readRDS("Datasets.rds") 

expr_genes_x_cells <- df[[1]]   # 14697 x 518 (numeric)
gene_info          <- df[[2]]   # 14697 x 2 (Gene, GeneSymbol)
cell_anno          <- df[[3]]   # 518 x 35 (annotations)

# Transposition 
expr_cells_x_genes <- t(expr_genes_x_cells)   
expr_cells_x_genes <- as.data.frame(expr_cells_x_genes)

# use GeneSymbol when available, otherwise fall back to Ensembl
symbols <- gene_info$GeneSymbol
ensembl <- gene_info$Gene

good_symbol <- !is.na(symbols) & symbols != ""
new_names <- ifelse(good_symbol, symbols, ensembl)

colnames(expr_cells_x_genes) <- new_names

# Merging annotations + expression into one dataframe
full_df <- cbind(cell_anno, expr_cells_x_genes)

# GMM cutoff for ID2 and TCF7
gmm_ID2 <- Mclust(full_df$ID2, G = 2)
summary(gmm_ID2)

gmm_TCF7 <- Mclust(full_df$TCF7, G = 2)
summary(gmm_TCF7)

gmm_cutoff <- function(gmm, x, grid_n = 5000) {
  x <- x[is.finite(x)]
  mu <- gmm$parameters$mean
  pi <- gmm$parameters$pro
  sigmasq <- gmm$parameters$variance$sigmasq
  
  if (length(sigmasq) == 1) sigmasq <- rep(sigmasq, length(mu))
  sd <- sqrt(sigmasq)
  
  ord <- order(mu)
  mu <- mu[ord]; sd <- sd[ord]; pi <- pi[ord]
  
  grid <- seq(min(x), max(x), length.out = grid_n)
  d1 <- pi[1] * dnorm(grid, mu[1], sd[1])
  d2 <- pi[2] * dnorm(grid, mu[2], sd[2])
  post2 <- d2 / (d1 + d2)
  
  grid[which.min(abs(post2 - 0.5))]
}

ID2_cutoff  <- gmm_cutoff(gmm_ID2,  full_df$ID2)
TCF7_cutoff <- gmm_cutoff(gmm_TCF7, full_df$TCF7)

# ID2 TCF7 GMM Cutoff Figure
p_GMM_cutoff <- ggplot(full_df, aes(x = TCF7, y = ID2)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_vline(xintercept = TCF7_cutoff, linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = ID2_cutoff, linetype = "dashed", linewidth = 1) +
  labs(
    x = "TCF7 expression",
    y = "ID2 expression",
    title = "ID2 vs TCF7 expression (all cells)"
  ) +
  theme_classic()

p_GMM_cutoff

ggsave(
  filename = "GMM_cutoff.svg",
  plot = p_GMM_cutoff,
  dpi = 300,
  width = 3.5,
  height = 3.5,
  units = "in"
)

# Classify between ID2+/ID2- and TCF7+/TCF7-
full_df <- full_df %>%
  mutate(
    ID2_status  = ifelse(ID2  >= ID2_cutoff,  "ID2+",  "ID2-"),
    TCF7_status = ifelse(TCF7 >= TCF7_cutoff, "TCF7+", "TCF7-"),
    ID2_TCF7_expression = paste0(ID2_status, TCF7_status)
  )

full_df$ID2_TCF7_expression <- factor(
  full_df$ID2_TCF7_expression,
  levels = c("ID2+TCF7+", "ID2+TCF7-", "ID2-TCF7+", "ID2-TCF7-")
)

saveRDS(full_df, "GMM_datasets.rds")


# Number of cells per subpop in each patients/ at each timepoints
cell_counts_long <- full_df %>%
  dplyr::count(Patient, TP, ID2_TCF7_expression, name = "n_cells")

write.csv(
  cell_counts_long,
  file = "cell_counts_byPatient_byTP_bySubpop.csv",
  row.names = FALSE
)

# Save session information for reproducibility
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


