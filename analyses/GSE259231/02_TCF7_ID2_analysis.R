# ============================================================
# TCF7/ID2 analysis of HBV-specific CD8+ T cells
# ============================================================
#
# Dataset:
#   Heim et al., "Attenuated effector T cells are linked to control
#   of chronic HBV infection", Nature Immunology (2024)
#   DOI: 10.1038/s41590-024-01928-4
#   GEO accession: GSE259231
#
# Purpose:
#   This script performs a secondary analysis of TCF7 and ID2 expression
#   in HBV polymerase-455-specific CD8+ T cells from endogenous-control
#   chronic HBV, acute HBV, and spontaneously resolved acute HBV.
#
# Method overview:
#   - Extract normalized TCF7 and ID2 expression from the RNA data layer
#   - Standardize each gene independently across all cells
#   - Define four expression states using a Z-score threshold of zero
#   - Project these states onto the authors' deposited UMAP
#   - Calculate cell-level proportions by clinical state and cluster
#   - Compare selected gene-expression patterns between the four states
#
# Dependency:
#   `01_load_prepare_HBV_data.R` must create the prepared Seurat object
#   named `hbv`, including `figure5_phase` and `figure5_cluster` metadata.
#
# Outputs:
#   UMAP, density, composition, and expression dot plots are written to
#   the current working directory.
# ============================================================
setwd("analyses/GSE259231")
source("01_load_prepare_HBV_data.R")


# ============================================================
# PART 1 — Calculate TCF7 and ID2 Z-scores
# ============================================================

# Use normalized RNA expression
DefaultAssay(hbv) <- "RNA"

expression_data <- FetchData(
  hbv,
  vars = c("TCF7", "ID2"),
  layer = "data"
)

# Standardize each gene independently across all cells
expression_z <- scale(expression_data)

# Store Z-scores in the Seurat metadata
hbv$TCF7_z <- as.numeric(expression_z[, "TCF7"])
hbv$ID2_z  <- as.numeric(expression_z[, "ID2"])

# Confirm that they were added
summary(hbv$TCF7_z)
summary(hbv$ID2_z)


library(Seurat)
library(ggplot2)

# ============================================================
# PART 2 — Define TCF7/ID2 populations and colors
# ============================================================

population <- rep("TCF7- ID2-", ncol(hbv))

# TCF7 above average, ID2 below average
population[
  hbv$TCF7_z > 0 &
    hbv$ID2_z < 0
] <- "TCF7+ ID2-"

# ID2 above average, TCF7 below average
population[
  hbv$ID2_z > 0 &
    hbv$TCF7_z < 0
] <- "TCF7- ID2+"

# Both above average
population[
  hbv$TCF7_z > 0 &
    hbv$ID2_z > 0
] <- "TCF7+ ID2+"

hbv$TCF7_ID2_population <- factor(
  population,
  levels = c(
    "TCF7- ID2-",
    "TCF7+ ID2-",
    "TCF7+ ID2+",
    "TCF7- ID2+"
  )
)

# Verify population sizes
table(hbv$TCF7_ID2_population)

# Overall percentages
round(
  prop.table(table(hbv$TCF7_ID2_population)) * 100,
  2
)

# Population colors used throughout the analysis
population_colors <- c(
  "TCF7- ID2-" = "#C4C4C4",
  "TCF7+ ID2-" = "#8C8BFF",
  "TCF7+ ID2+" = "#55E0B2",
  "TCF7- ID2+" = "#FCB468"
)


# ============================================================
# PART 3 — Visualize TCF7 and ID2 expression on the UMAP
# ============================================================

p_TCF7 <- FeaturePlot(
  hbv,
  features = "TCF7_z",
  reduction = "umap",
  order = TRUE,
  min.cutoff = -2,
  max.cutoff = 2
) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    name = "TCF7 Z-score"
  ) +
  ggtitle("TCF7 expression")

p_ID2 <- FeaturePlot(
  hbv,
  features = "ID2_z",
  reduction = "umap",
  order = TRUE,
  min.cutoff = -2,
  max.cutoff = 2
) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    name = "ID2 Z-score"
  ) +
  ggtitle("ID2 expression")

p_TCF7 | p_ID2


ggsave(
  "UMAP_TCF7_ID2_Zscore.png",
  p_TCF7 | p_ID2,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# PART 4 — Visualize TCF7+ID2+ cells and cluster density
# ============================================================

library(Seurat)
library(ggplot2)
library(MASS)

# A. Prepare UMAP coordinates

umap_density <- as.data.frame(
  Embeddings(hbv, reduction = "umap")
)

umap_density$Cluster <- factor(
  hbv$figure5_cluster,
  levels = paste0("C", 0:4)
)

# Identify double-positive cells directly from Z-scores
double_positive_names <- colnames(hbv)[
  hbv$TCF7_z > 0 &
    hbv$ID2_z > 0
]

double_positive_cells <- umap_density[
  rownames(umap_density) %in% double_positive_names,
  ,
  drop = FALSE
]


# B. Define one expanded KDE grid shared by every cluster


x_data_range <- range(umap_density$UMAP_1)
y_data_range <- range(umap_density$UMAP_2)

x_padding <- diff(x_data_range) * 0.15
y_padding <- diff(y_data_range) * 0.15

common_x_limits <- x_data_range + c(-x_padding, x_padding)
common_y_limits <- y_data_range + c(-y_padding, y_padding)

common_kde_limits <- c(
  common_x_limits[1],
  common_x_limits[2],
  common_y_limits[1],
  common_y_limits[2]
)


# C. Calculate manual KDE and extract closed contours


outer_polygon_list <- list()
contour_line_list <- list()

outer_index <- 1
line_index <- 1

for (current_cluster in levels(umap_density$Cluster)) {
  
  current_data <- umap_density[
    umap_density$Cluster == current_cluster,
    ,
    drop = FALSE
  ]
  
  # Calculate KDE on the common expanded grid
  current_kde <- MASS::kde2d(
    x = current_data$UMAP_1,
    y = current_data$UMAP_2,
    n = 400,
    lims = common_kde_limits
  )
  
  # Normalize each cluster density to a maximum of 1
  normalized_density <- current_kde$z / max(current_kde$z)
  
  # Outer filled region
  outer_contours <- contourLines(
    x = current_kde$x,
    y = current_kde$y,
    z = normalized_density,
    levels = 0.12
  )
  
  for (i in seq_along(outer_contours)) {
    outer_polygon_list[[outer_index]] <- data.frame(
      UMAP_1 = outer_contours[[i]]$x,
      UMAP_2 = outer_contours[[i]]$y,
      Cluster = current_cluster,
      Polygon = paste(current_cluster, "outer", i, sep = "_")
    )
    
    outer_index <- outer_index + 1
  }
  
  # Internal density contour lines
  internal_contours <- contourLines(
    x = current_kde$x,
    y = current_kde$y,
    z = normalized_density,
    levels = c(0.2, 0.4, 0.6, 0.8)
  )
  
  for (i in seq_along(internal_contours)) {
    contour_line_list[[line_index]] <- data.frame(
      UMAP_1 = internal_contours[[i]]$x,
      UMAP_2 = internal_contours[[i]]$y,
      Cluster = current_cluster,
      Density = internal_contours[[i]]$level,
      Path = paste(
        current_cluster,
        internal_contours[[i]]$level,
        i,
        sep = "_"
      )
    )
    
    line_index <- line_index + 1
  }
}

outer_polygon_data <- do.call(
  rbind,
  outer_polygon_list
)

contour_line_data <- do.call(
  rbind,
  contour_line_list
)

outer_polygon_data$Cluster <- factor(
  outer_polygon_data$Cluster,
  levels = paste0("C", 0:4)
)

contour_line_data$Cluster <- factor(
  contour_line_data$Cluster,
  levels = paste0("C", 0:4)
)


# D. Define cluster colors


cluster_colors <- c(
  C0 = "#8E79C6",
  C1 = "#3BA3EC",
  C2 = "#E89A6B",
  C3 = "#C8BE72",
  C4 = "#E06B75"
)


# E. Plot closed density regions


# Extract the default UMAP display limits from panel b (NEEDS THE SOURCE CODE)
p_all_build <- ggplot_build(p_all)

default_x_limits <- p_all_build$layout$panel_params[[1]]$x.range
default_y_limits <- p_all_build$layout$panel_params[[1]]$y.range

p_cluster_density_manual <- ggplot() +
  
  # Closed filled cluster regions
  geom_polygon(
    data = outer_polygon_data,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      group = Polygon,
      fill = Cluster
    ),
    alpha = 0.25,
    color = NA
  ) +
  
  # Internal density contours
  geom_path(
    data = contour_line_data,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      group = Path,
      color = Cluster
    ),
    linewidth = 0.7,
    alpha = 0.9,
    lineend = "round",
    linejoin = "round"
  ) +
  
  # TCF7+ ID2+ cells
  geom_point(
    data = double_positive_cells,
    aes(
      x = UMAP_1,
      y = UMAP_2
    ),
    shape = 21,
    fill = "#00C896",
    color = "#15382F",
    stroke = 0.25,
    size = 1.7,
    alpha = 0.45
  ) +
  
  scale_fill_manual(
    values = cluster_colors,
    name = "Cluster"
  ) +
  scale_color_manual(
    values = cluster_colors,
    name = "Cluster"
  ) +
  
  coord_equal(
    xlim = default_x_limits,
    ylim = default_y_limits,
    expand = FALSE,
    clip = "on"
  ) +
  
  labs(
    x = "UMAP_1",
    y = "UMAP_2",
    title = "TCF7+ ID2+ cells across transcriptional states",
    subtitle = "Filled regions and contours represent cluster density"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.ticks.length = grid::unit(2.5, "mm"),
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey35"),
    plot.margin = margin(12, 12, 12, 12)
  )

p_cluster_density_manual


# Save UMAP
ggsave(
  "UMAP_TCF7_ID2_populations_density.svg",
  p_cluster_density_manual,
  width = 8,
  height = 6
)


# ============================================================
# PART 5 — Calculate cluster proportions by clinical state
# ============================================================


# Count cells in each cluster and clinical group
cluster_prop <- as.data.frame(
  table(
    Phase = hbv$figure5_phase,
    Cluster = hbv$figure5_cluster
  )
)

# Calculate proportions within each clinical group
cluster_prop$Proportion <- with(
  cluster_prop,
  Freq / ave(Freq, Phase, FUN = sum)
)

cluster_colors <- c(
  C0 = "#A48AD3",
  C1 = "#1CC5FE",
  C2 = "#FBA27D",
  C3 = "#ECE6CA",
  C4 = "#FB7D80"
)


p_cluster_proportions <- ggplot(
  cluster_prop,
  aes(
    x = Phase,
    y = Proportion,
    fill = Cluster
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(
      label = ifelse(
        Proportion >= 0.04,
        paste0(round(Proportion * 100, 1), "%"),
        ""
      )
    ),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3.5
  ) +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Proportion of cells",
    fill = "Cluster"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "right"
  )

p_cluster_proportions

ggsave(
  "Cluster_p_depending_clinical_status.svg",
  p_cluster_proportions,
  width = 8,
  height = 6
)


# ============================================================
# PART 6 — Calculate TCF7/ID2 population proportions
# ============================================================

# 6A. Proportions within clinical states

phase_counts <- table(
  Phase = hbv$figure5_phase,
  Population = hbv$TCF7_ID2_population
)

# Cell counts
phase_counts

# Row percentages: each clinical group totals 100%
phase_percentages <- prop.table(
  phase_counts,
  margin = 1
) * 100

round(phase_percentages, 2)

# Prepare plotting data
phase_plot_data <- as.data.frame(phase_counts)

phase_plot_data$Proportion <- with(
  phase_plot_data,
  Freq / ave(Freq, Phase, FUN = sum)
)

# Percentage labels
phase_plot_data$Label <- ifelse(
  phase_plot_data$Proportion >= 0.025,
  paste0(round(phase_plot_data$Proportion * 100, 1), "%"),
  ""
)

p_phase_proportions <- ggplot(
  phase_plot_data,
  aes(
    x = Phase,
    y = Proportion,
    fill = Population
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 3.6,
    color = "black"
  ) +
  scale_fill_manual(
    values = population_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Proportion of cells",
    fill = "Population",
    title = "TCF7/ID2 populations by HBV disease phase",
    subtitle = "Z-score threshold = 0"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "right"
  )

p_phase_proportions

# Save clinical-phase proportions
ggsave(
  "TCF7_ID2_proportions_by_HBV_phase.svg",
  p_phase_proportions,
  width = 9,
  height = 6
)

# 6B. Proportions within clusters C0–C4

cluster_counts <- table(
  Cluster = hbv$figure5_cluster,
  Population = hbv$TCF7_ID2_population
)

# Cell counts
cluster_counts

# Row percentages: each cluster totals 100%
cluster_percentages <- prop.table(
  cluster_counts,
  margin = 1
) * 100

round(cluster_percentages, 2)

# Prepare plotting data
cluster_plot_data <- as.data.frame(cluster_counts)

cluster_plot_data$Proportion <- with(
  cluster_plot_data,
  Freq / ave(Freq, Cluster, FUN = sum)
)

# Percentage labels
cluster_plot_data$Label <- ifelse(
  cluster_plot_data$Proportion >= 0.025,
  paste0(round(cluster_plot_data$Proportion * 100, 1), "%"),
  ""
)

p_cluster_population_proportions <- ggplot(
  cluster_plot_data,
  aes(
    x = Cluster,
    y = Proportion,
    fill = Population
  )
) +
  geom_col(
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black"
  ) +
  scale_fill_manual(
    values = population_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Cluster",
    y = "Proportion of cells",
    fill = "Population",
    title = "TCF7/ID2 populations within Figure 5 clusters",
    subtitle = "Z-score threshold = 0"
  ) +
  theme_classic() +
  theme(
    legend.position = "right"
  )

p_cluster_population_proportions

# Save cluster proportions
ggsave(
  "TCF7_ID2_proportions_by_cluster.svg",
  p_cluster_population_proportions,
  width = 9,
  height = 6
)

# ============================================================
# PART 7 — Compare gene expression across TCF7/ID2 populations
# ============================================================

# Set the desired population order
hbv$TCF7_ID2_population <- factor(
  as.character(hbv$TCF7_ID2_population),
  levels = c(
    "TCF7+ ID2-",
    "TCF7+ ID2+",
    "TCF7- ID2+",
    "TCF7- ID2-"
  )
)

# Confirm population sizes
table(hbv$TCF7_ID2_population)

# Use the four populations as identities
Idents(hbv) <- "TCF7_ID2_population"

genes_population_comparison <- c("TCF7", "ID2", "TBX21", "PRDM1", "ZEB2", "RUNX3", "TOX", 
                                 "PDCD1", "LAG3","HMGB1",
                                 "IL7R", "GZMK", "SELL",
                                 "FOS", "FOSB", "JUN", "JUNB", 
                                 "CXCR4", "CD69", "DUSP1", "TSC22D3", "NFKBIA", 
                                 "KLF2", "KLF6", "E2F4", "CDC42", "KLF3"
)

# Retain genes available in the object
genes_population_present <- intersect(
  genes_population_comparison,
  rownames(hbv)
)

# Report missing genes, if any
setdiff(
  genes_population_comparison,
  genes_population_present
)

# Generate the dot plot
p_all_populations <- DotPlot(
  hbv,
  features = genes_population_present,
  group.by = "TCF7_ID2_population",
  scale = TRUE
) +
  RotatedAxis() +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  labs(
    x = NULL,
    y = "Population",
    color = "Average\nexpression",
    size = "Percent\nexpressed",
    title = "Gene expression across TCF7/ID2 populations",
    subtitle = "All cells; populations defined using a Z-score threshold of 0"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 9
    ),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold")
  )

p_all_populations

ggsave(
  "TCF7_ID2_DEG.svg",
  p_all_populations,
  width = 9,
  height = 6,
  dpi = 300
)

