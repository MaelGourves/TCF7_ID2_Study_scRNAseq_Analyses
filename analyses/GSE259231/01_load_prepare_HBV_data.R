# ============================================================
# Load and prepare the Heim et al. (2024) HBV scRNA-seq dataset
# ============================================================
#
# Publication:
# Heim et al., "Attenuated effector T cells are linked to control
# of chronic HBV infection", Nature Immunology (2024)
#
# DOI: 10.1038/s41590-024-01928-4
# GEO accession: GSE259231
#
# This script:
#   1. Downloads the processed Seurat object from GEO if necessary
#   2. Loads the double-gzipped RDS file
#   3. Uses the authors' deposited UMAP coordinates
#   4. Assigns the published C0-C4 cluster identities
#   5. Annotates the three clinical states used in the analysis
#
# This script does not repeat normalization, dimensionality reduction,
# integration, clustering, differential-expression testing, or SCENIC.
#
# Output:
#   A prepared Seurat object named `hbv`, used by
#   `02_TCF7_ID2_analysis.R`.
# ============================================================

options(stringsAsFactors = FALSE)

needed <- c("Seurat", "ggplot2", "patchwork")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install packages first: install.packages(c(",
       paste(sprintf("'%s'", missing), collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

dir.create("data", showWarnings = FALSE)
dir.create("results/figure5", recursive = TRUE, showWarnings = FALSE)

geo_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE259nnn/",
  "GSE259231/suppl/GSE259231_seurat_object_combined.rds.gz"
)
object_file <- "data/GSE259231_seurat_object_combined.rds.gz"
if (!file.exists(object_file)) {
  message("Downloading the 41 MB processed Seurat object from GEO...")
  download.file(geo_url, object_file, mode = "wb", method = "libcurl")
}

# Although GEO labels this as an RDS file, some submissions are actually an
# R workspace produced by save(). Try true RDS first, then load() the workspace
# and select its single Seurat object.
hbv <- tryCatch(
  readRDS(gzfile(object_file)),
  error = function(rds_error) {
    message("Single-gzip readRDS failed; testing GEO's double-gzip packaging...")
    double_gzip_result <- tryCatch({
      outer_connection <- gzfile(object_file, open = "rb")
      inner_connection <- gzcon(outer_connection)
      on.exit(close(inner_connection), add = TRUE)
      readRDS(inner_connection)
    }, error = identity)
    
    if (!inherits(double_gzip_result, "error")) {
      message("Loaded the double-gzipped RDS successfully.")
      return(double_gzip_result)
    }
    
    message("Double-gzip read failed; testing whether the file is an R workspace...")
    imported <- new.env(parent = emptyenv())
    object_names <- tryCatch(
      load(object_file, envir = imported),
      error = function(load_error) {
        stop(
          "The downloaded file is not a readable single/double-gzip RDS or workspace.\n",
          "single-gzip readRDS: ", conditionMessage(rds_error), "\n",
          "double-gzip readRDS: ", conditionMessage(double_gzip_result), "\n",
          "load: ", conditionMessage(load_error), "\n",
          "Delete ", object_file, " and rerun in case the download was corrupted."
        )
      }
    )
    seurat_names <- object_names[vapply(
      object_names,
      function(name) inherits(imported[[name]], "Seurat"),
      logical(1)
    )]
    if (length(seurat_names) != 1L) {
      stop(
        "Workspace objects: ", paste(object_names, collapse = ", "), ". ",
        "Expected exactly one Seurat object, found ", length(seurat_names), "."
      )
    }
    message("Loaded workspace object: ", seurat_names)
    imported[[seurat_names]]
  }
)
stopifnot(inherits(hbv, "Seurat"))
writeLines(capture.output(hbv), "results/figure5/object_summary.txt")
write.csv(
  data.frame(
    column = colnames(hbv[[]]),
    class = vapply(hbv[[]], function(x) class(x)[1], character(1))
  ),
  "results/figure5/metadata_columns.csv", row.names = FALSE
)

# Use deposited clusters and UMAP. Recomputing them would introduce seed- and
# version-dependent differences.
cluster_candidates <- c(
  "seurat_clusters_new", "seurat_clusters", "cluster", "clusters", "Cluster"
)
cluster_col <- cluster_candidates[cluster_candidates %in% colnames(hbv[[]])][1]
if (is.na(cluster_col)) stop("No cluster field found; inspect metadata_columns.csv.")

reduction_candidates <- c("umap", "harmony_umap", "umapharmony")
reduction <- reduction_candidates[reduction_candidates %in% Reductions(hbv)][1]
if (is.na(reduction)) {
  stop("No UMAP found. Available reductions: ", paste(Reductions(hbv), collapse = ", "))
}

clusters <- as.character(hbv[[cluster_col, drop = TRUE]])
hbv$figure5_cluster <- factor(
  clusters,
  levels = sort(unique(clusters)),
  labels = paste0("C", sort(unique(clusters)))
)

# Identify the clinical-phase field. If this fails, inspect metadata_columns.csv
# and set disease_col manually.
disease_candidates <- c(
  "diagnosis", "disease", "disease_phase", "Disease", "condition", "Condition",
  "clinical_phase", "infection", "group", "Group", "status"
)
disease_col <- NA_character_
for (candidate in intersect(disease_candidates, colnames(hbv[[]]))) {
  values <- unique(tolower(as.character(hbv[[candidate, drop = TRUE]])))
  if (sum(grepl("acute|resolved|chronic|endogenous", values)) >= 2) {
    disease_col <- candidate
    break
  }
}
if (is.na(disease_col)) {
  stop("Disease field not identified. Inspect metadata_columns.csv and set disease_col.")
}

phase_raw <- as.character(hbv[[disease_col, drop = TRUE]])
diagnosis_map <- c(
  chronic_treat_naive_pol = "Endogenous control/chronic HBV",
  acute_pol = "Acute HBV",
  acute_res_pol = "Acute resolved HBV"
)
phase <- unname(diagnosis_map[phase_raw])
if (anyNA(phase)) {
  stop("Unmapped disease labels: ", paste(unique(phase_raw[is.na(phase)]), collapse = ", "))
}
hbv$figure5_phase <- factor(
  phase,
  levels = c(
    "Endogenous control/chronic HBV",
    "Acute HBV",
    "Acute resolved HBV"
  )
)

message("Cluster field: ", cluster_col)
message("Disease field: ", disease_col)
message("Reduction: ", reduction)
print(table(hbv$figure5_phase, hbv$figure5_cluster))

# Figure 5a: overall cluster map and disease-phase split -------------------
cluster_colors <- c(
  C0 = "#A48AD3", C1 = "#1CC5FE", C2 = "#FBA27D",
  C3 = "#D9CE97", C4 = "#FB7D80"
)
cluster_colors <- cluster_colors[seq_len(length(levels(hbv$figure5_cluster)))]

p_all <- DimPlot(
  hbv, reduction = reduction, group.by = "figure5_cluster",
  label = TRUE, repel = TRUE, cols = cluster_colors, raster = FALSE
) + ggtitle("Clustering") + theme_classic()

p_phase <- DimPlot(
  hbv, reduction = reduction, group.by = "figure5_cluster",
  split.by = "figure5_phase", cols = cluster_colors,
  ncol = 3, raster = FALSE
) & theme_classic() & theme(legend.position = "none")

p5a <-  p_all / p_phase + plot_annotation(tag_levels = "a")


message("HBV Seurat object loaded and metadata prepared.")
message("Object name: hbv")

print(
  table(
    Clinical_state = hbv$figure5_phase,
    Cluster = hbv$figure5_cluster
  )
)
