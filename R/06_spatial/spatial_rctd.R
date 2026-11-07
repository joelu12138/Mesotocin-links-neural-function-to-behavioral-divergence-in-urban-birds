# Spatial cell-type annotation by RCTD (spacexr v2.2.1), full mode.
# Methods: cell identities inferred with RCTD in full mode using snRNA-seq
# expression + annotations as reference; annotation at the major cell-class
# level and, within annotated neurons, at the level of 13 neuronal subtypes.
# Each cellbin is assigned the identity with the highest estimated RCTD weight.
#
# Input contract (confirmed):
#   - query: cellbin spatial Seurat objects, per region (x habitat)
#   - regions {HPR, LSR, NCL, PVN}; group A1 = Rural, A2 = Urban
#   - reference (major class): sce.all$celltype  (from single_cell_diff.R)
#   - reference (neuronal subtype): subobj$celltype = Ex1..Ex10 / In1..In3
#     (from neuronal_subtype_diff.R)
# EDIT query_file() below to match your actual cellbin file names.

suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
  library(dplyr)
})

# ---- snRNA-seq references ----
load("sce_all_integrated.Rdata")   # produces sce.all  ($celltype = major classes)
load("neuronal_subtypes.Rdata")    # produces subobj   ($celltype = Ex#/In#)

regions <- c("HPR", "LSR", "NCL", "PVN")
groups  <- c("A1", "A2")           # A1 = Rural, A2 = Urban

query_file <- function(region, group) {
  sprintf("%s_%s.sub.Rdata", region, group)   # <-- EDIT to your cellbin file naming
}

load_obj <- function(f) {
  e <- new.env()
  load(f, envir = e)
  get(ls(e)[1], envir = e)
}

# Run RCTD (full mode) and assign each cellbin its max-weight identity.
# ref_celltype_col names the reference annotation column (major class or subtype).
annotate_rctd <- function(query_obj, ref_obj, ref_celltype_col = "celltype") {
  DefaultAssay(ref_obj) <- "RNA"
  sc_counts <- GetAssayData(ref_obj, assay = "RNA", slot = "counts")
  celltype <- as.factor(ref_obj@meta.data[[ref_celltype_col]])
  names(celltype) <- colnames(ref_obj)
  reference <- Reference(sc_counts, celltype)

  spatial_count <- query_obj@assays$Spatial$counts
  coords <- query_obj@images[[1]]@coordinates
  coords <- coords[colnames(spatial_count), , drop = FALSE]
  query <- SpatialRNA(coords, spatial_count)

  rctd <- create.RCTD(query, reference, max_cores = 1)
  rctd <- run.RCTD(rctd, doublet_mode = "full")
  wmat <- as.matrix(normalize_weights(rctd@results$weights))

  max_ct <- colnames(wmat)[apply(wmat, 1, which.max)]
  annot <- data.frame(spot = rownames(wmat), main_celltype = max_ct,
                      stringsAsFactors = FALSE)
  query_obj@meta.data$main_celltype <-
    annot$main_celltype[match(rownames(query_obj@meta.data), annot$spot)]
  query_obj@meta.data$main_celltype[is.na(query_obj@meta.data$main_celltype)] <- "Unknown"
  query_obj <- AddMetaData(query_obj, metadata = as.data.frame(wmat))
  query_obj
}

write_csv_bom <- function(d, path) {
  con <- file(path, open = "w", encoding = "UTF-8")
  writeLines("\ufeff", con, useBytes = TRUE, sep = "")
  write.csv(d, con, row.names = FALSE)
  close(con)
}

for (rg in regions) {
  for (g in groups) {
    q <- load_obj(query_file(rg, g))

    # Level 1: major cell classes (reference = sce.all)
    q <- annotate_rctd(q, sce.all, "celltype")
    write_csv_bom(as.data.frame(table(q$main_celltype)),
                  sprintf("%s_%s_major_class_counts.csv", rg, g))
    save(q, file = sprintf("%s_%s.annotated.Rdata", rg, g))

    # Level 2: within annotated neurons, 13 subtypes (reference = subobj)
    neurons <- subset(q, main_celltype == "Neuron")
    if (ncol(neurons) > 0) {
      neurons <- annotate_rctd(neurons, subobj, "celltype")  # main_celltype -> Ex#/In#
      write_csv_bom(as.data.frame(table(neurons$main_celltype)),
                    sprintf("%s_%s_neuron_subtype_counts.csv", rg, g))
      save(neurons, file = sprintf("%s_%s.neuron_subtype.Rdata", rg, g))
    }
    cat(sprintf("DONE %s_%s  major classes=%d  neurons=%d\n",
                rg, g, length(unique(q$main_celltype)), ncol(neurons)))
  }
}
cat("\nWrote *_annotated.Rdata (major class) and *_neuron_subtype.Rdata (13 subtypes).\n")