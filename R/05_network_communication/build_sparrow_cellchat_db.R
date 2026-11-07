###############################################################################
#  Build CellChatDB.treeSparrow v2
#  ================================
#  Translates CellChatDB.human → tree sparrow (Passer montanus) via
#  OrthoFinder ortholog mapping.
#
#  Key design decisions:
#    - 1-to-1 orthologs  → direct gene symbol replacement
#    - 1-to-many (human→sparrow paralogs) → expand into multiple interactions
#    - many-to-1 (multiple human → 1 sparrow) → deduplicate within row
#    - Unmapped genes → drop the interaction/complex
#    - Complex names are preserved; subunit genes are translated
#
#  Output:
#    - CellChatDB.treeSparrow.rda     (ready to use)
#    - CellChatDB.treeSparrow_conservative.rda  (1-to-1 only, fewer interactions)
#    - build_report.txt
#
#  Usage:
#    cellchat@DB <- CellChatDB.treeSparrow
#
#  Author : [Your Name]
#  Date   : 2026-08-07
###############################################################################

rm(list = ls())
library(CellChat)

# ============================================================================
# 0. Paths & settings
# ============================================================================

ORTHOLOG_ALL      <- "../ortholog_pipeline/sparrow_human_all_orthologs.tsv"
ORTHOLOG_ONE2ONE  <- "../ortholog_pipeline/sparrow_human_one2one.tsv"
OUTPUT_FULL       <- "CellChatDB.treeSparrow.rda"
OUTPUT_CONSERV    <- "CellChatDB.treeSparrow_conservative.rda"
REPORT_FILE       <- "build_report.txt"

sink(REPORT_FILE, split = TRUE)
cat("========================================\n")
cat("CellChatDB.treeSparrow — Build Report\n")
cat("========================================\n\n")

# ============================================================================
# 1. Load ortholog data
# ============================================================================

cat("--- 1. Loading ortholog data ---\n")

ortho_all <- read.table(ORTHOLOG_ALL, sep = "\t", header = TRUE,
                        comment.char = "", quote = "", stringsAsFactors = FALSE)
ortho_1to1 <- read.table(ORTHOLOG_ONE2ONE, sep = "\t", header = TRUE,
                         comment.char = "", quote = "", stringsAsFactors = FALSE)

cat(sprintf("  All ortholog pairs     : %d\n", nrow(ortho_all)))
cat(sprintf("  Strict 1-to-1 pairs    : %d\n", nrow(ortho_1to1)))
cat(sprintf("  Unique sparrow genes   : %d\n", length(unique(ortho_all$Sparrow_Gene))))
cat(sprintf("  Unique human genes     : %d\n", length(unique(ortho_all$Human_Symbol))))

tbl <- table(ortho_all$Relation)
cat("  Relation types:\n")
for (r in names(tbl)) cat(sprintf("    %-16s: %d\n", r, tbl[r]))

# ============================================================================
# 2. Build mapping tables
# ============================================================================

cat("\n--- 2. Building Human ↔ Sparrow maps ---\n")

# Full mapping (all ortholog types)
human2sparrow_all <- split(ortho_all$Sparrow_Gene, ortho_all$Human_Symbol)

# Conservative mapping (1-to-1 only)
human2sparrow_1to1 <- split(ortho_1to1$Sparrow_Gene, ortho_1to1$Human_Symbol)

# Reverse: sparrow → human
sparrow2human <- split(ortho_all$Human_Symbol, ortho_all$Sparrow_Gene)

cat(sprintf("  Full map  : %d human genes\n", length(human2sparrow_all)))
cat(sprintf("  1-to-1 map: %d human genes\n", length(human2sparrow_1to1)))

# ============================================================================
# 3. Human genes required by CellChatDB
# ============================================================================

cat("\n--- 3. Human genes in CellChatDB ---\n")

DB_HUMAN <- CellChatDB.human

collect_genes <- function(df, col) {
  genes <- unlist(strsplit(as.character(df[[col]]), ", "))
  genes <- trimws(genes)
  unique(genes[genes != ""])
}

int_lig_genes   <- collect_genes(DB_HUMAN$interaction, "ligand.symbol")
int_rec_genes   <- collect_genes(DB_HUMAN$interaction, "receptor.symbol")
complex_genes   <- unique(as.character(unlist(DB_HUMAN$complex[, 1:5])))
complex_genes   <- complex_genes[complex_genes != ""]
cof_genes       <- unique(as.character(unlist(DB_HUMAN$cofactor[, 1:16])))
cof_genes       <- cof_genes[cof_genes != ""]

all_needed <- unique(c(int_lig_genes, int_rec_genes, complex_genes, cof_genes))
cat(sprintf("  Total human genes needed: %d\n", length(all_needed)))

# Coverage
mapped_full  <- intersect(all_needed, names(human2sparrow_all))
mapped_1to1  <- intersect(all_needed, names(human2sparrow_1to1))
unmapped     <- setdiff(all_needed, names(human2sparrow_all))

cat(sprintf("  Full map coverage : %d / %d (%.1f%%)\n",
            length(mapped_full), length(all_needed),
            100 * length(mapped_full) / length(all_needed)))
cat(sprintf("  1-to-1 coverage   : %d / %d (%.1f%%)\n",
            length(mapped_1to1), length(all_needed),
            100 * length(mapped_1to1) / length(all_needed)))
cat(sprintf("  Unmapped genes    : %d\n", length(unmapped)))
cat("  First 40 unmapped:", paste(head(unmapped, 40), collapse = ", "), "\n")

# ============================================================================
# 4. Core translation functions
# ============================================================================

#' Translate one human gene → sparrow gene(s)
#' @param gene  Human gene symbol
#' @param map   Named list: human → sparrow gene vector
#' @return Character vector of sparrow genes (empty if no ortholog)
translate_one <- function(gene, map) {
  gene <- trimws(gene)
  if (gene == "" || !gene %in% names(map)) return(character(0))
  map[[gene]]
}

#' Translate comma-separated human gene list → all sparrow paralog combinations
#'
#' "TGFBR2, TGFBR1" where TGFBR2→[TgfbR2] and TGFBR1→[TgfbR1a, Tgfbr1b]
#' → c("TgfbR2, Tgfbr1a", "TgfbR2, Tgfbr1b")
#'
#' @param csv  Comma-separated human gene symbols
#' @param map  Human→sparrow mapping list
#' @return Character vector of sparrow gene combinations
translate_gene_set <- function(csv, map) {
  if (is.na(csv) || csv == "") return(character(0))
  human_genes <- trimws(unlist(strsplit(csv, ",")))
  human_genes <- human_genes[human_genes != ""]
  if (length(human_genes) == 0) return(character(0))

  sp_lists <- lapply(human_genes, translate_one, map = map)
  if (any(lengths(sp_lists) == 0)) return(character(0))

  combs <- expand.grid(sp_lists, stringsAsFactors = FALSE)
  apply(combs, 1, paste, collapse = ", ")
}

#' Combine ligand.symbol and receptor.symbol into an interaction_name
make_interaction_name <- function(lig_sym, rec_sym) {
  lig_part <- gsub(", ", "_", lig_sym)
  rec_part <- gsub(", ", "_", rec_sym)
  paste0(lig_part, "_", rec_part)
}

# ============================================================================
# 5. Translate complex table
# ============================================================================

translate_complex <- function(complex_human, map) {
  cat("\n--- Translating complex table ---\n")
  complex_sp <- complex_human
  valid <- rep(TRUE, nrow(complex_human))

  for (i in seq_len(nrow(complex_human))) {
    has_any <- FALSE
    for (j in 1:5) {
      hg <- complex_human[i, j]
      if (hg == "") next
      sg <- translate_one(hg, map)
      if (length(sg) == 0) {
        complex_sp[i, j] <- ""
      } else {
        complex_sp[i, j] <- sg[1]  # first paralog for complex definition
        has_any <- TRUE
      }
    }
    if (!has_any) valid[i] <- FALSE
  }

  complex_sp <- complex_sp[valid, , drop = FALSE]
  cat(sprintf("  %d → %d complexes (dropped %d)\n",
              nrow(complex_human), nrow(complex_sp), sum(!valid)))
  if (sum(!valid) > 0) {
    cat("  Dropped:", paste(rownames(complex_human)[!valid], collapse = ", "), "\n")
  }
  return(complex_sp)
}

# ============================================================================
# 6. Translate cofactor table
# ============================================================================

translate_cofactor <- function(cofactor_human, map) {
  cat("\n--- Translating cofactor table ---\n")
  cofactor_sp <- cofactor_human
  for (i in seq_len(nrow(cofactor_human))) {
    for (j in 1:16) {
      hg <- cofactor_human[i, j]
      if (hg == "") next
      sg <- translate_one(hg, map)
      cofactor_sp[i, j] <- if (length(sg) > 0) sg[1] else ""
    }
  }
  cat(sprintf("  %d cofactor rows\n", nrow(cofactor_sp)))
  return(cofactor_sp)
}

# ============================================================================
# 7. Translate interaction table
# ============================================================================

translate_interactions <- function(int_human, complex_human, complex_sparrow, map, version_tag) {
  cat("\n--- Translating interaction table ---\n")
  cat(sprintf("  Input: %d interactions\n", nrow(int_human)))

  # Build lookup: human complex name → sparrow complex name (same names, but only valid ones)
  valid_complex_names <- rownames(complex_sparrow)

  rows_list <- list()
  n_dropped_unmapped <- 0
  n_dropped_complex  <- 0
  n_total_input      <- nrow(int_human)

  for (i in seq_len(nrow(int_human))) {
    row_h <- int_human[i, ]

    # Translate ligand.symbol and receptor.symbol (these are the actual genes)
    lig_sym_sets <- translate_gene_set(row_h$ligand.symbol, map)
    rec_sym_sets <- translate_gene_set(row_h$receptor.symbol, map)

    if (length(lig_sym_sets) == 0 || length(rec_sym_sets) == 0) {
      n_dropped_unmapped <- n_dropped_unmapped + 1
      next
    }

    # Check if ligand/receptor references a complex that was dropped
    lig_ref <- row_h$ligand
    rec_ref <- row_h$receptor
    if (lig_ref %in% rownames(complex_human) && !lig_ref %in% valid_complex_names) {
      n_dropped_complex <- n_dropped_complex + 1
      next
    }
    if (rec_ref %in% rownames(complex_human) && !rec_ref %in% valid_complex_names) {
      n_dropped_complex <- n_dropped_complex + 1
      next
    }

    # For ligand/receptor fields:
    # - If a complex name → keep it
    # - If a single gene → translate it (for each paralog combination, use the
    #   matching ligand.symbol gene)
    # - If a multi-gene non-complex (rare) → translate each gene
    resolve_identifier <- function(ref, sym_set, complex_human_names) {
      if (ref %in% complex_human_names) {
        return(ref)  # complex name preserved
      }
      # Single gene reference — extract from sym_set
      sym_genes <- trimws(unlist(strsplit(sym_set, ", ")))
      if (length(sym_genes) == 1) {
        return(sym_genes)  # single sparrow gene
      }
      # multi-gene, not a known complex → create a name from the genes
      return(gsub(", ", "_", sym_set))
    }

    # Generate all paralog combinations
    for (ls in lig_sym_sets) {
      for (rs in rec_sym_sets) {
        row_new <- row_h

        row_new$ligand.symbol   <- ls
        row_new$receptor.symbol <- rs

        row_new$ligand   <- resolve_identifier(lig_ref, ls, rownames(complex_human))
        row_new$receptor <- resolve_identifier(rec_ref, rs, rownames(complex_human))

        row_new$interaction_name   <- make_interaction_name(ls, rs)
        row_new$interaction_name_2 <- paste0(gsub(", ", "_", ls),
                                             " - (", gsub(", ", "_", rs), ")")
        row_new$version <- version_tag

        rows_list[[length(rows_list) + 1]] <- row_new
      }
    }
  }

  int_sp <- do.call(rbind, rows_list)
  rownames(int_sp) <- NULL

  cat(sprintf("  Output: %d interactions\n", nrow(int_sp)))
  cat(sprintf("  Dropped (no ortholog)       : %d\n", n_dropped_unmapped))
  cat(sprintf("  Dropped (complex lost)      : %d\n", n_dropped_complex))
  cat(sprintf("  Kept / input                : %.1f%%\n",
              100 * (n_total_input - n_dropped_unmapped - n_dropped_complex) / n_total_input))
  return(int_sp)
}

# ============================================================================
# 8. Build geneInfo
# ============================================================================

build_geneinfo <- function(int_sp, complex_sp, cofactor_sp, sparrow2human, human_geneInfo) {
  cat("\n--- Building geneInfo ---\n")

  sp_genes <- unique(c(
    collect_genes(int_sp, "ligand.symbol"),
    collect_genes(int_sp, "receptor.symbol"),
    as.character(unlist(complex_sp[, 1:5])),
    as.character(unlist(cofactor_sp[, 1:16]))
  ))
  sp_genes <- sp_genes[sp_genes != ""]

  gi <- data.frame(
    EntryID.uniprot = NA_character_,
    Symbol          = sp_genes,
    Synonym         = NA_character_,
    Reviewed        = "unreviewed",
    EntryName       = NA_character_,
    ProteinName     = NA_character_,
    Keywords        = NA_character_,
    Location        = NA_character_,
    AntibodyName    = NA_character_,
    Human_Ortholog  = sapply(sp_genes, function(g) {
      if (g %in% names(sparrow2human)) paste(sparrow2human[[g]], collapse = "; ") else NA_character_
    }),
    stringsAsFactors = FALSE
  )

  # Inherit annotation for conserved gene symbols
  n_annot <- 0
  for (i in seq_len(nrow(gi))) {
    sg <- gi$Symbol[i]
    if (sg %in% human_geneInfo$Symbol) {
      hr <- human_geneInfo[human_geneInfo$Symbol == sg, ][1, ]
      gi$EntryID.uniprot[i] <- hr$EntryID.uniprot
      gi$Synonym[i]         <- hr$Synonym
      gi$EntryName[i]       <- hr$EntryName
      gi$ProteinName[i]     <- hr$ProteinName
      gi$Keywords[i]        <- hr$Keywords
      gi$Location[i]        <- hr$Location
      n_annot <- n_annot + 1
    }
  }

  cat(sprintf("  %d sparrow genes in DB (%d conserved, annotated from human)\n",
              nrow(gi), n_annot))
  return(gi)
}

# ============================================================================
# 9. Assemble & validate
# ============================================================================

build_db <- function(int_sp, complex_sp, cofactor_sp, geneinfo_sp, tag) {
  cat(sprintf("\n--- Assembling %s ---\n", tag))

  db <- list(
    interaction = int_sp,
    complex     = complex_sp,
    cofactor    = cofactor_sp,
    geneInfo    = geneinfo_sp
  )

  # Structural validation
  stopifnot(identical(sort(names(db)), sort(names(DB_HUMAN))))
  stopifnot(identical(colnames(db$interaction), colnames(DB_HUMAN$interaction)))
  stopifnot(identical(colnames(db$complex),     colnames(DB_HUMAN$complex)))
  stopifnot(identical(colnames(db$cofactor),    colnames(DB_HUMAN$cofactor)))
  stopifnot("Symbol" %in% colnames(db$geneInfo))

  cat("  Structure validation: PASS\n")
  cat(sprintf("  Interactions : %d\n", nrow(db$interaction)))
  cat(sprintf("  Complexes    : %d\n", nrow(db$complex)))
  cat(sprintf("  Cofactors    : %d\n", nrow(db$cofactor)))
  cat(sprintf("  geneInfo     : %d\n", nrow(db$geneInfo)))

  # Annotation breakdown
  cat("  Annotation breakdown:\n")
  for (ann in sort(unique(db$interaction$annotation))) {
    cat(sprintf("    %-22s: %d\n", ann, sum(db$interaction$annotation == ann)))
  }

  return(db)
}

# ============================================================================
# 10. Execute — Full version (all ortholog types)
# ============================================================================

cat("\n========================================\n")
cat("  BUILD 1: Full ortholog database\n")
cat("========================================\n")

complex_sp_full  <- translate_complex(DB_HUMAN$complex, human2sparrow_all)
cofactor_sp_full <- translate_cofactor(DB_HUMAN$cofactor, human2sparrow_all)
int_sp_full      <- translate_interactions(DB_HUMAN$interaction, DB_HUMAN$complex,
                                           complex_sp_full, human2sparrow_all,
                                           "CellChatDB.treeSparrow")
gi_full <- build_geneinfo(int_sp_full, complex_sp_full, cofactor_sp_full,
                          sparrow2human, DB_HUMAN$geneInfo)
CellChatDB.treeSparrow <- build_db(int_sp_full, complex_sp_full,
                                   cofactor_sp_full, gi_full,
                                   "CellChatDB.treeSparrow")

# ============================================================================
# 11. Execute — Conservative version (1-to-1 orthologs only)
# ============================================================================

cat("\n========================================\n")
cat("  BUILD 2: Conservative (1-to-1 only)\n")
cat("========================================\n")

complex_sp_cons  <- translate_complex(DB_HUMAN$complex, human2sparrow_1to1)
cofactor_sp_cons <- translate_cofactor(DB_HUMAN$cofactor, human2sparrow_1to1)
int_sp_cons      <- translate_interactions(DB_HUMAN$interaction, DB_HUMAN$complex,
                                           complex_sp_cons, human2sparrow_1to1,
                                           "CellChatDB.treeSparrow_conservative")
gi_cons <- build_geneinfo(int_sp_cons, complex_sp_cons, cofactor_sp_cons,
                          sparrow2human, DB_HUMAN$geneInfo)
CellChatDB.treeSparrow_conservative <- build_db(int_sp_cons, complex_sp_cons,
                                                cofactor_sp_cons, gi_cons,
                                                "CellChatDB.treeSparrow_conservative")

# ============================================================================
# 12. Comparison summary
# ============================================================================

cat("\n========================================\n")
cat("  COMPARISON: Human → Full → Conservative\n")
cat("========================================\n")

comp <- data.frame(
  Component   = c("Interactions", "Complexes", "Cofactors", "geneInfo genes"),
  Human       = c(nrow(DB_HUMAN$interaction), nrow(DB_HUMAN$complex),
                  nrow(DB_HUMAN$cofactor), nrow(DB_HUMAN$geneInfo)),
  Full        = c(nrow(int_sp_full), nrow(complex_sp_full),
                  nrow(cofactor_sp_full), nrow(gi_full)),
  Conservative = c(nrow(int_sp_cons), nrow(complex_sp_cons),
                   nrow(cofactor_sp_cons), nrow(gi_cons)),
  stringsAsFactors = FALSE
)
print(comp, row.names = FALSE)

# ============================================================================
# 13. Save
# ============================================================================

cat("\n--- Saving ---\n")
save(CellChatDB.treeSparrow, file = OUTPUT_FULL)
cat(sprintf("  %s (%.1f KB)\n", OUTPUT_FULL, file.info(OUTPUT_FULL)$size / 1024))

save(CellChatDB.treeSparrow_conservative, file = OUTPUT_CONSERV)
cat(sprintf("  %s (%.1f KB)\n", OUTPUT_CONSERV, file.info(OUTPUT_CONSERV)$size / 1024))

# ============================================================================
# 14. Functional test
# ============================================================================

cat("\n--- Functional test ---\n")

test_genes <- gi_full$Symbol[1:min(20, nrow(gi_full))]
test_data <- matrix(rpois(length(test_genes) * 10, 2), nrow = length(test_genes),
                    dimnames = list(test_genes, paste0("spot", 1:10)))
test_meta <- data.frame(labels = rep("Neuron", 10), samples = "test",
                        row.names = paste0("spot", 1:10))

# Test with full DB
cc <- createCellChat(object = test_data, meta = test_meta, group.by = "labels")
cc@DB <- CellChatDB.treeSparrow
cc <- subsetData(cc)
cat(sprintf("  Full DB:          subsetData OK (%d/%d genes)\n",
            nrow(cc@data), length(test_genes)))

# Test with conservative DB
cc2 <- createCellChat(object = test_data, meta = test_meta, group.by = "labels")
cc2@DB <- CellChatDB.treeSparrow_conservative
cc2 <- subsetData(cc2)
cat(sprintf("  Conservative DB:  subsetData OK (%d/%d genes)\n",
            nrow(cc2@data), length(test_genes)))

cat("\n========================================\n")
cat("Build complete. Both databases ready.\n")
cat("========================================\n")

sink()

cat("\nDone. See build_report.txt for details.\n")
