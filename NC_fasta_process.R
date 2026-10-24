# ==============================================================================
# HIV-1 NUCLEOCAPSID (NC) MUTATION ANALYSIS
# ==============================================================================
#
# Description:
#   Processes nucleotide sequences in FASTA format, translates codons including
#   ambiguous IUPAC nucleotide codes, and identifies amino-acid substitutions
#   across positions 1-55 of the HIV-1 nucleocapsid (NC) region.
#
#   Amino-acid substitutions are classified into:
#     1. NC resistance-associated mutations (np.NC-RAMs)
#     2. Subtype C-associated polymorphisms (pNC)
#     3. Other amino-acid substitutions
#
# Input:
#   A nucleotide FASTA file containing sequences aligned to the NC coding region.
#
# Output:
#   A CSV file containing:
#     - Sequence_ID
#     - np.NC-RAMs
#     - pNC
#     - other_mutations
#
# Requirements:
#   R packages:
#     - Biostrings (Bioconductor)
#     - tidyverse
#
# Notes:
#   If sequences begin exactly at NC position 1, use offset_to_pos1 = 0.
#   If sequences contain upstream Gag sequence, adjust offset_to_pos1 accordingly.
#
# ==============================================================================


# 1. CHECK REQUIRED PACKAGES ===================================================

required_packages <- c("Biostrings", "tidyverse")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following required package(s) are not installed: ",
    paste(missing_packages, collapse = ", "),
    ".\nPlease install them before running this script."
  )
}

library(Biostrings)
library(tidyverse)


# 2. DEFINE REFERENCE DATA =====================================================

Pos <- 1:55

HXB2 <- c(
  "M", "Q", "R", "G", "N", "F", "R", "N", "Q", "R",
  "K", "I", "V", "K", "C", "F", "N", "C", "G", "K",
  "E", "G", "H", "T", "A", "R", "N", "C", "R", "A",
  "P", "R", "K", "K", "G", "C", "W", "K", "C", "G",
  "K", "E", "G", "H", "Q", "M", "K", "D", "C", "T",
  "E", "R", "Q", "A", "N"
)

Subtype_C <- c(
  "M; I", "Q", "K; R", "S; G; N", "N", "F", "K; R", "G; S",
  "P; A; S", "K; R", "R; K", "I; N; T", "V; I; P", "K; R",
  "C", "F; Y", "N", "C", "G", "K; R", "E", "G", "H",
  "I; L; V", "A", "R; K", "N", "C", "R; K", "A", "P", "R",
  "K", "K; R", "G", "C", "W", "K", "C", "G", "K; R; Q",
  "E", "G", "H", "Q", "M", "K", "D; E", "C", "T; N; I",
  "E; G", "R", "Q", "A", "N"
)

RAMs <- c(
  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
  NA, "R", NA, "V", NA, NA, "S", NA, NA, NA,
  NA, NA, NA, NA, NA, NA, NA, NA, NA, "R",
  "G", NA, NA, NA, NA, "I; L", NA, "E", NA, "S",
  "G", NA, NA, NA, NA
)

reference_data <- tibble(
  Pos = Pos,
  HXB2 = HXB2,
  Subtype_C = Subtype_C,
  RAMs = RAMs
)


# 3. TRANSLATE AMBIGUOUS CODONS ================================================

# Cache translations so that each unique codon only needs to be processed once.
codon_cache <- new.env(hash = TRUE)

get_all_aa_translations <- function(codon_str) {

  codon_str <- toupper(codon_str)

  # Return cached result if this codon has already been processed.
  if (exists(codon_str, envir = codon_cache, inherits = FALSE)) {
    return(get(codon_str, envir = codon_cache, inherits = FALSE))
  }

  # Alignment gaps.
  if (grepl("-", codon_str, fixed = TRUE)) {
    assign(codon_str, "-", envir = codon_cache)
    return("-")
  }

  # Completely unknown or unsupported bases.
  if (grepl("N|\\?", codon_str, ignore.case = TRUE)) {
    assign(codon_str, "X", envir = codon_cache)
    return("X")
  }

  # Expand IUPAC ambiguous nucleotide codes.
  bases <- strsplit(codon_str, "")[[1]]

  expanded_bases <- lapply(bases, function(b) {
    if (b %in% names(IUPAC_CODE_MAP)) {
      strsplit(IUPAC_CODE_MAP[[b]], "")[[1]]
    } else {
      b
    }
  })

  grid <- expand.grid(expanded_bases, stringsAsFactors = FALSE)
  codons <- apply(grid, 1, paste, collapse = "")

  # Translate all possible codons.
  aas <- unique(
    vapply(
      codons,
      function(codon) {
        suppressWarnings(
          as.character(
            translate(DNAString(codon), if.fuzzy.codon = "solve")
          )
        )
      },
      character(1)
    )
  )

  # Save result in cache.
  assign(codon_str, aas, envir = codon_cache)

  aas
}


# 4. MAIN ANALYSIS FUNCTION ====================================================

process_nc_mutations <- function(
  fasta_file,
  ref_data = reference_data,
  offset_to_pos1 = 0
) {

  # Basic input checks.
  if (!file.exists(fasta_file)) {
    stop("FASTA file not found: ", fasta_file)
  }

  if (!is.numeric(offset_to_pos1) ||
      length(offset_to_pos1) != 1 ||
      is.na(offset_to_pos1) ||
      offset_to_pos1 < 0 ||
      offset_to_pos1 %% 1 != 0) {
    stop("offset_to_pos1 must be a single non-negative integer.")
  }

  # Read FASTA sequences.
  seqs <- readDNAStringSet(fasta_file, format = "fasta")

  total_seqs <- length(seqs)

  if (total_seqs == 0) {
    stop("No sequences were found in the FASTA file.")
  }

  # Provide identifiers if FASTA headers are absent.
  seq_ids <- names(seqs)

  if (is.null(seq_ids)) {
    seq_ids <- paste0("sequence_", seq_along(seqs))
  } else {
    missing_ids <- is.na(seq_ids) | seq_ids == ""
    seq_ids[missing_ids] <- paste0("sequence_", which(missing_ids))
  }

  message(
    "Successfully loaded ", total_seqs,
    " sequence(s). Starting processing..."
  )

  results_list <- vector("list", total_seqs)

  for (i in seq_along(seqs)) {

    # Progress indicator.
    if (i == 1 || i %% 10 == 0 || i == total_seqs) {
      message(
        "Processing sequence ", i,
        " of ", total_seqs, "..."
      )
    }

    seq_id <- seq_ids[i]
    seq_dna <- toupper(as.character(seqs[[i]]))

    np_rams_list <- character(0)
    pnc_list <- character(0)
    other_list <- character(0)

    for (j in seq_len(nrow(ref_data))) {

      pos <- ref_data$Pos[j]
      hxb2_aa <- ref_data$HXB2[j]

      start_idx <- offset_to_pos1 + (pos - 1) * 3 + 1
      end_idx <- start_idx + 2

      # Stop if the sequence does not extend to this position.
      if (end_idx > nchar(seq_dna)) {
        break
      }

      codon <- substr(seq_dna, start_idx, end_idx)

      translated_aas <- get_all_aa_translations(codon)

      for (aa in translated_aas) {

        # Ignore reference amino acids, gaps, unknown residues and stop codons.
        if (aa == hxb2_aa ||
            aa == "-" ||
            aa == "X" ||
            aa == "*") {
          next
        }

        mut_name <- paste0(hxb2_aa, pos, aa)
        is_classified <- FALSE

        # ----------------------------------------------------------------------
        # NC resistance-associated mutations
        # ----------------------------------------------------------------------

        ram_ref <- ref_data$RAMs[j]

        if (!is.na(ram_ref)) {

          ram_target <- trimws(
            unlist(strsplit(ram_ref, ";", fixed = TRUE))
          )

          if (aa %in% ram_target || mut_name %in% ram_target) {
            np_rams_list <- c(np_rams_list, mut_name)
            is_classified <- TRUE
          }
        }

        if (is_classified) {
          next
        }

        # ----------------------------------------------------------------------
        # Subtype C-associated polymorphisms
        # ----------------------------------------------------------------------

        subtype_c_ref <- ref_data$Subtype_C[j]

        if (!is.na(subtype_c_ref)) {

          sc_target <- trimws(
            unlist(strsplit(subtype_c_ref, ";", fixed = TRUE))
          )

          if (aa %in% sc_target || mut_name %in% sc_target) {
            pnc_list <- c(pnc_list, mut_name)
            is_classified <- TRUE
          }
        }

        if (is_classified) {
          next
        }

        # ----------------------------------------------------------------------
        # Other substitutions
        # ----------------------------------------------------------------------

        other_list <- c(other_list, mut_name)
      }
    }

    results_list[[i]] <- tibble(
      Sequence_ID = seq_id,
      `np.NC-RAMs` = paste(unique(np_rams_list), collapse = ", "),
      pNC = paste(unique(pnc_list), collapse = ", "),
      other_mutations = paste(unique(other_list), collapse = ", ")
    )
  }

  message("Processing complete. Compiling final table...")

  bind_rows(results_list)
}


# 5. RUN THE ANALYSIS ==========================================================

# Select the FASTA file interactively.
my_fasta_file <- file.choose()

# IMPORTANT:
#
# If the FASTA sequences start exactly at NC position 1:
#     offset_to_pos1 = 0
#
# If sequences contain upstream sequence before NC position 1,
# offset_to_pos1 must be changed to the corresponding nucleotide offset.
#
# For example, if NC begins after 377 amino acids of upstream Gag sequence,
# the offset is:
#
#     377 * 3 = 1131 nucleotides

final_results <- process_nc_mutations(
  fasta_file = my_fasta_file,
  ref_data = reference_data,
  offset_to_pos1 = 0
)


# 6. VIEW AND SAVE RESULTS =====================================================

# Display interactively when running in RStudio.
if (interactive()) {
  View(final_results)
}

# Save results to the current working directory.
write_csv(
  final_results,
  "NC_mutations_output.csv"
)

message(
  "Results saved to: ",
  normalizePath("NC_mutations_output.csv", mustWork = FALSE)
)
