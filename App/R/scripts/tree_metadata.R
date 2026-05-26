# ============================================================================
# tree_metadata.R — Assigning tree_id / tree_level / tree_parent_id to rows
# ============================================================================
#
# Tree-based correction algorithms need each contrast row in master_df to
# carry three columns:
#
#   tree_id          unique label per row (used to identify parent→child)
#   tree_level       integer level number (1 = root level)
#   tree_parent_id   tree_id of the row's parent at level (level - 1), or NA
#                    if this is a root-level row
#
# This helper takes a `tree_spec` list (one entry per level, in order) and
# the assembled master_df, and emits these columns. The `tree_spec` shape is:
#
#   tree_spec = list(
#     list(
#       name          = "L1 phase gate",
#       Test_Family   = "Genotype_Effect",
#       split_by      = "DMSO",                   # optional Split_By filter
#       pair_filter   = NULL,                     # optional Tested_Level regex
#       linkage_key   = NULL                      # L1 has no parent
#     ),
#     list(
#       name          = "L2 rescue",
#       Test_Family   = "Full_Interaction",
#       split_by      = NULL,
#       pair_filter   = NULL,                     # honor whatever Full_Interaction has
#       linkage_key   = list(parent_var = "Group",  child_var = "Group")
#     ),
#     list(
#       name          = "L3 specificity",
#       Test_Family   = "Drug_Effect",
#       split_by      = NULL,                     # both WT and KO contexts
#       linkage_key   = list(parent_var = "Drug",   child_var = "Tested_Level")
#     )
#   )
#
# linkage_key tells us which variable's identity threads from parent to child.
# For L1→L2, the obvious linkage is `Group` (phase): a phase gate is the
# parent of the rescue claims in THAT phase. For L2→L3, it's `Drug` (the drug
# being tested): a confirmed rescue at L2 is the parent of its specificity
# tests at L3.
#
# Rows in master_df that don't match any level's filter are returned with
# tree_level / tree_id / tree_parent_id = NA and a tree_status = "off_tree".
# Correction algorithms ignore those rows (they're computed but not gated).
# ============================================================================


# ---------------------------------------------------------------------------
# Apply a single level's filter to master_df. Returns row indices.
# ---------------------------------------------------------------------------
.match_level_rows <- function(master_df, level) {
  rows <- which(master_df$Test_Family == level$Test_Family)
  if (length(rows) == 0) return(integer(0))

  if (!is.null(level$split_by) && nzchar(level$split_by)) {
    # Split_By string may contain " + " for multi-variable splits. Match the
    # WHOLE Split_By column against the user's literal request.
    keep <- master_df$Split_By[rows] == level$split_by
    rows <- rows[keep]
  }

  if (!is.null(level$pair_filter) && nzchar(level$pair_filter)) {
    keep <- grepl(level$pair_filter,
                  master_df$Tested_Level[rows], perl = TRUE)
    rows <- rows[keep]
  }

  rows
}


# ---------------------------------------------------------------------------
# Public: assign tree metadata columns to master_df given a tree_spec.
# ---------------------------------------------------------------------------
assign_tree_metadata <- function(master_df, tree_spec) {
  n <- nrow(master_df)
  master_df$tree_id        <- NA_character_
  master_df$tree_level     <- NA_integer_
  master_df$tree_parent_id <- NA_character_

  if (is.null(tree_spec) || length(tree_spec) == 0) return(master_df)

  # First pass: assign tree_level + a unique tree_id per matched row.
  for (Lidx in seq_along(tree_spec)) {
    level <- tree_spec[[Lidx]]
    rows  <- .match_level_rows(master_df, level)
    if (length(rows) == 0) next

    # Build a deterministic id per row: "L{level}__{Tested_Level}__G{Group}__{Split_By}"
    new_ids <- sprintf("L%d__%s__G%s__%s",
                       Lidx,
                       master_df$Tested_Level[rows],
                       as.character(master_df$Group[rows]),
                       master_df$Split_By[rows])
    master_df$tree_id[rows]    <- new_ids
    master_df$tree_level[rows] <- Lidx
  }

  # Second pass: link each row at level L > 1 to its parent at level L-1
  # via the linkage_key.
  for (Lidx in seq_along(tree_spec)) {
    if (Lidx == 1L) next                       # L1 has no parent
    level   <- tree_spec[[Lidx]]
    linkage <- level$linkage_key
    if (is.null(linkage)) next

    child_var  <- linkage$child_var  %||% "Group"
    parent_var <- linkage$parent_var %||% child_var

    child_rows  <- which(master_df$tree_level == Lidx)
    parent_rows <- which(master_df$tree_level == (Lidx - 1L))
    if (length(child_rows) == 0 || length(parent_rows) == 0) next

    for (cr in child_rows) {
      # Extract the linkage value from the child row using child_var. The
      # child_var may be "Group" (use master_df$Group[cr]), or "Tested_Level"
      # (extract the drug name from a "Drug_X - DMSO" contrast string), or a
      # raw variable name (use master_df[[var]][cr] if such a column exists).
      child_value <- .extract_link_value(master_df, cr, child_var)

      # Walk parent rows, finding the one whose parent_var value matches.
      for (pr in parent_rows) {
        parent_value <- .extract_link_value(master_df, pr, parent_var)
        # Group is integer / numeric — coerce to character for comparison.
        if (identical(as.character(child_value), as.character(parent_value))) {
          master_df$tree_parent_id[cr] <- master_df$tree_id[pr]
          break
        }
      }
    }
  }

  master_df
}


# ---------------------------------------------------------------------------
# Pull a single linkage value out of one row of master_df, given a variable
# name. Handles the special cases:
#   "Group"           → master_df$Group[row]
#   "Tested_Level"    → first non-"DMSO"/reference token in the contrast
#                       (e.g. for "Drug_X - DMSO", returns "Drug_X")
#   other column      → master_df[[col]][row] if that column exists
# ---------------------------------------------------------------------------
.extract_link_value <- function(master_df, row, var_name) {
  if (var_name == "Group") return(as.character(master_df$Group[row]))
  if (var_name %in% colnames(master_df)) {
    return(as.character(master_df[[var_name]][row]))
  }
  # Fallback: parse "<focal> - <ref>" out of Tested_Level
  tl <- master_df$Tested_Level[row]
  parts <- strsplit(tl, " - ", fixed = TRUE)[[1]]
  if (length(parts) == 2L) {
    # The "treatment" side is usually first in trt.vs.ctrl naming. Return
    # both tokens stripped of spaces for downstream comparison.
    return(trimws(parts[1]))
  }
  return(trimws(tl))
}


# ---------------------------------------------------------------------------
# Build the DEFAULT rescue-paradigm tree spec from the run's variable names.
# Used when correction.json sidecar is missing — the rescue tree pre-populates
# (matches what the UI's "first open" will offer).
# ---------------------------------------------------------------------------
default_rescue_tree_spec <- function(var_names, ref_map) {
  if (length(var_names) < 2) return(NULL)
  # Convention: first variable is the gating variable (e.g. Genotype),
  # second variable is the screened variable (e.g. Drug).
  gate_var   <- var_names[1]
  screen_var <- var_names[2]
  ref_gate   <- ref_map[[gate_var]]    # e.g. "WT"
  ref_screen <- ref_map[[screen_var]]  # e.g. "DMSO"

  list(
    list(
      name        = paste0("L1 ", gate_var, " gate"),
      Test_Family = paste0(gate_var, "_Effect"),
      split_by    = ref_screen,         # gate in the reference-screen context
      pair_filter = NULL,
      linkage_key = NULL
    ),
    list(
      name        = "L2 rescue claim",
      Test_Family = "Full_Interaction",
      split_by    = NULL,
      pair_filter = NULL,
      linkage_key = list(parent_var = "Group", child_var = "Group")
    ),
    list(
      name        = paste0("L3 ", screen_var, " specificity"),
      Test_Family = paste0(screen_var, "_Effect"),
      split_by    = NULL,                # both genotype contexts
      pair_filter = NULL,
      linkage_key = list(parent_var = "Tested_Level", child_var = "Tested_Level")
    )
  )
}
