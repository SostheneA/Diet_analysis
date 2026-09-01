#------------------------------------------------------------------------------#
# PreyCategory.R - BOTTOM-UP PREY AGGREGATION AT A GIVEN STOMACH THRESHOLD
#
# Used by 3_taxonomic_groups_Final.R (pooled families) and
# 3PP_taxonomic_groups.R (per-predator family).
#
# rank_counts()  : stomachs and predators behind each name, per rank, on the
#                  target rows. Computed once; it does not depend on x.
# PreyCategory() : for one threshold x, walks each target row up from species
#                  to phylum until a name passes the thresholds; rows reaching
#                  no qualifying rank take their coarsest available name.
#                  Returns two plain vectors (level, name) aligned on dt, NA
#                  outside the target rows. Nothing is modified by reference.
#------------------------------------------------------------------------------#

rank_counts <- function(dt, ranks, target_rows = NULL) {
  idx <- if (is.null(target_rows)) seq_len(nrow(dt)) else which(target_rows)
  out <- list()
  for (rank in ranks) {
    v <- dt[[rank]]
    rows <- idx[!is.na(v[idx]) & v[idx] != ""]
    if (!length(rows)) next
    cn <- dt[rows, .(n_stomach  = uniqueN(unlist(list_stomach_id)),
                     n_predator = uniqueN(unlist(list_predator_id))),
             by = c(rank)]
    setnames(cn, rank, "name")
    out[[rank]] <- cn
  }
  out
}

PreyCategory <- function(dt,
                         target_ranks = c("phylum", "subphylum", "class", "subclass", "order",
                                          "infraorder", "family", "genus", "species"),
                         full_taxonomy,
                         target_rows = NULL,
                         n_stomach_threshold = 100,
                         n_predator_threshold = 2,
                         counts = NULL) {

  n <- nrow(dt)
  level <- rep(NA_character_, n)
  name  <- rep(NA_character_, n)

  all_idx <- if (is.null(target_rows)) seq_len(n) else which(target_rows)
  rem <- all_idx
  if (!length(rem)) return(list(level = level, name = name))

  valid_cols  <- intersect(full_taxonomy[full_taxonomy %in% target_ranks], names(dt))
  spec_to_gen <- rev(valid_cols)
  if (is.null(counts)) counts <- rank_counts(dt, spec_to_gen, target_rows)

  for (rank in spec_to_gen) {
    if (!length(rem)) break
    cn <- counts[[rank]]
    if (is.null(cn)) next
    valid <- cn[n_stomach >= n_stomach_threshold & n_predator >= n_predator_threshold, name]
    if (!length(valid)) next
    v   <- as.character(dt[[rank]])
    hit <- rem[!is.na(v[rem]) & v[rem] != "" & v[rem] %chin% valid]
    if (length(hit)) {
      level[hit] <- rank
      name[hit]  <- v[hit]
      rem <- setdiff(rem, hit)
    }
  }

  if (length(rem)) {
    m <- lapply(valid_cols, function(cc) { x <- as.character(dt[[cc]][rem]); x[x == ""] <- NA_character_; x })
    level[rem] <- valid_cols[1]
    name[rem]  <- do.call(fcoalesce, m)
  }

  list(level = level, name = name)
}
