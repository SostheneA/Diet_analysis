# =============================================================================
# 0_set_prey_family.R — Bascule PREY_FAMILY ("_1" / "_2") dans tous les scripts
# =============================================================================
# Usage : mettre approach <- 1 ou 2, puis sourcer ce script depuis le dossier
# du projet. Il parcourt tous les .R du dossier (sous-dossiers inclus, donc
# R_helpers/ aussi) et remplace toute affectation
#     PREY_FAMILY <- "_x"   ou   PREY_FAMILY = "_x"
# par la valeur demandée, en préservant indentation et commentaires de fin de
# ligne. Seules les lignes d'affectation sont touchées; les usages en aval
# (paste0(..., PREY_FAMILY, ...)) restent intacts par construction.
# =============================================================================

approach <- 1    # <<< 1 ou 2

stopifnot(approach %in% c(1, 2))
target <- paste0('"_', approach, '"')

# Motif : debut de ligne (indentation permise), PREY_FAMILY, <- ou =, "_1" ou "_2"
pat <- '^(\\s*(?:if \\(!exists\\("PREY_FAMILY"\\)\\) )?PREY_FAMILY\\s*(<-|=)\\s*)"_(1|2|PP)"'

files <- list.files(".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("0_set_prey_family\\.R$", files)]   # ne pas se modifier soi-meme

n_changed <- 0L
for (f in files) {
  lines <- readLines(f, warn = FALSE)
  hits  <- grep(pat, lines)
  if (length(hits) == 0) next

  new_lines <- lines
  new_lines[hits] <- sub(pat, paste0("\\1", target), lines[hits])

  if (!identical(new_lines, lines)) {
    writeLines(new_lines, f)
    n_changed <- n_changed + 1L
    cat(sprintf("[modifie] %s  (%d ligne%s)\n", f, length(hits),
                ifelse(length(hits) > 1, "s", "")))
  } else {
    cat(sprintf("[deja ok] %s\n", f))
  }
}

if (n_changed == 0) {
  cat("Aucun fichier a modifier : tout est deja sur", target, "\n")
} else {
  cat(sprintf("\nTermine : %d fichier(s) bascule(s) sur PREY_FAMILY %s\n",
              n_changed, target))
}

# Verification finale — liste toutes les affectations restantes
cat("\n--- Etat actuel des affectations PREY_FAMILY ---\n")
for (f in files) {
  lines <- readLines(f, warn = FALSE)
  hits  <- grep("PREY_FAMILY\\s*(<-|=)\\s*\"", lines)
  for (h in hits) cat(sprintf("%s : l.%d : %s\n", f, h, trimws(lines[h])))
}
