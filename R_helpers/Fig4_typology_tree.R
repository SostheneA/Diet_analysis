# =============================================================================
# Fig4_typology_tree.R - FIGURE 4: DECISION TREE OF THE DIAGNOSTIC TYPOLOGY
# -----------------------------------------------------------------------------
# Draws the tree as an SVG built from plain strings (no plotting package), then
# rasterises it to PNG with the rsvg package when it is installed.
#
#   source("R_helpers/Fig4_typology_tree.R")
#   draw_typology_tree()                                   # defaults
#   draw_typology_tree(out_prefix = "Output_Figures/Fig4_typology_tree", width = 3300)
#
# Text, chips and colours are in typology_tree_config(); pass an edited copy:
#   cfg <- typology_tree_config()
#   cfg$nodes$R2$sub <- "turnover + breadth change"
#   draw_typology_tree(cfg)
# Box positions are in the GEOMETRY block. The first version of this figure was
# produced by Typology-tree2.py (kept beside this file); the drawing is the same.
# =============================================================================

# ---------------------------------------------------------------- config
typology_tree_config <- function() {
  list(
    title            = "A decision tree for the diagnostic typology",
    subtitle         = "",
    gate             = "Cell eligible at this resolution?",
    gate_sub         = "predator x size-class cell: N_MIN = 5 stomachs with prey in both periods",
    inconclusive     = "Inconclusive",
    inconclusive_sub = "coverage flag — outside the families",
    q1               = "Q1 · Prey composition changed?",
    q1_sub           = "PERMANOVA on Bray–Curtis, trawl set as replicate",
    branch_no        = "composition unchanged  (p ≥ 0.10)",
    branch_yes       = "composition changed  (p < 0.05)",
    branch_trend     = "trend",
    emerging         = "Emerging Shift",
    emerging_sub1    = "comp. trend only (0.05 ≤ p < 0.10)",
    emerging_sub2    = "→ folded into Stability",
    left_header      = "Composition UNCHANGED — same prey identities",
    right_header     = "Composition CHANGED — prey turnover",
    footer           = "",
    # one entry per family: box fill / border / text, and the "sig" chip colours
    palette = list(
      stable = list(fill = "#1F3B4D", stroke = "#142A38", text = "#FFFFFF", chip_fill = "#CFDBE2", chip_text = "#1F3B4D"),
      subst  = list(fill = "#E2742A", stroke = "#B95A1C", text = "#3A1D06", chip_fill = "#F6D8BF", chip_text = "#8A3D12"),
      func   = list(fill = "#9CB7C5", stroke = "#5E8090", text = "#16242B", chip_fill = "#DCE7EC", chip_text = "#1F3B4D"),
      reorg  = list(fill = "#9D2B25", stroke = "#6F1712", text = "#FFFFFF", chip_fill = "#ECC4C0", chip_text = "#7E2018"),
      incon  = list(fill = "#E4E7EA", stroke = "#9AA6B0", text = "#4B5563")
    ),
    # the eight leaves = the 2^3 signal states; chips = c(Comp, H', Bs)
    nodes = list(
      L1 = list(label = "Stability",         sub = "no turnover, no functional change",
                family = "stable", chips = c("ns",  "ns",  "ns")),
      L2 = list(label = "Functional change", sub = "same prey list; diversity re-weighted",
                family = "func",   chips = c("ns",  "sig", "ns")),
      L3 = list(label = "Functional change", sub = "same prey list; niche breadth shifts",
                family = "func",   chips = c("ns",  "ns",  "sig")),
      L4 = list(label = "Functional change", sub = "same prey list; diversity and breadth shift",
                family = "func",   chips = c("ns",  "sig", "sig")),
      R1 = list(label = "Substitution",      sub = "like-for-like turnover; diversity and breadth hold",
                family = "subst",  chips = c("sig", "ns",  "ns")),
      R2 = list(label = "Reorganisation",    sub = "turnover + diversity change",
                family = "reorg",  chips = c("sig", "sig", "ns")),
      R3 = list(label = "Reorganisation",    sub = "turnover + breadth change",
                family = "reorg",  chips = c("sig", "ns",  "sig")),
      R4 = list(label = "Reorganisation",    sub = "turnover + diversity and breadth change",
                family = "reorg",  chips = c("sig", "sig", "sig"))
    ),
    legend = list(
      list(family = "stable", x = 62,  text = "Stability (incl. Emerging)"),
      list(family = "subst",  x = 290, text = "Substitution — like-for-like turnover"),
      list(family = "func",   x = 645, text = "Functional change (same prey)"),
      list(family = "reorg",  x = 905, text = "Reorganisation — multi-axis")
    )
  )
}

# ---------------------------------------------------------------- geometry
.TT_CANVAS <- c(1210, 672)
.TT_LEAF   <- c(w = 215, h = 92)
.TT_SLOTS  <- list(                      # top-left corner of each leaf box
  L1 = c(72, 356),  L2 = c(298, 356), L3 = c(72, 468),  L4 = c(298, 468),
  R1 = c(697, 356), R2 = c(923, 356), R3 = c(697, 468), R4 = c(923, 468)
)

# ---------------------------------------------------------------- helpers
.tt_esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;",  s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}

.tt_text <- function(x, y, txt, size, fill, weight = NULL, anchor = "middle", style = NULL) {
  sprintf('<text x="%s" y="%s" text-anchor="%s" font-size="%s"%s%s fill="%s">%s</text>',
          x, y, anchor, size,
          if (is.null(weight)) "" else sprintf(' font-weight="%s"', weight),
          if (is.null(style))  "" else sprintf(' font-style="%s"', style),
          fill, .tt_esc(txt))
}

.tt_rect <- function(x, y, w, h, fill, stroke, sw = 1.5, rx = 8) {
  sprintf('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s" stroke="%s" stroke-width="%s"/>',
          x, y, w, h, rx, fill, stroke, sw)
}

.tt_arrow <- function(x1, y1, x2, y2, col = "#5a6b7b") {
  sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="1.6" marker-end="url(#arrow)"/>',
          x1, y1, x2, y2, col)
}

.tt_wrap <- function(text, n = 34) {
  if (!nzchar(text)) return(character(0))
  words <- strsplit(text, " ")[[1]]; lines <- character(0); cur <- ""
  for (w in words) {
    if (nchar(cur) + nchar(w) + 1 <= n) cur <- trimws(paste(cur, w)) else { lines <- c(lines, cur); cur <- w }
  }
  if (nzchar(cur)) lines <- c(lines, cur)
  head(lines, 2)
}

.tt_chip <- function(x, y, w, label, sig, fam) {
  if (sig) { fill <- fam$chip_fill; stroke <- fam$stroke; tcol <- if (is.null(fam$chip_text)) fam$text else fam$chip_text }
  else     { fill <- "#ECEFF1"; stroke <- "#B0BAC2"; tcol <- "#607080" }
  paste0(sprintf('<rect x="%s" y="%s" width="%s" height="18" rx="4" fill="%s" stroke="%s"/>', x, y, w, fill, stroke),
         sprintf('<text x="%.1f" y="%s" text-anchor="middle" fill="%s">%s</text>', x + w / 2, y + 12, tcol, .tt_esc(label)))
}

.tt_leaf <- function(slot, node, palette) {
  xy <- .TT_SLOTS[[slot]]; x <- xy[1]; y <- xy[2]
  w <- .TT_LEAF[["w"]]; h <- .TT_LEAF[["h"]]
  fam <- palette[[node$family]]
  lines <- node$label; n <- length(lines)
  out <- .tt_rect(x, y, w, h, fam$fill, fam$stroke, 1.5)
  cx <- x + w / 2
  ly <- y + if (n == 1) 32 else 26
  for (i in seq_along(lines))
    out <- c(out, .tt_text(cx, ly + (i - 1) * 16, lines[i], 14, fam$text, weight = 700))
  sub_y <- ly + n * 16 + 2
  subl  <- .tt_wrap(node$sub)
  for (j in seq_along(subl))
    out <- c(out, .tt_text(cx, sub_y + (j - 1) * 11, subl[j], 9, fam$text))
  ch <- node$chips; cy <- y + 64
  out <- c(out, '<g font-size="9" font-weight="700">',
           .tt_chip(x + 8,   cy, 64, paste("Comp", ch[1]), ch[1] == "sig", fam),
           .tt_chip(x + 78,  cy, 64, paste("H'",   ch[2]), ch[2] == "sig", fam),
           .tt_chip(x + 148, cy, 59, paste("Bs",   ch[3]), ch[3] == "sig", fam),
           '</g>')
  paste(out, collapse = "")
}

# ---------------------------------------------------------------- builder
typology_tree_svg <- function(config = typology_tree_config()) {
  c0  <- config; pal <- c0$palette
  W <- .TT_CANVAS[1]; H <- .TT_CANVAS[2]
  inc <- pal$incon; em <- pal$stable
  s <- c(
    sprintf('<svg viewBox="0 0 %s %s" xmlns="http://www.w3.org/2000/svg" font-family="Helvetica, Arial, sans-serif">', W, H),
    '<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="7.5" refY="3" orient="auto"><path d="M0,0 L8,3 L0,6 Z" fill="#5a6b7b"/></marker></defs>',
    sprintf('<rect x="0" y="0" width="%s" height="%s" fill="#ffffff"/>', W, H),

    # title
    .tt_text(605, 32, c0$title, 22, "#1f2d3d", weight = 700),
    .tt_text(605, 54, c0$subtitle, 13, "#5a6b7b"),

    # eligibility gate and the Inconclusive box
    .tt_rect(420, 78, 370, 46, "#fff", "#34495e", 1.5),
    .tt_text(605, 98,  c0$gate, 12.5, "#2c3e50", weight = 700),
    .tt_text(605, 114, c0$gate_sub, 10, "#5a6b7b"),
    .tt_rect(910, 80, 200, 44, inc$fill, inc$stroke, 1.3),
    .tt_text(1010, 100, c0$inconclusive, 12.5, inc$text, weight = 700),
    .tt_text(1010, 115, c0$inconclusive_sub, 9.5, "#6b7785"),
    .tt_arrow(790, 102, 906, 102),
    .tt_text(846, 96, "no", 10, "#5a6b7b"),
    .tt_arrow(605, 124, 605, 148),
    .tt_text(617, 140, "yes", 10, "#5a6b7b", anchor = "start"),

    # Q1
    .tt_rect(440, 150, 330, 56, "#fff", "#34495e", 1.6),
    .tt_text(605, 172, c0$q1, 13, "#2c3e50", weight = 700),
    .tt_text(605, 190, c0$q1_sub, 10.5, "#5a6b7b"),

    # branches
    .tt_arrow(470, 206, 292, 318, "#2E7D32"),
    .tt_text(300, 240, c0$branch_no, 10, "#2E7D32", weight = 700),
    .tt_arrow(740, 206, 917, 318, "#C0392B"),
    .tt_text(910, 240, c0$branch_yes, 10, "#C0392B", weight = 700),
    .tt_arrow(605, 206, 605, 248),
    .tt_text(612, 230, c0$branch_trend, 9.5, "#5a6b7b", anchor = "start"),

    # Emerging Shift
    .tt_rect(540, 250, 130, 58, em$fill, em$stroke, 1.3),
    .tt_text(605, 272, c0$emerging, 12.5, em$text, weight = 700),
    .tt_text(605, 287, c0$emerging_sub1, 8.5, em$chip_fill),
    .tt_text(605, 299, c0$emerging_sub2, 8.5, em$chip_fill),

    # the two containers
    .tt_rect(60, 320, 465, 256, "#F5F8F5", "#2E7D32", 1.2, rx = 10),
    .tt_text(292, 342, c0$left_header, 12, "#2E7D32", weight = 700),
    .tt_rect(685, 320, 465, 256, "#FBF4F3", "#C0392B", 1.2, rx = 10),
    .tt_text(917, 342, c0$right_header, 12, "#C0392B", weight = 700)
  )

  # leaves
  for (slot in names(.TT_SLOTS)) s <- c(s, .tt_leaf(slot, c0$nodes[[slot]], pal))

  # legend (colours taken from the palette, so they stay in sync)
  s <- c(s, '<line x1="60" y1="592" x2="1150" y2="592" stroke="#d6dde3" stroke-width="1"/>',
         '<g font-size="11">')
  for (lg in c0$legend) {
    fam <- pal[[lg$family]]
    s <- c(s, sprintf('<rect x="%s" y="602" width="14" height="14" rx="3" fill="%s" stroke="%s"/>', lg$x, fam$fill, fam$stroke),
           .tt_text(lg$x + 20, 613, lg$text, 11, "#37474F", anchor = "start"))
  }
  s <- c(s, sprintf('<rect x="62" y="626" width="14" height="14" rx="3" fill="%s" stroke="%s"/>', inc$fill, inc$stroke),
         .tt_text(82, 637, "Inconclusive", 11, "#37474F", anchor = "start"),
         .tt_text(300, 637, "Chips = the three signals (Comp / H' / Bs):  sig = significant change,  ns = no change.",
                  11, "#5a6b7b", anchor = "start"),
         '</g>',
         .tt_text(62, 660, c0$footer, 11, "#5a6b7b", anchor = "start", style = "italic"),
         '</svg>')
  paste(s, collapse = "\n")
}

# ---------------------------------------------------------------- writer
draw_typology_tree <- function(config = typology_tree_config(),
                               out_prefix = "Output_Figures/Fig4_typology_tree",
                               width = 3300) {
  svg <- typology_tree_svg(config)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  svg_path <- paste0(out_prefix, ".svg")
  writeLines(svg, svg_path, useBytes = TRUE)
  written <- svg_path
  if (requireNamespace("rsvg", quietly = TRUE)) {
    png_path <- paste0(out_prefix, ".png")
    rsvg::rsvg_png(charToRaw(svg), png_path, width = as.integer(width))
    written <- c(written, png_path)
  } else {
    message("Package rsvg not installed: SVG only (install.packages(\"rsvg\") for the PNG).")
  }
  message("written: ", paste(written, collapse = ", "))
  invisible(svg)
}
