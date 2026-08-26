"""
Typology-tree2.py
=================
Build the ten-level diagnostic-typology decision tree as SVG (+ optional PNG/PDF).
Everything you are likely to change — titles, labels, sub-labels, the three chips
(Comp / H' / Bs = "sig"/"ns"), and the colours — lives in a single config dict
returned by get_default_config(). The geometry (box positions) is in the GEOMETRY
block near the bottom; edit it only if you want to move boxes.

----------------------------------------------------------------------
USE FROM R (reticulate)
----------------------------------------------------------------------
    library(reticulate)
    # py_install("cairosvg")           # once, only needed for PNG/PDF
    source_python("typology_tree.py")  # exposes build() and get_default_config()

    # 1) defaults:
    build(out_prefix = "typology_tree", width = 1600)

    # 2) customise from R (R list <-> Python dict):
    cfg <- get_default_config()
    cfg$title              <- "Diagnostic typology (sGSL, 4T)"
    cfg$nodes$R1$family    <- "ghost"          # change a colour family
    cfg$nodes$R2$label     <- list("Partial", "Diet Shift")  # multi-line label
    cfg$nodes$R2$chips     <- list("sig","sig","ns")          # Comp / H' / Bs

    # swap Partial <-> Structural if your Table 2 defines them the other way:
    cfg$nodes$R2$chips <- list("sig","ns","sig"); cfg$nodes$R2$sub <- "turnover + breadth change"
    cfg$nodes$R3$chips <- list("sig","sig","ns"); cfg$nodes$R3$sub <- "turnover + diversity change"

    build(config = cfg, out_prefix = "fig_typology", width = 1600,
          make_png = TRUE, make_pdf = TRUE)
----------------------------------------------------------------------

Returns the SVG string. Writes <out_prefix>.svg always; .png/.pdf if cairosvg is
installed and requested.
"""
from html import escape


# ===================== EDITABLE CONFIG =====================
def get_default_config():
    return {
        "title": "A decision tree for the diagnostic typology",
        "subtitle": "",
        "gate": "Stratum eligible at this resolution?",
        "gate_sub": "at least N_MIN = 5 valid stomachs in both periods",
        "inconclusive": "Inconclusive",
        "inconclusive_sub": "coverage flag \u2014 outside the families",
        "q1": "Q1 \u00b7 Prey composition changed?",
        "q1_sub": "PERMANOVA \u2014 Bray\u2013Curtis (biomass and occurrence)",
        "branch_no": "composition unchanged  (p > 0.10)",
        "branch_yes": "composition changed  (p \u2264 0.05)",
        "branch_trend": "trend",
        "emerging": "Emerging Shift",
        "emerging_sub1": "comp. trend only (0.05<p\u22640.10)",
        "emerging_sub2": "\u2192 folded into Stability",
        "left_header": "Composition UNCHANGED \u2014 same prey identities",
        "right_header": "Composition CHANGED \u2014 prey turnover",
        "footer": "",
        # family -> colours. chip_fill = fill of a "sig" chip; chip_text optional.
        # ============================================================
        # COLOURS — change a family here and you change every box of
        # that family.  Keys per family:
        #   fill      = box background
        #   stroke    = box border
        #   text      = title + sub-text inside the box
        #   chip_fill = background of a "sig" chip (Comp/H'/Bs)
        #   chip_text = text colour of a "sig" chip (optional; falls
        #               back to `text` if omitted)
        # Which key feeds which leaves:
        #   stable  -> Stability  (+ the "Emerging Shift" box)
        #   func*   -> the 3 Functional-change leaves
        #   subst   -> Substitution
        #   reorg*  -> the 3 Reorganisation leaves
        #   incon   -> Inconclusive
        # Contrast rule: dark fill -> light `text` (#FFFFFF) + dark
        # `chip_text`; light fill -> dark `text`.
        # NB: two things live OUTSIDE this dict and must be kept in
        # sync by hand — the legend swatches (`leg = [...]` in
        # build_svg) and the "Emerging Shift" sub-text colour.
        # ============================================================
        "palette": {
            "stable": {"fill": "#1F3B4D", "stroke": "#142A38", "text": "#FFFFFF", "chip_fill": "#CFDBE2", "chip_text": "#1F3B4D"},
            "subst":  {"fill": "#E2742A", "stroke": "#B95A1C", "text": "#3A1D06", "chip_fill": "#F6D8BF", "chip_text": "#8A3D12"},
            "func":   {"fill": "#9CB7C5", "stroke": "#5E8090", "text": "#16242B", "chip_fill": "#DCE7EC", "chip_text": "#1F3B4D"},
            "func_a": {"fill": "#9CB7C5", "stroke": "#5E8090", "text": "#16242B", "chip_fill": "#DCE7EC", "chip_text": "#1F3B4D"},
            "func_b": {"fill": "#9CB7C5", "stroke": "#5E8090", "text": "#16242B", "chip_fill": "#DCE7EC", "chip_text": "#1F3B4D"},
            "func_c": {"fill": "#9CB7C5", "stroke": "#5E8090", "text": "#16242B", "chip_fill": "#DCE7EC", "chip_text": "#1F3B4D"},
            "reorg":  {"fill": "#9D2B25", "stroke": "#6F1712", "text": "#FFFFFF", "chip_fill": "#ECC4C0", "chip_text": "#7E2018"},
            "reorg_a":{"fill": "#9D2B25", "stroke": "#6F1712", "text": "#FFFFFF", "chip_fill": "#ECC4C0", "chip_text": "#7E2018"},
            "reorg_b":{"fill": "#9D2B25", "stroke": "#6F1712", "text": "#FFFFFF", "chip_fill": "#ECC4C0", "chip_text": "#7E2018"},
            "reorg_c":{"fill": "#9D2B25", "stroke": "#6F1712", "text": "#FFFFFF", "chip_fill": "#ECC4C0", "chip_text": "#7E2018"},
            "incon":  {"fill": "#E4E7EA", "stroke": "#9AA6B0", "text": "#4B5563"},
        },
        # the eight leaves = the 2^3 signal states. chips = (Comp, H', Bs).
        "nodes": {
            "L1": {"label": ["Stability"],         "sub": "no turnover, no functional change",
                   "family": "stable", "chips": ["ns", "ns", "ns"]},
            "L2": {"label": ["Functional change"], "sub": "same prey list; diversity re-weighted",
                   "family": "func_a", "chips": ["ns", "sig", "ns"]},
            "L3": {"label": ["Functional change"], "sub": "same prey list; niche breadth shifts",
                   "family": "func_b", "chips": ["ns", "ns", "sig"]},
            "L4": {"label": ["Functional change"], "sub": "same prey list; diversity and breadth shift",
                   "family": "func_c", "chips": ["ns", "sig", "sig"]},
            "R1": {"label": ["Substitution"],      "sub": "like-for-like turnover; diversity and breadth hold",
                   "family": "subst",  "chips": ["sig", "ns", "ns"]},
            "R2": {"label": ["Reorganisation"],    "sub": "turnover + diversity change",
                   "family": "reorg_a","chips": ["sig", "sig", "ns"]},
            "R3": {"label": ["Reorganisation"],    "sub": "turnover + breadth change",
                   "family": "reorg_b","chips": ["sig", "ns", "sig"]},
            "R4": {"label": ["Reorganisation"],    "sub": "turnover + diversity and breadth change",
                   "family": "reorg_c","chips": ["sig", "sig", "sig"]},
        },
    }


# ===================== GEOMETRY (edit to move boxes) =====================
CANVAS = (1210, 672)
LEAF_W, LEAF_H = 215, 92
SLOTS = {  # top-left corner of each leaf box
    "L1": (72, 356), "L2": (298, 356), "L3": (72, 468), "L4": (298, 468),
    "R1": (697, 356), "R2": (923, 356), "R3": (697, 468), "R4": (923, 468),
}


# ===================== SVG BUILDER =====================
def _esc(s):
    return escape(str(s), quote=False)


def _chip(x, y, w, label, sig, fam):
    if sig:
        fill, stroke = fam["chip_fill"], fam["stroke"]
        tcol = fam.get("chip_text", fam["text"])
    else:
        fill, stroke, tcol = "#ECEFF1", "#B0BAC2", "#607080"
    cx = x + w / 2
    return (f'<rect x="{x}" y="{y}" width="{w}" height="18" rx="4" fill="{fill}" stroke="{stroke}"/>'
            f'<text x="{cx:.1f}" y="{y + 12}" text-anchor="middle" fill="{tcol}">{_esc(label)}</text>')


def _leaf(slot, node, palette):
    x, y = SLOTS[slot]
    fam = palette[node["family"]]
    lines = node["label"]
    n = len(lines)
    parts = [f'<rect x="{x}" y="{y}" width="{LEAF_W}" height="{LEAF_H}" rx="8" '
             f'fill="{fam["fill"]}" stroke="{fam["stroke"]}" stroke-width="1.5"/>']
    cx = x + LEAF_W / 2
    ly = y + (32 if n == 1 else 26)
    for i, ln in enumerate(lines):
        parts.append(f'<text x="{cx}" y="{ly + i * 16}" text-anchor="middle" font-size="14" '
                     f'font-weight="700" fill="{fam["text"]}">{_esc(ln)}</text>')
    sub_y = ly + n * 16 + 2
    for j, sl in enumerate(_wrap(node.get("sub", ""), 34)):
        parts.append(f'<text x="{cx}" y="{sub_y + j * 11}" text-anchor="middle" font-size="9" '
                     f'fill="{fam["text"]}">{_esc(sl)}</text>')
    # chips
    comp, hp, bs = node["chips"][0], node["chips"][1], node["chips"][2]
    cy = y + 64
    parts.append('<g font-size="9" font-weight="700">')
    parts.append(_chip(x + 8,   cy, 64, f"Comp {comp}", comp == "sig", fam))
    parts.append(_chip(x + 78,  cy, 64, f"H' {hp}",     hp == "sig",   fam))
    parts.append(_chip(x + 148, cy, 59, f"Bs {bs}",     bs == "sig",   fam))
    parts.append('</g>')
    return "".join(parts)


def _wrap(text, n):
    if not text:
        return []
    words, lines, cur = text.split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 <= n:
            cur = (cur + " " + w).strip()
        else:
            lines.append(cur); cur = w
    if cur:
        lines.append(cur)
    return lines[:2]


def build_svg(config=None):
    c = config or get_default_config()
    pal = c["palette"]
    W, H = CANVAS
    s = [f'<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" '
         'font-family="Helvetica, Arial, sans-serif">',
         '<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="7.5" refY="3" '
         'orient="auto"><path d="M0,0 L8,3 L0,6 Z" fill="#5a6b7b"/></marker></defs>',
         f'<rect x="0" y="0" width="{W}" height="{H}" fill="#ffffff"/>']

    # title
    s.append(f'<text x="605" y="32" text-anchor="middle" font-size="22" font-weight="700" '
             f'fill="#1f2d3d">{_esc(c["title"])}</text>')
    s.append(f'<text x="605" y="54" text-anchor="middle" font-size="13" '
             f'fill="#5a6b7b">{_esc(c["subtitle"])}</text>')

    # gate + inconclusive
    s.append('<rect x="420" y="78" width="370" height="46" rx="8" fill="#fff" stroke="#34495e" stroke-width="1.5"/>')
    s.append(f'<text x="605" y="98" text-anchor="middle" font-size="12.5" font-weight="700" fill="#2c3e50">{_esc(c["gate"])}</text>')
    s.append(f'<text x="605" y="114" text-anchor="middle" font-size="11" fill="#5a6b7b">{_esc(c["gate_sub"])}</text>')
    inc = pal["incon"]
    s.append(f'<rect x="910" y="80" width="200" height="44" rx="8" fill="{inc["fill"]}" stroke="{inc["stroke"]}" stroke-width="1.3"/>')
    s.append(f'<text x="1010" y="100" text-anchor="middle" font-size="12.5" font-weight="700" fill="{inc["text"]}">{_esc(c["inconclusive"])}</text>')
    s.append(f'<text x="1010" y="115" text-anchor="middle" font-size="9.5" fill="#6b7785">{_esc(c["inconclusive_sub"])}</text>')
    s.append('<line x1="790" y1="102" x2="906" y2="102" stroke="#5a6b7b" stroke-width="1.6" marker-end="url(#arrow)"/>')
    s.append('<text x="846" y="96" text-anchor="middle" font-size="10" fill="#5a6b7b">no</text>')
    s.append('<line x1="605" y1="124" x2="605" y2="148" stroke="#5a6b7b" stroke-width="1.6" marker-end="url(#arrow)"/>')
    s.append('<text x="617" y="140" font-size="10" fill="#5a6b7b">yes</text>')

    # Q1
    s.append('<rect x="440" y="150" width="330" height="56" rx="8" fill="#fff" stroke="#34495e" stroke-width="1.6"/>')
    s.append(f'<text x="605" y="172" text-anchor="middle" font-size="13" font-weight="700" fill="#2c3e50">{_esc(c["q1"])}</text>')
    s.append(f'<text x="605" y="190" text-anchor="middle" font-size="10.5" fill="#5a6b7b">{_esc(c["q1_sub"])}</text>')

    # branches
    s.append('<line x1="470" y1="206" x2="292" y2="318" stroke="#2E7D32" stroke-width="1.6" marker-end="url(#arrow)"/>')
    s.append(f'<text x="300" y="240" text-anchor="middle" font-size="10" font-weight="700" fill="#2E7D32">{_esc(c["branch_no"])}</text>')
    s.append('<line x1="740" y1="206" x2="917" y2="318" stroke="#C0392B" stroke-width="1.6" marker-end="url(#arrow)"/>')
    s.append(f'<text x="910" y="240" text-anchor="middle" font-size="10" font-weight="700" fill="#C0392B">{_esc(c["branch_yes"])}</text>')
    s.append('<line x1="605" y1="206" x2="605" y2="248" stroke="#5a6b7b" stroke-width="1.6" marker-end="url(#arrow)"/>')
    s.append(f'<text x="612" y="230" font-size="9.5" fill="#5a6b7b">{_esc(c["branch_trend"])}</text>')

    # emerging
    em = pal["stable"]
    s.append(f'<rect x="540" y="250" width="130" height="58" rx="8" fill="{em["fill"]}" stroke="{em["stroke"]}" stroke-width="1.3"/>')
    s.append(f'<text x="605" y="272" text-anchor="middle" font-size="12.5" font-weight="700" fill="{em["text"]}">{_esc(c["emerging"])}</text>')
    # "Emerging Shift" sub-text colour is hard-coded (#DCE7EC): keep it
    # readable on the `stable` fill — light text on a dark fill, or vice versa.
    s.append(f'<text x="605" y="287" text-anchor="middle" font-size="8.5" fill="#DCE7EC">{_esc(c["emerging_sub1"])}</text>')
    s.append(f'<text x="605" y="299" text-anchor="middle" font-size="8.5" fill="#DCE7EC">{_esc(c["emerging_sub2"])}</text>')

    # containers
    s.append('<rect x="60" y="320" width="465" height="256" rx="10" fill="#F5F8F5" stroke="#2E7D32" stroke-width="1.2"/>')
    s.append(f'<text x="292" y="342" text-anchor="middle" font-size="12" font-weight="700" fill="#2E7D32">{_esc(c["left_header"])}</text>')
    s.append('<rect x="685" y="320" width="465" height="256" rx="10" fill="#FBF4F3" stroke="#C0392B" stroke-width="1.2"/>')
    s.append(f'<text x="917" y="342" text-anchor="middle" font-size="12" font-weight="700" fill="#C0392B">{_esc(c["right_header"])}</text>')

    # leaves
    for slot in ("L1", "L2", "L3", "L4", "R1", "R2", "R3", "R4"):
        s.append(_leaf(slot, c["nodes"][slot], pal))

    # legend
    s.append('<line x1="60" y1="592" x2="1150" y2="592" stroke="#d6dde3" stroke-width="1"/>')
    # legend swatches are HARD-CODED here — keep these hex values in
    # sync with the `palette` dict above whenever you recolour a family.
    leg = [("#1F3B4D", "#142A38", 62,  "Stability (incl. Emerging)"),
           ("#E2742A", "#B95A1C", 290, "Substitution \u2014 like-for-like turnover"),
           ("#9CB7C5", "#5E8090", 645, "Functional change (same prey)"),
           ("#9D2B25", "#6F1712", 905, "Reorganisation \u2014 multi-axis")]
    s.append('<g font-size="11">')
    for fill, stroke, x, txt in leg:
        s.append(f'<rect x="{x}" y="602" width="14" height="14" rx="3" fill="{fill}" stroke="{stroke}"/>'
                 f'<text x="{x + 20}" y="613" fill="#37474F">{_esc(txt)}</text>')
    s.append(f'<rect x="62" y="626" width="14" height="14" rx="3" fill="{inc["fill"]}" stroke="{inc["stroke"]}"/>'
             '<text x="82" y="637" fill="#37474F">Inconclusive</text>')
    s.append('<text x="300" y="637" fill="#5a6b7b">Chips = the three signals (Comp / H\' / Bs):  '
             'sig = significant change,  ns = no change.</text>')
    s.append('</g>')
    s.append(f'<text x="62" y="660" font-size="11" font-style="italic" fill="#5a6b7b">{_esc(c["footer"])}</text>')

    s.append('</svg>')
    return "".join(s)


def build(config=None, out_prefix="typology_tree", width=1600,
          make_png=True, make_pdf=True):
    """Write <out_prefix>.svg (+ .png/.pdf if cairosvg available). Returns the SVG string."""
    svg = build_svg(config)
    svg_path = f"{out_prefix}.svg"
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg)
    written = [svg_path]
    if make_png or make_pdf:
        try:
            import cairosvg
            if make_png:
                cairosvg.svg2png(bytestring=svg.encode("utf-8"),
                                 write_to=f"{out_prefix}.png", output_width=int(width))
                written.append(f"{out_prefix}.png")
            if make_pdf:
                cairosvg.svg2pdf(bytestring=svg.encode("utf-8"), write_to=f"{out_prefix}.pdf")
                written.append(f"{out_prefix}.pdf")
        except ImportError:
            print("cairosvg not installed -> SVG only. Install with: pip install cairosvg")
    print("written:", ", ".join(written))
    return svg


if __name__ == "__main__":
    build(out_prefix="Typology-tree2")
