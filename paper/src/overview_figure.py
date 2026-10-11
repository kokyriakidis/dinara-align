"""The overview figure: how a pair reaches an algorithm, drawn once and rendered to SVG and to TikZ."""

WIDTH, HEIGHT = 760, 318

LANES = [(14, "Unit costs"), (110, "Gap-affine and two-piece costs"), (250, "Batches of short pairs")]

# Boxes: name -> (x, y, width, height, title, note).
BOXES = {
    "dt": (10, 24, 160, 54, "Diagonal transition", "near-identical pairs"),
    "two": (200, 24, 160, 54, "Two-ended transition", "while cheaper than a band"),
    "band": (390, 24, 170, 54, "Band doubling, tiles", "seed heuristic on long pairs"),
    "trace": (590, 24, 160, 54, "Tile-by-tile traceback", "within the memory limit"),
    "wave": (10, 120, 200, 48, "Bidirectional wavefront", "WFA, BiWFA"),
    "kept": (560, 120, 190, 48, "Traceback through fronts", "split past the memory limit"),
    "sweep": (330, 184, 200, 48, "Full difference sweep", "8-bit lanes, same alignment"),
    "lanes": (10, 260, 240, 48, "One pair per vector lane", "8-bit, band certified per lane"),
    "again": (330, 260, 240, 48, "Uncertified pairs again", "a wider band; 16 bits if a byte saturates"),
}

# Arrows: polyline points and an optional label (x, y, text lines, anchor).
ARROWS = [
    ([(170, 51), (200, 51)], None),
    ([(360, 51), (390, 51)], None),
    ([(560, 51), (590, 51)], None),
    ([(210, 144), (560, 144)], (385, 137, ["the wavefront finishes first"], "middle")),
    ([(110, 168), (110, 208), (330, 208)], (122, 192, ["projected work above the", "sweep's, after 1/24 of it"], "start")),
    ([(530, 208), (655, 208), (655, 168)], None),
    ([(250, 284), (330, 284)], (290, 277, ["not certified"], "middle")),
]

NOTES = [(10, 96, "A pair settled by diagonal transition is traced back from its own wavefronts.")]


def svg():
    """The figure as inline SVG, coloured by the page's theme tokens."""
    out = [f'<svg class="chart overview" viewBox="0 0 {WIDTH} {HEIGHT}" role="img" '
           'aria-label="How dinara-align assigns a pair to an algorithm">',
           '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
           'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="var(--ink-soft)"/></marker></defs>']
    for y, title in LANES:
        out.append(f'<text x="10" y="{y}" class="lane">{title.upper()}</text>')
    for x, y, w, h, title, note in BOXES.values():
        out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" class="box"/>')
        out.append(f'<text x="{x + w / 2}" y="{y + h / 2 - 4}" text-anchor="middle" class="box-title">{title}</text>')
        out.append(f'<text x="{x + w / 2}" y="{y + h / 2 + 13}" text-anchor="middle" class="box-note">{note}</text>')
    for points, label in ARROWS:
        pts = " ".join(f"{x},{y}" for x, y in points)
        out.append(f'<polyline points="{pts}" fill="none" class="flow" marker-end="url(#arrow)"/>')
        if label:
            x, y, lines, anchor = label
            for k, line in enumerate(lines):
                out.append(f'<text x="{x}" y="{y + 13 * k}" text-anchor="{anchor}" class="flow-label">{line}</text>')
    for x, y, text in NOTES:
        out.append(f'<text x="{x}" y="{y}" class="flow-label">{text}</text>')
    out.append("</svg>")
    return "\n".join(out)


def tikz(scale=0.0224):
    """The figure as TikZ, the SVG's coordinates scaled to centimetres, y pointing up."""
    def P(x, y):
        return f"({x * scale:.3f},{-y * scale:.3f})"
    out = []
    for y, title in LANES:
        out.append(f"\\node[anchor=west,font=\\scriptsize\\scshape] at {P(10, y - 4)} {{{title}}};")
    for x, y, w, h, title, note in BOXES.values():
        out.append(f"\\draw[rounded corners=2pt,fill=gray!8,draw=gray!60] {P(x, y)} rectangle {P(x + w, y + h)};")
        out.append(f"\\node[align=center,font=\\scriptsize] at {P(x + w / 2, y + h / 2)} "
                   f"{{\\textbf{{{title}}}\\\\{{\\color{{gray!70!black}}{note}}}}};")
    for points, label in ARROWS:
        out.append("\\draw[->,>=stealth,gray!70!black] " + " -- ".join(P(x, y) for x, y in points) + ";")
        if label:
            x, y, lines, anchor = label
            where = {"middle": "", "start": "anchor=west,"}[anchor]
            text = "\\\\".join(lines)
            out.append(f"\\node[{where}align={'left' if anchor == 'start' else 'center'},font=\\tiny] at "
                       f"{P(x, y + 6 * (len(lines) - 1) - 3)} {{{text}}};")
    for x, y, text in NOTES:
        out.append(f"\\node[anchor=west,font=\\tiny] at {P(x, y - 3)} {{{text}}};")
    return "\n".join(out)
