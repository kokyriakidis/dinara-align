"""The overview figure: how a pair reaches an algorithm, drawn once and rendered to SVG and to TikZ."""

WIDTH, HEIGHT = 760, 338

LANES = [(14, "Unit costs"), (116, "Gap-affine and two-piece costs"), (272, "Batches of short pairs")]

# Boxes: name -> (x, y, width, height, title, note).
BOXES = {
    "dt": (10, 24, 160, 54, "Diagonal transition", "near-identical pairs"),
    "two": (200, 24, 160, 54, "Two-ended transition", "while cheaper than a band"),
    "band": (390, 24, 170, 54, "Band doubling, tiles", "seed heuristic on long pairs"),
    "trace": (590, 24, 160, 54, "Tile-by-tile traceback", "within the memory limit"),
    "one": (10, 130, 180, 48, "One search, fronts kept", "from the side the traceback reads"),
    "both": (290, 130, 180, 48, "Two searches", "as BiWFA"),
    "kept": (570, 130, 180, 48, "Traceback through fronts", "split past the memory limit"),
    "sweep": (290, 206, 180, 44, "Full difference sweep", "own flags, same alignment"),
    "lanes": (10, 282, 240, 48, "One pair per vector lane", "8-bit, band certified per lane"),
    "again": (330, 282, 240, 48, "Uncertified pairs again", "wider band; 16-bit scores on saturation"),
}

# Arrows: polyline points and an optional label (x, y, text lines, anchor).
ARROWS = [
    ([(170, 51), (200, 51)], None),
    ([(360, 51), (390, 51)], (375, 18, ["after 1/24 of the cheapest band"], "middle")),
    ([(560, 51), (590, 51)], None),
    ([(190, 154), (290, 154)], (240, 147, ["steps more than two"], "middle")),
    ([(470, 154), (570, 154)], (520, 147, ["searches meet"], "middle")),
    ([(100, 130), (100, 123), (660, 123), (660, 130)], (380, 118, ["reaches the end alone: traced at once"], "middle")),
    ([(380, 178), (380, 206)], (392, 189, ["projected work above the", "sweep's, after 1/24 of it"], "start")),
    ([(250, 306), (330, 306)], (290, 299, ["not certified"], "middle")),
]

NOTES = [(10, 96, "A pair aligned by diagonal transition is traced back from its own wavefronts.")]


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
