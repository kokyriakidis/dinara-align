"""The Methods additions and the supplement, written once and rendered to HTML and to LaTeX.
Blocks: ("p", head, text) a paragraph with an optional bold lead; ("ul", [items]); ("h", level, title, id);
("table", caption, header, rows)."""
import re

DT_ITEMS = [
 "<b>Fronts eight at a time.</b> Fronts are padded with unreached diagonals and stepped eight diagonals at a time with branch-free validity masks, the step in a function of its own. A step costs about 1.3 ns, down from 5. The two-ended distance checks the fronts' overlap with vector compares and takes about half the steps.",
 "<b>Budgets.</b> A step costs about 1.3 ns and a band column about 6 ns, so the diagonal transition gets a step budget proportional to the band it would replace. Cutting the alignment budget from 12 to 5 steps a column took 10 kbp pairs at 5% from 241 to 144 µs.",
 "<b>One front or two.</b> The switch between one front and two from both ends was fitted on 1,300 pairs of 1 to 30 kbp. The CIGAR is written straight from the fronts.",
 "<b>Short pairs stay.</b> A pair of up to 2,048 columns had its band aimed wide over the whole matrix, and 42% of such pairs went to a band the diagonal transition would have beaten. Such a band is now charged 30,000 steps more. ont-1k became 25% faster on the Skylake-X, 42.3 to 31.8 ms over 1,221 reads, which put dinara-align ahead of WFA on short reads.",
 "<b>Short noisy pairs.</b> A short noisy pair gives up the diagonal search for the whole matrix: a 1 kbp pair at 15% fell from 19.2 to 11.8 µs. A short pair's first bound is 1.7 times its projection, at most 128 past it.",
]

METHODS_UNIT = [
 ("p", "Trusting a projection.", "The probe projects the distance from its first edits. Real reads gather their errors at the ends, so that projection can be hundreds of times too high: a 30 kbp SARS-CoV-2 pair at distance 49 was projected at 14,448 and swept whole. dinara-align therefore treats a projection with suspicion:"),
 ("ul", [
   "The diagonal transition gives up on its projection only after spending a twenty-fourth of what the cheapest band would cost.",
   "A projection aims the band only when the one-front projection from the first edits and the search's projection from many more agree within a third.",
   "An untrusted projection starts the first bound 256 past what is already known, the floor or the heuristic at the origin, and every retry at most doubles, as A*PA2's band doubling grows.",
   "A later checkpoint never sets the bound below a quarter past the last bound, and a failed round grows from the bound it ended on, not the one it started on.",
 ]),
 ("p", None, "These rules took ont-500k-genvar from 4.33 s, with 13 of 15 reads finished, to 754 ms, the largest single step in its history."),
 ("p", "Seeded retries.", "With seeds, a retry aims a quarter of the projected climb above the origin's bound past the estimate; half a climb had failed 90% of the way across a 527 kbp read. A seeded round that dies within its first eighth of the columns does not project at all: its margin over the origin's bound grows fourfold instead, as A*PA2's does. An 884 kbp read had died at column 1,280 and projected 125,000 for a distance of 75,000. Only the first round re-aims at checkpoints; a lowered later round died most of the way across, four times over on one 897 kbp read."),
 ("p", "Tracing the band back.", "The traceback walks the band right to left a tile at a time, from the tile edges the sweep kept."),
 ("ul", [
   "Each tile is traced by a forward wavefront search from its recorded left edge, kept to the diagonals that can still reach the traced cell. It steps eight diagonals at a time, with gathered slides on AVX-512. On 100 kbp pairs at 15% the tile searches take 1.4 ms a pair on the Skylake-X, against 3.55 ms before.",
   "A tile the search gives up on is recomputed over a window of rows above the traced cell, starting at four words and doubling, rather than over the band's whole height. The window's top reads +1, so its scores are real paths' scores. A 600 kbp nanopore pair at 13% divergence aligned in 617 ms, where it had taken 1.38 s, 840 ms of it recomputing 900-word tiles.",
   "Recompute buffers are reused and left uninitialised, and tile edges are sized once and written through pointers. The edge loads had been a quarter of the traceback's time; ont-500k became 2 to 3% faster.",
   "The CIGAR is written by reversing moves sixteen at a time and splitting diagonal runs by comparing eight bases at a time, within 1 to 3% of not writing one.",
 ]),
 ("p", "Diagonal transition, engineered.", "Beyond the slides, the diagonal transition is shaped by measured costs:"),
 ("ul", DT_ITEMS),
]

SEED_ITEMS = [
 "<b>Inexact matches by halves.</b> A seed matched within one edit matches one of its 8-base halves exactly. Each half indexes a table of the seeds holding it, and every window of the query is looked up as both halves. An entry holds its seed beside its code, so one load gives both.",
 "<b>Only seeds that can still reach the end.</b> Within a bucket the seeds go in order, and only a range of them can still chain to the end from a given row, so a lookup tries that range alone. Local pruning's diagonal transition updates its fronts in place, with an unreachable front on either side in place of three bounds tests per diagonal. With the branch-free test, this took the inexact setup on genvar from 5.02 to 2.90 s and the exact setup from 1.69 to 1.21 s, and whole alignments of genvar from 11.77 to 10.21 s, on the Skylake-X.",
]

METHODS_MODES = [
 ("h", 3, "Other modes", "modes"),
 ("p", "Extension.", "An extension runs one wavefront search from its anchor, the match reward folded into costs, and takes the best stop from each front's furthest anti-diagonal. It stops once a proven bound shows that no later cost can win. The prefix up to the best stop is then aligned globally. With an end bonus, the same search keeps the best point on the query's last row and goes on only while such an alignment could still come within the bonus. An extension with an end bonus took 186 µs where it had taken 1.02 ms; it now takes 158 µs, against KSW2's 1.23 ms."),
 ("p", "Free ends.", "Under gap costs with free ends, the forward search starts at cost zero on every diagonal a free start lies on. The span comes first: one search from the free side finds the end on the highest diagonal, and a narrow search back finds the start. The span is then aligned globally. A read at a reference's start, an overlap or a suffix aligns 1.4 to 3.5 times faster on x86; a read anywhere in the reference aligns 0 to 14% slower. Keeping the search outside an optional wrapper while it runs, so that its fields stay in registers, took an infix at (4, 6, 2) from 3.92 to 3.00 ms."),
 ("p", "Infix and prefix at unit costs.", "An infix at unit costs runs the bit-parallel sweep with Edlib's cutoff, the bound doubling from 64. Words entering the band restart from +1, and the last word steps a column at a time. On the M2 a 10 kbp read in 100 kbp went from 21 to 13 ms, and 100 kbp in 1 Mbp from 2 to 0.75 s. A prefix search sweeps only the band its bound allows:"),
 ("ul", [
   "the band's top moves down with the diagonal;",
   "no column past the read's length plus the best end so far is swept;",
   "a try whose band has emptied stops;",
   "at an eighth, a quarter and a half of the columns, the try projects the distance and gives itself up when it is doomed.",
 ]),
 ("p", None, "Giving up doomed tries took its distance from 23.7 to 20.5 µs on the Skylake-X; its alignment takes 43 µs, against WFA2-lib's 47 µs and Edlib's 114 µs."),
]

TABLE_BAND = ("p", "A band from the score.", "A global alignment under a table of one match and one mismatch score stores only the diagonals its score allows. In the folded costs of Eq. (1), every gapped letter costs at least \\(2e+a\\) of the total \\(C=a(n+m)-2S\\), and a path straying \\(t\\) diagonals from both the main and the end diagonal has at least \\(2t\\) gapped letters. So only diagonals within \\(C/(2(2e+a))\\) of them can hold an optimal path. A 1 kbp pair at 5% stores 137 diagonals instead of 1,000 and aligns in 0.21 ms instead of 0.36 ms. Under any other table, the score from the sweep bounds the band the same way.")

SUPPLEMENT = [
 ("h", 2, "Supplementary material", "supplement"),
 ("p", None, "This supplement gives the measured effect of each change behind the results, as the commit that made it records it, and the work that was tried and dropped. Times are as measured at the time of the change, not at the published commit."),
 ("h", 3, "S1. Measured effect of each change", "s1"),
 ("table", "Table S1. Changes and their measured effects. Machines: the i9-7900X (Skylake-X, AVX-512) and the Apple M2 (NEON).",
  ["Commit", "Change", "Measured effect", "Machine"],
  [
   ["a6101e5", "Inexact seeds (r = 2), matched by halves", "100 kbp 15.6 → 9.5 ms; 300 kbp 105 → 47 ms; 1 Mbp 1.36 → 0.41 s", "M2"],
   ["7435edd", "Leaner seed setup: branch-free test, reachable seeds only, fronts in place", "inexact setup on genvar 5.02 → 2.90 s; genvar 11.77 → 10.21 s; ont-500k 8.04 → 7.11 s", "Skylake-X"],
   ["28b9e67", "Inexact-seed cut-off by vector width (20% on AVX-512)", "ont-500k (50 reads) 8.73 → 8.03 s; genvar (48) 12.03 → 11.73 s; 100 kbp at 6% 114 → 64 ms", "Skylake-X"],
   ["708f151", "Seeded retry a quarter of the climb past the estimate", "genvar 9.86 → 9.66 s; ont-500k 6.81 → 6.30 s; 100 kbp at 15% 4% faster", "Skylake-X"],
   ["ff26cea", "Seed gate at 86 kbp on AVX-512; exact-seed filter", "30 kbp at 5% 1.48 ms → 898 µs; SARS-CoV-2 373 → 283 µs; ont-50k 897 → 800 µs", "Skylake-X"],
   ["ff26cea", "Exact-seed filter alone", "matching on 320 mid-length reads 87 → 32 ms; 30 kbp at 5% 15% faster", "M2"],
   ["45697e0", "AVX-512 gathers in the diagonal transition", "ont-1k 26.8 → 25.0 µs; SARS-CoV-2 335 → 306 µs; 100 kbp at 1% 1.70 → 1.45 ms, at 2% 5.60 → 4.54 ms", "Skylake-X"],
   ["2c7c46e", "Short pairs stay in the diagonal transition", "ont-1k 42.3 → 31.8 ms over 1,221 reads", "Skylake-X"],
   ["e0ea104", "Alignment step budget 12 → 5 a column", "10 kbp at 5% 241 → 144 µs", "not recorded"],
   ["a593ade", "Short noisy pairs give up the diagonal search", "1 kbp at 15% 19.2 → 11.8 µs", "Skylake-X"],
   ["beb7a8c", "Tiles traced by a forward search", "ont-50k 3% and ont-500k 8% faster", "M2"],
   ["5a0e514", "Tile search eight diagonals at a time, gathered", "100 kbp at 15% tile searches 3.55 → 1.4 ms a pair", "Skylake-X"],
   ["7f1a813", "Windowed recompute of a failed tile", "600 kbp at 13% 1.38 s → 617 ms", "not recorded"],
   ["3becb4b", "Recompute buffers reused, uninitialised", "ont-1k 23.4 → 22.6 µs; ont-10k 151 → 142 µs", "Skylake-X"],
   ["e2125b0", "Tile edges sized once, written through pointers", "ont-500k 2 to 3%, ont-10k about 2% faster", "both"],
   ["38837f2", "Fronts sixteen offsets a store; unchecked reads", "ont-1k 21.4 → 20.7 µs; SARS-CoV-2 213 → 199 µs", "Skylake-X"],
   ["34cb623", "Pointers in hot lookups", "genvar 10.16 → 9.90 s; ont-500k 7.04 → 6.82 s", "Skylake-X"],
   ["f29a945", "Separate seeded and unseeded band copies", "unseeded pairs back within 2% (had slowed 3 to 7%)", "Skylake-X"],
   ["aca56c7", "AVX-512 gathers in the gap-affine step", "ont-1k 149 → 129 µs; ont-10k 3.52 → 2.94 ms; 10 kbp at 15% 36.4 → 24.5 ms", "Skylake-X"],
   ["53e15d4", "Gather mask from the group's own lanes; one reduction after the loop", "ont-1k 129 → 111 µs; ont-10k 2.34 → 1.84 ms; 10 kbp at 15% 17.9 → 14.0 ms", "Skylake-X"],
   ["ced0b8d", "Ring rows padded against 4 KB aliasing", "ont-10k 1.86 → 1.68 ms; 100 kbp at 5% 349 → 305 ms", "Skylake-X"],
   ["17fd781", "Slides straight from the step", "100 kbp at 5% 315 → 287 ms; ont-10k 1.80 → 1.67 ms", "M2"],
   ["01c1e9b", "Slide columns through scalar loads and stores", "100 kbp at 5% 287 → 271 ms; ont-1k 91 → 86 µs", "M2"],
   ["72c3696", "Kept fronts in fixed blocks, never moved", "ont-10k 1.48 → 1.37 ms (memmove had been 7.5% of samples)", "M2"],
   ["b7e5232", "Gap-affine traceback through the wavefront's own fronts", "SARS-CoV-2 1.34 ms → 638 µs; ont-10k 14.6 → 4.05 ms", "M2"],
   ["e6a435f", "End bonus weighed in the extension's own search", "with an end bonus 1.02 ms → 186 µs (KSW2 1.23 ms)", "Skylake-X"],
   ["8ba2a6a", "Free ends: span first, start at cost zero on free diagonals", "start, overlap, suffix 1.4 to 3.5× faster; anywhere 0 to 14% slower", "Skylake-X"],
   ["f0bb161", "Unit-cost infix: Ukkonen cutoff, restarted words", "10 kbp in 100 kbp 21 → 13 ms; 100 kbp in 1 Mbp 2 → 0.75 s", "M2"],
   ["2f8b935", "Prefix search gives up a doomed try", "distance 23.7 → 20.5 µs; alignment 45.8 → 42.9 µs", "Skylake-X"],
   ["2b038c8", "Global band bounded by the score (one match, one mismatch score)", "1 kbp at 5%: 137 of 1,000 diagonals; 0.36 → 0.21 ms", "not recorded"],
   ["756f3e1", "Table sweeps vectorised, any table", "1 kbp global 13.3 ms → 550 µs; local 133 → 2.89 ms", "Skylake-X"],
   ["433a8b5", "Local under gap costs: lanes along the shorter sequence, extension back", "47 vs 55 µs short pairs; 1.93 vs 3.22 ms 1 kbp in 10 kbp (SSW)", "Skylake-X"],
   ["cf87b2b", "Local sweep in 16-bit lanes", "3.4 ms, from 4.6", "Skylake-X"],
   ["a8355aa", "Batch letters transposed in registers", "affine 0.59 → 0.43 µs a pair; unit 0.55 → 0.40 µs", "not recorded"],
   ["20425a7", "Batch first pass in saturating bytes", "affine 0.43 → 0.40 µs a pair", "Skylake-X"],
   ["64f0cc3", "Second-pass bands grouped by exact width", "1 kbp at 10% scores 18 → 15 µs; 46 → 38 µs at (4, 6, 2)", "not recorded"],
   ["0f7b5af", "A lane group stops once every lane has passed its cap", "1 kbp at 10% aligns 33 → 22 µs, scores 18 → 12 µs", "Skylake-X"],
   ["792ba9f", "Pairs with N take seeds; a seed with N goes uncounted", "500 kbp reads 2.5 to 3× faster (M2), 1.4 to 1.7× (Skylake-X)", "both"],
  ]),
 ("h", 3, "S2. Robustness", "s2"),
 ("ul", [
   "Layers that spill into a staircase, and a bounded pruning slide: poly-A of 32,000 bases against 16,000 took 17 s and now takes 0.15 s.",
   "Seeds given up past 64 matches a seed on average: a 200 kbp tandem repeat took 73 s and 760 MB to set up, and now 0.045 s and 34 MB.",
   "A front that reached nothing checks nothing: a case that took 1.5 s now takes 13 ms.",
 ]),
 ("h", 3, "S3. Tried and dropped", "s3"),
 ("ul", [
   "<b>Threading within a pair.</b> Striped bands over eight threads aligned 1 Mbp pairs at 15% in 575 ms against 1,122 ms on one. It was built and then removed: the library runs on its caller's thread, and spreading pairs over threads gains more.",
   "<b>A third bit-parallel group on AVX-512.</b> It changed nothing within 1.5%.",
   "<b>One sixteen-lane group</b> instead of two of eight: as good on long reads, up to 2% slower on short divergent pairs.",
   "<b>A fully branch-free seed scan.</b> It would save 1 to 2% for about 25 million cycles a pair spent.",
   "<b>A Bloom gate for inexact keys.</b> It halved the candidates, but alignments gained only 1.5 to 2%.",
   "<b>Shorter lookahead, or no local pruning.</b> Without pruning, genvar was 21% and ont-500k 72% slower.",
   "<b>A* on the diagonal transition</b> and <b>an upper bound by beam search.</b> Each needed several tuned parts for small gains.",
   "<b>Traceback while sweeping.</b> The tile edges take 0.6 to 13 MB, against 10 to 17 MB of sequences and planes, so there was little memory to save.",
 ]),
]

# ---------------------------------------------------------------------------------------------- rendering
def html_text(t):
    t = t.replace("dinara-align", '<span class="sys">dinara-align</span>')
    for unit in ("µs", "ms", "kbp", "Mbp", "MB", "KB", "bp", "s"):
        t = re.sub(r"(\d) " + unit + r"\b", r"\1&nbsp;" + unit, t)
    return t

def tex_text(t):
    t = t.replace("<b>", "\\textbf{").replace("</b>", "}")
    t = t.replace("dinara-align", "\\dinara{}")
    t = t.replace("%", "\\%").replace("→", "$\\to$").replace("×", "$\\times$").replace("≤", "$\\le$")
    t = t.replace("µs", "\\textmu s").replace("*", "*")
    t = re.sub(r"\\\((.*?)\\\)", lambda m: "$" + m.group(1) + "$", t)
    t = t.replace("Eq. (1)", "Eq.~(1)")
    t = re.sub(r"(\d),(\d{3})", r"\1{,}\2", t)
    return t

def to_html(blocks, sec_num=None):
    out = []
    for b in blocks:
        if b[0] == "p":
            head = f'<span class="para-head">{b[1]}</span> ' if b[1] else ""
            out.append(f"<p>{head}{html_text(b[2])}</p>")
        elif b[0] == "ul":
            out.append("<ul>\n" + "\n".join(f"  <li>{html_text(i)}</li>" for i in b[1]) + "\n</ul>")
        elif b[0] == "h":
            tag = "h2" if b[1] == 2 else "h3"
            num = f'<span class="num">{sec_num}</span>' if (sec_num and b[1] == 3) else ""
            out.append(f'<{tag} id="{b[3]}">{num}{b[2]}</{tag}>')
        elif b[0] == "table":
            rows = "\n".join("    <tr>" + "".join(f'<td class="wrap">{html_text(c)}</td>' if k in (1, 2) else f"<td>{html_text(c)}</td>" for k, c in enumerate(r)) + "</tr>" for r in b[3])
            head = "".join(f"<th style=\"text-align:left\">{h}</th>" for h in b[2])
            cap = b[1].split(". ", 1)
            out.append(f'</div>\n<div class="table-wrap">\n<table class="supp">\n  <caption><b>{cap[0]}.</b> {html_text(cap[1])}</caption>\n  <thead><tr>{head}</tr></thead>\n  <tbody>\n{rows}\n  </tbody>\n</table>\n</div>\n<div class="col">')
    return "\n".join(out)

def to_tex(blocks):
    out = []
    for b in blocks:
        if b[0] == "p":
            head = f"\\paragraph{{{b[1]}}} " if b[1] else ""
            out.append(head + tex_text(b[2]) + "\n")
        elif b[0] == "ul":
            out.append("\\begin{itemize}\n" + "\n".join(f"  \\item {tex_text(i)}" for i in b[1]) + "\n\\end{itemize}\n")
        elif b[0] == "h":
            cmd = {2: "\\section", 3: "\\subsection"}[b[1]]
            title = re.sub(r"^S\d\. ", "", b[2])
            out.append(f"{cmd}{{{title}}}\\label{{sec:{b[3]}}}\n")
        elif b[0] == "table":
            rows = "\n".join(" & ".join(tex_text(c) for c in r) + "\\\\" for r in b[3])
            cap = tex_text(b[1].split(". ", 1)[1])
            out.append("{\\footnotesize\n\\setlength{\\LTcapwidth}{\\linewidth}\n\\begin{longtable}{@{}l>{\\raggedright\\arraybackslash}p{4.6cm}>{\\raggedright\\arraybackslash}p{6.2cm}l@{}}\n"
                       f"\\caption{{{cap}}}\\label{{tab:changes}}\\\\\n\\toprule\nCommit & Change & Measured effect & Machine\\\\\n\\midrule\n\\endhead\n{rows}\n\\bottomrule\n\\end{{longtable}}\n}}\n")
    return "\n".join(out)
