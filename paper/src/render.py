"""Renders the paper and its supplement to HTML (for the artifact) and LaTeX (for submission).

    python3 paper/src/render.py      # writes paper/paper.html, paper/supplement.html,
                                     # paper/dinara-align.tex, paper/supplement.tex
"""
import html as htmlmod
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# Where the supplement is published, for the main page's link to it.
SUPPLEMENT_URL = "https://claude.ai/artifact/DXUYvmZnYRMz8X5jSCEsBE"
OUT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import main_text as M
import supplement as S

# Full references for the HTML list: authors, title, where (ending in the year), DOI.
FULL = {
 "Gotoh1982": ["Gotoh,O.", "An improved algorithm for matching biological sequences.", "J. Mol. Biol., 162, 705–708.", "10.1016/0022-2836(82)90398-9"],
 "Ukkonen1985": ["Ukkonen,E.", "Algorithms for approximate string matching.", "Information and Control, 64, 100–118.", "10.1016/S0019-9958(85)80046-2"],
 "Myers1999": ["Myers,G.", "A fast bit-vector algorithm for approximate string matching based on dynamic programming.", "J. ACM, 46, 395–415.", "10.1145/316542.316550"],
 "Sosic2017": ["Šošić,M. and Šikić,M.", "Edlib: a C/C++ library for fast, exact sequence alignment using edit distance.", "Bioinformatics, 33, 1394–1395.", "10.1093/bioinformatics/btw753"],
 "Myers1986": ["Myers,E.W.", "An O(ND) difference algorithm and its variations.", "Algorithmica, 1, 251–266.", "10.1007/BF01840446"],
 "MarcoSola2021": ["Marco-Sola,S., Moure,J.C., Moreto,M. and Espinosa,A.", "Fast gap-affine pairwise alignment using the wavefront algorithm.", "Bioinformatics, 37, 456–463.", "10.1093/bioinformatics/btaa777"],
 "MarcoSola2023": ["Marco-Sola,S., Eizenga,J.M., Guarracino,A., Paten,B., Garrison,E. and Moreto,M.", "Optimal gap-affine alignment in O(s) space.", "Bioinformatics, 39, btad074.", "10.1093/bioinformatics/btad074"],
 "GrootKoerkamp2024a": ["Groot Koerkamp,R. and Ivanov,P.", "Exact global alignment using A* with chaining seed heuristic and match pruning.", "Bioinformatics, 40, btae032.", "10.1093/bioinformatics/btae032"],
 "GrootKoerkamp2024b": ["Groot Koerkamp,R.", "A*PA2: up to 19× faster exact global alignment.", "In 24th International Workshop on Algorithms in Bioinformatics (WABI 2024), LIPIcs 312, 17:1–17:25.", "10.4230/LIPIcs.WABI.2024.17"],
 "Wozniak1997": ["Wozniak,A.", "Using video-oriented instructions to speed up sequence comparison.", "Comput. Appl. Biosci., 13, 145–150.", "10.1093/bioinformatics/13.2.145"],
 "Farrar2007": ["Farrar,M.", "Striped Smith–Waterman speeds database searches six times over other SIMD implementations.", "Bioinformatics, 23, 156–161.", "10.1093/bioinformatics/btl582"],
 "Zhao2013": ["Zhao,M., Lee,W.-P., Garrison,E.P. and Marth,G.T.", "SSW library: an SIMD Smith–Waterman C/C++ library for use in genomic applications.", "PLoS One, 8, e82138.", "10.1371/journal.pone.0082138"],
 "Suzuki2018": ["Suzuki,H. and Kasahara,M.", "Introducing difference recurrence relations for faster semi-global alignment of long sequences.", "BMC Bioinformatics, 19, 45.", "10.1186/s12859-018-2014-8"],
 "Li2018": ["Li,H.", "Minimap2: pairwise alignment for nucleotide sequences.", "Bioinformatics, 34, 3094–3100.", "10.1093/bioinformatics/bty191"],
 "Daily2016": ["Daily,J.", "Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments.", "BMC Bioinformatics, 17, 81.", "10.1186/s12859-016-0930-z"],
 "Liu2023": ["Liu,D. and Steinegger,M.", "Block Aligner: an adaptive SIMD-accelerated aligner for sequences and position-specific scoring matrices.", "Bioinformatics, 39, btad487.", "10.1093/bioinformatics/btad487"],
 "Rognes2011": ["Rognes,T.", "Faster Smith–Waterman database searches with inter-sequence SIMD parallelisation.", "BMC Bioinformatics, 12, 221.", "10.1186/1471-2105-12-221"],
 "Rahn2018": ["Rahn,R., Budach,S., Costanza,P., Ehrhardt,M., Hancox,J. and Reinert,K.", "Generic accelerated sequence alignment in SeqAn using vectorization and multi-threading.", "Bioinformatics, 34, 3437–3445.", "10.1093/bioinformatics/bty380"],
 "Kallenborn2026": ["Kallenborn,F., Dabbaghie,F., Steinegger,M. and Schmidt,B.", "Accelign: a GPU-based library for accelerating pairwise sequence alignment.", "BMC Bioinformatics, 27, 137.", "10.1186/s12859-026-06521-0"],
 "Modular": ["Modular Inc.", "Mojo programming language.", "https://docs.modular.com/mojo/", None],
 "Lattner2021": ["Lattner,C., Amini,M., Bondhugula,U. et al.", "MLIR: scaling compiler infrastructure for domain specific computation.", "In 2021 IEEE/ACM International Symposium on Code Generation and Optimization (CGO), pp. 2–14.", "10.1109/CGO51591.2021.9370308"],
 "Eizenga2022": ["Eizenga,J.M. and Paten,B.", "Improving the time and space complexity of the WFA algorithm and generalizing its scoring.", "bioRxiv.", "10.1101/2022.01.12.476087"],
 "hyalite": ["Ferguson,J.", "hyalite: exact, SIMD-accelerated pairwise and database sequence alignment in pure Rust, version 0.4.0.", "https://github.com/Psy-Fer/hyalite", None],
 "Gibrat2018": ["Gibrat,J.-F.", "A short note on dynamic programming in a band.", "BMC Bioinformatics, 19, 226.", "10.1186/s12859-018-2228-9"],
 "Fujiki2020": ["Fujiki,D., Wu,S., Ozog,N., Goliya,K., Blaauw,D., Narayanasamy,S. and Das,R.", "SeedEx: a genome sequencing accelerator for optimal alignments in subminimal space.", "In 2020 53rd Annual IEEE/ACM International Symposium on Microarchitecture (MICRO), pp. 937–950.", "10.1109/MICRO50266.2020.00080"],
 "Gao2021": ["Gao,Y., Liu,Y., Ma,Y., Liu,B., Wang,Y. and Xing,Y.", "abPOA: an SIMD-based C library for fast partial order alignment using adaptive band.", "Bioinformatics, 37, 2209–2211.", "10.1093/bioinformatics/btaa963"],
}

# ------------------------------------------------------------------------------------------- numbering
def number_items(body, tables):
    """Tables and figures numbered in order of first appearance, and section numbers."""
    tabs, figs = {}, {}
    for b in body:
        if b[0] == "table" and b[1] not in tabs:
            tabs[b[1]] = tables[b[1]].get("label", str(len(tabs) + 1))
        if b[0] == "figure":
            figs[b[1]] = str(len(figs) + 1)
    return tabs, figs

# ------------------------------------------------------------------------------------------- inline
UNIT = r"(µs|ms|ns|s|kbp|Mbp|bp|MB|GB|KB|GHz)"

def split_math(text):
    """Pieces of text and maths, the maths kept as written."""
    parts = re.split(r"(\\\(.*?\\\)|\\\[.*?\\\])", text)
    return [(p, p.startswith("\\(") or p.startswith("\\[")) for p in parts if p]

def cite_html(keys, narrative, refs):
    out = []
    for k in keys:
        name, year = refs[k]
        out.append(f'<a class="cite" href="#ref-{k}">{name} ({year})</a>' if narrative else f'<a class="cite" href="#ref-{k}">{name}, {year}</a>')
    return ", ".join(out) if narrative else "(" + "; ".join(out) + ")"

def inline_html(text, tabs, figs, refs):
    out = []
    for piece, is_math in split_math(text):
        if is_math:
            out.append(piece); continue
        t = htmlmod.escape(piece, quote=False)
        t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
        t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
        t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", t)
        t = re.sub(r"\[@!([\w]+)\]", lambda m: cite_html([m.group(1)], True, refs), t)
        t = re.sub(r"\[(@[\w]+(?:;\s*@[\w]+)*)\]", lambda m: cite_html([k.strip()[1:] for k in m.group(1).split(";")], False, refs), t)
        t = re.sub(r"\[#tab:(\w+)\]", lambda m: f'<a href="#tab-{m.group(1)}">Table {tabs[m.group(1)]}</a>', t)
        t = re.sub(r"\[#fig:(\w+)\]", lambda m: f'<a href="#fig-{m.group(1)}">Fig. {figs[m.group(1)]}</a>', t)
        t = re.sub(r"\^([\w/]+)", r"<sup>\1</sup>", t)
        t = re.sub(r"(\d) " + UNIT + r"\b", r"\1&nbsp;\2", t)
        t = re.sub(r"(https://[^\s<]+[^\s<.,;)])", r'<a href="\1">\1</a>', t)
        out.append(t)
    return "".join(out)

def tex_escape(t):
    t = t.replace("\\", "\\textbackslash{}")
    for a, b in (("&", "\\&"), ("%", "\\%"), ("#", "\\#"), ("_", "\\_"), ("$", "\\$")):
        t = t.replace(a, b)
    return t

def inline_tex(text, bibkeys):
    out = []
    for piece, is_math in split_math(text):
        if is_math:
            out.append(piece); continue
        codes = []
        def keep_code(m):
            codes.append(m.group(1)); return f"\x00{len(codes) - 1}\x00"
        t = re.sub(r"`([^`]+)`", keep_code, piece)
        refs = []
        def keep_ref(m):
            refs.append(m.group(0)); return f"\x01{len(refs) - 1}\x01"
        t = re.sub(r"\[@!?[\w]+(?:;\s*@[\w]+)*\]|\[#(?:tab|fig):\w+\]", keep_ref, t)
        urls = []
        def keep_url(m):
            urls.append(m.group(1)); return f"\x02{len(urls) - 1}\x02"
        t = re.sub(r"(https://[^\s]+[^\s.,;)])", keep_url, t)
        t = tex_escape(t)
        t = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", t)
        t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", t)
        t = re.sub(r"\^([\w/]+)", r"$^{\1}$", t)
        for a, b in (("µs", "\\textmu s"), ("×", "$\\times$"), ("≤", "$\\le$"), ("→", "$\\to$"), ("≈", "$\\approx$"),
                     ("“", "``"), ("”", "''"), ("’", "'"), ("A*PA", "A*PA")):
            t = t.replace(a, b)
        t = re.sub(r"(\d),(\d{3})", r"\1{,}\2", t)
        t = re.sub(r"(\d) " + UNIT.replace("µs", "\\\\textmu s") + r"(?![a-zA-Z])", r"\1~\2", t)
        def put_ref(m):
            r = refs[int(m.group(1))]
            if r.startswith("[#tab:"):
                return "Table~\\ref{tab:" + r[6:-1] + "}"
            if r.startswith("[#fig:"):
                return "Fig.~\\ref{fig:" + r[6:-1] + "}"
            if r.startswith("[@!"):
                return "\\citet{" + bibkeys[r[3:-1]] + "}"
            keys = [k.strip()[1:] for k in r[1:-1].split(";")]
            return "\\citep{" + ",".join(bibkeys[k] for k in keys) + "}"
        t = re.sub(r"\x01(\d+)\x01", put_ref, t)
        t = re.sub(r"\x00(\d+)\x00", lambda m: "\\texttt{" + tex_escape(codes[int(m.group(1))]) + "}", t)
        t = re.sub(r"\x02(\d+)\x02", lambda m: "\\url{" + urls[int(m.group(1))] + "}", t)
        out.append(t)
    return "".join(out)

# ------------------------------------------------------------------------------------------- HTML
def section_numbers(body, prefix=""):
    sec = sub = 0
    nums = []
    for b in body:
        if b[0] == "sec":
            sec += 1; sub = 0; nums.append(f"{prefix}{sec}")
        elif b[0] == "sub":
            sub += 1; nums.append(f"{prefix}{sec}.{sub}")
        else:
            nums.append(None)
    return nums

def html_table(key, t, tabs, refs, figs):
    wrap = set(t.get("wrap", []))
    head = "".join(f'<th{" class=\"l\"" if t["align"][i] == "l" else ""}>{inline_html(h, tabs, figs, refs)}</th>' for i, h in enumerate(t["header"]))
    rows = []
    for r in t["rows"]:
        cells = []
        for i, c in enumerate(r):
            cls = []
            if t["align"][i] == "l": cls.append("l")
            if i in wrap: cls.append("wrap")
            if c.startswith("**"): cls.append("best")
            attr = f' class="{" ".join(cls)}"' if cls else ""
            cells.append(f"<td{attr}>{inline_html(c, tabs, figs, refs)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return (f'<div class="table-wrap" id="tab-{key}"><table>\n<caption><b>Table {tabs[key]}.</b> {inline_html(t["caption"], tabs, figs, refs)}</caption>\n'
            f'<thead><tr>{head}</tr></thead>\n<tbody>\n' + "\n".join(rows) + "\n</tbody></table></div>")

def html_figure(key, caption, num, tabs, figs, refs):
    cap = f'<figcaption><b>Fig. {num}.</b> {inline_html(caption, tabs, figs, refs)}</figcaption>'
    assert key == "sweeps", key
    return f'''<figure id="fig-{key}">
  <div class="panel"><div class="panel-label">a</div>
    <div class="fig-box chart-wrap"><svg class="chart" id="chart-divergence" viewBox="0 0 760 360" role="img" aria-label="Time per alignment by divergence"></svg></div>
    <div class="legend" id="legend-divergence"></div></div>
  <div class="panel"><div class="panel-label">b</div>
    <div class="controls" role="group" aria-label="Divergence">
      <button type="button" id="len-5" aria-pressed="false">5% divergence</button>
      <button type="button" id="len-15" aria-pressed="true">15% divergence</button></div>
    <div class="fig-box chart-wrap"><svg class="chart" id="chart-length" viewBox="0 0 760 360" role="img" aria-label="Time per alignment by length"></svg></div>
    <div class="legend" id="legend-length"></div></div>
  {cap}
</figure>'''

def html_body(body, tables, figures, refs, prefix=""):
    tabs, figs = number_items(body, tables)
    nums = section_numbers(body, prefix)
    out, open_col = [], False
    def col():
        nonlocal open_col
        if not open_col:
            out.append('<div class="col">'); open_col = True
    def end_col():
        nonlocal open_col
        if open_col:
            out.append("</div>"); open_col = False
    for b, n in zip(body, nums):
        if b[0] == "sec":
            end_col(); col(); out.append(f'<h2 id="s{n}"><span class="num">{n}</span>{inline_html(b[1], tabs, figs, refs)}</h2>')
        elif b[0] == "sub":
            col(); out.append(f'<h3 id="s{n}"><span class="num">{n}</span>{inline_html(b[1], tabs, figs, refs)}</h3>')
        elif b[0] == "p":
            col()
            if b[1].startswith("\\["):
                out.append(f'<div class="math-block">{b[1]}</div>')
            else:
                out.append(f"<p>{inline_html(b[1], tabs, figs, refs)}</p>")
        elif b[0] in ("list", "bullets"):
            col(); tag = "ol" if b[0] == "list" else "ul"
            out.append(f"<{tag}>" + "".join(f"<li>{inline_html(i, tabs, figs, refs)}</li>" for i in b[1]) + f"</{tag}>")
        elif b[0] == "table":
            end_col(); out.append(html_table(b[1], tables[b[1]], tabs, refs, figs))
        elif b[0] == "figure":
            end_col(); out.append(html_figure(b[1], figures[b[1]], figs[b[1]], tabs, figs, refs))
    end_col()
    return "\n".join(out)

def cited_keys(body, tables):
    text = json.dumps(body) + json.dumps(tables)
    return sorted(set(re.findall(r"@!?(\w+)", text)) & set(FULL))

def reference_list_html(keys, refs):
    def sort_key(k):
        return (FULL[k][0].lower(), refs[k][1])
    items = []
    for k in sorted(keys, key=sort_key):
        a, t, w, doi = FULL[k]
        d = f' <a class="doi" href="https://doi.org/{doi}">doi:{doi}</a>' if doi else ""
        items.append(f'<li id="ref-{k}">{htmlmod.escape(a)} ({refs[k][1]}) {htmlmod.escape(t)} <i>{inline_html(w, {}, {}, refs)}</i>{d}</li>')
    return "<ul class=\"refs\">\n" + "\n".join(items) + "\n</ul>"

def page_assets():
    """The page's stylesheet and chart script, kept in assets/ beside this file."""
    style = open(os.path.join(HERE, "assets", "style.css")).read()
    import charts_data
    charts = open(os.path.join(HERE, "assets", "charts.js")).read().replace("/*DATA*/", charts_data.js())
    return style, charts

def author_line_html():
    """The authors, each from affiliation 1 and a corresponding author."""
    return " and ".join(f"{name}<sup>1,*</sup>" for name, _ in M.AUTHORS)


def author_line_tex():
    """`author_line_html` in LaTeX."""
    return " and ".join(f"{name}$^{{1,*}}$" for name, _ in M.AUTHORS)


def render_html_main():
    style, charts = page_assets()
    abstract = "\n".join(f'<p><b>{h}:</b> {inline_html(t, {}, {}, M.REFS)}</p>' for h, t in M.ABSTRACT)
    body = html_body(M.BODY, M.TABLES, M.FIGURES, M.REFS)
    back = "\n".join(f'<h4>{h}</h4><p>{inline_html(t, {}, {}, M.REFS)}</p>' for h, t in M.BACK)
    refs = reference_list_html(cited_keys(M.BODY + M.ABSTRACT, M.TABLES), M.REFS)
    title = inline_html(M.TITLE, {}, {}, M.REFS).replace("dinara-align:", '<span class="sys">dinara-align</span>:', 1)
    return f'''<title>dinara-align paper</title>
<style>
{style}
</style>
<main>
<header class="col masthead">
  <div class="journal">Bioinformatics · Original Paper · Sequence analysis <span class="draft">Draft</span></div>
  <h1>{title}</h1>
  <div class="authors">{author_line_html()}</div>
  <div class="affil"><sup>1</sup>{M.AFFILIATION}</div>
  <div class="affil">*To whom correspondence should be addressed: {", ".join(email for _, email in M.AUTHORS)}</div>
  <div class="supp-link">Supplementary material: <a href="{SUPPLEMENT_URL}">companion page</a></div>
</header>
<section class="col abstract" id="abstract">
<h4>Abstract</h4>
{abstract}
</section>
{body}
<section class="col back">
{back}
</section>
<section class="col" id="references">
<h2>References</h2>
{refs}
</section>
</main>
<script>
{charts}
</script>
<script>
window.MathJax = {{ tex: {{ inlineMath: [["\\\\(", "\\\\)"]], displayMath: [["\\\\[", "\\\\]"]] }}, svg: {{ fontCache: "global" }} }};
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>
'''

def render_html_supplement():
    style, _ = page_assets()
    body = html_body(S.BODY, S.TABLES, {}, M.REFS, prefix="S")
    return f'''<title>dinara-align supplement</title>
<style>
{style}
</style>
<main>
<header class="col masthead">
  <div class="journal">Bioinformatics · Supplementary material <span class="draft">Draft</span></div>
  <h1>{inline_html(S.TITLE, {}, {}, M.REFS)}</h1>
  <div class="authors">{" and ".join(name for name, _ in M.AUTHORS)}</div>
</header>
{body}
</main>
<script>
window.MathJax = {{ tex: {{ inlineMath: [["\\\\(", "\\\\)"]], displayMath: [["\\\\[", "\\\\]"]] }}, svg: {{ fontCache: "global" }} }};
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>
'''

# ------------------------------------------------------------------------------------------- LaTeX
def tex_table(key, t, wide, bibkeys):
    wrap = set(t.get("wrap", []))
    cols = []
    for i, a in enumerate(t["align"]):
        if i in wrap:
            width = {True: "5.2cm", False: "3.1cm"}[wide]
            cols.append(f">{{\\raggedright\\arraybackslash}}p{{{width}}}")
        else:
            cols.append(a)
    env = "table*" if wide else "table"
    # A table marked `fit` is scaled to the text width when its columns would not fit.
    fit = ("\\resizebox{\\textwidth}{!}{", "}") if t.get("fit") else ("", "")
    body = "\n".join(" & ".join(inline_tex(c, bibkeys) for c in r) + " \\\\" for r in t["rows"])
    head = " & ".join(inline_tex(h, bibkeys) for h in t["header"]) + " \\\\"
    return (f"\\begin{{{env}}}[t]\n\\caption{{{inline_tex(t['caption'], bibkeys)}}}\\label{{tab:{key}}}\n\\centering\\footnotesize\n"
            f"\\setlength{{\\tabcolsep}}{{4pt}}\n{fit[0]}\\begin{{tabular}}{{@{{}}{''.join(cols)}@{{}}}}\n\\toprule\n{head}\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}{fit[1]}\n\\end{{{env}}}\n")

def tikz_line_plot(series, xs, xlog, width, height, ymin, ymax, xticks, xlabel, styles, rotate=False):
    """A TikZ line plot of `series` (name -> list of values or None) on a log y axis."""
    def X(x):
        v = math.log10(x) if xlog else x
        lo = math.log10(xs[0]) if xlog else xs[0]
        hi = math.log10(xs[-1]) if xlog else xs[-1]
        return (v - lo) / (hi - lo) * width
    def Y(y):
        return (math.log10(y) - math.log10(ymin)) / (math.log10(ymax) - math.log10(ymin)) * height
    lines = [f"\\draw[->] (0,0) -- ({width + 0.3:.2f},0);", f"\\draw[->] (0,0) -- (0,{height + 0.3:.2f});"]
    p = math.ceil(math.log10(ymin))
    while p <= math.floor(math.log10(ymax)):
        v = 10 ** p
        lab = (f"{v * 1000:g}\\,\\textmu s" if v < 1 else f"{v:g}\\,ms" if v < 1000 else f"{v / 1000:g}\\,s")
        lines.append(f"\\draw[gray!25] (0,{Y(v):.3f}) -- ({width:.2f},{Y(v):.3f}); \\node[left,font=\\scriptsize] at (0,{Y(v):.3f}) {{{lab}}};")
        p += 1
    for x, lab in xticks:
        node = "anchor=north east,rotate=40,inner sep=1pt" if rotate else "below"
        lines.append(f"\\draw ({X(x):.3f},0) -- ({X(x):.3f},-0.06) node[{node},font=\\scriptsize] {{{lab}}};")
    if xlabel:
        lines.append(f"\\node[font=\\small] at ({width / 2:.2f},-0.75) {{{xlabel}}};")
    for name, values in series.items():
        pts = [(X(x), Y(v[0] if isinstance(v, list) else v), (isinstance(v, list) and len(v) > 1)) for x, v in zip(xs, values) if v]
        coords = " ".join(f"({a:.3f},{b:.3f})" for a, b, _ in pts)
        lines.append(f"\\draw[{styles[name]}] plot coordinates {{{coords}}};")
        for a, b, part in pts:
            lines.append(f"\\draw[{styles[name]},solid,fill={'white' if part else '.'}] ({a:.3f},{b:.3f}) circle (1.1pt);")
    return "\n".join(lines)

def tex_figure(key, caption, bibkeys):
    import charts_data as D
    styles = {"dinara-align": "dinara", "A*PA2-full": "apafull", "A*PA2-simple": "apasimple", "A*PA": "apa",
              "Edlib": "edlib", "BiWFA": "biwfa", "WFA": "wfa"}
    legend = ("\\begin{tikzpicture}[font=\\scriptsize]\n" + "\n".join(
        f"\\draw[{s}] ({k * 2.9:.2f},0) -- ({k * 2.9 + 0.5:.2f},0) node[right,black] {{{n}}};" for k, (n, s) in enumerate(styles.items())) + "\n\\end{tikzpicture}")
    assert key == "sweeps", key
    a = tikz_line_plot(D.DIVERGENCE, list(range(16)), False, 7.0, 4.2, 0.01, 1000, [(d, f"{d}\\%") for d in range(0, 16, 3)],
                       "divergence (\\%), 100\\,kbp pairs", styles)
    b = tikz_line_plot(D.LENGTH15, [3, 10, 30, 100, 300, 1000], True, 7.0, 4.2, 0.01, 100000,
                       [(x, (f"{x}\\,kbp" if x < 1000 else "1\\,Mbp")) for x in [3, 10, 30, 100, 300, 1000]], "length, 15\\% divergence", styles)
    return (f"\\begin{{figure*}}[t]\n\\centering\n\\begin{{minipage}}{{0.49\\textwidth}}\\centering\\textbf{{a}}\\par\\resizebox{{\\linewidth}}{{!}}{{\\begin{{tikzpicture}}\n{a}\n\\end{{tikzpicture}}}}\\end{{minipage}}\\hfill"
            f"\\begin{{minipage}}{{0.49\\textwidth}}\\centering\\textbf{{b}}\\par\\resizebox{{\\linewidth}}{{!}}{{\\begin{{tikzpicture}}\n{b}\n\\end{{tikzpicture}}}}\\end{{minipage}}\n"
            f"\\par\\medskip\\resizebox{{\\linewidth}}{{!}}{{{legend}}}\n\\caption{{{inline_tex(caption, bibkeys)}}}\\label{{fig:{key}}}\n\\end{{figure*}}\n")

TEX_PREAMBLE = r"""\usepackage[a4paper,margin=1.8cm,columnsep=0.7cm]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{array}
\usepackage{longtable}
\usepackage{graphicx}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\usepackage{xcolor}
\usepackage{tikz}
\usepackage[round]{natbib}
\usepackage[hidelinks]{hyperref}
\usepackage{url}
\newcommand{\dinara}{dinara-align}
\tikzset{
  dinara/.style={very thick, black},
  apafull/.style={thick, blue!80!black},
  apasimple/.style={thick, cyan!70!black, dashed},
  apa/.style={thick, violet, dotted},
  edlib/.style={thick, orange!85!black},
  biwfa/.style={thick, red!75!black, dashed},
  wfa/.style={thick, green!50!black},
}
"""

def tex_body(body, tables, figures, bibkeys, wide_tables=()):
    out = []
    for b in body:
        if b[0] == "sec":
            out.append(f"\\section{{{inline_tex(b[1], bibkeys)}}}\n")
        elif b[0] == "sub":
            out.append(f"\\subsection{{{inline_tex(b[1], bibkeys)}}}\n")
        elif b[0] == "p":
            if b[1].startswith("\\["):
                # Two-column pages are narrow: a pair of equations set side by side goes on two lines.
                eq = b[1][2:-2].strip().replace(r", \qquad", r",\\")
                out.append("\\begin{equation}\\begin{gathered}\n" + eq + "\n\\end{gathered}\\end{equation}\n")
            else:
                out.append(inline_tex(b[1], bibkeys) + "\n")
        elif b[0] in ("list", "bullets"):
            env = "enumerate" if b[0] == "list" else "itemize"
            out.append(f"\\begin{{{env}}}\n" + "\n".join(f"  \\item {inline_tex(i, bibkeys)}" for i in b[1]) + f"\n\\end{{{env}}}\n")
        elif b[0] == "table":
            t = tables[b[1]]
            if t.get("long"):
                out.append(tex_longtable(b[1], t, bibkeys))
            else:
                out.append(tex_table(b[1], t, b[1] in wide_tables, bibkeys))
        elif b[0] == "figure":
            out.append(tex_figure(b[1], figures[b[1]], bibkeys))
    return "\n".join(out)

def tex_longtable(key, t, bibkeys):
    wrap = set(t.get("wrap", []))
    cols = "".join(f">{{\\raggedright\\arraybackslash}}p{{{'5.0cm' if i == 1 else '6.6cm'}}}" if i in wrap else a for i, a in enumerate(t["align"]))
    head = " & ".join(inline_tex(h, bibkeys) for h in t["header"]) + " \\\\"
    rows = "\n".join(" & ".join(inline_tex(c, bibkeys) for c in r) + " \\\\" for r in t["rows"])
    return (f"{{\\footnotesize\n\\setlength{{\\LTcapwidth}}{{\\linewidth}}\n\\begin{{longtable}}{{@{{}}{cols}@{{}}}}\n"
            f"\\caption{{{inline_tex(t['caption'], bibkeys)}}}\\label{{tab:{key}}}\\\\\n\\toprule\n{head}\n\\midrule\n\\endfirsthead\n"
            f"\\caption[]{{(continued)}}\\\\\n\\toprule\n{head}\n\\midrule\n\\endhead\n{rows}\n\\bottomrule\n\\end{{longtable}}\n}}\n")

def render_tex_main():
    bk = M.BIBKEYS
    abstract = "\n".join(f"\\noindent\\textbf{{{h}:}} {inline_tex(t, bk)}\\par" for h, t in M.ABSTRACT)
    body = tex_body(M.BODY, M.TABLES, M.FIGURES, bk, wide_tables=("engines", "real", "affine", "batch", "ablation"))
    back = "\n".join(f"\\section*{{{h}}}\n{inline_tex(t, bk)}\n" for h, t in M.BACK)
    return f"""% Generated by paper/src/render.py from main_text.py; edit the source, not this file.
\\documentclass[10pt,twocolumn]{{article}}
{TEX_PREAMBLE}
\\begin{{document}}
\\twocolumn[
\\begin{{@twocolumnfalse}}
{{\\small\\textit{{Bioinformatics}}, Original Paper, Sequence analysis --- draft}}\\par\\medskip
{{\\LARGE\\bfseries {inline_tex(M.TITLE, bk)}\\par}}\\medskip
{{\\large {author_line_tex()}}}\\par
{{\\small $^{{1}}${tex_escape(M.AFFILIATION)}}}\\par
{{\\small $^{{*}}$To whom correspondence should be addressed: {", ".join(email for _, email in M.AUTHORS)}}}\\par\\bigskip
\\fbox{{\\begin{{minipage}}{{0.97\\textwidth}}\\small
{abstract}
\\end{{minipage}}}}
\\bigskip
\\end{{@twocolumnfalse}}
]
{body}
{back}
\\bibliographystyle{{plainnat}}
\\bibliography{{references}}
\\end{{document}}
"""

def render_tex_supplement():
    bk = M.BIBKEYS
    for t in S.TABLES.values():
        if t["label"] == "S1":
            t["long"] = True
    body = tex_body(S.BODY, S.TABLES, {}, bk)
    return f"""% Generated by paper/src/render.py from supplement.py; edit the source, not this file.
\\documentclass[10pt]{{article}}
{TEX_PREAMBLE}
\\renewcommand{{\\thesection}}{{S\\arabic{{section}}}}
\\renewcommand{{\\thetable}}{{S\\arabic{{table}}}}
\\begin{{document}}
{{\\Large\\bfseries {inline_tex(S.TITLE, bk)}\\par}}\\medskip
{{\\large {" and ".join(name for name, _ in M.AUTHORS)}}}\\par\\bigskip
{body}
\\end{{document}}
"""

if __name__ == "__main__":
    open(os.path.join(OUT, "paper.html"), "w").write(render_html_main())
    open(os.path.join(OUT, "supplement.html"), "w").write(render_html_supplement())
    open(os.path.join(OUT, "dinara-align.tex"), "w").write(render_tex_main())
    open(os.path.join(OUT, "supplement.tex"), "w").write(render_tex_supplement())
    print("rendered")
