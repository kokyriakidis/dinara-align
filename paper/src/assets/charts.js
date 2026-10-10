/* Charts: one log-scale line chart, drawn to scale from the published tables (benchmarks/README.md, ebe290f). */
const TOOLS = [
  ["dinara-align", "--c-dinara", 3.2, ""], ["A*PA2-full", "--c-apafull", 1.8, ""], ["A*PA2-simple", "--c-apasimple", 1.8, "6 4"],
  ["A*PA", "--c-apa", 1.8, "2 3"], ["Edlib", "--c-edlib", 1.8, ""], ["BiWFA", "--c-biwfa", 1.8, "6 4"], ["WFA", "--c-wfa", 1.8, ""],
];
/*DATA*/

const css = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const fmt = ms => ms >= 1000 ? (ms / 1000).toPrecision(3).replace(/\.?0+$/, "") + " s" : ms >= 1 ? (+ms.toPrecision(3)) + " ms" : (+(ms * 1000).toPrecision(3)) + " µs";
const NS = "http://www.w3.org/2000/svg";
function el(tag, attrs, parent) { const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }

function lineChart(svg, cfg) {
  // cfg: xs (numbers), xLog, xLabel, xTick(x) label, series {name: [[ms, partial?] | null]}, styles {name: [token, width, dash]},
  // hidden Set, yMin, yMax (ms, powers of ten), tooltipTitle(x)
  const W = 760, H = +svg.getAttribute("viewBox").split(" ")[3], m = { l: 64, r: 34, t: 16, b: 52 };
  svg.replaceChildren();
  const xv = x => cfg.xLog ? Math.log10(x) : x;
  const x0 = xv(cfg.xs[0]), x1 = xv(cfg.xs[cfg.xs.length - 1]);
  const X = x => m.l + (xv(x) - x0) / (x1 - x0) * (W - m.l - m.r);
  const ly0 = Math.log10(cfg.yMin), ly1 = Math.log10(cfg.yMax);
  const Y = v => H - m.b - (Math.log10(v) - ly0) / (ly1 - ly0) * (H - m.t - m.b);
  for (let p = Math.ceil(ly0); p <= Math.floor(ly1); p++) {
    const v = Math.pow(10, p), y = Y(v);
    el("line", { x1: m.l, x2: W - m.r, y1: y, y2: y, class: "grid" }, svg);
    el("text", { x: m.l - 8, y: y + 4, "text-anchor": "end" }, svg).textContent = fmt(v);
  }
  cfg.xs.forEach((x, i) => {
    if (cfg.xTickEvery && i % cfg.xTickEvery) return;
    el("line", { x1: X(x), x2: X(x), y1: H - m.b, y2: H - m.b + 5, class: "axis" }, svg);
    el("text", { x: X(x), y: H - m.b + 18, "text-anchor": "middle" }, svg).textContent = cfg.xTick(x, i);
  });
  el("line", { x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: "axis" }, svg);
  el("line", { x1: m.l, x2: m.l, y1: m.t, y2: H - m.b, class: "axis" }, svg);
  el("text", { x: (m.l + W - m.r) / 2, y: H - 10, "text-anchor": "middle", class: "axis-title" }, svg).textContent = cfg.xLabel;
  el("text", { x: 14, y: (m.t + H - m.b) / 2, "text-anchor": "middle", class: "axis-title", transform: `rotate(-90 14 ${(m.t + H - m.b) / 2})` }, svg).textContent = "time per alignment";
  // Draw the highlighted series last, so it sits on top.
  const names = Object.keys(cfg.series).sort((a, b) => (a === cfg.top) - (b === cfg.top));
  for (const name of names) {
    if (cfg.hidden.has(name)) continue;
    const [token, width, dash] = cfg.styles[name], color = css(token);
    const pts = cfg.series[name].map((v, i) => v ? [X(cfg.xs[i]), Y(v[0]), v[1]] : null).filter(Boolean);
    el("polyline", { points: pts.map(p => p[0] + "," + p[1]).join(" "), fill: "none", stroke: color, "stroke-width": width,
      "stroke-dasharray": dash || "none", "stroke-linejoin": "round" }, svg);
    pts.forEach(p => el("circle", { cx: p[0], cy: p[1], r: name === cfg.top ? 3.6 : 2.8, fill: p[2] ? css("--paper") : color, stroke: color, "stroke-width": 1.5 }, svg));
  }
  // Hover: a crosshair at the nearest x and a tooltip with every shown tool's value, the fastest in bold.
  const box = svg.parentElement; let tip = box.querySelector(".tip");
  if (!tip) { tip = document.createElement("div"); tip.className = "tip"; tip.hidden = true; box.appendChild(tip); }
  const hover = el("line", { y1: m.t, y2: H - m.b, class: "hover-line", visibility: "hidden" }, svg);
  const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent" }, svg);
  function show(evt) {
    const r = svg.getBoundingClientRect(), sx = (evt.clientX - r.left) / r.width * W;
    let best = 0; cfg.xs.forEach((x, i) => { if (Math.abs(X(x) - sx) < Math.abs(X(cfg.xs[best]) - sx)) best = i; });
    const xp = X(cfg.xs[best]); hover.setAttribute("x1", xp); hover.setAttribute("x2", xp); hover.setAttribute("visibility", "visible");
    const rows = Object.keys(cfg.series).filter(n => !cfg.hidden.has(n) && cfg.series[n][best]).map(n => [n, cfg.series[n][best]]);
    rows.sort((a, b) => a[1][0] - b[1][0]);
    tip.replaceChildren();
    const head = document.createElement("div"); head.style.fontWeight = "700"; head.style.marginBottom = "0.2rem"; head.textContent = cfg.tooltipTitle(cfg.xs[best], best); tip.appendChild(head);
    rows.forEach(([n, v], k) => { const row = document.createElement("div"); row.className = "row" + (k === 0 ? " best" : "");
      const a = document.createElement("span"); a.textContent = n; a.style.color = css(cfg.styles[n][0]);
      const b = document.createElement("span"); b.textContent = fmt(v[0]) + (v[1] ? " (some)" : ""); row.append(a, b); tip.appendChild(row); });
    tip.hidden = false;
    const bx = box.getBoundingClientRect(), left = evt.clientX - bx.left + box.scrollLeft + 14;
    tip.style.left = Math.min(left, box.scrollWidth - tip.offsetWidth - 4) + "px"; tip.style.top = (evt.clientY - bx.top + 12) + "px";
  }
  hit.addEventListener("pointermove", show);
  hit.addEventListener("pointerleave", () => { tip.hidden = true; hover.setAttribute("visibility", "hidden"); });
}

function legend(container, names, styles, hidden, redraw) {
  container.replaceChildren();
  names.forEach(n => { const b = document.createElement("button"); b.type = "button"; b.setAttribute("aria-pressed", String(!hidden.has(n)));
    const s = document.createElement("span"); s.className = "swatch"; s.style.color = css(styles[n][0]);
    if (styles[n][2]) s.style.borderTopStyle = "dashed";
    b.append(s, n); b.addEventListener("click", () => { hidden.has(n) ? hidden.delete(n) : hidden.add(n); b.setAttribute("aria-pressed", String(!hidden.has(n))); redraw(); });
    container.appendChild(b); });
}

const STYLES = Object.fromEntries(TOOLS.map(([n, t, w, d]) => [n, [t, w, d]]));
const toPairs = s => Object.fromEntries(Object.entries(s).map(([k, v]) => [k, v.map(x => x == null ? null : Array.isArray(x) ? x : [x])]));
const divHidden = new Set(), lenHidden = new Set();
let lenRate = 15;
function drawAll() {
  lineChart(document.getElementById("chart-divergence"), { xs: DIVERGENCE.x, xLog: false, xLabel: "divergence (%)", xTickEvery: 1,
    xTick: x => x + "%", series: toPairs(DIVERGENCE.series), styles: STYLES, hidden: divHidden, yMin: 0.01, yMax: 1000, top: "dinara-align",
    tooltipTitle: x => "100 kbp at " + x + "%" });
  lineChart(document.getElementById("chart-length"), { xs: LENGTH_X, xLog: true, xLabel: "length (kbp), " + lenRate + "% divergence",
    xTick: x => x >= 1000 ? "1 Mbp" : x + " kbp", series: LENGTH[lenRate], styles: STYLES, hidden: lenHidden, yMin: 0.01, yMax: 100000, top: "dinara-align",
    tooltipTitle: x => (x >= 1000 ? "1 Mbp" : x + " kbp") + " at " + lenRate + "%" });
}
function drawLegends() {
  legend(document.getElementById("legend-divergence"), Object.keys(DIVERGENCE.series), STYLES, divHidden, drawAll);
  legend(document.getElementById("legend-length"), Object.keys(LENGTH[15]), STYLES, lenHidden, drawAll);
}
document.getElementById("len-5").addEventListener("click", () => { lenRate = 5; setRate(); });
document.getElementById("len-15").addEventListener("click", () => { lenRate = 15; setRate(); });
function setRate() { document.getElementById("len-5").setAttribute("aria-pressed", String(lenRate === 5));
  document.getElementById("len-15").setAttribute("aria-pressed", String(lenRate === 15)); drawAll(); }
drawLegends(); drawAll();
// Theme changes swap the token colours the charts read.
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { drawLegends(); drawAll(); });
new MutationObserver(() => { drawLegends(); drawAll(); }).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
