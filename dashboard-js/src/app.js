const RESULTS_ROOT = "/results/ANDV_trees_aln-hyphy";
const TABLE_ROOT = `${RESULTS_ROOT}/dashboard_tables`;
const METHODS = ["FEL", "MEME", "BUSTED", "aBSREL", "RELAX", "MSS", "CFEL", "GARD"];
const state = {
  tab: "overview",
  pThreshold: 0.1,
  selectedRun: null,
  search: "",
  treeSearch: "",
  treeFilter: "all",
  absrelShowLabels: false,
  showTestedCodons: true,
  showAlignmentQuality: true,
  showDomainTrack: true,
  treeViewBoxes: {},
};

const tableFiles = {
  analysis: "analysis_summary.tsv",
  qc: "qc_summary.tsv",
  dropped: "dropped_sequences.tsv",
  fel: "fel_sites.tsv",
  meme: "meme_sites.tsv",
  absrel: "absrel_branches.tsv",
  relax: "relax_results.tsv",
  mss: "mss_results.tsv",
  warnings: "warnings.tsv",
  quality: "alignment_quality_by_site.tsv",
};

const tableExports = new Map();

const app = document.querySelector("#app");
const threshold = document.querySelector("#p-threshold");
const thresholdLabel = document.querySelector("#p-threshold-label");

threshold.addEventListener("input", () => {
  state.pThreshold = Number(threshold.value);
  thresholdLabel.textContent = state.pThreshold.toFixed(2);
  render();
});

function parseTsv(text) {
  const lines = text.trim().split(/\r?\n/);
  if (!lines.length || !lines[0]) return [];
  const headers = lines.shift().split("\t");
  return lines.map(line => {
    const values = line.split("\t");
    return Object.fromEntries(headers.map((header, i) => [header, coerce(values[i] ?? "")]));
  });
}

function coerce(value) {
  if (value === "True") return true;
  if (value === "False") return false;
  if (value === "") return "";
  const number = Number(value);
  return Number.isFinite(number) && String(value).match(/^-?\d+(\.\d+)?(e-?\d+)?$/i) ? number : value;
}

async function loadTable(file) {
  const response = await fetch(`${TABLE_ROOT}/${file}`);
  if (!response.ok) throw new Error(`Could not load ${file}`);
  return parseTsv(await response.text());
}

async function loadData() {
  const entries = await Promise.all(
    Object.entries(tableFiles).map(async ([key, file]) => [key, await loadTable(file)])
  );
  const data = Object.fromEntries(entries);
  data.genes = buildGeneRows(data);
  data.trees = await loadTrees(data.genes);
  state.selectedRun = data.genes[0]?.gene_id ?? null;
  return data;
}

async function loadTrees(genes) {
  const entries = await Promise.all(genes.map(async gene => {
    const path = `${RESULTS_ROOT}/inputs/${gene.segment}_${gene.label_set}.hyphy_ready.treefile`;
    try {
      const response = await fetch(path);
      return [gene.gene_id, response.ok ? await response.text() : ""];
    } catch {
      return [gene.gene_id, ""];
    }
  }));
  return Object.fromEntries(entries);
}

function buildGeneRows(data) {
  const runs = new Map();
  for (const row of data.analysis) {
    const gene_id = `${row.segment}_${row.label_set}`;
    if (!runs.has(gene_id)) {
      runs.set(gene_id, {
        gene_id,
        segment: row.segment,
        label_set: row.label_set,
        annotation: `${row.segment} segment, ${row.label_set} labeled tree`,
        n_sequences: row.n_sequences || "",
        codons: row.codons || "",
        fel_sites: 0,
        meme_sites: 0,
        cfel_sites: 0,
        busted_q: "",
        absrel_branches: 0,
        relax_k: "",
        relax_p: "",
        mss_models: "",
        test_branches: "",
        background_branches: "",
        stop_codons: 0,
        warning_flags: "",
      });
    }
    const gene = runs.get(gene_id);
    if (row.method === "MEME") gene.meme_sites = Number(row.significant_count || 0);
    if (row.method === "FEL") gene.fel_sites = Number(row.significant_count || 0);
    if (row.method === "aBSREL") gene.absrel_branches = Number(row.significant_count || 0);
    if (row.method === "RELAX") {
      gene.relax_k = row.k;
      gene.relax_p = row.p_value;
    }
    if (row.method === "MSS") gene.mss_models = row.significant_count;
    if (!gene.n_sequences && row.n_sequences) gene.n_sequences = row.n_sequences;
    if (!gene.codons && row.codons) gene.codons = row.codons;
  }
  for (const row of data.relax) {
    const gene = runs.get(`${row.segment}_${row.label_set}`);
    if (gene) {
      gene.test_branches = row.test_branches;
      gene.background_branches = row.reference_branches;
    }
  }
  for (const row of data.qc) {
    const gene = runs.get(`${row.segment}_${row.label_set}`);
    if (gene) {
      gene.stop_codons = Number(row.dropped_sequences || 0);
      gene.qc_status = row.status;
    }
  }
  for (const gene of runs.values()) {
    const warnings = getWarnings(gene, data.warnings);
    gene.warning_flags = warnings.map(w => w.warning).join("; ");
    gene.evidence_tier = evidenceTier(gene, warnings);
  }
  return [...runs.values()];
}

function evidenceTier(gene, warnings) {
  const support = [
    gene.fel_sites > 0,
    gene.meme_sites > 0,
    gene.cfel_sites > 0,
    gene.absrel_branches > 0,
    Number(gene.relax_p) <= 0.05,
  ].filter(Boolean).length;
  const severeWarning = warnings.some(w => w.severity === "high");
  if (severeWarning && support > 0) return "Caution";
  if (support >= 4) return "Very high";
  if (support >= 2) return "High";
  if (support === 1) return "Moderate";
  return "Low";
}

function getWarnings(gene, warnings) {
  return warnings.filter(row => row.segment === gene.segment && row.label_set === gene.label_set);
}

function badge(value) {
  const cls = String(value).toLowerCase().replace(/\s+/g, "-");
  const severity = cls === "high" ? "high-warning" : cls;
  return `<span class="badge ${severity}">${value}</span>`;
}

function fmt(value, digits = 3) {
  if (value === "" || value === null || value === undefined) return "";
  if (typeof value === "number") return Number.isInteger(value) ? value : value.toPrecision(digits);
  return value;
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;");
}

function fileSafe(value) {
  return String(value).replace(/[^a-z0-9._-]+/gi, "_").replace(/^_+|_+$/g, "");
}

function serializeRows(rows, columns) {
  const keys = columns.map(col => col.key);
  const labels = columns.map(col => col.label);
  const body = rows.map(row => keys.map(key => String(row[key] ?? "").replaceAll("\t", " ").replaceAll("\n", " ")).join("\t"));
  return [labels.join("\t"), ...body].join("\n") + "\n";
}

function downloadBlob(text, filename, type) {
  const blob = new Blob([text], { type });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function downloadSvg(svgId, filename) {
  const svg = document.getElementById(svgId);
  if (!svg) return;
  const clone = svg.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.querySelectorAll("title").forEach(title => title.remove());
  const source = `<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(clone)}\n`;
  downloadBlob(source, filename, "image/svg+xml;charset=utf-8");
}

function exportToolbar({ tableId = "", svgId = "", filename = "selectionscope_export" } = {}) {
  const buttons = [];
  if (tableId) buttons.push(`<button class="export-button" data-export-table="${tableId}" data-filename="${fileSafe(filename)}.tsv">Download table TSV</button>`);
  if (svgId) buttons.push(`<button class="export-button" data-export-svg="${svgId}" data-filename="${fileSafe(filename)}.svg">Download figure SVG</button>`);
  return buttons.length ? `<div class="export-toolbar">${buttons.join("")}</div>` : "";
}

function table(rows, columns, exportName = "") {
  if (exportName) tableExports.set(exportName, { rows, columns });
  if (!rows.length) return `${exportToolbar({ tableId: exportName, filename: exportName })}<p class="muted">No rows to display.</p>`;
  const head = columns.map(col => `<th>${col.label}</th>`).join("");
  const body = rows.map(row => (
    `<tr>${columns.map(col => `<td>${col.render ? col.render(row[col.key], row) : escapeHtml(fmt(row[col.key]))}</td>`).join("")}</tr>`
  )).join("");
  return `${exportToolbar({ tableId: exportName, filename: exportName })}
    <div class="table-wrap"><table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>`;
}

function cards(items) {
  return `<div class="grid cards">${items.map(item => `
    <div class="card">
      <div class="label">${item.label}</div>
      <div class="value">${fmt(item.value)}</div>
    </div>
  `).join("")}</div>`;
}

function tabs() {
  const labels = [
    ["overview", "Overview"],
    ["qc", "Input QC"],
    ["genes", "Gene-Level"],
    ["sites", "Site-Level"],
    ["branches", "Branch-Level"],
    ["input-trees", "Input Trees"],
    ["browser", "Gene Browser"],
    ["warnings", "Warnings"],
    ["export", "Export"],
  ];
  return `<nav class="tabs">${labels.map(([id, label]) => (
    `<button class="tab ${state.tab === id ? "active" : ""}" data-tab="${id}">${label}</button>`
  )).join("")}</nav>`;
}

function statusHeatmap(data, id = "method-status-heatmap") {
  tableExports.set("method_status_heatmap", {
    rows: data.analysis,
    columns: [
      { key: "segment", label: "Segment" },
      { key: "label_set", label: "Tree" },
      { key: "method", label: "Method" },
      { key: "status", label: "Status" },
      { key: "significant_count", label: "Signal count" },
      { key: "p_value", label: "p" },
      { key: "k", label: "K" },
      { key: "interpretation", label: "Interpretation" },
    ],
  });
  const colors = {
    signal: "#2d6cdf",
    pass: "#2f8f5b",
    not_run: "#8792a2",
    missing: "#8792a2",
    fail: "#c84630",
    warning: "#c97924",
  };
  const width = 900;
  const rowHeight = 38;
  const height = 76 + data.genes.length * rowHeight;
  const methodWidth = 96;
  const cells = data.genes.map((gene, rowIndex) => {
    const y = 62 + rowIndex * rowHeight;
    const methodCells = METHODS.map((method, methodIndex) => {
      const row = data.analysis.find(item => item.segment === gene.segment && item.label_set === gene.label_set && item.method === method);
      const status = row?.status ?? "missing";
      const significant = Number(row?.significant_count || 0) > 0;
      const label = status === "pass" && significant ? "signal" : status;
      const x = 132 + methodIndex * methodWidth;
      return `<g>
        <rect x="${x}" y="${y}" width="86" height="26" rx="3" fill="${colors[label] || colors.missing}"/>
        <text x="${x + 43}" y="${y + 17}" text-anchor="middle" font-size="11" font-weight="700" fill="#ffffff">${label}</text>
      </g>`;
    }).join("");
    return `<g>
      <text x="24" y="${y + 17}" font-size="12" font-weight="700" fill="#17202a">${gene.gene_id}</text>
      ${methodCells}
    </g>`;
  }).join("");
  const header = METHODS.map((method, index) => (
    `<text x="${132 + index * methodWidth + 43}" y="52" text-anchor="middle" font-size="12" font-weight="700" fill="#374151">${method}</text>`
  )).join("");
  return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="24" y="24" font-size="15" font-weight="700" fill="#17202a">Method status and selection signal</text>
    <text x="24" y="52" font-size="12" font-weight="700" fill="#374151">Run</text>
    ${header}
    ${cells}
  </svg>`;
}

function barChart(rows, labelKey, valueKey, color = "#2d6cdf", id = "bar-chart", title = "Bar chart", yLabel = "Count") {
  const width = 760;
  const height = 300;
  const max = Math.max(1, ...rows.map(row => Number(row[valueKey] || 0)));
  const band = width / Math.max(1, rows.length);
  const ticks = [0, 0.25, 0.5, 0.75, 1].map(fraction => {
    const y = 235 - fraction * 185;
    const label = Math.round(max * fraction);
    return `<g>
      <line x1="58" y1="${y}" x2="720" y2="${y}" stroke="#e5eaf1" />
      <text x="50" y="${y + 4}" text-anchor="end" font-size="11" fill="#4b5563">${label}</text>
    </g>`;
  }).join("");
  const bars = rows.map((row, i) => {
    const value = Number(row[valueKey] || 0);
    const barHeight = (value / max) * 185;
    const x = i * band + 70;
    const y = 235 - barHeight;
    return `<g>
      <rect x="${x}" y="${y}" width="${Math.max(20, band - 42)}" height="${barHeight}" fill="${color}" rx="2"></rect>
      <text x="${x + Math.max(20, band - 42) / 2}" y="258" text-anchor="middle" font-size="11" fill="#17202a">${escapeHtml(row[labelKey])}</text>
      <text x="${x + Math.max(20, band - 42) / 2}" y="${Math.max(38, y - 7)}" text-anchor="middle" font-size="11" fill="#17202a">${value}</text>
    </g>`;
  }).join("");
  return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="60" y="24" font-size="15" font-weight="700" fill="#17202a">${escapeHtml(title)}</text>
    <text x="16" y="150" transform="rotate(-90 16 150)" text-anchor="middle" font-size="12" fill="#374151">${escapeHtml(yLabel)}</text>
    ${ticks}
    <line x1="58" y1="235" x2="720" y2="235" stroke="#17202a"/>
    <line x1="58" y1="50" x2="58" y2="235" stroke="#17202a"/>
    ${bars}
  </svg>`;
}

function scatterPlot(rows, { id, title, xKey, yKey, xLabel, yLabel, colorKey = null, sizeKey = null, showLabels = true }) {
  const width = 760;
  const height = 330;
  const plot = { left: 68, right: 720, top: 48, bottom: 255 };
  const xValues = rows.map(row => Number(row[xKey])).filter(Number.isFinite);
  const yValues = rows.map(row => Number(row[yKey])).filter(Number.isFinite);
  const xMin = Math.min(0, ...xValues);
  const xMax = Math.max(1, ...xValues);
  const yMin = Math.min(0, ...yValues);
  const yMax = Math.max(1, ...yValues);
  const xScale = value => plot.left + ((Number(value) - xMin) / Math.max(1e-9, xMax - xMin)) * (plot.right - plot.left);
  const yScale = value => plot.bottom - ((Number(value) - yMin) / Math.max(1e-9, yMax - yMin)) * (plot.bottom - plot.top);
  const xTicks = [0, 0.25, 0.5, 0.75, 1].map(fraction => xMin + fraction * (xMax - xMin));
  const yTicks = [0, 0.25, 0.5, 0.75, 1].map(fraction => yMin + fraction * (yMax - yMin));
  const grid = [
    ...xTicks.map(tick => `<g><line x1="${xScale(tick)}" y1="${plot.top}" x2="${xScale(tick)}" y2="${plot.bottom}" stroke="#e5eaf1"/><text x="${xScale(tick)}" y="${plot.bottom + 18}" text-anchor="middle" font-size="11" fill="#4b5563">${fmt(tick, 2)}</text></g>`),
    ...yTicks.map(tick => `<g><line x1="${plot.left}" y1="${yScale(tick)}" x2="${plot.right}" y2="${yScale(tick)}" stroke="#e5eaf1"/><text x="${plot.left - 10}" y="${yScale(tick) + 4}" text-anchor="end" font-size="11" fill="#4b5563">${fmt(tick, 2)}</text></g>`),
  ].join("");
  const points = rows.map(row => {
    const x = xScale(row[xKey]);
    const y = yScale(row[yKey]);
    const color = colorKey ? row[colorKey] : "#2d6cdf";
    const radius = sizeKey ? Math.max(5, Math.min(15, Math.sqrt(Number(row[sizeKey] || 1)) / 4)) : 7;
    const labelClass = showLabels ? "point-label always-visible" : "point-label hover-visible";
    return `<g class="scatter-point" tabindex="0">
      <circle cx="${x}" cy="${y}" r="${radius}" fill="${color}" stroke="#17202a" stroke-width="0.8" opacity="0.82">
        <title>${escapeHtml(row.label || "")}</title>
      </circle>
      <text class="${labelClass}" x="${x}" y="${y - radius - 6}" text-anchor="middle" font-size="11">${escapeHtml(row.label || "")}</text>
    </g>`;
  }).join("");
  return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="${plot.left}" y="24" font-size="15" font-weight="700" fill="#17202a">${escapeHtml(title)}</text>
    ${grid}
    <line x1="${plot.left}" y1="${plot.bottom}" x2="${plot.right}" y2="${plot.bottom}" stroke="#17202a"/>
    <line x1="${plot.left}" y1="${plot.top}" x2="${plot.left}" y2="${plot.bottom}" stroke="#17202a"/>
    <text x="${(plot.left + plot.right) / 2}" y="310" text-anchor="middle" font-size="12" fill="#374151">${escapeHtml(xLabel)}</text>
    <text x="18" y="${(plot.top + plot.bottom) / 2}" transform="rotate(-90 18 ${(plot.top + plot.bottom) / 2})" text-anchor="middle" font-size="12" fill="#374151">${escapeHtml(yLabel)}</text>
    ${points}
  </svg>`;
}

function getTicks(geneLength) {
  if (geneLength <= 300) return { major: 50, minor: 10 };
  if (geneLength <= 1000) return { major: 100, minor: 25 };
  return { major: 250, minor: 50 };
}

function tickValues(length, step) {
  const values = [1];
  for (let value = step; value < length; value += step) values.push(value);
  if (!values.includes(length)) values.push(length);
  return values;
}

function bubblePlot(data, id = "selection-burden-bubble") {
  const width = 900;
  const height = 300;
  const maxMeme = Math.max(1, ...data.genes.map(g => g.meme_sites));
  const maxAbsrel = Math.max(1, ...data.genes.map(g => g.absrel_branches));
  const palette = ["#138a8a", "#2d6cdf", "#c97924", "#7c3aed", "#2f8f5b", "#c84630", "#64748b", "#0f766e", "#b45309"];
  const points = data.genes.map((gene, index) => {
    const x = 70 + (gene.absrel_branches / maxAbsrel) * 560;
    const y = 230 - (gene.meme_sites / maxMeme) * 180;
    const r = 8 + Math.min(20, Math.sqrt(gene.absrel_branches) * 5);
    const color = palette[index % palette.length];
    return `<g class="bubble" data-run="${gene.gene_id}">
      <circle cx="${x}" cy="${y}" r="${r}" fill="${color}" stroke="#17202a" stroke-width="1.1" opacity="0.82">
        <title>${escapeHtml(gene.gene_id)}; MEME ${gene.meme_sites}; aBSREL ${gene.absrel_branches}</title>
      </circle>
      <text x="${x}" y="${y + 4}" text-anchor="middle" font-size="11" font-weight="700" fill="#ffffff">${index + 1}</text>
    </g>`;
  }).join("");
  const legend = data.genes.map((gene, index) => {
    const y = 54 + index * 22;
    const color = palette[index % palette.length];
    return `<g>
      <circle cx="705" cy="${y - 4}" r="6" fill="${color}" stroke="#17202a" stroke-width="0.8"/>
      <text x="718" y="${y}" font-size="11" fill="#17202a">${index + 1}. ${escapeHtml(gene.gene_id)}</text>
    </g>`;
  }).join("");
  return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="60" y="25" font-size="15" font-weight="700" fill="#17202a">Selection burden by run</text>
    <line x1="60" y1="50" x2="650" y2="50" stroke="#e5eaf1"/>
    <line x1="60" y1="110" x2="650" y2="110" stroke="#e5eaf1"/>
    <line x1="60" y1="170" x2="650" y2="170" stroke="#e5eaf1"/>
    <line x1="60" y1="230" x2="650" y2="230" stroke="#8792a2"/>
    <line x1="60" y1="40" x2="60" y2="230" stroke="#8792a2"/>
    <text x="355" y="280" text-anchor="middle" font-size="12" fill="#374151">aBSREL selected branches</text>
    <text x="18" y="145" transform="rotate(-90 18 145)" text-anchor="middle" font-size="12" fill="#374151">MEME candidate sites</text>
    <text x="700" y="30" font-size="12" font-weight="700" fill="#374151">Runs</text>
    ${points}
    ${legend}
  </svg>`;
}

function lollipop(data, gene, { method = "MEME", id = `${method.toLowerCase()}-lollipop-${fileSafe(gene.gene_id)}` } = {}) {
  const isFel = method === "FEL";
  const sites = (isFel ? data.fel : data.meme)
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id && Number(row.p_value) <= state.pThreshold)
    .sort((a, b) => a.codon - b.codon);
  const width = 900;
  const height = 360;
  const codons = Number(gene.codons || 1);
  const quality = data.quality.filter(row => `${row.segment}_${row.label_set}` === gene.gene_id);
  const { major, minor } = getTicks(codons);
  const x = codon => 70 + ((Number(codon) - 1) / Math.max(1, codons - 1)) * 780;
  const plotTop = 76;
  const plotBottom = 205;
  const plotHeight = plotBottom - plotTop;
  const maxSignal = Math.max(3, ...sites.map(site => Number(site.neg_log10_p || 0)).filter(Number.isFinite));
  const yMax = Math.ceil(maxSignal * 1.18 * 10) / 10;
  const zeroY = isFel ? plotTop + plotHeight / 2 : plotBottom;
  const y = value => {
    const numeric = Number(value || 0);
    if (isFel) {
      const clamped = Math.max(-yMax, Math.min(yMax, numeric));
      return zeroY - (clamped / yMax) * (plotHeight / 2);
    }
    return plotBottom - (Math.min(yMax, Math.max(0, numeric)) / yMax) * plotHeight;
  };
  const signedFelSignal = site => {
    if (!isFel) return Number(site.neg_log10_p || 0);
    if (site.direction === "purifying") return -Number(site.neg_log10_p || 0);
    if (site.direction === "diversifying") return Number(site.neg_log10_p || 0);
    return 0;
  };
  const effectKey = isFel ? "omega" : "omega_plus";
  const maxEffect = Math.max(1, ...sites.map(site => Number(site[effectKey] || 0)));
  const siteColor = site => {
    if (!isFel) return "#2d6cdf";
    if (site.direction === "diversifying") return "#2d6cdf";
    if (site.direction === "purifying") return "#2f8f5b";
    return "#8792a2";
  };
  const siteStroke = site => {
    if (!isFel) return "#174ea6";
    if (site.direction === "diversifying") return "#174ea6";
    if (site.direction === "purifying") return "#1d6b42";
    return "#64748b";
  };
  const thresholdSpecs = [
    { p: 0.1, label: "p = 0.10" },
    { p: 0.05, label: "p = 0.05" },
    { p: 0.01, label: "p = 0.01" },
  ];
  const thresholdLines = (isFel
    ? thresholdSpecs.flatMap(item => [
      { ...item, value: -Math.log10(item.p), label: `+ ${item.label}` },
      { ...item, value: Math.log10(item.p), label: `- ${item.label}` },
    ])
    : thresholdSpecs.map(item => ({ ...item, value: -Math.log10(item.p) }))
  ).map(item => {
    const yy = y(item.value);
    return `<g>
      <line x1="70" y1="${yy}" x2="850" y2="${yy}" stroke="#c97924" stroke-width="1" stroke-dasharray="5 4"/>
      <text x="855" y="${yy + 4}" font-size="11" fill="#8a4b12">${item.label}</text>
    </g>`;
  }).join("");
  const minorTicks = tickValues(codons, minor).map(tick => `<line x1="${x(tick)}" y1="${plotBottom}" x2="${x(tick)}" y2="${plotBottom + 6}" stroke="#a8b2c1"/>`).join("");
  const majorTicks = tickValues(codons, major).map(tick => `<g>
    <line x1="${x(tick)}" y1="${plotBottom}" x2="${x(tick)}" y2="${plotBottom + 11}" stroke="#17202a"/>
    <text x="${x(tick)}" y="232" text-anchor="middle" font-size="11" fill="#17202a">${tick}</text>
  </g>`).join("");
  const yTickValues = isFel
    ? [-yMax, -yMax * 0.5, 0, yMax * 0.5, yMax]
    : [0, yMax * 0.25, yMax * 0.5, yMax * 0.75, yMax];
  const yTicks = Array.from(new Set(yTickValues)).map(tick => `<g>
    <line x1="64" y1="${y(tick)}" x2="70" y2="${y(tick)}" stroke="#17202a"/>
    <text x="58" y="${y(tick) + 4}" text-anchor="end" font-size="11" fill="#17202a">${fmt(tick, 2)}</text>
  </g>`).join("");
  const testedRug = state.showTestedCodons
    ? tickValues(codons, Math.max(1, Math.floor(codons / 160))).map(tick => `<line class="tested-codon-rug" x1="${x(tick)}" y1="250" x2="${x(tick)}" y2="258" stroke="#9ca3af" stroke-width="0.7"/>`).join("")
    : "";
  const bins = 120;
  const qualityBins = state.showAlignmentQuality ? Array.from({ length: bins }, (_, index) => {
    const start = Math.floor(index * codons / bins) + 1;
    const end = Math.floor((index + 1) * codons / bins);
    const rows = quality.filter(row => row.codon >= start && row.codon <= end);
    const gap = rows.length ? rows.reduce((sum, row) => sum + Number(row.gap_fraction || 0), 0) / rows.length : 0;
    const entropy = rows.length ? rows.reduce((sum, row) => sum + Number(row.entropy || 0), 0) / rows.length : 0;
    const intensity = Math.max(gap, Math.min(1, entropy / 2));
    const color = intensity > 0.2 ? "#c84630" : intensity > 0.05 ? "#c97924" : "#2f8f5b";
    return `<rect x="${x(start)}" y="284" width="${Math.max(1, x(end) - x(start))}" height="14" fill="${color}" opacity="${0.25 + intensity * 0.65}">
      <title>Codons ${start}-${end}; mean gap ${fmt(gap, 2)}; mean entropy ${fmt(entropy, 2)}</title>
    </rect>`;
  }).join("") : "";
  const topLabels = [...sites].sort((a, b) => Number(b.neg_log10_p || 0) - Number(a.neg_log10_p || 0)).slice(0, 5);
  const sticks = sites.map(site => {
    const xx = x(site.codon);
    const signal = signedFelSignal(site);
    const yy = y(signal);
    const radius = 4 + Math.min(8, (Number(site[effectKey] || 0) / maxEffect) * 8);
    const borderline = Number(site.p_value) > 0.05;
    const shouldLabel = topLabels.includes(site);
    const tooltip = isFel
      ? `Codon ${site.codon}; ${site.direction}; p=${fmt(site.p_value, 3)}; omega ${fmt(site.omega, 3)}; alpha ${fmt(site.alpha, 3)}; beta ${fmt(site.beta, 3)}`
      : `Codon ${site.codon}; p=${fmt(site.p_value, 3)}; omega+ ${fmt(site.omega_plus, 3)}; branches ${fmt(site.branches_under_selection)}`;
    return `<g>
      <line x1="${xx}" y1="${zeroY}" x2="${xx}" y2="${yy}" stroke="${siteColor(site)}" stroke-width="1.5" opacity="${borderline ? 0.55 : 1}" />
      <circle cx="${xx}" cy="${yy}" r="${radius}" fill="${siteColor(site)}" stroke="${siteStroke(site)}" stroke-width="1.2" opacity="${borderline ? 0.62 : 0.95}">
        <title>${escapeHtml(tooltip)}</title>
      </circle>
      ${shouldLabel ? `<text x="${xx}" y="${signal < 0 ? Math.min(plotBottom - 4, yy + radius + 13) : Math.max(58, yy - radius - 8)}" text-anchor="middle" font-size="10" fill="#17202a">${site.codon}</text>` : ""}
    </g>`;
  }).join("");
  const legend = isFel ? `<g>
    <circle cx="650" cy="30" r="5" fill="#2d6cdf" stroke="#174ea6"/>
    <text x="660" y="34" font-size="11" fill="#374151">Diversifying</text>
    <circle cx="740" cy="30" r="5" fill="#2f8f5b" stroke="#1d6b42"/>
    <text x="750" y="34" font-size="11" fill="#374151">Purifying</text>
  </g>` : "";
  return `<svg id="${id}" class="lollipop publication-svg" viewBox="0 0 ${width} ${height}">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="44" y="30" font-size="15" font-weight="700" fill="#17202a">${escapeHtml(gene.gene_id)} ${method} site evidence</text>
    <text x="44" y="50" font-size="12" fill="#4b5563">Sites shown at p <= ${state.pThreshold.toFixed(2)}</text>
    ${isFel ? `<text x="260" y="50" font-size="12" fill="#4b5563">Positive = diversifying; negative = purifying</text>` : ""}
    ${legend}
    <rect x="70" y="${plotTop}" width="780" height="${plotHeight}" fill="#fbfcfe" stroke="#e5eaf1"/>
    ${thresholdLines}
    ${yTicks}
    <line x1="70" y1="${zeroY}" x2="850" y2="${zeroY}" stroke="#17202a" stroke-width="1.5"/>
    <line x1="70" y1="${plotTop}" x2="70" y2="${plotBottom}" stroke="#17202a" stroke-width="1.5"/>
    ${minorTicks}
    ${majorTicks}
    ${sticks}
    <text x="460" y="246" text-anchor="middle" font-size="12" fill="#374151">Codon position</text>
    <text x="18" y="135" transform="rotate(-90 18 135)" text-anchor="middle" font-size="12" fill="#374151">${isFel ? "signed -log10 p-value" : "-log10 p-value"}</text>
    ${state.showTestedCodons ? `<text class="tested-codon-rug" x="70" y="268" font-size="11" fill="#374151">All tested codons</text>` : ""}
    ${testedRug}
    ${state.showAlignmentQuality ? `<text x="70" y="279" font-size="11" fill="#374151">Alignment quality: green clean, orange/red higher gap or entropy</text>` : ""}
    ${qualityBins}
    ${state.showDomainTrack ? `<rect class="domain-track-placeholder" x="70" y="318" width="120" height="16" fill="#eef2f7" stroke="#cbd5e1"/>
    <text class="domain-track-placeholder" x="198" y="330" font-size="11" fill="#64748b">Domains/GARD tracks unavailable for this run</text>` : ""}
  </svg>`;
}

function parseNewick(text) {
  let index = 0;
  let nodeId = 0;
  function skipWhitespace() {
    while (/\s/.test(text[index] || "")) index += 1;
  }
  function readLabel() {
    skipWhitespace();
    let label = "";
    while (index < text.length && ![":", ",", ")", "(", ";"].includes(text[index])) {
      label += text[index];
      index += 1;
    }
    return label.trim();
  }
  function readLength() {
    skipWhitespace();
    if (text[index] !== ":") return 0;
    index += 1;
    let value = "";
    while (index < text.length && ![",", ")", ";"].includes(text[index])) {
      value += text[index];
      index += 1;
    }
    const number = Number(value);
    return Number.isFinite(number) ? number : 0;
  }
  function makeNode(rawLabel = "") {
    const isTest = rawLabel.includes("{Foreground}");
    const cleanName = rawLabel.replaceAll("{Foreground}", "") || `internal_${nodeId}`;
    return { id: `n${nodeId++}`, name: cleanName, rawLabel, isTest, length: 0, children: [] };
  }
  function parseNode() {
    skipWhitespace();
    let node;
    if (text[index] === "(") {
      index += 1;
      node = makeNode();
      while (index < text.length) {
        node.children.push(parseNode());
        skipWhitespace();
        if (text[index] === ",") {
          index += 1;
          continue;
        }
        if (text[index] === ")") {
          index += 1;
          break;
        }
      }
      const rawLabel = readLabel();
      if (rawLabel) {
        node.rawLabel = rawLabel;
        node.isTest = rawLabel.includes("{Foreground}");
        node.name = rawLabel.replaceAll("{Foreground}", "") || node.name;
      }
      node.length = readLength();
    } else {
      const rawLabel = readLabel();
      node = makeNode(rawLabel);
      node.length = readLength();
    }
    return node;
  }
  return parseNode();
}

function treeLayout(root) {
  const leaves = [];
  const lengths = [];
  function walk(node, depth = 0) {
    node.depth = depth;
    lengths.push(node.length || 0);
    if (!node.children.length) {
      node.y = leaves.length;
      leaves.push(node);
    } else {
      node.children.forEach(child => walk(child, depth + (child.length || 0)));
      node.y = node.children.reduce((sum, child) => sum + child.y, 0) / node.children.length;
    }
  }
  walk(root, 0);
  return { leaves, maxDepth: Math.max(...leaves.map(leaf => leaf.depth), 1), lengths };
}

function viewBoxString(box) {
  return `${box.x} ${box.y} ${box.width} ${box.height}`;
}

function parseViewBox(value) {
  const [x, y, width, height] = String(value).split(/\s+/).map(Number);
  return { x, y, width, height };
}

function defaultTreeViewBox(width, height) {
  return { x: 0, y: 0, width, height };
}

function niceScaleValue(maxDepth, targetFraction = 0.15) {
  const raw = Math.max(maxDepth * targetFraction, Number.EPSILON);
  const power = 10 ** Math.floor(Math.log10(raw));
  const normalized = raw / power;
  const nice = normalized >= 5 ? 5 : normalized >= 2 ? 2 : 1;
  return nice * power;
}

function treeScaleBar({ xStart, yStart, pixelsPerDepth, maxDepth }) {
  const value = niceScaleValue(maxDepth);
  const width = value * pixelsPerDepth;
  const label = value >= 0.1 ? fmt(value, 2) : fmt(value, 1);
  return `<g class="tree-scale-bar">
    <line x1="${xStart}" y1="${yStart}" x2="${xStart + width}" y2="${yStart}" stroke="#17202a" stroke-width="2"/>
    <line x1="${xStart}" y1="${yStart - 5}" x2="${xStart}" y2="${yStart + 5}" stroke="#17202a" stroke-width="1.4"/>
    <line x1="${xStart + width}" y1="${yStart - 5}" x2="${xStart + width}" y2="${yStart + 5}" stroke="#17202a" stroke-width="1.4"/>
    <text x="${xStart + width / 2}" y="${yStart + 18}" text-anchor="middle" font-size="11" font-weight="700" fill="#17202a">${label} substitutions/site</text>
  </g>`;
}

function treePanel(data, gene, id = `tree-${fileSafe(gene.gene_id)}`) {
  const treeText = data.trees[gene.gene_id];
  if (!treeText) return `<p class="muted">No tree file available for ${gene.gene_id}.</p>`;
  const root = parseNewick(treeText);
  const { leaves, maxDepth, lengths } = treeLayout(root);
  const selected = new Set(data.absrel.filter(row => `${row.segment}_${row.label_set}` === gene.gene_id && row.status === "selected").map(row => row.branch));
  const longThreshold = [...lengths].sort((a, b) => a - b)[Math.floor(lengths.length * 0.95)] || 1;
  const width = 960;
  const rowHeight = gene.segment === "S" ? 9 : 13;
  const bottomPadding = 72;
  const height = Math.max(520, 92 + leaves.length * rowHeight + bottomPadding);
  const viewportHeight = Math.min(height, 1120);
  const baseViewBox = defaultTreeViewBox(width, viewportHeight);
  const storedViewBox = state.treeViewBoxes[id];
  const viewBox = storedViewBox && storedViewBox.height <= viewportHeight * 1.5 ? storedViewBox : baseViewBox;
  const x = depth => 58 + (depth / maxDepth) * 650;
  const y = node => 58 + node.y * rowHeight;
  const scaleBar = treeScaleBar({ xStart: 58, yStart: viewportHeight - 44, pixelsPerDepth: 650 / maxDepth, maxDepth });
  const branches = [];
  function draw(node) {
    if (!node.children.length) return;
    const yy = y(node);
    const childYs = node.children.map(child => y(child));
    branches.push(`<line x1="${x(node.depth)}" y1="${Math.min(...childYs)}" x2="${x(node.depth)}" y2="${Math.max(...childYs)}" stroke="#94a3b8" stroke-width="1"/>`);
    for (const child of node.children) {
      const isSelected = selected.has(child.name);
      const isLong = (child.length || 0) >= longThreshold;
      const color = isSelected ? "#c84630" : child.isTest ? "#2d6cdf" : "#8792a2";
      const strokeWidth = isSelected ? 3.2 : child.isTest ? 2.1 : 1.2;
      branches.push(`<line x1="${x(node.depth)}" y1="${y(child)}" x2="${x(child.depth)}" y2="${y(child)}" stroke="${isLong ? "#c97924" : color}" stroke-width="${strokeWidth}" ${isLong ? 'stroke-dasharray="5 3"' : ""}>
        <title>${escapeHtml(child.name)}; ${child.isTest ? "RELAX Test" : "Background"}; length ${fmt(child.length, 3)}${isSelected ? "; aBSREL selected" : ""}</title>
      </line>`);
      draw(child);
    }
  }
  draw(root);
  const leafLabels = leaves.map(leaf => {
    const showLabel = selected.has(leaf.name) || leaf.isTest || leaves.length <= 90;
    return showLabel ? `<text x="${x(leaf.depth) + 5}" y="${y(leaf) + 3}" font-size="${leaves.length > 120 ? 7 : 9}" fill="#334155">${escapeHtml(leaf.name)}</text>` : "";
  }).join("");
  return `<div class="tree-controls" data-tree-controls="${id}">
    <button class="export-button" data-tree-zoom="${id}" data-zoom-factor="0.75">Zoom in</button>
    <button class="export-button" data-tree-zoom="${id}" data-zoom-factor="1.333333">Zoom out</button>
    <button class="export-button" data-tree-reset="${id}" data-base-viewbox="${viewBoxString(baseViewBox)}">Reset view</button>
    <span class="tree-help">Mouse wheel zooms. Drag the tree to pan.</span>
  </div>
  <div class="tree-scroll">
  <svg id="${id}" class="tree-svg linked-tree-svg publication-svg zoomable-tree" viewBox="${viewBoxString(viewBox)}" width="${width}" height="${viewportHeight}" data-base-viewbox="${viewBoxString(baseViewBox)}" data-zoomable-tree="true" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="24" y="24" font-size="15" font-weight="700" fill="#17202a">${escapeHtml(gene.gene_id)} RELAX/aBSREL tree view</text>
    <text x="24" y="43" font-size="11" fill="#4b5563">Blue = RELAX Test branches; gray = Background; red = aBSREL selected; dashed orange = long branch</text>
    ${branches.join("")}
    ${leafLabels}
    ${scaleBar}
  </svg>
  </div>`;
}

function walkTree(node, visitor, parent = null) {
  visitor(node, parent);
  node.children.forEach(child => walkTree(child, visitor, node));
}

function annotatedBranchRows(treeText) {
  if (!treeText) return [];
  const root = parseNewick(treeText);
  const rows = [];
  walkTree(root, (node, parent) => {
    if (!parent || !node.isTest) return;
    rows.push({
      branch: node.name,
      raw_label: node.rawLabel,
      branch_length: node.length,
      node_type: node.children.length ? "internal" : "terminal",
      descendants: countLeaves(node),
    });
  });
  return rows;
}

function totalBranchCount(treeText) {
  if (!treeText) return 0;
  let count = 0;
  walkTree(parseNewick(treeText), (node, parent) => {
    if (parent) count += 1;
  });
  return count;
}

function treeFilePath(gene) {
  return `${RESULTS_ROOT}/inputs/${gene.segment}_${gene.label_set}.hyphy_ready.treefile`;
}

function matchesTreeSearch(row, query) {
  if (!query) return true;
  const needle = query.toLowerCase();
  return [row.branch, row.raw_label, row.node_type].some(value => String(value ?? "").toLowerCase().includes(needle));
}

function filterAnnotatedBranches(rows) {
  const typeFiltered = rows.filter(row => {
    if (state.treeFilter === "terminal") return row.node_type === "terminal";
    if (state.treeFilter === "internal") return row.node_type === "internal";
    return true;
  });
  return typeFiltered.filter(row => matchesTreeSearch(row, state.treeSearch));
}

function countLeaves(node) {
  if (!node.children.length) return 1;
  return node.children.reduce((sum, child) => sum + countLeaves(child), 0);
}

function inputTreePanel(data, gene, visibleAnnotations, id = `input-tree-${fileSafe(gene.gene_id)}`) {
  const treeText = data.trees[gene.gene_id];
  if (!treeText) return `<p class="muted">No input tree file available for ${gene.gene_id}.</p>`;
  const root = parseNewick(treeText);
  const { leaves, maxDepth, lengths } = treeLayout(root);
  const longThreshold = [...lengths].sort((a, b) => a - b)[Math.floor(lengths.length * 0.95)] || 1;
  const query = state.treeSearch.trim().toLowerCase();
  const visibleBranches = new Set(visibleAnnotations.map(row => row.branch));
  const width = 980;
  const rowHeight = leaves.length > 140 ? 7 : leaves.length > 90 ? 9 : 12;
  const bottomPadding = 72;
  const height = Math.max(520, 96 + leaves.length * rowHeight + bottomPadding);
  const viewportHeight = Math.min(height, 960);
  const baseViewBox = defaultTreeViewBox(width, viewportHeight);
  const storedViewBox = state.treeViewBoxes[id];
  const viewBox = storedViewBox && storedViewBox.height <= viewportHeight * 1.5 ? storedViewBox : baseViewBox;
  const x = depth => 58 + (depth / maxDepth) * 670;
  const y = node => 62 + node.y * rowHeight;
  const scaleBar = treeScaleBar({ xStart: 58, yStart: viewportHeight - 44, pixelsPerDepth: 670 / maxDepth, maxDepth });
  const branches = [];
  function draw(node) {
    if (!node.children.length) return;
    const childYs = node.children.map(child => y(child));
    branches.push(`<line x1="${x(node.depth)}" y1="${Math.min(...childYs)}" x2="${x(node.depth)}" y2="${Math.max(...childYs)}" stroke="#cbd5e1" stroke-width="1"/>`);
    for (const child of node.children) {
      const isLong = (child.length || 0) >= longThreshold;
      const isVisibleAnnotation = child.isTest && visibleBranches.has(child.name);
      const isSearchHit = query && (child.name.toLowerCase().includes(query) || child.rawLabel.toLowerCase().includes(query));
      const color = isVisibleAnnotation ? "#2d6cdf" : child.isTest ? "#7aa7f7" : "#94a3b8";
      const strokeWidth = isVisibleAnnotation ? 3.2 : child.isTest ? 2.1 : 1.15;
      const opacity = child.isTest && !isVisibleAnnotation ? 0.38 : 1;
      branches.push(`<line x1="${x(node.depth)}" y1="${y(child)}" x2="${x(child.depth)}" y2="${y(child)}" stroke="${isLong && !child.isTest ? "#c97924" : color}" stroke-width="${strokeWidth}" opacity="${opacity}" ${isLong && !child.isTest ? 'stroke-dasharray="5 3"' : ""}>
        <title>${escapeHtml(child.name)}; ${child.isTest ? "annotated {Foreground}" : "unannotated"}; length ${fmt(child.length, 3)}</title>
      </line>`);
      if (child.isTest) {
        branches.push(`<circle cx="${x(child.depth)}" cy="${y(child)}" r="${isVisibleAnnotation ? 4 : 2.8}" fill="${isVisibleAnnotation ? "#2d6cdf" : "#7aa7f7"}" opacity="${opacity}" stroke="#ffffff" stroke-width="1">
          <title>${escapeHtml(child.rawLabel)} includes {Foreground}</title>
        </circle>`);
      }
      if (isSearchHit) {
        branches.push(`<circle cx="${x(child.depth)}" cy="${y(child)}" r="7" fill="none" stroke="#2f8f5b" stroke-width="2">
          <title>${escapeHtml(child.name)} matches search</title>
        </circle>`);
      }
      draw(child);
    }
  }
  draw(root);
  const leafLabels = leaves.map(leaf => {
    const showLabel = leaf.isTest || leaves.length <= 100;
    return showLabel ? `<text x="${x(leaf.depth) + 6}" y="${y(leaf) + 3}" font-size="${leaves.length > 120 ? 7 : 9}" fill="${leaf.isTest ? "#174ea6" : "#334155"}">${escapeHtml(leaf.name)}</text>` : "";
  }).join("");
  return `<div class="tree-controls" data-tree-controls="${id}">
    <button class="export-button" data-tree-zoom="${id}" data-zoom-factor="0.75">Zoom in</button>
    <button class="export-button" data-tree-zoom="${id}" data-zoom-factor="1.333333">Zoom out</button>
    <button class="export-button" data-tree-reset="${id}" data-base-viewbox="${viewBoxString(baseViewBox)}">Reset view</button>
    <span class="tree-help">Drag to pan. Wheel or buttons zoom. Blue branches are visible {Foreground} annotations.</span>
  </div>
  <div class="tree-scroll">
  <svg id="${id}" class="tree-svg input-tree-svg publication-svg zoomable-tree" viewBox="${viewBoxString(viewBox)}" width="${width}" height="${viewportHeight}" data-base-viewbox="${viewBoxString(baseViewBox)}" data-zoomable-tree="true" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="24" y="24" font-size="15" font-weight="700" fill="#17202a">${escapeHtml(gene.gene_id)} input Newick annotations</text>
    <text x="24" y="43" font-size="11" fill="#4b5563">Blue = {Foreground} annotation; gray = unannotated; dashed orange = long unannotated branch</text>
    ${branches.join("")}
    ${leafLabels}
    ${scaleBar}
  </svg>
  </div>`;
}

function annotationList(rows, gene) {
  if (!rows.length) return `<div class="empty-list">No annotated branches match the current filters.</div>`;
  const items = rows.map(row => `
    <article class="branch-row">
      <div class="branch-row-head">
        <strong class="branch-name">${escapeHtml(row.branch)}</strong>
        <span class="branch-type">${row.node_type === "terminal" ? "Tip" : "Internal"}</span>
      </div>
      <div class="branch-meta">
        <span>${fmt(row.descendants)} descendant ${Number(row.descendants) === 1 ? "tip" : "tips"}</span>
        <span>Length ${fmt(row.branch_length)}</span>
      </div>
      <code>${escapeHtml(row.raw_label)}</code>
    </article>
  `).join("");
  tableExports.set(`${gene.gene_id}_annotated_branches`, {
    rows,
    columns: [
      { key: "branch", label: "Branch" },
      { key: "node_type", label: "Type" },
      { key: "descendants", label: "Descendant tips" },
      { key: "branch_length", label: "Branch length" },
      { key: "raw_label", label: "Raw Newick label" },
    ],
  });
  return `<div class="branch-list">${items}</div>`;
}

function relaxEffectPlot(data) {
  const rows = data.relax.map(row => ({
    label: `${row.segment}_${row.label_set}`,
    log2k: Number(row.k) > 0 ? Math.log2(Number(row.k)) : 0,
    neglogp: Number(row.neg_log10_p || 0),
    color: Number(row.k) > 1 ? "#2d6cdf" : Number(row.k) < 1 ? "#c97924" : "#8792a2",
    codons: data.genes.find(gene => gene.gene_id === `${row.segment}_${row.label_set}`)?.codons || 1,
  }));
  return scatterPlot(rows, {
    id: "relax-effect-plot",
    title: "RELAX effect plot",
    xKey: "log2k",
    yKey: "neglogp",
    xLabel: "log2(K): relaxation left, intensification right",
    yLabel: "-log10 p-value",
    colorKey: "color",
    sizeKey: "codons",
  });
}

function relaxReliabilityPlot(data) {
  const rows = data.relax.map(row => ({
    label: `${row.segment}_${row.label_set}`,
    test: Number(row.test_branches || 0),
    neglogp: Number(row.neg_log10_p || 0),
    color: Number(row.k) > 1 ? "#2d6cdf" : Number(row.k) < 1 ? "#c97924" : "#8792a2",
    codons: data.genes.find(gene => gene.gene_id === `${row.segment}_${row.label_set}`)?.codons || 1,
  }));
  return scatterPlot(rows, {
    id: "relax-reliability-plot",
    title: "RELAX reliability by Test branch count",
    xKey: "test",
    yKey: "neglogp",
    xLabel: "Number of RELAX Test branches",
    yLabel: "-log10 p-value",
    colorKey: "color",
    sizeKey: "codons",
  });
}

function absrelEvidencePlot(data, gene = null, id = "absrel-evidence-plot") {
  const source = gene ? data.absrel.filter(row => `${row.segment}_${row.label_set}` === gene.gene_id) : data.absrel;
  const rows = source.map(row => ({
    label: `${row.segment}_${row.label_set} ${row.branch}`,
    branchLength: Number(row.branch_length || 0),
    neglogq: Number(row.neg_log10_corrected_p || 0),
    color: row.status === "selected" ? "#c84630" : "#c97924",
  }));
  return scatterPlot(rows, {
    id,
    title: gene ? `${gene.gene_id} aBSREL branch evidence` : "aBSREL branch evidence",
    xKey: "branchLength",
    yKey: "neglogq",
    xLabel: "Branch length",
    yLabel: "-log10 corrected p-value",
    colorKey: "color",
    showLabels: state.absrelShowLabels,
  });
}

function applyTreeView(svg, box) {
  svg.setAttribute("viewBox", viewBoxString(box));
  state.treeViewBoxes[svg.id] = box;
}

function zoomTree(svg, factor, clientX = null, clientY = null) {
  const box = parseViewBox(svg.getAttribute("viewBox"));
  const rect = svg.getBoundingClientRect();
  const px = clientX === null ? rect.width / 2 : clientX - rect.left;
  const py = clientY === null ? rect.height / 2 : clientY - rect.top;
  const anchorX = box.x + (px / rect.width) * box.width;
  const anchorY = box.y + (py / rect.height) * box.height;
  const next = {
    x: anchorX - (anchorX - box.x) * factor,
    y: anchorY - (anchorY - box.y) * factor,
    width: box.width * factor,
    height: box.height * factor,
  };
  applyTreeView(svg, next);
}

function attachTreeExploration() {
  app.querySelectorAll("svg[data-zoomable-tree]").forEach(svg => {
    svg.addEventListener("wheel", event => {
      event.preventDefault();
      zoomTree(svg, event.deltaY > 0 ? 1.15 : 0.87, event.clientX, event.clientY);
    }, { passive: false });

    let dragging = false;
    let last = null;
    svg.addEventListener("pointerdown", event => {
      dragging = true;
      last = { x: event.clientX, y: event.clientY };
      svg.setPointerCapture(event.pointerId);
      svg.classList.add("dragging");
    });
    svg.addEventListener("pointermove", event => {
      if (!dragging || !last) return;
      const box = parseViewBox(svg.getAttribute("viewBox"));
      const rect = svg.getBoundingClientRect();
      const dx = event.clientX - last.x;
      const dy = event.clientY - last.y;
      last = { x: event.clientX, y: event.clientY };
      applyTreeView(svg, {
        x: box.x - dx * (box.width / rect.width),
        y: box.y - dy * (box.height / rect.height),
        width: box.width,
        height: box.height,
      });
    });
    svg.addEventListener("pointerup", event => {
      dragging = false;
      last = null;
      svg.releasePointerCapture(event.pointerId);
      svg.classList.remove("dragging");
    });
    svg.addEventListener("pointerleave", () => {
      dragging = false;
      last = null;
      svg.classList.remove("dragging");
    });
  });

  app.querySelectorAll("[data-tree-zoom]").forEach(button => {
    button.addEventListener("click", () => {
      const svg = document.getElementById(button.dataset.treeZoom);
      if (svg) zoomTree(svg, Number(button.dataset.zoomFactor));
    });
  });
  app.querySelectorAll("[data-tree-reset]").forEach(button => {
    button.addEventListener("click", () => {
      const svg = document.getElementById(button.dataset.treeReset);
      if (!svg) return;
      const box = parseViewBox(button.dataset.baseViewbox);
      applyTreeView(svg, box);
    });
  });
}

function interpretation(data, gene) {
  const warnings = getWarnings(gene, data.warnings);
  const support = [];
  if (gene.meme_sites > 0) support.push(`${gene.meme_sites} MEME episodic candidate sites`);
  if (gene.fel_sites > 0) support.push(`${gene.fel_sites} FEL pervasive candidate sites`);
  if (gene.absrel_branches > 0) support.push(`${gene.absrel_branches} aBSREL selected branch`);
  if (Number(gene.relax_p) <= 0.05) support.push(`RELAX p=${fmt(gene.relax_p)}`);
  if (Number(gene.mss_models) > 0) support.push(`${gene.mss_models} MSS-GA models evaluated`);
  const relaxText = gene.relax_k === "" ? "" : `RELAX estimates K=${fmt(gene.relax_k)}, interpreted as ${Number(gene.relax_k) > 1 ? "intensified" : Number(gene.relax_k) < 1 ? "relaxed" : "neutral"} selection for this branch set.`;
  return `<div class="interpretation">
    <p>${gene.gene_id} is currently classified as <strong>${gene.evidence_tier}</strong>.</p>
    <p>${support.length ? `Supporting signals: ${support.join(", ")}.` : "No significant selection signal is detected in the methods currently present."}</p>
    <p>${relaxText}</p>
    <p>${warnings.length ? `${warnings.length} warning(s) should be reviewed before biological interpretation.` : "No dashboard warnings were generated for this run."}</p>
  </div>`;
}

function overview(data) {
  const completed = data.analysis.filter(row => row.status === "pass").length;
  const notRun = data.analysis.filter(row => row.status === "not_run").length;
  const memeSites = data.meme.filter(row => Number(row.p_value) <= state.pThreshold).length;
  const felSites = data.fel.filter(row => row.direction === "diversifying" && Number(row.p_value) <= state.pThreshold).length;
  const selectedBranches = data.absrel.filter(row => row.status === "selected").length;
  const mssModels = data.mss.reduce((sum, row) => sum + Number(row.model_count || 0), 0);
  const highWarnings = data.warnings.filter(row => row.severity === "high").length;
  return `<section class="page grid">
    ${cards([
      { label: "Runs", value: data.genes.length },
      { label: "Completed method runs", value: completed },
      { label: "Methods not run", value: notRun },
      { label: "FEL sites", value: felSites },
      { label: "MEME sites", value: memeSites },
      { label: "aBSREL branches", value: selectedBranches },
      { label: "MSS models", value: mssModels },
      { label: "High warnings", value: highWarnings },
    ])}
    <div class="grid two-col">
      <section class="panel">
        <h2>Method Status Heatmap</h2>
        ${exportToolbar({ tableId: "method_status_heatmap", svgId: "method-status-heatmap", filename: "method_status_heatmap" })}
        ${statusHeatmap(data)}
      </section>
      <section class="panel">
        <h2>Selection Burden Bubble Plot</h2>
        ${exportToolbar({ svgId: "selection-burden-bubble", filename: "selection_burden_bubble_plot" })}
        ${bubblePlot(data)}
      </section>
    </div>
  </section>`;
}

function qcPage(data) {
  return `<section class="page grid">
    <div class="grid two-col">
      <section class="panel">
        <h2>QC Table</h2>
        ${table(data.qc, [
          { key: "segment", label: "Segment" },
          { key: "label_set", label: "Tree" },
          { key: "kept_sequences", label: "Sequences" },
          { key: "codons", label: "Codons" },
          { key: "dropped_sequences", label: "Dropped" },
          { key: "status", label: "Status", render: value => badge(value) },
          { key: "notes", label: "Notes" },
        ], "qc_metrics_visible")}
      </section>
      <section class="panel">
        <h2>Branch Label Coverage</h2>
        ${exportToolbar({ svgId: "branch-label-coverage", filename: "branch_label_coverage" })}
        ${barChart(data.relax.map(row => ({ run: `${row.segment}_${row.label_set}`, test: row.test_branches })), "run", "test", "#138a8a", "branch-label-coverage", "RELAX test branch coverage", "Test branches")}
      </section>
    </div>
    <section class="panel">
      <h2>Dropped Sequences</h2>
      ${table(data.dropped, [
        { key: "segment", label: "Segment" },
        { key: "label_set", label: "Tree" },
        { key: "id", label: "Sequence" },
        { key: "reasons", label: "Reason" },
        { key: "stop_codons", label: "Stop codons" },
      ], "dropped_sequences_visible")}
    </section>
  </section>`;
}

function genesPage(data) {
  return `<section class="page grid">
    <section class="panel">
      <div class="toolbar">
        <input class="search" id="gene-search" placeholder="Search runs" value="${state.search}">
      </div>
      ${table(data.genes.filter(g => g.gene_id.toLowerCase().includes(state.search.toLowerCase())), [
        { key: "gene_id", label: "Run" },
        { key: "n_sequences", label: "Sequences" },
        { key: "codons", label: "Codons" },
        { key: "fel_sites", label: "FEL" },
        { key: "meme_sites", label: "MEME" },
        { key: "absrel_branches", label: "aBSREL" },
        { key: "relax_k", label: "RELAX K" },
        { key: "relax_p", label: "RELAX p" },
        { key: "mss_models", label: "MSS models" },
        { key: "evidence_tier", label: "Evidence", render: value => badge(value) },
      ], "gene_summary_visible")}
    </section>
  </section>`;
}

function sitesPage(data) {
  const gene = data.genes.find(row => row.gene_id === state.selectedRun) ?? data.genes[0];
  const sites = data.meme
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id && Number(row.p_value) <= state.pThreshold)
    .map(site => {
      const quality = data.quality.find(row => `${row.segment}_${row.label_set}` === gene.gene_id && row.codon === site.codon) || {};
      const warning = Number(quality.gap_fraction || 0) > 0.1 ? "gappy" : Number(quality.entropy || 0) > 1.5 ? "high entropy" : "";
      return { ...site, gap_fraction: quality.gap_fraction ?? "", entropy: quality.entropy ?? "", non_gap_sequences: quality.non_gap_sequences ?? "", warning };
    });
  const felSites = data.fel
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id && Number(row.p_value) <= state.pThreshold)
    .map(site => {
      const quality = data.quality.find(row => `${row.segment}_${row.label_set}` === gene.gene_id && row.codon === site.codon) || {};
      const warning = Number(quality.gap_fraction || 0) > 0.1 ? "gappy" : Number(quality.entropy || 0) > 1.5 ? "high entropy" : "";
      return { ...site, gap_fraction: quality.gap_fraction ?? "", entropy: quality.entropy ?? "", non_gap_sequences: quality.non_gap_sequences ?? "", warning };
    });
  return `<section class="page grid">
    <section class="panel">
      <h2>Integrated Site-Level Browser</h2>
      ${runPicker(data, gene)}
      <label class="checkbox-control">
        <input type="checkbox" id="quality-track-toggle" ${state.showAlignmentQuality ? "checked" : ""}>
        Show alignment quality track
      </label>
      <label class="checkbox-control">
        <input type="checkbox" id="tested-codons-toggle" ${state.showTestedCodons ? "checked" : ""}>
        Show all tested codons
      </label>
      <label class="checkbox-control">
        <input type="checkbox" id="domain-track-toggle" ${state.showDomainTrack ? "checked" : ""}>
        Show Domains/GARD track
      </label>
      <h3>MEME Episodic Site Evidence</h3>
      ${exportToolbar({ svgId: `meme-lollipop-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_meme_lollipop` })}
      ${lollipop(data, gene, { method: "MEME", id: `meme-lollipop-${fileSafe(gene.gene_id)}` })}
      <h3>FEL Pervasive Site Evidence</h3>
      ${exportToolbar({ svgId: `fel-lollipop-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_fel_lollipop` })}
      ${lollipop(data, gene, { method: "FEL", id: `fel-lollipop-${fileSafe(gene.gene_id)}` })}
    </section>
    <section class="panel">
      <h2>FEL Site Table</h2>
      ${table(felSites, [
        { key: "codon", label: "Codon" },
        { key: "p_value", label: "p" },
        { key: "alpha", label: "alpha" },
        { key: "beta", label: "beta" },
        { key: "omega", label: "omega" },
        { key: "direction", label: "Direction" },
        { key: "gap_fraction", label: "Gap fraction" },
        { key: "entropy", label: "Entropy" },
        { key: "non_gap_sequences", label: "Non-gap seqs" },
        { key: "warning", label: "Warning" },
        { key: "lrt", label: "LRT" },
        { key: "status", label: "Status", render: value => badge(value) },
      ], `${gene.gene_id}_fel_sites_visible`)}
    </section>
    <section class="panel">
      <h2>MEME Site Table</h2>
      ${table(sites, [
        { key: "codon", label: "Codon" },
        { key: "p_value", label: "p" },
        { key: "omega_plus", label: "omega+" },
        { key: "branch_fraction", label: "Branch fraction" },
        { key: "branches_under_selection", label: "Branches" },
        { key: "gap_fraction", label: "Gap fraction" },
        { key: "entropy", label: "Entropy" },
        { key: "non_gap_sequences", label: "Non-gap seqs" },
        { key: "warning", label: "Warning" },
        { key: "lrt", label: "LRT" },
      ], `${gene.gene_id}_meme_sites_visible`)}
    </section>
  </section>`;
}

function branchesPage(data) {
  const gene = data.genes.find(row => row.gene_id === state.selectedRun) ?? data.genes[0];
  return `<section class="page grid">
    <section class="panel">
      <h2>Linked RELAX/aBSREL Tree</h2>
      ${runPicker(data, gene)}
      ${exportToolbar({ svgId: `tree-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_relax_absrel_tree` })}
      ${treePanel(data, gene)}
    </section>
    <div class="grid two-col">
      <section class="panel">
        <h2>RELAX Effect Plot</h2>
        ${exportToolbar({ svgId: "relax-effect-plot", filename: "relax_effect_plot" })}
        ${relaxEffectPlot(data)}
      </section>
      <section class="panel">
        <h2>RELAX Reliability Plot</h2>
        ${exportToolbar({ svgId: "relax-reliability-plot", filename: "relax_reliability_plot" })}
        ${relaxReliabilityPlot(data)}
      </section>
    </div>
    <section class="panel">
      <h2>aBSREL Branch Evidence Plot</h2>
      <label class="checkbox-control">
        <input type="checkbox" id="absrel-label-toggle" ${state.absrelShowLabels ? "checked" : ""}>
        Show branch labels
      </label>
      ${exportToolbar({ svgId: "absrel-evidence-plot", filename: "absrel_branch_evidence_plot" })}
      ${absrelEvidencePlot(data)}
    </section>
    <section class="panel">
      <h2>aBSREL Branch-Level Selection</h2>
      ${table(data.absrel, [
        { key: "segment", label: "Segment" },
        { key: "label_set", label: "Tree" },
        { key: "branch", label: "Branch" },
        { key: "corrected_p_value", label: "Corrected p" },
        { key: "lrt", label: "LRT" },
        { key: "branch_length", label: "Branch length" },
        { key: "status", label: "Status", render: value => badge(value) },
      ], "absrel_branches_visible")}
    </section>
    <section class="panel">
      <h2>RELAX Context</h2>
      ${table(data.relax, [
        { key: "segment", label: "Segment" },
        { key: "label_set", label: "Tree" },
        { key: "test_branches", label: "Test branches" },
        { key: "reference_branches", label: "Reference branches" },
        { key: "k", label: "K" },
        { key: "p_value", label: "p" },
        { key: "interpretation", label: "Interpretation" },
      ], "relax_results_visible")}
      <h2>MSS Context</h2>
      ${table(data.mss, [
        { key: "segment", label: "Segment" },
        { key: "label_set", label: "Tree" },
        { key: "files", label: "Files" },
        { key: "model_count", label: "Models" },
        { key: "best_ic", label: "Best IC" },
        { key: "classes", label: "Classes" },
        { key: "interpretation", label: "Interpretation" },
      ], "mss_results_visible")}
    </section>
  </section>`;
}

function inputTreesPage(data) {
  const gene = data.genes.find(row => row.gene_id === state.selectedRun) ?? data.genes[0];
  const treeText = data.trees[gene.gene_id] || "";
  const annotated = annotatedBranchRows(treeText);
  const visibleAnnotated = filterAnnotatedBranches(annotated);
  const totalBranches = totalBranchCount(treeText);
  const terminalAnnotated = annotated.filter(row => row.node_type === "terminal").length;
  const internalAnnotated = annotated.filter(row => row.node_type === "internal").length;
  return `<section class="page grid">
    <section class="panel tree-explorer">
      <div class="tree-explorer-header">
        <div>
          <h2>Input Newick Explorer</h2>
          <p class="muted">Inspect the exact HyPhy-ready tree and the branches labeled with {Foreground}.</p>
        </div>
        ${runPicker(data, gene)}
      </div>
      <div class="tree-file-path">${escapeHtml(treeFilePath(gene))}</div>
      <div class="tree-explorer-toolbar">
        <label class="field-label" for="tree-search">Find annotated branch</label>
        <input class="search tree-search" id="tree-search" placeholder="Search branch labels" value="${escapeHtml(state.treeSearch)}">
        <div class="segmented" role="group" aria-label="Annotated branch filter">
          ${[
            ["all", "All"],
            ["terminal", "Tips"],
            ["internal", "Internal"],
          ].map(([value, label]) => `<button class="segment-button ${state.treeFilter === value ? "active" : ""}" data-tree-filter="${value}">${label}</button>`).join("")}
        </div>
      </div>
      ${cards([
        { label: "Branches", value: totalBranches },
        { label: "Annotated", value: annotated.length },
        { label: "Matching", value: visibleAnnotated.length },
        { label: "Annotated tips", value: terminalAnnotated },
        { label: "Annotated internal", value: internalAnnotated },
      ])}
    </section>
    <div class="grid tree-workspace">
      <section class="panel">
        <div class="panel-title-row">
          <h2>Annotated Newick Tree</h2>
          <div class="tree-legend" aria-label="Tree legend">
            <span><i class="legend-swatch annotated"></i>{Foreground}</span>
            <span><i class="legend-swatch search-hit"></i>Search hit</span>
            <span><i class="legend-swatch long-branch"></i>Long branch</span>
          </div>
        </div>
        ${exportToolbar({ svgId: `input-tree-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_input_tree_annotations` })}
        ${inputTreePanel(data, gene, visibleAnnotated)}
      </section>
      <aside class="panel annotation-sidebar">
        <div class="panel-title-row">
          <h2>Annotated Branches</h2>
          ${exportToolbar({ tableId: `${gene.gene_id}_annotated_branches`, filename: `${gene.gene_id}_annotated_branches` })}
        </div>
        ${annotationList(visibleAnnotated, gene)}
      </aside>
    </div>
    <section class="panel">
      <details class="raw-newick-details">
        <summary>Raw Newick</summary>
        <pre class="raw-newick">${escapeHtml(treeText)}</pre>
      </details>
    </section>
  </section>`;
}

function browserPage(data) {
  const gene = data.genes.find(row => row.gene_id === state.selectedRun) ?? data.genes[0];
  const warnings = getWarnings(gene, data.warnings);
  return `<section class="page grid two-col">
    <section class="panel">
      <h2>Integrated Gene Browser</h2>
      ${runPicker(data, gene)}
      ${cards([
        { label: "MEME sites", value: gene.meme_sites },
        { label: "FEL sites", value: gene.fel_sites },
        { label: "RELAX K", value: gene.relax_k || "NA" },
        { label: "MSS models", value: gene.mss_models || "NA" },
        { label: "aBSREL branches", value: gene.absrel_branches },
        { label: "Test branches", value: gene.test_branches || "NA" },
      ])}
      ${exportToolbar({ svgId: `meme-lollipop-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_integrated_lollipop` })}
      ${lollipop(data, gene)}
      ${exportToolbar({ svgId: `tree-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_integrated_tree` })}
      ${treePanel(data, gene)}
      ${exportToolbar({ svgId: `absrel-evidence-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_absrel_evidence` })}
      ${absrelEvidencePlot(data, gene, `absrel-evidence-${fileSafe(gene.gene_id)}`)}
      <h3>Evidence Summary</h3>
      ${table([gene], [
        { key: "gene_id", label: "Run" },
        { key: "fel_sites", label: "FEL" },
        { key: "meme_sites", label: "MEME" },
        { key: "absrel_branches", label: "aBSREL" },
        { key: "relax_k", label: "RELAX K" },
        { key: "relax_p", label: "RELAX p" },
        { key: "mss_models", label: "MSS models" },
        { key: "evidence_tier", label: "Tier", render: value => badge(value) },
      ], `${gene.gene_id}_evidence_summary`)}
    </section>
    <aside class="panel">
      <h2>Interpretation</h2>
      ${interpretation(data, gene)}
      <h2>Warnings</h2>
      ${warningList(warnings)}
    </aside>
  </section>`;
}

function runPicker(data, gene) {
  return `<div class="toolbar">
    <select class="select" id="run-picker">
      ${data.genes.map(row => `<option value="${row.gene_id}" ${row.gene_id === gene.gene_id ? "selected" : ""}>${row.gene_id}</option>`).join("")}
    </select>
  </div>`;
}

function warningList(rows) {
  if (!rows.length) return `<p class="muted">No warnings.</p>`;
  return `<div class="warning-list">${rows.map(row => `
    <div class="warning-item ${row.severity}">
      <strong>${row.method}: ${row.warning}</strong>
      <p>${row.suggested_action}</p>
    </div>
  `).join("")}</div>`;
}

function warningsPage(data) {
  return `<section class="page grid">
    <section class="panel">
      <h2>Warnings and Robustness</h2>
      ${table(data.warnings, [
        { key: "severity", label: "Severity", render: value => badge(value) },
        { key: "segment", label: "Segment" },
        { key: "label_set", label: "Tree" },
        { key: "method", label: "Method" },
        { key: "warning", label: "Warning" },
        { key: "suggested_action", label: "Suggested action" },
      ], "warnings_visible")}
    </section>
  </section>`;
}

function exportPage(data) {
  return `<section class="page grid">
    <section class="panel">
      <h2>Export</h2>
      <p class="muted">Download the normalized tables used by this dashboard.</p>
      <div class="toolbar">
        ${Object.values(tableFiles).map(file => `<a class="tab" href="${TABLE_ROOT}/${file}" download>${file}</a>`).join("")}
      </div>
      <h3>Draft Results Text</h3>
      <p class="interpretation">Across the current hantavirus analyses, FEL identified ${data.fel.filter(row => row.direction === "diversifying").length} candidate pervasive diversifying sites at p <= 0.1, MEME identified ${data.meme.length} candidate episodic sites at p <= 0.1, aBSREL identified ${data.absrel.filter(row => row.status === "selected").length} selected branches after correction, RELAX found ${data.relax.filter(row => row.significant === true).length} significant branch-set shifts at p <= 0.05, and MSS-GA evaluated ${data.mss.reduce((sum, row) => sum + Number(row.model_count || 0), 0)} synonymous-rate class models. Results should be interpreted with the duplicate-sequence and RELAX convergence warnings shown in the warning panel.</p>
    </section>
  </section>`;
}

let DATA = null;

function render() {
  if (!DATA) return;
  tableExports.clear();
  const pages = {
    overview,
    qc: qcPage,
    genes: genesPage,
    sites: sitesPage,
    branches: branchesPage,
    "input-trees": inputTreesPage,
    browser: browserPage,
    warnings: warningsPage,
    export: exportPage,
  };
  app.innerHTML = `${tabs()}${pages[state.tab](DATA)}`;
  app.querySelectorAll(".tab[data-tab]").forEach(button => {
    button.addEventListener("click", () => {
      state.tab = button.dataset.tab;
      render();
    });
  });
  app.querySelector("#run-picker")?.addEventListener("change", event => {
    state.selectedRun = event.target.value;
    render();
  });
  app.querySelector("#gene-search")?.addEventListener("input", event => {
    state.search = event.target.value;
    render();
  });
  app.querySelector("#tree-search")?.addEventListener("input", event => {
    state.treeSearch = event.target.value;
    render();
  });
  app.querySelectorAll("[data-tree-filter]").forEach(button => {
    button.addEventListener("click", () => {
      state.treeFilter = button.dataset.treeFilter;
      render();
    });
  });
  app.querySelector("#absrel-label-toggle")?.addEventListener("change", event => {
    state.absrelShowLabels = event.target.checked;
    render();
  });
  app.querySelector("#quality-track-toggle")?.addEventListener("change", event => {
    state.showAlignmentQuality = event.target.checked;
    render();
  });
  app.querySelector("#tested-codons-toggle")?.addEventListener("change", event => {
    state.showTestedCodons = event.target.checked;
    render();
  });
  app.querySelector("#domain-track-toggle")?.addEventListener("change", event => {
    state.showDomainTrack = event.target.checked;
    render();
  });
  app.querySelectorAll("[data-export-table]").forEach(button => {
    button.addEventListener("click", () => {
      const exportData = tableExports.get(button.dataset.exportTable);
      if (!exportData) return;
      downloadBlob(
        serializeRows(exportData.rows, exportData.columns),
        button.dataset.filename || `${button.dataset.exportTable}.tsv`,
        "text/tab-separated-values;charset=utf-8"
      );
    });
  });
  app.querySelectorAll("[data-export-svg]").forEach(button => {
    button.addEventListener("click", () => {
      downloadSvg(button.dataset.exportSvg, button.dataset.filename || `${button.dataset.exportSvg}.svg`);
    });
  });
  attachTreeExploration();
}

loadData()
  .then(data => {
    DATA = data;
    render();
  })
  .catch(error => {
    app.innerHTML = `<section class="loading">Could not load dashboard data: ${error.message}</section>`;
  });
