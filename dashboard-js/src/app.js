const RESULTS_ROOT = "/results/ANDV_trees_aln-hyphy";
const TABLE_ROOT = `${RESULTS_ROOT}/dashboard_tables`;
const METHODS = ["FEL", "MEME", "BUSTED", "aBSREL", "RELAX", "MSS", "CFEL", "GARD"];
const state = {
  tab: "overview",
  pThreshold: 0.1,
  siteMetric: "q",
  selectedRun: null,
  search: "",
  treeSearch: "",
  treeFilter: "all",
  absrelShowLabels: false,
  showTestedCodons: true,
  showAlignmentQuality: true,
  showDomainTrack: true,
  showRelaxNonSignificant: false,
  treeViewBoxes: {},
  tableSorts: {},
};

const tableFiles = {
  analysis: "analysis_summary.tsv",
  qc: "qc_summary.tsv",
  dropped: "dropped_sequences.tsv",
  fel: "fel_sites.tsv",
  meme: "meme_sites.tsv",
  memeBranchEbf: "meme_branch_ebf.tsv",
  absrel: "absrel_branches.tsv",
  relax: "relax_results.tsv",
  mss: "mss_results.tsv",
  mssModels: "mss_models.tsv",
  mssParameters: "mss_parameter_frequency.tsv",
  warnings: "warnings.tsv",
  quality: "alignment_quality_by_site.tsv",
};

const tableExports = new Map();

const CODON_TO_AA = {
  TTT: "Phe", TTC: "Phe", TTA: "Leu", TTG: "Leu",
  TCT: "Ser", TCC: "Ser", TCA: "Ser", TCG: "Ser",
  TAT: "Tyr", TAC: "Tyr", TAA: "Stop", TAG: "Stop",
  TGT: "Cys", TGC: "Cys", TGA: "Stop", TGG: "Trp",
  CTT: "Leu", CTC: "Leu", CTA: "Leu", CTG: "Leu",
  CCT: "Pro", CCC: "Pro", CCA: "Pro", CCG: "Pro",
  CAT: "His", CAC: "His", CAA: "Gln", CAG: "Gln",
  CGT: "Arg", CGC: "Arg", CGA: "Arg", CGG: "Arg",
  ATT: "Ile", ATC: "Ile", ATA: "Ile", ATG: "Met",
  ACT: "Thr", ACC: "Thr", ACA: "Thr", ACG: "Thr",
  AAT: "Asn", AAC: "Asn", AAA: "Lys", AAG: "Lys",
  AGT: "Ser", AGC: "Ser", AGA: "Arg", AGG: "Arg",
  GTT: "Val", GTC: "Val", GTA: "Val", GTG: "Val",
  GCT: "Ala", GCC: "Ala", GCA: "Ala", GCG: "Ala",
  GAT: "Asp", GAC: "Asp", GAA: "Glu", GAG: "Glu",
  GGT: "Gly", GGC: "Gly", GGA: "Gly", GGG: "Gly",
};

const AA_COLORS = {
  Ala: "#2d6cdf",
  Arg: "#7c3aed",
  Cys: "#c97924",
  Gly: "#138a8a",
  Ile: "#2f8f5b",
  Leu: "#64748b",
  Pro: "#b45309",
  Ser: "#c84630",
  Thr: "#0f766e",
  Val: "#5b7cfa",
};

const AA_FULL_NAMES = {
  Ala: "Alanine",
  Arg: "Arginine",
  Asn: "Asparagine",
  Asp: "Aspartic acid",
  Cys: "Cysteine",
  Gln: "Glutamine",
  Glu: "Glutamic acid",
  Gly: "Glycine",
  His: "Histidine",
  Ile: "Isoleucine",
  Leu: "Leucine",
  Lys: "Lysine",
  Met: "Methionine",
  Phe: "Phenylalanine",
  Pro: "Proline",
  Ser: "Serine",
  Thr: "Threonine",
  Trp: "Tryptophan",
  Tyr: "Tyrosine",
  Val: "Valine",
};

const AA_TO_CODONS = Object.entries(CODON_TO_AA).reduce((acc, [codon, aa]) => {
  if (aa === "Stop") return acc;
  if (!acc[aa]) acc[aa] = [];
  acc[aa].push(codon);
  return acc;
}, {});
Object.values(AA_TO_CODONS).forEach(codons => codons.sort());

const app = document.querySelector("#app");
const threshold = document.querySelector("#p-threshold");
const thresholdTitle = document.querySelector("#p-threshold-title");
const thresholdLabel = document.querySelector("#p-threshold-label");
const siteMetricButtons = document.querySelectorAll("[data-site-metric]");

function setLoadingMessage(message) {
  app.innerHTML = `<section class="loading">${escapeHtml(message)}</section>`;
}

function updateSiteMetricControls() {
  const metric = siteMetricMeta();
  thresholdTitle.textContent = `Site ${metric.longLabel} threshold`;
  thresholdLabel.textContent = state.pThreshold.toFixed(2);
  siteMetricButtons.forEach(button => {
    button.classList.toggle("active", button.dataset.siteMetric === state.siteMetric);
    button.setAttribute("aria-pressed", String(button.dataset.siteMetric === state.siteMetric));
  });
}

threshold.addEventListener("input", () => {
  state.pThreshold = Number(threshold.value);
  updateSiteMetricControls();
  render();
});

siteMetricButtons.forEach(button => {
  button.addEventListener("click", () => {
    state.siteMetric = button.dataset.siteMetric === "p" ? "p" : "q";
    updateSiteMetricControls();
    render();
  });
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
  if (value === "inf" || value === "Infinity") return Infinity;
  if (value === "-inf" || value === "-Infinity") return -Infinity;
  if (value === "") return "";
  const number = Number(value);
  return Number.isFinite(number) && String(value).match(/^-?\d+(\.\d+)?(e-?\d+)?$/i) ? number : value;
}

async function loadTable(file) {
  const text = await loadText(`${TABLE_ROOT}/${file}`, file);
  return parseTsv(text);
}

async function loadText(url, label = url) {
  if (typeof fetch === "function") {
    const response = await fetch(url);
    if (!response.ok) throw new Error(`Could not load ${label} (${response.status})`);
    return response.text();
  }
  return new Promise((resolve, reject) => {
    const request = new XMLHttpRequest();
    request.open("GET", url, true);
    request.onload = () => {
      if (request.status >= 200 && request.status < 300) {
        resolve(request.responseText);
      } else {
        reject(new Error(`Could not load ${label} (${request.status})`));
      }
    };
    request.onerror = () => reject(new Error(`Could not load ${label}`));
    request.send();
  });
}

async function loadData() {
  const entries = [];
  const files = Object.entries(tableFiles);
  for (const [index, [key, file]] of files.entries()) {
    setLoadingMessage(`Loading dashboard table ${index + 1} of ${files.length}: ${file}`);
    entries.push([key, await loadTable(file)]);
  }
  setLoadingMessage("Building dashboard views...");
  const data = Object.fromEntries(entries);
  data.genes = buildGeneRows(data);
  setLoadingMessage("Loading labeled input trees...");
  data.trees = await loadTrees(data.genes);
  state.selectedRun = data.genes[0]?.gene_id ?? null;
  return data;
}

async function loadTrees(genes) {
  const entries = await Promise.all(genes.map(async gene => {
    const path = `${RESULTS_ROOT}/inputs/${gene.segment}_${gene.label_set}.hyphy_ready.treefile`;
    try {
      return [gene.gene_id, await loadText(path, `${gene.gene_id} tree`)];
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
        relax_lrt: "",
        relax_p: "",
        relax_q: "",
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
      gene.relax_lrt = row.lrt;
      gene.relax_p = row.p_value;
      gene.relax_q = row.q_value;
    }
    if (row.method === "MSS") gene.mss_models = row.tested || row.significant_count;
    if (!gene.n_sequences && row.n_sequences) gene.n_sequences = row.n_sequences;
    if (!gene.codons && row.codons) gene.codons = row.codons;
  }
  for (const row of data.relax) {
    const gene = runs.get(`${row.segment}_${row.label_set}`);
    if (gene) {
      gene.test_branches = row.test_branches;
      gene.background_branches = row.reference_branches;
      gene.relax_lrt = row.lrt;
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
    Number(gene.relax_q || gene.relax_p) <= 0.05,
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

function lrtCell(value) {
  const number = Number(value);
  if (!Number.isFinite(number)) return escapeHtml(fmt(value));
  const cls = number < 0 ? "negative-lrt" : "positive-lrt";
  const label = number < 0 ? `${fmt(number, 4)} (negative)` : fmt(number, 4);
  return `<span class="lrt-cell ${cls}">${escapeHtml(label)}</span>`;
}

function passesRelaxContextThreshold(row) {
  const q = Number(row.q_value);
  const p = Number(row.p_value);
  return (Number.isFinite(q) && q <= 0.1) || (Number.isFinite(p) && p <= 0.1);
}

function felSelectionInterpretation(row) {
  if (row.selection_interpretation) return row.selection_interpretation;
  if (row.direction === "diversifying") return "Positive selection (dN > dS)";
  if (row.direction === "purifying") return "Negative selection (dN < dS)";
  if (row.direction === "neutral") return "Neutral (dN = dS)";
  return "Unknown";
}

function felSelectionClass(row) {
  const call = row.selection_call || row.direction;
  if (call === "positive_selection" || call === "diversifying") return "positive-selection";
  if (call === "negative_selection" || call === "purifying") return "negative-selection";
  if (call === "neutral") return "neutral-selection";
  return "unknown-selection";
}

function felSummaryRows(rows, totalSites = "") {
  const categories = [
    {
      direction: "purifying",
      label: "Purifying selection",
      interpretation: "Negative selection: dN < dS",
    },
    {
      direction: "diversifying",
      label: "Positive selection",
      interpretation: "Diversifying selection: dN > dS",
    },
    {
      direction: "neutral",
      label: "Neutral",
      interpretation: "No rate difference: dN = dS",
    },
    {
      direction: "unknown",
      label: "Unknown",
      interpretation: "FEL rate estimates could not be classified",
    },
  ];
  return categories
    .map(category => {
      const categoryRows = rows.filter(row => row.direction === category.direction);
      return {
        selection_type: category.label,
        total_sites: totalSites,
        uncorrected_p_lt_0_1: categoryRows.filter(row => Number(row.p_value) < 0.1).length,
        bh_q_lte_0_1: categoryRows.filter(row => Number(row.q_value) <= 0.1).length,
        strongest_p: minFinite(categoryRows.map(row => row.p_value)),
        strongest_q: minFinite(categoryRows.map(row => row.q_value)),
        interpretation: category.interpretation,
      };
    })
    .filter(row => row.uncorrected_p_lt_0_1 > 0 || row.bh_q_lte_0_1 > 0);
}

function minFinite(values) {
  const numeric = values.map(Number).filter(Number.isFinite);
  return numeric.length ? Math.min(...numeric) : "";
}

function fmt(value, digits = 3) {
  if (value === "" || value === null || value === undefined) return "";
  if (typeof value === "number") return Number.isInteger(value) ? value : value.toPrecision(digits);
  return value;
}

function siteMetricMeta() {
  return state.siteMetric === "p"
    ? { valueKey: "p_value", signalKey: "neg_log10_p", shortLabel: "p", longLabel: "raw p-value" }
    : { valueKey: "q_value", signalKey: "neg_log10_q", shortLabel: "q", longLabel: "BH FDR q-value" };
}

function siteQValue(row) {
  const q = Number(row.q_value);
  if (Number.isFinite(q)) return q;
  const p = Number(row.p_value);
  return Number.isFinite(p) ? p : Infinity;
}

function siteMetricValue(row) {
  const metric = siteMetricMeta();
  const value = Number(row[metric.valueKey]);
  if (Number.isFinite(value)) return value;
  return state.siteMetric === "q" ? siteQValue(row) : Infinity;
}

function siteSignal(row) {
  const metric = siteMetricMeta();
  const signal = Number(row[metric.signalKey]);
  if (Number.isFinite(signal)) return signal;
  if (state.siteMetric === "q") {
    const pSignal = Number(row.neg_log10_p);
    return Number.isFinite(pSignal) ? pSignal : 0;
  }
  return 0;
}

function passesSiteThreshold(row) {
  return siteMetricValue(row) <= state.pThreshold;
}

function siteValueColumns() {
  return state.siteMetric === "p"
    ? [
      { key: "p_value", label: "p" },
      { key: "q_value", label: "BH q" },
    ]
    : [
      { key: "q_value", label: "BH q" },
      { key: "p_value", label: "p" },
    ];
}

function finiteNumber(value, fallback = 0) {
  const number = Number(value);
  return Number.isFinite(number) ? number : fallback;
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

function figureSvgSource(svgId) {
  const svg = document.getElementById(svgId);
  if (!svg) return null;
  const clone = svg.cloneNode(true);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
  style.textContent = `
    text { font-family: Inter, Arial, Helvetica, sans-serif; letter-spacing: 0; }
    line, path, rect, circle { vector-effect: non-scaling-stroke; }
    .point-label { fill: #0f172a; font-weight: 750; paint-order: stroke; stroke: #ffffff; stroke-linejoin: round; stroke-width: 4px; }
    .hover-visible { opacity: 1; }
    .tested-codon-rug, .domain-track-placeholder { display: inline; }
  `;
  clone.insertBefore(style, clone.firstChild);
  clone.querySelectorAll("title").forEach(title => title.remove());
  return {
    source: `<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(clone)}\n`,
    viewBox: clone.getAttribute("viewBox"),
    width: Number(clone.getAttribute("width")) || 0,
    height: Number(clone.getAttribute("height")) || 0,
  };
}

function downloadSvg(svgId, filename) {
  const figure = figureSvgSource(svgId);
  if (!figure) return;
  const { source } = figure;
  downloadBlob(source, filename, "image/svg+xml;charset=utf-8");
}

function figureDimensions(figure) {
  if (figure.viewBox) {
    const values = figure.viewBox.split(/\s+/).map(Number);
    if (values.length === 4 && values.every(Number.isFinite)) {
      return { width: values[2], height: values[3] };
    }
  }
  return {
    width: figure.width || 1200,
    height: figure.height || 800,
  };
}

function downloadPng(svgId, filename, scale = 4) {
  const figure = figureSvgSource(svgId);
  if (!figure) return;
  const { width, height } = figureDimensions(figure);
  const image = new Image();
  const svgBlob = new Blob([figure.source], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(svgBlob);
  image.onload = () => {
    const canvas = document.createElement("canvas");
    canvas.width = Math.ceil(width * scale);
    canvas.height = Math.ceil(height * scale);
    const context = canvas.getContext("2d");
    if (!context) {
      URL.revokeObjectURL(url);
      return;
    }
    context.fillStyle = "#ffffff";
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.drawImage(image, 0, 0, canvas.width, canvas.height);
    URL.revokeObjectURL(url);
    canvas.toBlob(blob => {
      if (!blob) return;
      const pngUrl = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = pngUrl;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(pngUrl);
    }, "image/png");
  };
  image.onerror = () => URL.revokeObjectURL(url);
  image.src = url;
}

function exportToolbar({ tableId = "", svgId = "", filename = "selectionscope_export" } = {}) {
  const buttons = [];
  if (tableId) buttons.push(`<button class="export-button" data-export-table="${tableId}" data-filename="${fileSafe(filename)}.tsv">Download table TSV</button>`);
  if (svgId) {
    buttons.push(`<button class="export-button" data-export-svg="${svgId}" data-filename="${fileSafe(filename)}.svg">Download figure SVG</button>`);
    buttons.push(`<button class="export-button" data-export-png="${svgId}" data-filename="${fileSafe(filename)}_4x.png">Download hi-res PNG</button>`);
  }
  return buttons.length ? `<div class="export-toolbar">${buttons.join("")}</div>` : "";
}

function compareTableValues(a, b) {
  const aNumber = Number(a);
  const bNumber = Number(b);
  const aNumeric = a !== "" && a !== null && a !== undefined && Number.isFinite(aNumber);
  const bNumeric = b !== "" && b !== null && b !== undefined && Number.isFinite(bNumber);
  if (aNumeric && bNumeric) return aNumber - bNumber;
  return String(a ?? "").localeCompare(String(b ?? ""), undefined, { numeric: true, sensitivity: "base" });
}

function sortedRows(rows, sortSpec) {
  if (!sortSpec?.key) return rows;
  const direction = sortSpec.direction === "desc" ? -1 : 1;
  return [...rows].sort((a, b) => {
    const result = compareTableValues(a[sortSpec.key], b[sortSpec.key]);
    return result * direction;
  });
}

function table(rows, columns, exportName = "", options = {}) {
  const sortable = Boolean(exportName && options.sortable !== false);
  const sortSpec = sortable ? state.tableSorts[exportName] : null;
  const displayRows = sortable ? sortedRows(rows, sortSpec) : rows;
  if (exportName) tableExports.set(exportName, { rows: displayRows, columns });
  if (!displayRows.length) return `${exportToolbar({ tableId: exportName, filename: exportName })}<p class="muted">No rows to display.</p>`;
  const head = columns.map(col => {
    if (!sortable) return `<th>${escapeHtml(col.label)}</th>`;
    const isActive = sortSpec?.key === col.key;
    const direction = isActive ? sortSpec.direction : "";
    const ariaSort = isActive ? (direction === "desc" ? "descending" : "ascending") : "none";
    const stateClass = isActive ? ` is-active is-${direction}` : "";
    const stateText = isActive ? ` sorted ${ariaSort}` : "";
    return `<th aria-sort="${ariaSort}"><button class="sortable-heading${stateClass}" data-sort-table="${exportName}" data-sort-key="${col.key}" aria-label="Sort by ${escapeHtml(col.label)}${stateText}"><span class="sort-label">${escapeHtml(col.label)}</span><span class="sort-indicator" aria-hidden="true"></span></button></th>`;
  }).join("");
  const body = displayRows.map(row => (
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
    ["genes", "Gene-Level Summary"],
    ["sites", "Site-Level"],
    ["branches", "Branch-Level"],
    ["mss", "MSS Charts"],
    ["input-trees", "Input Trees"],
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
    const rawLabel = max * fraction;
    const label = max <= 1.25 ? rawLabel.toFixed(2) : Math.round(rawLabel);
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
    const barColor = row.color || color;
    return `<g>
      <rect x="${x}" y="${y}" width="${Math.max(20, band - 42)}" height="${barHeight}" fill="${barColor}" rx="2"></rect>
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

function branchLabelCoverageChart(rows, id = "branch-label-coverage") {
  const width = 760;
  const rowHeight = 26;
  const height = 112 + rows.length * rowHeight;
  const plot = { left: 112, right: 710, top: 56, bottom: height - 44 };
  const max = Math.max(1, ...rows.map(row => Number(row.test_branches || 0) + Number(row.reference_branches || 0)));
  const x = value => plot.left + (Number(value || 0) / max) * (plot.right - plot.left);
  const ticks = [0, 0.25, 0.5, 0.75, 1].map(fraction => {
    const value = Math.round(max * fraction);
    const xx = x(value);
    return `<g>
      <line x1="${xx}" y1="${plot.top - 6}" x2="${xx}" y2="${plot.bottom}" stroke="#e5eaf1"/>
      <text x="${xx}" y="${plot.bottom + 18}" text-anchor="middle" font-size="10.5" fill="#64748b">${value}</text>
    </g>`;
  }).join("");
  const bars = rows.map((row, index) => {
    const y = plot.top + index * rowHeight;
    const test = Number(row.test_branches || 0);
    const reference = Number(row.reference_branches || 0);
    const testWidth = Math.max(0, x(test) - plot.left);
    const referenceWidth = Math.max(0, x(test + reference) - x(test));
    const label = `${row.segment}-${row.label_set}`;
    return `<g>
      <text x="${plot.left - 12}" y="${y + 13}" text-anchor="end" font-size="11" fill="#334155">${escapeHtml(label)}</text>
      <rect x="${plot.left}" y="${y + 3}" width="${testWidth}" height="14" rx="2" fill="#138a8a">
        <title>${escapeHtml(label)}; Test branches ${test}; Reference branches ${reference}</title>
      </rect>
      <rect x="${x(test)}" y="${y + 3}" width="${referenceWidth}" height="14" rx="2" fill="#d9e2ec">
        <title>${escapeHtml(label)}; Test branches ${test}; Reference branches ${reference}</title>
      </rect>
      <text x="${Math.min(plot.right + 6, x(test + reference) + 6)}" y="${y + 14}" font-size="10.5" fill="#475569">${test}/${reference}</text>
    </g>`;
  }).join("");
  return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="${plot.left}" y="24" font-size="15" font-weight="700" fill="#17202a">RELAX branch label coverage</text>
    <g class="chart-legend">
      <rect x="${plot.left}" y="34" width="10" height="10" rx="2" fill="#138a8a"/>
      <text x="${plot.left + 15}" y="43" font-size="11" fill="#475569">Test</text>
      <rect x="${plot.left + 68}" y="34" width="10" height="10" rx="2" fill="#d9e2ec"/>
      <text x="${plot.left + 83}" y="43" font-size="11" fill="#475569">Reference</text>
      <text x="${plot.right - 80}" y="43" font-size="11" fill="#475569">labels: Test/Reference</text>
    </g>
    ${ticks}
    <line x1="${plot.left}" y1="${plot.bottom}" x2="${plot.right}" y2="${plot.bottom}" stroke="#8792a2"/>
    ${bars}
  </svg>`;
}

function scatterPlot(rows, { id, title, xKey, yKey, xLabel, yLabel, colorKey = null, sizeKey = null, showLabels = true, className = "" }) {
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
  return `<svg id="${id}" class="chart ${className} publication-svg" viewBox="0 0 ${width} ${height}" role="img">
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

function linePlot(rows, { id, title, xKey, yKey, xLabel, yLabel, color = "#2d6cdf", labelKey = null }) {
  const width = 820;
  const height = 340;
  const plot = { left: 72, right: 780, top: 50, bottom: 265 };
  const xValues = rows.map(row => Number(row[xKey])).filter(Number.isFinite);
  const yValues = rows.map(row => Number(row[yKey])).filter(Number.isFinite);
  if (!xValues.length || !yValues.length) {
    return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}" role="img">
      <rect width="${width}" height="${height}" fill="#ffffff"/>
      <text x="${plot.left}" y="28" font-size="15" font-weight="700" fill="#17202a">${escapeHtml(title)}</text>
      <text x="${plot.left}" y="70" font-size="12" fill="#637083">No plottable rows.</text>
    </svg>`;
  }
  const xMin = Math.min(...xValues);
  const xMax = Math.max(...xValues);
  const yMin = Math.min(0, ...yValues);
  const yMax = Math.max(1, ...yValues);
  const xScale = value => plot.left + ((Number(value) - xMin) / Math.max(1e-9, xMax - xMin)) * (plot.right - plot.left);
  const yScale = value => plot.bottom - ((Number(value) - yMin) / Math.max(1e-9, yMax - yMin)) * (plot.bottom - plot.top);
  const ticks = [0, 0.25, 0.5, 0.75, 1];
  const grid = ticks.map(fraction => {
    const yTick = yMin + fraction * (yMax - yMin);
    const yy = yScale(yTick);
    return `<g>
      <line x1="${plot.left}" y1="${yy}" x2="${plot.right}" y2="${yy}" stroke="#e5eaf1"/>
      <text x="${plot.left - 10}" y="${yy + 4}" text-anchor="end" font-size="11" fill="#4b5563">${fmt(yTick, 2)}</text>
    </g>`;
  }).join("");
  const ordered = [...rows]
    .filter(row => Number.isFinite(Number(row[xKey])) && Number.isFinite(Number(row[yKey])))
    .sort((a, b) => Number(a[xKey]) - Number(b[xKey]));
  const path = ordered.map((row, index) => `${index ? "L" : "M"} ${xScale(row[xKey])} ${yScale(row[yKey])}`).join(" ");
  const points = ordered.filter((_, index) => index % Math.max(1, Math.floor(ordered.length / 80)) === 0).map(row => (
    `<circle cx="${xScale(row[xKey])}" cy="${yScale(row[yKey])}" r="3" fill="${color}">
      <title>${escapeHtml(labelKey ? row[labelKey] : `${xKey} ${fmt(row[xKey])}; ${yKey} ${fmt(row[yKey])}`)}</title>
    </circle>`
  )).join("");
  return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="${plot.left}" y="28" font-size="15" font-weight="700" fill="#17202a">${escapeHtml(title)}</text>
    ${grid}
    <line x1="${plot.left}" y1="${plot.bottom}" x2="${plot.right}" y2="${plot.bottom}" stroke="#17202a"/>
    <line x1="${plot.left}" y1="${plot.top}" x2="${plot.left}" y2="${plot.bottom}" stroke="#17202a"/>
    <path d="${path}" fill="none" stroke="${color}" stroke-width="2.2"/>
    ${points}
    <text x="${(plot.left + plot.right) / 2}" y="318" text-anchor="middle" font-size="12" fill="#374151">${escapeHtml(xLabel)}</text>
    <text x="18" y="${(plot.top + plot.bottom) / 2}" transform="rotate(-90 18 ${(plot.top + plot.bottom) / 2})" text-anchor="middle" font-size="12" fill="#374151">${escapeHtml(yLabel)}</text>
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

function parseMssParameter(parameter) {
  const match = String(parameter || "").match(/^alpha_([ACGT]{3})_([ACGT]{3})$/);
  if (!match) {
    return {
      codon1: "",
      codon2: "",
      codon_pair: String(parameter || "").replace("alpha_", ""),
      amino_acid: "Unknown",
      synonymous: false,
    };
  }
  const codon1 = match[1];
  const codon2 = match[2];
  const aa1 = CODON_TO_AA[codon1] || "Unknown";
  const aa2 = CODON_TO_AA[codon2] || "Unknown";
  return {
    codon1,
    codon2,
    codon_pair: `${codon1} <-> ${codon2}`,
    amino_acid: aa1 === aa2 ? aa1 : `${aa1}/${aa2}`,
    synonymous: aa1 === aa2 && aa1 !== "Unknown" && aa1 !== "Stop",
  };
}

function aminoAcidName(value) {
  return String(value || "")
    .split("/")
    .map(part => AA_FULL_NAMES[part] || part || "Unknown")
    .join("/");
}

function decorateMssParameter(row) {
  const parsed = parseMssParameter(row.parameter);
  return {
    ...row,
    ...parsed,
    amino_acid_name: aminoAcidName(parsed.amino_acid),
    parameter_label: String(row.parameter || "").replace("alpha_", ""),
    active_fraction: Number(row.active_fraction || 0),
  };
}

function mssRunSummaryPlot(data, id = "mss-run-summary-plot") {
  const rows = data.mss.map(row => ({
    label: `${row.segment}_${row.label_set}`,
    models: Number(row.model_count || 0),
    active: Number(row.median_active_parameters || 0),
    color: row.segment === "L" ? "#2d6cdf" : row.segment === "M" ? "#138a8a" : "#c97924",
  }));
  const width = 900;
  const height = 340;
  const plot = { left: 70, right: 640, top: 52, bottom: 262 };
  const xMax = Math.max(1, ...rows.map(row => row.active));
  const yMax = Math.max(1, ...rows.map(row => row.models));
  const xScale = value => plot.left + (Number(value) / xMax) * (plot.right - plot.left);
  const yScale = value => plot.bottom - (Number(value) / yMax) * (plot.bottom - plot.top);
  const ticks = [0, 0.25, 0.5, 0.75, 1];
  const grid = [
    ...ticks.map(fraction => {
      const xx = plot.left + fraction * (plot.right - plot.left);
      const label = Math.round(xMax * fraction);
      return `<g><line x1="${xx}" y1="${plot.top}" x2="${xx}" y2="${plot.bottom}" stroke="#e5eaf1"/><text x="${xx}" y="${plot.bottom + 18}" text-anchor="middle" font-size="11" fill="#4b5563">${label}</text></g>`;
    }),
    ...ticks.map(fraction => {
      const yy = plot.bottom - fraction * (plot.bottom - plot.top);
      const label = Math.round(yMax * fraction);
      return `<g><line x1="${plot.left}" y1="${yy}" x2="${plot.right}" y2="${yy}" stroke="#e5eaf1"/><text x="${plot.left - 10}" y="${yy + 4}" text-anchor="end" font-size="11" fill="#4b5563">${label}</text></g>`;
    }),
  ].join("");
  const points = rows.map((row, index) => `<g>
    <circle cx="${xScale(row.active)}" cy="${yScale(row.models)}" r="7" fill="${row.color}" stroke="#17202a" stroke-width="0.9">
      <title>${escapeHtml(row.label)}; median active ${row.active}; models ${row.models}</title>
    </circle>
    <text x="${xScale(row.active)}" y="${yScale(row.models) + 4}" text-anchor="middle" font-size="10" font-weight="800" fill="#ffffff">${index + 1}</text>
  </g>`).join("");
  const legend = rows.map((row, index) => `<g>
    <circle cx="690" cy="${62 + index * 22}" r="5.5" fill="${row.color}" stroke="#17202a" stroke-width="0.8"/>
    <text x="704" y="${66 + index * 22}" font-size="11" fill="#17202a">${index + 1}. ${escapeHtml(row.label)}</text>
  </g>`).join("");
  return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="${plot.left}" y="28" font-size="15" font-weight="700" fill="#17202a">MSS-GA search breadth by run</text>
    ${grid}
    <line x1="${plot.left}" y1="${plot.bottom}" x2="${plot.right}" y2="${plot.bottom}" stroke="#17202a"/>
    <line x1="${plot.left}" y1="${plot.top}" x2="${plot.left}" y2="${plot.bottom}" stroke="#17202a"/>
    ${points}
    <text x="${(plot.left + plot.right) / 2}" y="318" text-anchor="middle" font-size="12" fill="#374151">Median active synonymous parameters per model</text>
    <text x="18" y="${(plot.top + plot.bottom) / 2}" transform="rotate(-90 18 ${(plot.top + plot.bottom) / 2})" text-anchor="middle" font-size="12" fill="#374151">Models evaluated</text>
    <text x="686" y="36" font-size="12" font-weight="700" fill="#374151">Runs</text>
    ${legend}
  </svg>`;
}

function mssModelFrontierPlot(data, gene, id = `mss-frontier-${fileSafe(gene.gene_id)}`) {
  const rows = data.mssModels
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id)
    .filter(row => Number(row.rank) <= 250)
    .map(row => ({
      ...row,
      label: `Rank ${row.rank}; delta IC ${fmt(row.delta_ic)}; active ${row.active_parameters}`,
    }));
  return linePlot(rows, {
    id,
    title: `${gene.gene_id} MSS-GA model frontier`,
    xKey: "rank",
    yKey: "delta_ic",
    xLabel: "Model rank by information criterion",
    yLabel: "Delta IC from best model",
    color: "#138a8a",
    labelKey: "label",
  });
}

function mssComplexityPlot(data, gene, id = `mss-complexity-${fileSafe(gene.gene_id)}`) {
  const rows = data.mssModels
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id)
    .filter(row => Number(row.rank) <= 500)
    .map(row => ({
      label: `Rank ${row.rank}`,
      active: Number(row.active_parameters || 0),
      delta: Number(row.delta_ic || 0),
      color: Number(row.delta_ic || 0) <= 2 ? "#2f8f5b" : Number(row.delta_ic || 0) <= 10 ? "#c97924" : "#8792a2",
      rank: Number(row.rank || 0),
    }));
  return scatterPlot(rows, {
    id,
    title: `${gene.gene_id} MSS-GA complexity/evidence`,
    xKey: "active",
    yKey: "delta",
    xLabel: "Active synonymous parameters",
    yLabel: "Delta IC from best model",
    colorKey: "color",
    sizeKey: "rank",
    showLabels: false,
  });
}

function mssParameterPlot(data, gene, id = `mss-parameters-${fileSafe(gene.gene_id)}`) {
  const rows = data.mssParameters
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id)
    .sort((a, b) => Number(b.active_fraction || 0) - Number(a.active_fraction || 0))
    .slice(0, 14)
    .map(decorateMssParameter);
  const width = 820;
  const rowHeight = 25;
  const height = 88 + rows.length * rowHeight;
  const plot = { left: 178, right: 760, top: 54, bottom: height - 38 };
  const x = value => plot.left + Math.max(0, Math.min(1, Number(value || 0))) * (plot.right - plot.left);
  const ticks = [0, 0.25, 0.5, 0.75, 1].map(value => {
    const xx = x(value);
    return `<g>
      <line x1="${xx}" y1="${plot.top - 6}" x2="${xx}" y2="${plot.bottom}" stroke="#e5eaf1"/>
      <text x="${xx}" y="${plot.bottom + 18}" text-anchor="middle" font-size="10.5" fill="#64748b">${value.toFixed(2)}</text>
    </g>`;
  }).join("");
  const points = rows.map((row, index) => {
    const y = plot.top + index * rowHeight;
    const xx = x(row.active_fraction);
    const color = AA_COLORS[row.amino_acid] || "#2d6cdf";
    const label = `${row.codon_pair} (${row.amino_acid_name})`;
    return `<g>
      <text x="${plot.left - 12}" y="${y + 5}" text-anchor="end" font-size="11" fill="#334155">${escapeHtml(label)}</text>
      <line x1="${plot.left}" y1="${y}" x2="${xx}" y2="${y}" stroke="${color}" stroke-width="2.2"/>
      <circle cx="${xx}" cy="${y}" r="5" fill="${color}" stroke="#17202a" stroke-width="0.8">
        <title>${escapeHtml(row.parameter)}; ${escapeHtml(row.amino_acid_name)}; active fraction ${fmt(row.active_fraction, 3)}</title>
      </circle>
      <text x="${Math.min(plot.right + 6, xx + 8)}" y="${y + 4}" font-size="10.5" fill="#475569">${fmt(row.active_fraction, 3)}</text>
    </g>`;
  }).join("");
  return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="${plot.left}" y="24" font-size="15" font-weight="700" fill="#17202a">${gene.gene_id} top MSS codon-pair parameters</text>
    <text x="${plot.left}" y="40" font-size="11" fill="#64748b">Each point is a synonymous codon-pair parameter; x-axis is active fraction across fitted MSS-GA models.</text>
    ${ticks}
    <line x1="${plot.left}" y1="${plot.bottom}" x2="${plot.right}" y2="${plot.bottom}" stroke="#8792a2"/>
    ${points}
  </svg>`;
}

function mssHeatmapRows(data, gene, limit = 40) {
  return data.mssParameters
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id)
    .map(decorateMssParameter)
    .filter(row => row.synonymous)
    .sort((a, b) => Number(b.active_fraction || 0) - Number(a.active_fraction || 0))
    .slice(0, limit);
}

function mssHeatmapInterpretation(data, gene) {
  const rows = mssHeatmapRows(data, gene, 40);
  if (!rows.length) {
    return `<p class="interpretation">No synonymous codon-pair MSS parameters are available for this run.</p>`;
  }
  const top = rows[0];
  const byAa = new Map();
  for (const row of rows) {
    if (!byAa.has(row.amino_acid)) byAa.set(row.amino_acid, []);
    byAa.get(row.amino_acid).push(row);
  }
  const familySummaries = [...byAa.entries()]
    .sort((a, b) => Math.max(...b[1].map(row => Number(row.active_fraction || 0))) - Math.max(...a[1].map(row => Number(row.active_fraction || 0))))
    .slice(0, 3)
    .map(([aa, familyRows]) => {
      const best = familyRows[0];
      return `${aminoAcidName(aa)}: ${best.codon_pair} (${fmt(best.active_fraction, 3)})`;
    });
  const recurringFamilies = [...byAa.entries()]
    .filter(([, familyRows]) => familyRows.length >= 2)
    .map(([aa]) => aa);
  return `<div class="mss-interpretation">
    <p><strong>Main signal:</strong> ${escapeHtml(gene.gene_id)} is led by ${escapeHtml(top.codon_pair)} in the ${escapeHtml(top.amino_acid_name)} family, active in ${fmt(top.active_fraction, 3)} of MSS-GA models.</p>
    <p><strong>Family concentration:</strong> top synonymous-pair support spans ${byAa.size} amino-acid families. Strongest families are ${escapeHtml(familySummaries.join("; "))}.</p>
    <p><strong>How to read the heatmaps:</strong> darker upper-triangle cells mark codon pairs that recur more often across fitted MSS models. Lower triangles are hidden because ${escapeHtml("ACA-ACG")} and ${escapeHtml("ACG-ACA")} represent the same symmetric pair.</p>
    ${recurringFamilies.length ? `<p><strong>Repeated within-family structure:</strong> ${escapeHtml(recurringFamilies.join(", "))} each contain multiple supported codon pairs, suggesting distributed synonymous-rate heterogeneity within those families.</p>` : ""}
  </div>`;
}

function interpolateColor(start, end, t) {
  const parse = hex => [1, 3, 5].map(index => parseInt(hex.slice(index, index + 2), 16));
  const [r1, g1, b1] = parse(start);
  const [r2, g2, b2] = parse(end);
  const clamped = Math.max(0, Math.min(1, Number(t || 0)));
  const channel = (a, b) => Math.round(a + (b - a) * clamped);
  return `rgb(${channel(r1, r2)}, ${channel(g1, g2)}, ${channel(b1, b2)})`;
}

function mssFamilyHeatmaps(data, gene, id = `mss-family-heatmaps-${fileSafe(gene.gene_id)}`) {
  const rows = mssHeatmapRows(data, gene, 40);
  const byAa = new Map();
  rows.forEach(row => {
    if (!byAa.has(row.amino_acid)) byAa.set(row.amino_acid, []);
    byAa.get(row.amino_acid).push(row);
  });
  const families = [...byAa.entries()]
    .sort((a, b) => Math.max(...b[1].map(row => row.active_fraction)) - Math.max(...a[1].map(row => row.active_fraction)))
    .slice(0, 8);
  const width = 980;
  const panelWidth = 455;
  const panelHeight = 238;
  const columns = 2;
  const rowsOfPanels = Math.max(1, Math.ceil(families.length / columns));
  const height = 104 + rowsOfPanels * panelHeight;
  const scaleMin = 0.5;
  const scaleMax = 1;
  const color = value => {
    const t = (Number(value || 0) - scaleMin) / (scaleMax - scaleMin);
    return interpolateColor("#eaf4fb", "#084081", t);
  };
  const legendX = width - 300;
  const legendY = 31;
  const legendWidth = 190;
  const legendStops = Array.from({ length: 38 }, (_, index) => {
    const t = index / 37;
    return `<rect x="${legendX + t * legendWidth}" y="${legendY}" width="${legendWidth / 37 + 1}" height="10" fill="${interpolateColor("#eaf4fb", "#084081", t)}"/>`;
  }).join("");
  const panels = families.map(([aa, familyRows], panelIndex) => {
    const col = panelIndex % columns;
    const rowIndex = Math.floor(panelIndex / columns);
    const originX = 30 + col * panelWidth;
    const originY = 82 + rowIndex * panelHeight;
    const codons = AA_TO_CODONS[aa] || [...new Set(familyRows.flatMap(row => [row.codon1, row.codon2]))].sort();
    const cell = Math.min(34, Math.floor(186 / Math.max(1, codons.length)));
    const matrixX = originX + 96;
    const matrixY = originY + 46;
    const lookup = new Map();
    familyRows.forEach(row => {
      lookup.set(`${row.codon1}|${row.codon2}`, row.active_fraction);
      lookup.set(`${row.codon2}|${row.codon1}`, row.active_fraction);
    });
    const cells = codons.flatMap((codon1, yIndex) => codons.map((codon2, xIndex) => {
      if (xIndex <= yIndex) {
        return "";
      }
      const value = codon1 === codon2 ? "" : lookup.get(`${codon1}|${codon2}`);
      const fill = value === undefined || value === "" ? "#f7f9fb" : color(value);
      const textFill = Number(value || 0) >= 0.78 ? "#ffffff" : "#17202a";
      const tooltip = `${aa} ${codon1} <-> ${codon2}${value === undefined || value === "" ? "; no top parameter" : `; active fraction ${fmt(value, 3)}`}`;
      return `<g>
        <rect x="${matrixX + xIndex * cell}" y="${matrixY + yIndex * cell}" width="${cell - 1.5}" height="${cell - 1.5}" rx="2" fill="${fill}" stroke="#ffffff">
          <title>${escapeHtml(tooltip)}</title>
        </rect>
        ${value === undefined || value === "" || cell < 26 ? "" : `<text x="${matrixX + xIndex * cell + cell / 2}" y="${matrixY + yIndex * cell + cell / 2 + 3}" text-anchor="middle" font-size="8.8" font-weight="700" fill="${textFill}">${fmt(value, 2)}</text>`}
      </g>`;
    })).join("");
    const xLabels = codons.map((codon, index) => `<text x="${matrixX + index * cell + cell / 2}" y="${matrixY - 9}" text-anchor="middle" font-family="ui-monospace, SFMono-Regular, Menlo, monospace" font-size="10.5" fill="#334155">${codon}</text>`).join("");
    const yLabels = codons.map((codon, index) => `<text x="${matrixX - 10}" y="${matrixY + index * cell + cell / 2 + 4}" text-anchor="end" font-family="ui-monospace, SFMono-Regular, Menlo, monospace" font-size="10.5" fill="#334155">${codon}</text>`).join("");
    const familyColor = AA_COLORS[aa] || "#2d6cdf";
    const matrixWidth = codons.length * cell;
    const matrixHeight = codons.length * cell;
    return `<g>
      <circle cx="${originX + 18}" cy="${originY + 20}" r="5" fill="${familyColor}"/>
      <text x="${originX + 30}" y="${originY + 24}" font-size="14" font-weight="760" fill="#17202a">${escapeHtml(aminoAcidName(aa))} family</text>
      <text x="${matrixX + matrixWidth / 2}" y="${matrixY + matrixHeight + 24}" text-anchor="middle" font-size="10.5" fill="#64748b">codon 2</text>
      <text x="${matrixX - 58}" y="${matrixY + matrixHeight / 2}" text-anchor="middle" font-size="10.5" fill="#64748b" transform="rotate(-90 ${matrixX - 58} ${matrixY + matrixHeight / 2})">codon 1</text>
      ${xLabels}
      ${yLabels}
      ${cells}
    </g>`;
  }).join("");
  return `<svg id="${id}" class="chart publication-svg" viewBox="0 0 ${width} ${height}" role="img">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="30" y="28" font-size="17" font-weight="760" fill="#17202a">${gene.gene_id} synonymous-family MSS heatmaps</text>
    <text x="30" y="48" font-size="11.5" fill="#64748b">Upper triangle only; cell labels show active fraction for recurrent synonymous codon-pair parameters.</text>
    ${legendStops}
    <rect x="${legendX}" y="${legendY}" width="${legendWidth}" height="10" fill="none" stroke="#cbd5e1"/>
    <text x="${legendX}" y="${legendY + 25}" font-size="10.5" fill="#64748b">${scaleMin.toFixed(2)}</text>
    <text x="${legendX + legendWidth}" y="${legendY + 25}" text-anchor="end" font-size="10.5" fill="#64748b">${scaleMax.toFixed(2)}</text>
    <text x="${legendX + legendWidth + 12}" y="${legendY + 9}" font-size="10.5" fill="#64748b">active fraction</text>
    ${panels || `<text x="28" y="62" font-size="12" fill="#64748b">No synonymous codon-pair parameters available.</text>`}
  </svg>`;
}

function mssInterpretationPanel(data) {
  const rows = data.mss.filter(row => Number(row.model_count || 0) > 0);
  if (!rows.length) {
    return `<section class="panel">
      <h2>MSS Interpretation</h2>
      <p class="interpretation">No MSS-GA results are available in the current dashboard tables.</p>
    </section>`;
  }
  const bySegment = new Map();
  for (const row of rows) {
    if (!bySegment.has(row.segment)) bySegment.set(row.segment, []);
    bySegment.get(row.segment).push(row);
  }
  const segmentText = [...bySegment.entries()]
    .sort(([a], [b]) => a.localeCompare(b))
    .map(([segment, segmentRows]) => {
      const freqs = segmentRows.map(row => Number(row.top_parameter_frequency || 0)).filter(Number.isFinite);
      const active = segmentRows.map(row => Number(row.median_active_parameters || 0)).filter(Number.isFinite);
      const meanFreq = freqs.length ? freqs.reduce((sum, value) => sum + value, 0) / freqs.length : "";
      const meanActive = active.length ? active.reduce((sum, value) => sum + value, 0) / active.length : "";
      return `${segment}: mean top support ${fmt(meanFreq, 3)}, mean median active parameters ${fmt(meanActive, 3)}`;
    })
    .join("; ");
  const recurrent = new Map();
  for (const row of rows) {
    if (!row.top_parameter) continue;
    recurrent.set(row.top_parameter, (recurrent.get(row.top_parameter) || 0) + 1);
  }
  const recurrentText = [...recurrent.entries()]
    .filter(([, count]) => count > 1)
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
    .map(([parameter, count]) => `${parameter.replace("alpha_", "")} (${count} runs)`)
    .join(", ");
  return `<section class="panel">
    <h2>MSS Interpretation</h2>
    <div class="mss-interpretation">
      <p><strong>What MSS is saying:</strong> MSS-GA is a model-search summary for synonymous-rate heterogeneity. The key quantity shown here is model-inclusion frequency: how often a synonymous codon-pair parameter appears among fitted MSS models.</p>
      <p><strong>What it is not:</strong> these active fractions are not raw p-values, BH q-values, or formal multiple-testing discoveries. They should be interpreted as recurrent model support, not as site-level or branch-level significance calls.</p>
      <p><strong>Current narrative:</strong> MSS supports segment-level synonymous-rate heterogeneity across the ANDV analyses. ${escapeHtml(segmentText)}. Top parameter identity changes across tree label sets, so the robust story is broader synonymous-rate structure rather than one universal codon-pair effect.</p>
      ${recurrentText ? `<p><strong>Recurring top parameters:</strong> ${escapeHtml(recurrentText)} recur as the top parameter in more than one run; most other top parameters are segment/tree-specific.</p>` : ""}
      <p><strong>How to use this with selection results:</strong> treat MSS as context for synonymous-rate model structure. It can motivate cautious interpretation of dS-sensitive selection tests, but it does not itself identify adaptive amino-acid evolution.</p>
    </div>
  </section>`;
}

function mssPage(data) {
  const gene = data.genes.find(row => row.gene_id === state.selectedRun) ?? data.genes[0];
  const runMss = data.mss.filter(row => `${row.segment}_${row.label_set}` === gene.gene_id);
  const runModels = data.mssModels
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id)
    .filter(row => Number(row.rank) <= 100);
  const runParameters = data.mssParameters
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id)
    .sort((a, b) => Number(b.active_fraction || 0) - Number(a.active_fraction || 0))
    .slice(0, 30)
    .map(decorateMssParameter);
  return `<section class="page grid">
    <section class="panel">
      <div class="panel-title-row">
        <div>
          <h2>MSS-GA Synonymous-Rate Model Charts</h2>
          <p class="muted">Use these plots to separate model search breadth, near-optimal model complexity, and repeatedly selected synonymous substitution parameters.</p>
        </div>
        ${runPicker(data, gene)}
      </div>
    </section>
    ${mssInterpretationPanel(data)}
    <section class="panel">
      <h2>MSS Context</h2>
      ${table(data.mss, [
        { key: "segment", label: "Segment" },
        { key: "label_set", label: "Tree" },
        { key: "files", label: "Files" },
        { key: "model_count", label: "Models" },
        { key: "best_ic", label: "Best IC" },
        { key: "median_active_parameters", label: "Median active params" },
        { key: "top_parameter", label: "Top parameter" },
        { key: "top_parameter_frequency", label: "Top frequency" },
        { key: "interpretation", label: "Interpretation" },
      ], "mss_context_visible", { sortable: true })}
    </section>
    <div class="grid two-col">
      <section class="panel">
        <h2>Across-Run MSS Summary</h2>
        ${exportToolbar({ svgId: "mss-run-summary-plot", filename: "mss_run_summary_plot" })}
        ${mssRunSummaryPlot(data)}
      </section>
      <section class="panel">
        <h2>Model Frontier</h2>
        ${exportToolbar({ svgId: `mss-frontier-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_mss_model_frontier` })}
        ${mssModelFrontierPlot(data, gene)}
      </section>
    </div>
    <div class="grid two-col">
      <section class="panel">
        <h2>Complexity Versus Fit</h2>
        ${exportToolbar({ svgId: `mss-complexity-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_mss_complexity` })}
        ${mssComplexityPlot(data, gene)}
      </section>
      <section class="panel">
        <h2>Recurrent Synonymous Parameters</h2>
        ${exportToolbar({ svgId: `mss-parameters-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_mss_parameter_frequency` })}
        ${mssParameterPlot(data, gene)}
      </section>
    </div>
    <section class="panel">
      <h2>Synonymous Codon-Pair Heatmaps</h2>
      ${mssHeatmapInterpretation(data, gene)}
      ${exportToolbar({ svgId: `mss-family-heatmaps-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_mss_synonymous_family_heatmaps` })}
      ${mssFamilyHeatmaps(data, gene)}
    </section>
    <section class="panel">
      <h2>MSS Run Summary</h2>
      ${table(runMss, [
        { key: "segment", label: "Segment" },
        { key: "label_set", label: "Tree" },
        { key: "files", label: "Files" },
        { key: "model_count", label: "Models" },
        { key: "best_ic", label: "Best IC" },
        { key: "min_active_parameters", label: "Min active" },
        { key: "median_active_parameters", label: "Median active" },
        { key: "max_active_parameters", label: "Max active" },
        { key: "top_parameter", label: "Top parameter" },
        { key: "top_parameter_frequency", label: "Top frequency" },
      ], `${gene.gene_id}_mss_summary_visible`)}
    </section>
    <div class="grid two-col">
      <section class="panel">
        <h2>Top MSS Models</h2>
        ${table(runModels, [
          { key: "rank", label: "Rank" },
          { key: "ic", label: "IC" },
          { key: "delta_ic", label: "Delta IC" },
          { key: "log_likelihood", label: "Log likelihood" },
          { key: "active_parameters", label: "Active params" },
          { key: "class_count", label: "Classes" },
          { key: "class_rates", label: "Class rates" },
        ], `${gene.gene_id}_mss_top_models_visible`)}
      </section>
      <section class="panel">
        <h2>Top MSS Parameters</h2>
        ${table(runParameters, [
          { key: "parameter", label: "Parameter" },
          { key: "codon_pair", label: "Codon pair" },
          { key: "amino_acid_name", label: "Amino acid family" },
          { key: "active_models", label: "Active models" },
          { key: "model_count", label: "Models" },
          { key: "active_fraction", label: "Active fraction" },
        ], `${gene.gene_id}_mss_parameters_visible`)}
      </section>
    </div>
  </section>`;
}

function lollipop(data, gene, { method = "MEME", id = `${method.toLowerCase()}-lollipop-${fileSafe(gene.gene_id)}` } = {}) {
  const isFel = method === "FEL";
  const metric = siteMetricMeta();
  const sites = (isFel ? data.fel : data.meme)
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id && (isFel || (row.method || "MEME") === method) && passesSiteThreshold(row))
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
  const maxSignal = Math.max(3, ...sites.map(siteSignal).filter(Number.isFinite));
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
    if (!isFel) return siteSignal(site);
    if (site.direction === "purifying") return -siteSignal(site);
    if (site.direction === "diversifying") return siteSignal(site);
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
    { p: 0.1, label: `${metric.shortLabel} = 0.10` },
    { p: 0.05, label: `${metric.shortLabel} = 0.05` },
    { p: 0.01, label: `${metric.shortLabel} = 0.01` },
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
  const topLabels = [...sites].sort((a, b) => siteSignal(b) - siteSignal(a)).slice(0, 5);
  const sticks = sites.map(site => {
    const xx = x(site.codon);
    const signal = signedFelSignal(site);
    const yy = y(signal);
    const radius = 4 + Math.min(8, (Number(site[effectKey] || 0) / maxEffect) * 8);
    const borderline = siteMetricValue(site) > 0.05;
    const shouldLabel = topLabels.includes(site);
    const tooltip = isFel
      ? `Codon ${site.codon}; ${site.direction}; q=${fmt(site.q_value, 3)}; p=${fmt(site.p_value, 3)}; omega ${fmt(site.omega, 3)}; alpha ${fmt(site.alpha, 3)}; beta ${fmt(site.beta, 3)}`
      : `Codon ${site.codon}; q=${fmt(site.q_value, 3)}; p=${fmt(site.p_value, 3)}; omega+ ${fmt(site.omega_plus, 3)}; branches ${fmt(site.branches_under_selection)}`;
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
    <text x="44" y="50" font-size="12" fill="#4b5563">Sites shown at ${metric.longLabel} <= ${state.pThreshold.toFixed(2)}</text>
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
    <text x="18" y="135" transform="rotate(-90 18 135)" text-anchor="middle" font-size="12" fill="#374151">${isFel ? `signed -log10 ${metric.shortLabel}-value` : `-log10 ${metric.shortLabel}-value`}</text>
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
  const viewportHeight = height;
  const baseViewBox = defaultTreeViewBox(width, viewportHeight);
  const storedViewBox = state.treeViewBoxes[id];
  const viewBox = storedViewBox && storedViewBox.height <= viewportHeight * 1.5 && storedViewBox.width <= width * 1.5 ? storedViewBox : baseViewBox;
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
        <title>${escapeHtml(child.name)}; ${child.isTest ? "RELAX Test" : "Background"}; length ${fmt(child.length, 3)}${isSelected ? "; foreground aBSREL selected" : ""}</title>
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
    <text x="24" y="24" font-size="15" font-weight="700" fill="#17202a">${escapeHtml(gene.gene_id)} RELAX/foreground aBSREL tree view</text>
    <text x="24" y="43" font-size="11" fill="#4b5563">Blue = RELAX Test branches; gray = Background; red = foreground aBSREL selected; dashed orange = long branch</text>
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
  const viewportHeight = height;
  const baseViewBox = defaultTreeViewBox(width, viewportHeight);
  const storedViewBox = state.treeViewBoxes[id];
  const viewBox = storedViewBox && storedViewBox.height <= viewportHeight * 1.5 && storedViewBox.width <= width * 1.5 ? storedViewBox : baseViewBox;
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
    neglogq: Number(row.neg_log10_q || row.neg_log10_p || 0),
    color: Number(row.k) > 1 ? "#2d6cdf" : Number(row.k) < 1 ? "#c97924" : "#8792a2",
    codons: data.genes.find(gene => gene.gene_id === `${row.segment}_${row.label_set}`)?.codons || 1,
  }));
  return scatterPlot(rows, {
    id: "relax-effect-plot",
    title: "RELAX effect plot",
    xKey: "log2k",
    yKey: "neglogq",
    xLabel: "log2(K): relaxation left, intensification right",
    yLabel: "-log10 BH FDR q-value",
    colorKey: "color",
    sizeKey: "codons",
  });
}

function relaxReliabilityPlot(data) {
  const rows = data.relax.map(row => ({
    label: `${row.segment}_${row.label_set}`,
    test: Number(row.test_branches || 0),
    neglogq: Number(row.neg_log10_q || row.neg_log10_p || 0),
    color: Number(row.k) > 1 ? "#2d6cdf" : Number(row.k) < 1 ? "#c97924" : "#8792a2",
    codons: data.genes.find(gene => gene.gene_id === `${row.segment}_${row.label_set}`)?.codons || 1,
  }));
  return scatterPlot(rows, {
    id: "relax-reliability-plot",
    title: "RELAX reliability by Test branch count",
    xKey: "test",
    yKey: "neglogq",
    xLabel: "Number of RELAX Test branches",
    yLabel: "-log10 BH FDR q-value",
    colorKey: "color",
    sizeKey: "codons",
  });
}

function absrelEvidencePlot(data, gene = null, id = "absrel-evidence-plot") {
  const source = gene ? data.absrel.filter(row => `${row.segment}_${row.label_set}` === gene.gene_id) : data.absrel;
  const rows = source.map(row => ({
    label: `${row.segment}_${row.label_set} ${row.branch}`,
    branchLength: Number(row.branch_length || 0),
    neglogq: Number(row.neg_log10_q || row.neg_log10_corrected_p || 0),
    color: row.status === "selected" ? "#c84630" : "#c97924",
  }));
  return scatterPlot(rows, {
    id,
    title: gene ? `${gene.gene_id} foreground aBSREL branch evidence` : "Foreground aBSREL branch evidence",
    xKey: "branchLength",
    yKey: "neglogq",
    xLabel: "Branch length",
    yLabel: "-log10 BH q-value",
    colorKey: "color",
    showLabels: state.absrelShowLabels,
    className: "absrel-evidence-chart",
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
  if (gene.meme_sites > 0) support.push(`${gene.meme_sites} foreground MEME episodic candidate sites`);
  if (gene.fel_sites > 0) support.push(`${gene.fel_sites} foreground FEL pervasive candidate sites`);
  if (gene.absrel_branches > 0) support.push(`${gene.absrel_branches} foreground aBSREL selected branch`);
  if (Number(gene.relax_q || gene.relax_p) <= 0.05) support.push(`RELAX q=${fmt(gene.relax_q || gene.relax_p)}`);
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
  const metric = siteMetricMeta();
  const completed = data.analysis.filter(row => row.status === "pass").length;
  const notRun = data.analysis.filter(row => row.status === "not_run").length;
  const memeSites = data.meme.filter(passesSiteThreshold).length;
  const felSites = data.fel.filter(row => row.direction === "diversifying" && passesSiteThreshold(row)).length;
  const selectedBranches = data.absrel.filter(row => row.status === "selected").length;
  const mssModels = data.mss.reduce((sum, row) => sum + Number(row.model_count || 0), 0);
  const highWarnings = data.warnings.filter(row => row.severity === "high").length;
  return `<section class="page grid">
    ${cards([
      { label: "Runs", value: data.genes.length },
      { label: "Completed method runs", value: completed },
      { label: "Methods not run", value: notRun },
      { label: `Foreground FEL sites (${metric.shortLabel} <= ${state.pThreshold.toFixed(2)})`, value: felSites },
      { label: `Foreground MEME sites (${metric.shortLabel} <= ${state.pThreshold.toFixed(2)})`, value: memeSites },
      { label: "Foreground aBSREL branches", value: selectedBranches },
      { label: "MSS models", value: mssModels },
      { label: "High warnings", value: highWarnings },
    ])}
    <section class="panel">
      <h2>Method Status Heatmap</h2>
      ${exportToolbar({ tableId: "method_status_heatmap", svgId: "method-status-heatmap", filename: "method_status_heatmap" })}
      ${statusHeatmap(data)}
    </section>
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
        ${branchLabelCoverageChart(data.relax)}
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
        { key: "n_sequences", label: "Seqs" },
        { key: "codons", label: "Codons" },
        { key: "fel_sites", label: "FEL sites" },
        { key: "meme_sites", label: "MEME sites" },
        { key: "absrel_branches", label: "aBSREL branches" },
        { key: "relax_k", label: "RELAX K" },
        { key: "relax_lrt", label: "RELAX LRT", render: value => lrtCell(value) },
        { key: "relax_p", label: "RELAX p" },
        { key: "relax_q", label: "RELAX q" },
        { key: "mss_models", label: "MSS models" },
        { key: "evidence_tier", label: "Evidence", render: value => badge(value) },
      ], "gene_summary_visible")}
    </section>
  </section>`;
}

function sitesPage(data) {
  const gene = data.genes.find(row => row.gene_id === state.selectedRun) ?? data.genes[0];
  const metric = siteMetricMeta();
  const annotateQuality = site => {
    const quality = data.quality.find(row => `${row.segment}_${row.label_set}` === gene.gene_id && row.codon === site.codon) || {};
    const warning = Number(quality.gap_fraction || 0) > 0.1 ? "gappy" : Number(quality.entropy || 0) > 1.5 ? "high entropy" : "";
    return { ...site, gap_fraction: quality.gap_fraction ?? "", entropy: quality.entropy ?? "", non_gap_sequences: quality.non_gap_sequences ?? "", warning };
  };
  const sites = data.meme
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id && passesSiteThreshold(row))
    .map(annotateQuality);
  const allFelRawSites = data.fel
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id);
  const felSummarySource = data.analysis.find(row => row.segment === gene.segment && row.label_set === gene.label_set && row.method === "FEL");
  const felSummary = felSummaryRows(allFelRawSites, felSummarySource?.tested || gene.codons || "");
  const felSites = data.fel
    .filter(row => `${row.segment}_${row.label_set}` === gene.gene_id && passesSiteThreshold(row))
    .map(annotateQuality);
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
      <h3>MEME Foreground-Branch Episodic Site Evidence</h3>
      ${exportToolbar({ svgId: `meme-lollipop-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_meme_lollipop` })}
      ${lollipop(data, gene, { method: "MEME", id: `meme-lollipop-${fileSafe(gene.gene_id)}` })}
      <h3>Foreground FEL Pervasive Site Evidence</h3>
      ${exportToolbar({ svgId: `fel-lollipop-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_fel_lollipop` })}
      ${lollipop(data, gene, { method: "FEL", id: `fel-lollipop-${fileSafe(gene.gene_id)}` })}
      <div class="subpanel">
        <h3>Uncorrected FEL Summary</h3>
        <p class="muted">Counts use raw FEL p-values without multiple-testing correction. The BH q column is shown for comparison.</p>
        ${table(felSummary, [
          { key: "selection_type", label: "Selection type" },
          { key: "total_sites", label: "Total FEL sites tested" },
          { key: "uncorrected_p_lt_0_1", label: "Raw p < 0.10 sites" },
          { key: "bh_q_lte_0_1", label: "BH q <= 0.10 sites" },
          { key: "strongest_p", label: "Min raw p" },
          { key: "strongest_q", label: "Min BH q" },
          { key: "interpretation", label: "Interpretation" },
        ], `${gene.gene_id}_fel_uncorrected_summary`)}
      </div>
    </section>
    <section class="panel">
      <h2>Foreground FEL Site Table</h2>
      ${table(felSites, [
        { key: "codon", label: "Codon" },
        ...siteValueColumns(),
        { key: "alpha", label: "alpha" },
        { key: "beta", label: "beta" },
        { key: "omega", label: "omega" },
        {
          key: "selection_interpretation",
          label: "Selection call",
          render: (_, row) => `<span class="selection-call ${felSelectionClass(row)}">${escapeHtml(felSelectionInterpretation(row))}</span>`,
        },
        { key: "gap_fraction", label: "Gap fraction" },
        { key: "entropy", label: "Entropy" },
        { key: "non_gap_sequences", label: "Non-gap seqs" },
        { key: "warning", label: "Warning" },
        { key: "lrt", label: "LRT" },
        { key: "fdr_status", label: "FDR", render: value => badge(value) },
        { key: "status", label: "Status", render: value => badge(value) },
      ], `${gene.gene_id}_fel_sites_visible`, { sortable: true })}
    </section>
    <section class="panel">
      <h2>MEME Foreground Site Table</h2>
      ${table(sites, [
        { key: "codon", label: "Codon" },
        ...siteValueColumns(),
        { key: "omega_plus", label: "omega+" },
        { key: "branch_fraction", label: "Branch fraction" },
        { key: "branches_under_selection", label: "Branches" },
        { key: "gap_fraction", label: "Gap fraction" },
        { key: "entropy", label: "Entropy" },
        { key: "non_gap_sequences", label: "Non-gap seqs" },
        { key: "warning", label: "Warning" },
        { key: "lrt", label: "LRT" },
        { key: "fdr_status", label: "FDR", render: value => badge(value) },
      ], `${gene.gene_id}_meme_sites_visible`)}
    </section>
  </section>`;
}

function branchesPage(data) {
  const gene = data.genes.find(row => row.gene_id === state.selectedRun) ?? data.genes[0];
  const relaxRows = state.showRelaxNonSignificant ? data.relax : data.relax.filter(passesRelaxContextThreshold);
  return `<section class="page grid">
    <section class="panel">
      <h2>Linked RELAX/Foreground aBSREL Tree</h2>
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
      <h2>Foreground aBSREL Branch Evidence Plot</h2>
      <label class="checkbox-control">
        <input type="checkbox" id="absrel-label-toggle" ${state.absrelShowLabels ? "checked" : ""}>
        Show branch labels
      </label>
      ${exportToolbar({ svgId: "absrel-evidence-plot", filename: "absrel_branch_evidence_plot" })}
      ${absrelEvidencePlot(data)}
    </section>
    <section class="panel">
      <h2>Foreground aBSREL Branch-Level Selection</h2>
      <p class="muted">aBSREL reports a HyPhy-corrected branch p-value; the dashboard also reports a per-run BH q-value computed from raw branch p-values.</p>
      ${table(data.absrel, [
        { key: "segment", label: "Segment" },
        { key: "label_set", label: "Tree" },
        { key: "branch", label: "Branch" },
        { key: "corrected_p_value", label: "HyPhy corrected p" },
        { key: "q_value", label: "Raw-p BH q" },
        { key: "lrt", label: "LRT" },
        { key: "branch_length", label: "Branch length" },
        { key: "status", label: "Status", render: value => badge(value) },
      ], "absrel_branches_visible")}
    </section>
    <section class="panel">
      <h2>RELAX Context</h2>
      <label class="checkbox-control">
        <input type="checkbox" id="relax-nonsig-toggle" ${state.showRelaxNonSignificant ? "checked" : ""}>
        Include RELAX rows with p and q > 0.10
      </label>
      ${table(relaxRows, [
        { key: "segment", label: "Segment" },
        { key: "label_set", label: "Tree" },
        { key: "test_branches", label: "Test branches" },
        { key: "reference_branches", label: "Reference branches" },
        { key: "k", label: "K" },
        { key: "lrt", label: "LRT", render: value => lrtCell(value) },
        { key: "p_value", label: "p" },
        { key: "q_value", label: "BH q" },
        { key: "interpretation", label: "Interpretation" },
      ], "relax_results_visible")}
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
        { label: "Foreground MEME sites", value: gene.meme_sites },
        { label: "Foreground FEL sites", value: gene.fel_sites },
        { label: "RELAX K", value: gene.relax_k || "NA" },
        { label: "MSS models", value: gene.mss_models || "NA" },
        { label: "Foreground aBSREL branches", value: gene.absrel_branches },
        { label: "Test branches", value: gene.test_branches || "NA" },
      ])}
      ${exportToolbar({ svgId: `meme-lollipop-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_integrated_lollipop` })}
      ${lollipop(data, gene)}
      ${exportToolbar({ svgId: `tree-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_integrated_tree` })}
      ${treePanel(data, gene)}
      ${exportToolbar({ svgId: `absrel-evidence-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_absrel_evidence` })}
      ${absrelEvidencePlot(data, gene, `absrel-evidence-${fileSafe(gene.gene_id)}`)}
      ${exportToolbar({ svgId: `mss-frontier-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_mss_frontier` })}
      ${mssModelFrontierPlot(data, gene)}
      ${exportToolbar({ svgId: `mss-parameters-${fileSafe(gene.gene_id)}`, filename: `${gene.gene_id}_mss_parameters` })}
      ${mssParameterPlot(data, gene)}
      <h3>Evidence Summary</h3>
      ${table([gene], [
        { key: "gene_id", label: "Run" },
        { key: "fel_sites", label: "FEL" },
        { key: "meme_sites", label: "MEME sites" },
        { key: "absrel_branches", label: "aBSREL branches" },
        { key: "relax_k", label: "RELAX K" },
        { key: "relax_lrt", label: "RELAX LRT", render: value => lrtCell(value) },
        { key: "relax_p", label: "RELAX p" },
        { key: "relax_q", label: "RELAX q" },
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
  const metric = siteMetricMeta();
  return `<section class="page grid">
    <section class="panel">
      <h2>Export</h2>
      <p class="muted">Download the normalized tables used by this dashboard.</p>
      <div class="toolbar">
        ${Object.values(tableFiles).map(file => `<a class="tab" href="${TABLE_ROOT}/${file}" download>${file}</a>`).join("")}
      </div>
      <h3>Draft Results Text</h3>
      <p class="interpretation">Across the current hantavirus analyses, foreground-branch FEL identified ${data.fel.filter(row => row.direction === "diversifying" && passesSiteThreshold(row)).length} candidate pervasive diversifying sites at ${metric.longLabel} <= ${state.pThreshold.toFixed(2)}, foreground-branch MEME identified ${data.meme.filter(passesSiteThreshold).length} candidate episodic sites, foreground aBSREL identified ${data.absrel.filter(row => row.status === "selected").length} selected branches at per-run BH FDR q <= 0.05, RELAX found ${data.relax.filter(row => row.significant === true).length} significant branch-set shifts at per-run BH FDR q <= 0.05, and MSS-GA evaluated ${data.mss.reduce((sum, row) => sum + Number(row.model_count || 0), 0)} synonymous-rate class models.</p>
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
    mss: mssPage,
    "input-trees": inputTreesPage,
    warnings: warningsPage,
    export: exportPage,
  };
  if (!pages[state.tab]) state.tab = "overview";
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
  app.querySelectorAll("[data-sort-table]").forEach(button => {
    button.addEventListener("click", () => {
      const tableName = button.dataset.sortTable;
      const key = button.dataset.sortKey;
      const current = state.tableSorts[tableName];
      state.tableSorts[tableName] = {
        key,
        direction: current?.key === key && current.direction === "asc" ? "desc" : "asc",
      };
      render();
    });
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
  app.querySelector("#relax-nonsig-toggle")?.addEventListener("change", event => {
    state.showRelaxNonSignificant = event.target.checked;
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
  app.querySelectorAll("[data-export-png]").forEach(button => {
    button.addEventListener("click", () => {
      downloadPng(button.dataset.exportPng, button.dataset.filename || `${button.dataset.exportPng}_4x.png`);
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
