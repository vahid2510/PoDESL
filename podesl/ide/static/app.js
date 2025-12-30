const KEYWORDS = [
  "PROBLEM",
  "GIVEN",
  "SOLVE",
  "EQUATIONS",
  "REPORT",
  "print",
  "export",
  "plot",
  "help",
  "linspace",
  "zeros",
  "ones",
  "array",
  "diag",
];

const NEW_TEMPLATE = `PROBLEM "Solid : Truss2D : Static"
GIVEN
  nodes = [[0,0],[1,0],[1,1]]
  elems = [[0,1],[1,2]]
  E = 210e9
  A = 1.0e-4
  loads = [[2, 0.0, -1000.0]]
  fix = [[0,"both",0.0],[1,"uy",0.0]]
REPORT
  print U
`;

const state = {
  examples: [],
  filteredExamples: [],
  problems: [],
  identifiers: [],
  tabs: [],
  activeTabId: null,
  helpMatches: [],
  helpDetail: "",
  history: [],
  notebookCells: [],
  mode: "editor",
  jobName: "custom_job",
  exportPath: "build/ide_outputs/model.inp",
  theme: "light",
};

const els = {
  exampleList: document.getElementById("exampleList"),
  exampleSearch: document.getElementById("exampleSearch"),
  refreshExamples: document.getElementById("refreshExamples"),
  serverStatus: document.getElementById("serverStatus"),
  helpSearch: document.getElementById("helpSearch"),
  helpMatches: document.getElementById("helpMatches"),
  helpDetails: document.getElementById("helpDetails"),
  loadHelp: document.getElementById("loadHelp"),
  consoleOutput: document.getElementById("consoleOutput"),
  resultView: document.getElementById("resultView"),
  clearResult: document.getElementById("clearResult"),
  menuBar: document.getElementById("menuBar"),
  statusBar: document.getElementById("statusBar"),
  tabBar: document.getElementById("tabBar"),
  chartX: document.getElementById("chartX"),
  chartY: document.getElementById("chartY"),
  plotChart: document.getElementById("plotChart"),
  reportNotes: document.getElementById("reportNotes"),
  reportMetrics: document.getElementById("reportMetrics"),
  generateReport: document.getElementById("generateReport"),
  historyList: document.getElementById("historyList"),
  clearHistory: document.getElementById("clearHistory"),
  notebookPane: document.getElementById("notebookPane"),
  notebookCells: document.getElementById("notebookCells"),
  addNotebookCell: document.getElementById("addNotebookCell"),
  runNotebookAll: document.getElementById("runNotebookAll"),
  toggleNotebookMode: document.getElementById("toggleNotebookMode"),
  modalOverlay: document.getElementById("modalOverlay"),
  modalTitle: document.getElementById("modalTitle"),
  modalBody: document.getElementById("modalBody"),
  closeModal: document.getElementById("closeModal"),
};

const editor = ace.edit("editor");
editor.setTheme("ace/theme/xcode");
editor.session.setMode("ace/mode/python");
editor.setShowPrintMargin(false);
editor.setOptions({
  enableBasicAutocompletion: true,
  enableLiveAutocompletion: true,
  enableSnippets: true,
  tabSize: 2,
  useSoftTabs: true,
});

const languageTools = ace.require("ace/ext/language_tools");

languageTools.setCompleters([
  {
    getCompletions: (editorInstance, session, pos, prefix, callback) => {
      const completions = buildCompletions(prefix, session, pos);
      callback(null, completions);
    },
  },
]);

editor.session.on("change", () => {
  syncActiveTabContent();
  setDirty(true, true);
  scheduleIdentifierRefresh();
});

editor.session.setUseWrapMode(true);

editor.commands.addCommand({
  name: "runCode",
  bindKey: { win: "Ctrl-Enter", mac: "Command-Enter" },
  exec: () => runSource(),
});

editor.commands.addCommand({
  name: "formatCode",
  bindKey: { win: "Ctrl-Shift-F", mac: "Command-Shift-F" },
  exec: () => formatSource(),
});

editor.commands.addCommand({
  name: "saveFile",
  bindKey: { win: "Ctrl-S", mac: "Command-S" },
  exec: () => saveFile(),
});

els.refreshExamples.addEventListener("click", () => loadExamples());
els.exampleSearch.addEventListener("input", () => renderExamples());
els.helpSearch.addEventListener("input", () => debounceHelp());
els.loadHelp.addEventListener("click", () => fetchHelp(els.helpSearch.value.trim()));
els.clearResult.addEventListener("click", () => {
  els.consoleOutput.textContent = "";
  els.resultView.textContent = "";
});
els.plotChart.addEventListener("click", () => plotChart());
els.generateReport.addEventListener("click", () => generateReport());
els.clearHistory?.addEventListener("click", () => clearHistory());
els.historyList?.addEventListener("click", (event) => {
  const shareBtn = event.target.closest("[data-share]");
  if (shareBtn) {
    const id = shareBtn.getAttribute("data-share");
    shareHistoryEntry(id);
  }
});
els.closeModal?.addEventListener("click", () => closeModal());
els.modalOverlay?.addEventListener("click", (event) => {
  if (event.target === els.modalOverlay) {
    closeModal();
  }
});
els.menuBar?.addEventListener("click", (event) => {
  const btn = event.target.closest("[data-menu]");
  if (!btn) return;
  event.preventDefault();
  const action = btn.getAttribute("data-menu");
  handleMenuAction(action);
});
els.addNotebookCell?.addEventListener("click", () => addNotebookCell());
els.runNotebookAll?.addEventListener("click", () => runNotebook());
els.toggleNotebookMode?.addEventListener("click", () => setMode("editor"));
els.tabBar?.addEventListener("click", (event) => {
  const closeBtn = event.target.closest("[data-close-tab]");
  if (closeBtn) {
    closeTab(closeBtn.getAttribute("data-close-tab"));
    return;
  }
  const tabBtn = event.target.closest("[data-tab-id]");
  if (tabBtn) {
    setActiveTab(tabBtn.getAttribute("data-tab-id"));
  }
});


async function requestJSON(url, options = {}) {
  const response = await fetch(url, {
    headers: { "Content-Type": "application/json" },
    ...options,
  });
  const data = await response.json().catch(() => ({}));
  if (!response.ok) {
    const error = new Error(data.error || "Request failed");
    error.details = data;
    throw error;
  }
  return data;
}

async function loadExamples() {
  els.serverStatus.textContent = "Loading...";
  try {
    const data = await requestJSON("/api/examples");
    state.examples = data.examples || [];
    renderExamples();
    els.serverStatus.textContent = `Examples: ${state.examples.length}`;
  } catch (error) {
    console.error(error);
    els.serverStatus.textContent = "Failed";
  }
}

function renderExamples() {
  const filter = (els.exampleSearch.value || "").trim().toLowerCase();
  const entries = state.examples.filter((item) =>
    item.path.toLowerCase().includes(filter)
  );
  state.filteredExamples = entries;
  els.exampleList.innerHTML = "";
  const activePath = getActiveTab()?.path || null;
  entries.slice(0, 120).forEach((item) => {
    const li = document.createElement("li");
    const label = document.createElement("span");
    label.className = "file-path";
    label.textContent = item.path;
    const badge = document.createElement("span");
    const kind = item.kind || (item.path.toLowerCase().endsWith(".dsl") ? "dsl" : "doc");
    badge.className = `tag ${kind === "doc" ? "doc" : "dsl"}`;
    badge.textContent = kind === "doc" ? "DOC" : "DSL";
    li.appendChild(label);
    li.appendChild(badge);
    if (activePath && activePath === item.path) {
      li.classList.add("active");
    }
    li.addEventListener("click", () => openExample(item.path));
    els.exampleList.appendChild(li);
  });
}

async function openExample(path) {
  try {
    const existing = state.tabs.find((tab) => tab.path === path);
    if (existing) {
      setActiveTab(existing.id);
      setConsole(`Activated ${path}`);
      return;
    }
    const data = await requestJSON(`/api/file?path=${encodeURIComponent(path)}`);
    createTab({
      path,
      name: path.split("/").pop(),
      content: data.content,
    });
    setDirty(false, true);
    renderExamples();
    setConsole(`Loaded ${path}`);
  } catch (error) {
    setConsole(error.message || "Failed to open file", true);
  }
}

function setDirty(flag, skipContentUpdate = false) {
  const tab = getActiveTab();
  if (!tab) {
    return;
  }
  tab.dirty = flag;
  if (!skipContentUpdate) {
    tab.content = editor.getValue();
  }
  const label = tab.name || "untitled.dsl";
  if (els.statusBar) {
    els.statusBar.textContent = flag ? `Unsaved — ${label}` : `Saved — ${label}`;
  }
  updateTabBar();
}

async function runSource() {
  const tab = getActiveTab();
  const payload = {
    source: editor.getValue(),
    filename: tab?.path || tab?.name || "<workspace>",
  };
  setConsole("Running...");
  try {
    const data = await requestJSON("/api/run", {
      method: "POST",
      body: JSON.stringify(payload),
    });
    const combined = `${data.stdout || ""}${data.stderr || ""}`;
    setConsole(combined || "✔ Execution completed");
    renderResult(data.result);
    clearDiagnostics();
    loadHistory();
  } catch (error) {
    setConsole(error.message || "Execution failed", true);
    if (error.details && error.details.diagnostics) {
      applyDiagnostics(error.details.diagnostics);
    }
    loadHistory();
  }
}

function applyDiagnostics(diags = []) {
  const annotations = diags.map((diag) => ({
    row: Math.max(0, (diag.line || 1) - 1),
    column: Math.max(0, diag.column || 0),
    text: diag.message || "Error",
    type: "error",
  }));
  editor.session.setAnnotations(annotations);
}

function clearDiagnostics() {
  editor.session.clearAnnotations();
}

async function formatSource() {
  const payload = { source: editor.getValue() };
  try {
    const data = await requestJSON("/api/format", {
      method: "POST",
      body: JSON.stringify(payload),
    });
    editor.setValue(data.text, -1);
    setConsole("Document formatted");
    setDirty(true, true);
  } catch (error) {
    setConsole(error.message || "Formatting failed", true);
  }
}

async function saveFile(forcePrompt = false) {
  const tab = getActiveTab();
  if (!tab) {
    return;
  }
  let targetPath = tab.path;
  if (!targetPath || forcePrompt) {
    targetPath = window.prompt("Save path (relative to workspace):", "examples/custom.dsl");
    if (!targetPath) return;
    tab.path = targetPath;
    tab.name = targetPath.split(/[\\/]/).pop();
  }
  try {
    await requestJSON("/api/save", {
      method: "POST",
      body: JSON.stringify({ path: targetPath, source: editor.getValue() }),
    });
    setDirty(false, true);
    loadExamples();
    setConsole(`Saved ${targetPath}`);
  } catch (error) {
    setConsole(error.message || "Save failed", true);
  }
}

async function exportAbaqus() {
  const metadata = {};
  if (state.jobName) {
    metadata.job_name = state.jobName;
  }
  const tab = getActiveTab();
  const payload = {
    source: editor.getValue(),
    filename: tab?.path || tab?.name || "<workspace>",
    output_path: state.exportPath || null,
    metadata,
  };
  setConsole("Exporting to Abaqus...");
  try {
    const data = await requestJSON("/api/export", {
      method: "POST",
      body: JSON.stringify(payload),
    });
    setConsole(
      `Abaqus deck ready: ${data.path}\nNodes: ${data.node_count}, Elements: ${data.element_count}`
    );
  } catch (error) {
    setConsole(error.message || "Export failed", true);
  }
}

function renderResult(result) {
  const payload = result === undefined ? null : result;
  try {
    els.resultView.textContent = JSON.stringify(payload, null, 2);
  } catch (error) {
    els.resultView.textContent = String(payload);
  }
}

function setConsole(message, isError = false) {
  els.consoleOutput.textContent = message;
  els.consoleOutput.style.color = isError ? "#ff9e9e" : "#c5c7f7";
}

function copyResult() {
  const text = `${els.consoleOutput.textContent}\n${els.resultView.textContent}`;
  navigator.clipboard.writeText(text).then(
    () => setConsole("Copied result to clipboard"),
    () => setConsole("Copy failed", true)
  );
}

function getActiveTab() {
  return state.tabs.find((tab) => tab.id === state.activeTabId) || null;
}

function saveActiveTabContent() {
  const tab = getActiveTab();
  if (!tab) {
    return;
  }
  tab.content = editor.getValue();
}

function setActiveTab(id, initializing = false) {
  if (!initializing) {
    saveActiveTabContent();
  }
  state.activeTabId = id;
  const tab = getActiveTab();
  if (!tab) {
    return;
  }
  editor.setValue(tab.content, -1);
  setDirty(tab.dirty, true);
  updateTabBar();
}

function createTab({ path = null, name, content = "" }) {
  saveActiveTabContent();
  const id = `tab-${Date.now()}-${Math.random().toString(36).slice(2, 8)}`;
  const tab = {
    id,
    path,
    name: name || (path ? path.split("/").pop() : "untitled.dsl"),
    content,
    dirty: false,
  };
  state.tabs.push(tab);
  setActiveTab(id, true);
  editor.setValue(content, -1);
  setDirty(false, true);
  updateTabBar();
}

function updateTabBar() {
  if (!els.tabBar) {
    return;
  }
  els.tabBar.innerHTML = "";
  state.tabs.forEach((tab) => {
    const wrapper = document.createElement("div");
    wrapper.className = `tab${tab.id === state.activeTabId ? " active" : ""}`;
    wrapper.setAttribute("data-tab-id", tab.id);
    const label = document.createElement("span");
    label.textContent = tab.name;
    wrapper.appendChild(label);
    if (tab.dirty) {
      const dirty = document.createElement("span");
      dirty.className = "dirty-indicator";
      dirty.textContent = "*";
      wrapper.appendChild(dirty);
    }
    const close = document.createElement("button");
    close.className = "close";
    close.textContent = "×";
    close.type = "button";
    close.setAttribute("data-close-tab", tab.id);
    wrapper.appendChild(close);
    els.tabBar.appendChild(wrapper);
  });
}

function syncActiveTabContent() {
  const tab = getActiveTab();
  if (!tab) {
    return;
  }
  tab.content = editor.getValue();
}

function closeTab(id) {
  const idx = state.tabs.findIndex((tab) => tab.id === id);
  if (idx === -1) {
    return;
  }
  const closingActive = state.tabs[idx].id === state.activeTabId;
  state.tabs.splice(idx, 1);
  if (!state.tabs.length) {
    createTab({ content: NEW_TEMPLATE });
    return;
  }
  if (closingActive) {
    const next = state.tabs[Math.max(0, idx - 1)];
    state.activeTabId = next.id;
    editor.setValue(next.content, -1);
    setDirty(next.dirty, true);
  }
  updateTabBar();
}

let untitledCounter = 1;

function newFile() {
  const tab = getActiveTab();
  if (tab && tab.dirty) {
    const confirmed = window.confirm("Discard current changes and start a new DSL file?");
    if (!confirmed) {
      return;
    }
  }
  createTab({
    name: `untitled-${untitledCounter++}.dsl`,
    content: NEW_TEMPLATE,
  });
  setDirty(true, true);
  setConsole("Loaded new DSL template");
}

function insertTemplatePayload({ type, span, load, E, param }) {
  let template = "";
  if (type === "truss2d") {
    template = `PROBLEM "Solid : Truss2D : Static"
GIVEN
  nodes = [[0,0],[${span},0],[${span},${span / 2}]]
  elems = [[0,1],[1,2]]
  E = ${E}
  A = ${param}
  loads = [[2, 0.0, ${load}]]
  fix = [[0,"both",0.0],[1,"uy",0.0]]
REPORT
  print U
  print element_forces
`;
  } else if (type === "frame2d") {
    template = `PROBLEM "Solid : Frame2D : Static"
GIVEN
  nodes = [[0,0],[${span},0]]
  elems = [[0,1]]
  E = ${E}
  A = ${param}
  I = ${param * 1e-4}
  loads = [[1,"uy",${load}]]
  supports = [[0,1,1,1]]
REPORT
  print U
`;
  } else {
    template = `PROBLEM "Thermal : Heat2D : Steady"
GIVEN
  nodes = [[0,0],[${span},0],[${span},${span}],[0,${span}]]
  elems = [[0,1,2,3]]
  k = ${param}
  heat = ${Math.abs(load)}
  dirichlet = [[0,"temp",100.0],[1,"temp",50.0]]
REPORT
  print T
`;
  }
  createTab({
    name: `${type}-${untitledCounter++}.dsl`,
    content: template,
  });
  setDirty(true, true);
  setConsole("Inserted starter template");
}

function openModal(title, builder) {
  if (!els.modalOverlay) {
    return;
  }
  els.modalTitle.textContent = title;
  const body = els.modalBody;
  body.innerHTML = "";
  builder(body);
  els.modalOverlay.classList.remove("hidden");
}

function closeModal() {
  els.modalOverlay?.classList.add("hidden");
  if (els.modalBody) {
    els.modalBody.innerHTML = "";
  }
}

function showSettingsModal() {
  openModal("Job Settings", (body) => {
    const jobLabel = document.createElement("label");
    jobLabel.textContent = "Job name";
    const jobInput = document.createElement("input");
    jobInput.type = "text";
    jobInput.value = state.jobName;
    jobLabel.appendChild(jobInput);
    const pathLabel = document.createElement("label");
    pathLabel.textContent = "Target path";
    const pathInput = document.createElement("input");
    pathInput.type = "text";
    pathInput.value = state.exportPath;
    pathLabel.appendChild(pathInput);
    const saveBtn = document.createElement("button");
    saveBtn.textContent = "Save";
    saveBtn.addEventListener("click", () => {
      state.jobName = jobInput.value.trim() || "job";
      state.exportPath = pathInput.value.trim() || "build/ide_outputs/model.inp";
      closeModal();
      setConsole("Job settings updated.");
    });
    body.appendChild(jobLabel);
    body.appendChild(pathLabel);
    body.appendChild(saveBtn);
  });
}

function showTemplateModal() {
  openModal("Template Builder", (body) => {
    const typeLabel = document.createElement("label");
    typeLabel.textContent = "Problem type";
    const typeSelect = document.createElement("select");
    ["truss2d", "frame2d", "heat2d"].forEach((value) => {
      const opt = document.createElement("option");
      opt.value = value;
      opt.textContent =
        value === "truss2d"
          ? "Truss 2D Static"
          : value === "frame2d"
          ? "Frame 2D Static"
          : "Heat 2D Steady";
      typeSelect.appendChild(opt);
    });
    typeLabel.appendChild(typeSelect);
    const spanLabel = document.createElement("label");
    spanLabel.textContent = "Span (m)";
    const spanInput = document.createElement("input");
    spanInput.type = "number";
    spanInput.step = "0.1";
    spanInput.value = "1.0";
    spanLabel.appendChild(spanInput);
    const loadLabel = document.createElement("label");
    loadLabel.textContent = "Load (N)";
    const loadInput = document.createElement("input");
    loadInput.type = "number";
    loadInput.value = "-1000";
    loadLabel.appendChild(loadInput);
    const eLabel = document.createElement("label");
    eLabel.textContent = "E (Pa)";
    const eInput = document.createElement("input");
    eInput.type = "number";
    eInput.value = "210000000000";
    eLabel.appendChild(eInput);
    const paramLabel = document.createElement("label");
    paramLabel.textContent = "Area / Conductivity";
    const paramInput = document.createElement("input");
    paramInput.type = "number";
    paramInput.step = "0.0001";
    paramInput.value = "0.0001";
    paramLabel.appendChild(paramInput);
    const insertBtn = document.createElement("button");
    insertBtn.textContent = "Insert Template";
    insertBtn.addEventListener("click", () => {
      insertTemplatePayload({
        type: typeSelect.value,
        span: parseFloat(spanInput.value) || 1.0,
        load: parseFloat(loadInput.value) || -1000.0,
        E: parseFloat(eInput.value) || 210e9,
        param: parseFloat(paramInput.value) || 1e-4,
      });
      closeModal();
    });
    body.appendChild(typeLabel);
    body.appendChild(spanLabel);
    body.appendChild(loadLabel);
    body.appendChild(eLabel);
    body.appendChild(paramLabel);
    body.appendChild(insertBtn);
  });
}

function handleMenuAction(action) {
  switch (action) {
    case "new-file":
      newFile();
      break;
    case "save-file":
      saveFile();
      break;
    case "save-as":
      saveFile(true);
      break;
    case "run":
      runSource();
      break;
    case "run-notebook":
      runNotebook();
      break;
    case "export":
      exportAbaqus();
      break;
    case "copy":
      copyResult();
      break;
    case "chart":
      plotChart();
      break;
    case "report":
      generateReport();
      break;
    case "format":
      formatSource();
      break;
    case "open-template":
      showTemplateModal();
      break;
    case "open-settings":
      showSettingsModal();
      break;
    case "open-notebook":
      setMode("notebook");
      break;
    case "open-editor":
      setMode("editor");
      break;
    default:
      break;
  }
}

function setMode(mode) {
  state.mode = mode;
  const isNotebook = mode === "notebook";
  els.notebookPane?.classList.toggle("hidden", !isNotebook);
  const editorContainer = document.getElementById("editor");
  if (editorContainer) {
    editorContainer.style.display = isNotebook ? "none" : "block";
  }
  if (isNotebook && (!state.notebookCells || !state.notebookCells.length)) {
    addNotebookCell('PROBLEM "Notebook Model"\n');
    addNotebookCell("GIVEN\n  # add DSL inputs\n");
  }
  renderNotebook();
}

function addNotebookCell(content = "") {
  state.notebookCells = state.notebookCells || [];
  const cell = {
    id: `cell-${Date.now()}-${Math.random().toString(36).slice(2, 6)}`,
    content,
    output: "",
  };
  state.notebookCells.push(cell);
  renderNotebook();
}

function renderNotebook() {
  if (!els.notebookCells) {
    return;
  }
  const container = els.notebookCells;
  container.innerHTML = "";
  (state.notebookCells || []).forEach((cell) => {
    const wrapper = document.createElement("div");
    wrapper.className = "notebook-cell";
    const textarea = document.createElement("textarea");
    textarea.value = cell.content;
    textarea.addEventListener("input", () => {
      cell.content = textarea.value;
    });
    wrapper.appendChild(textarea);
    const toolbar = document.createElement("div");
    toolbar.className = "cell-toolbar";
    const runBtn = document.createElement("button");
    runBtn.textContent = "Run Cell";
    runBtn.addEventListener("click", () => runNotebook(cell.id));
    toolbar.appendChild(runBtn);
    const deleteBtn = document.createElement("button");
    deleteBtn.textContent = "Delete";
    deleteBtn.addEventListener("click", () => removeNotebookCell(cell.id));
    toolbar.appendChild(deleteBtn);
    wrapper.appendChild(toolbar);
    const output = document.createElement("div");
    output.className = "cell-output";
    output.textContent = cell.output || "No output";
    wrapper.appendChild(output);
    container.appendChild(wrapper);
  });
  els.notebookPane?.classList.toggle("hidden", state.mode !== "notebook");
}

function removeNotebookCell(id) {
  state.notebookCells = (state.notebookCells || []).filter((cell) => cell.id !== id);
  renderNotebook();
}

function assembleNotebookSource() {
  return (state.notebookCells || []).map((cell) => cell.content || "").join("\n\n");
}

async function runNotebook(targetCellId = null) {
  if (!state.notebookCells || !state.notebookCells.length) {
    setConsole("Notebook is empty.", true);
    return;
  }
  const source = assembleNotebookSource();
  try {
    const data = await requestJSON("/api/run", {
      method: "POST",
      body: JSON.stringify({ source, filename: "notebook.dsl" }),
    });
    const combined = `${data.stdout || ""}${data.stderr || ""}`;
    const outputText = combined || JSON.stringify(data.result, null, 2);
    if (targetCellId) {
      const cell = state.notebookCells.find((c) => c.id === targetCellId);
      if (cell) {
        cell.output = outputText;
      }
    } else {
      state.notebookCells.forEach((cell) => {
        cell.output = outputText;
      });
    }
    renderNotebook();
    setConsole("Notebook executed.");
    loadHistory();
  } catch (error) {
    setConsole(error.message || "Notebook run failed", true);
    loadHistory();
  }
}

const chartCtx = document.getElementById("chartCanvas").getContext("2d");
let chartInstance = null;

async function plotChart() {
  const xExpr = (els.chartX.value || "").trim();
  const yExpr = (els.chartY.value || "").trim();
  if (!xExpr || !yExpr) {
    setConsole("Chart expressions are required.", true);
    return;
  }
  try {
    const data = await requestJSON("/api/chart", {
      method: "POST",
      body: JSON.stringify({ x_expr: xExpr, y_expr: yExpr }),
    });
    renderChart(data.x, data.y, yExpr);
    setConsole("Chart generated.");
  } catch (error) {
    setConsole(error.message || "Chart request failed", true);
  }
}

function renderChart(xValues, yValues, label = "Series") {
  if (chartInstance) {
    chartInstance.destroy();
  }
  chartInstance = new Chart(chartCtx, {
    type: "line",
    data: {
      labels: xValues,
      datasets: [
        {
          label,
          data: yValues,
          fill: false,
          tension: 0.15,
          borderColor: "#4fb3ff",
          pointRadius: 1.5,
        },
      ],
    },
    options: {
      responsive: true,
      plugins: {
        legend: { display: true },
      },
      scales: {
        x: { title: { display: true, text: els.chartX.value || "x" } },
        y: { title: { display: true, text: els.chartY.value || "y" } },
      },
    },
  });
}

async function generateReport() {
  const metrics = parseMetrics(els.reportMetrics.value || "");
  try {
    const summary = await requestJSON("/api/report", {
      method: "POST",
      body: JSON.stringify({ fields: metrics, notes: els.reportNotes.value || "" }),
    });
    downloadReport(summary);
    setConsole("Report downloaded.");
  } catch (error) {
    setConsole(error.message || "Report generation failed", true);
  }
}

function parseMetrics(text) {
  return text
    .split("\n")
    .map((line) => line.trim())
    .filter(Boolean)
    .map((line) => {
      const [labelPart, ...rest] = line.split("=");
      if (rest.length === 0) {
        return { label: labelPart.trim(), expr: labelPart.trim() };
      }
      return { label: labelPart.trim(), expr: rest.join("=").trim() };
    });
}

function downloadReport(summary) {
  const lines = [];
  lines.push("# PoDESL Report");
  const title = summary.problem?.title || "Unknown problem";
  lines.push(`**Problem:** ${title}`);
  const now = new Date().toISOString();
  lines.push(`**Generated:** ${now}`);
  lines.push("");
  lines.push("## Metrics");
  if (summary.metrics && summary.metrics.length) {
    summary.metrics.forEach((metric) => {
      const value =
        typeof metric.value === "object"
          ? "`" + JSON.stringify(metric.value) + "`"
          : metric.value;
      lines.push(`- **${metric.label || metric.expr}** = ${value}`);
    });
  } else {
    lines.push("_No metrics requested._");
  }
  lines.push("");
  lines.push("## Notes");
  lines.push(summary.notes || "_No notes provided._");
  lines.push("");
  lines.push("## Raw JSON");
  lines.push("```json");
  lines.push(JSON.stringify(summary, null, 2));
  lines.push("```");
  const blob = new Blob([lines.join("\n")], { type: "text/markdown" });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = `podesl-report-${Date.now()}.md`;
  document.body.appendChild(anchor);
  anchor.click();
  document.body.removeChild(anchor);
  URL.revokeObjectURL(url);
}

async function loadHistory() {
  try {
    const data = await requestJSON("/api/history");
    state.history = data.history || [];
    renderHistory();
  } catch (error) {
    console.warn("Unable to load history", error);
  }
}

function renderHistory() {
  if (!els.historyList) {
    return;
  }
  els.historyList.innerHTML = "";
  state.history.forEach((entry) => {
    const li = document.createElement("li");
    li.className = `history-item ${entry.status || "ok"}`;
    const meta = document.createElement("div");
    meta.className = "history-meta";
    const info = document.createElement("span");
    const when = new Date((entry.timestamp || Date.now()) * 1000);
    info.textContent = `${when.toLocaleTimeString()} • ${entry.path || "untitled"}`;
    const status = document.createElement("span");
    status.textContent = entry.status?.toUpperCase() || "OK";
    meta.appendChild(info);
    meta.appendChild(status);
    li.appendChild(meta);
    const message = document.createElement("div");
    message.textContent = entry.message || "";
    li.appendChild(message);
    const actions = document.createElement("div");
    actions.className = "history-actions";
    const share = document.createElement("button");
    share.textContent = "Copy JSON";
    share.type = "button";
    share.setAttribute("data-share", entry.id);
    actions.appendChild(share);
    li.appendChild(actions);
    els.historyList.appendChild(li);
  });
}

async function clearHistory() {
  try {
    await requestJSON("/api/history", {
      method: "POST",
      body: JSON.stringify({ action: "clear" }),
    });
    state.history = [];
    renderHistory();
    setConsole("History cleared.");
  } catch (error) {
    setConsole(error.message || "Failed to clear history", true);
  }
}

function shareHistoryEntry(id) {
  const entry = state.history.find((item) => item.id === id);
  if (!entry) {
    return;
  }
  navigator.clipboard
    .writeText(JSON.stringify(entry, null, 2))
    .then(() => setConsole("History entry copied."), () => setConsole("Copy failed", true));
}

function applyTheme() {
  document.body.classList.add("theme-light");
}

function buildCompletions(prefix, session, pos) {
  const token = (prefix || "").toLowerCase();
  const dataset = new Set([...KEYWORDS, ...state.problems, ...state.identifiers]);
  const completions = [];
  dataset.forEach((word) => {
    if (!token || word.toLowerCase().startsWith(token)) {
      const isProblem = word.includes(":");
      completions.push({
        caption: word,
        value: isProblem ? `"${word}"` : word,
        meta: isProblem ? "PROBLEM" : "DSL",
        score: isProblem ? 1000 : 500,
      });
    }
  });
  return completions.slice(0, 60);
}

let identifierTimer = null;
function scheduleIdentifierRefresh() {
  clearTimeout(identifierTimer);
  identifierTimer = setTimeout(() => {
    const source = editor.getValue();
    const regex = /^\s*([A-Za-z_]\w*)\s*=/gm;
    const ids = new Set();
    let match;
    while ((match = regex.exec(source))) {
      ids.add(match[1]);
    }
    state.identifiers = Array.from(ids);
  }, 150);
}

let helpTimer = null;
function debounceHelp() {
  clearTimeout(helpTimer);
  helpTimer = setTimeout(() => fetchHelp(els.helpSearch.value.trim()), 200);
}

async function fetchHelp(term) {
  try {
    const data = await requestJSON(`/api/help?query=${encodeURIComponent(term || "")}`);
    state.helpMatches = data.matches || [];
    state.helpDetail = data.detail || "";
    renderHelp();
  } catch (error) {
    setConsole("Unable to load help", true);
  }
}

function renderHelp() {
  els.helpMatches.innerHTML = "";
  state.helpMatches.slice(0, 30).forEach((item) => {
    const btn = document.createElement("button");
    btn.textContent = item;
    btn.addEventListener("click", () => fetchHelp(item));
    els.helpMatches.appendChild(btn);
  });
  els.helpDetails.textContent = state.helpDetail || "Select a help topic.";
}

async function loadProblems() {
  try {
    const data = await requestJSON("/api/problems");
    state.problems = data.problems || [];
  } catch (error) {
    console.warn("Unable to load problem catalog", error);
  }
}

async function refreshStatus() {
  try {
    const data = await requestJSON("/api/status");
    els.serverStatus.textContent = `Ready · ${new Date(data.time * 1000).toLocaleTimeString()}`;
  } catch (error) {
    els.serverStatus.textContent = "Server offline?";
  } finally {
    setTimeout(refreshStatus, 15000);
  }
}

function setupBeforeUnload() {
  window.addEventListener("beforeunload", (event) => {
    if (!state.tabs.some((tab) => tab.dirty)) return;
    event.preventDefault();
    event.returnValue = "";
  });
}

function init() {
  applyTheme();
  if (!state.tabs.length) {
    createTab({
      name: `untitled-${untitledCounter++}.dsl`,
      content: NEW_TEMPLATE,
    });
    setDirty(false, true);
  }
  setMode(state.mode);
  loadExamples();
  loadProblems();
  fetchHelp("");
  loadHistory();
  refreshStatus();
  scheduleIdentifierRefresh();
  setupBeforeUnload();
}

init();
