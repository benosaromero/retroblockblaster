const coverageForm = document.getElementById("coverage-form");
const coverageResult = document.getElementById("coverage-result");
const coverageParams = new URLSearchParams(window.location.search);
const parentSmilesParam = coverageParams.get("parent_smiles");

if (parentSmilesParam) {
  document.getElementById("coverage-parent-smiles").value = parentSmilesParam;
}

function linesToSmiles(value) {
  return value
    .split("\n")
    .map((line) => line.trim().split(/\s+/)[0])
    .filter(Boolean);
}

function toggleCoverageOptional(sectionId, button) {
  const section = document.getElementById(sectionId);

  if (!section) {
    return;
  }

  section.classList.toggle("is-open");
  button.classList.toggle("is-active", section.classList.contains("is-open"));

  if (section.classList.contains("is-open")) {
    button.textContent = button.textContent.replace("Add", "Hide");
  } else {
    button.textContent = button.textContent.replace("Hide", "Add");
  }
}

function openRetroForCoverage() {
  const jobId = coverageParams.get("job_id");

  if (jobId) {
    window.open(`/synthons/${encodeURIComponent(jobId)}`, "_blank", "noopener");
    return;
  }

  const parentSmiles = document.getElementById("coverage-parent-smiles").value.trim();
  const url = `/?prefill_smiles=${encodeURIComponent(parentSmiles)}`;
  window.open(url, "_blank", "noopener");
}

function summaryTable(summary) {
  return `
    <table class="synthon-results-table">
      <tbody>
        <tr><th>Total close analogs</th><td>${summary.total_close_analogs || 0}</td></tr>
        <tr><th>Buyable</th><td>${summary.buyable || 0}</td></tr>
        <tr><th>Makeable</th><td>${summary.makeable || 0}</td></tr>
        <tr><th>Both</th><td>${summary.both || 0}</td></tr>
        <tr><th>Buyable only</th><td>${summary.buyable_only || 0}</td></tr>
        <tr><th>Makeable only</th><td>${summary.makeable_only || 0}</td></tr>
        <tr><th>Neither</th><td>${summary.neither || 0}</td></tr>
      </tbody>
    </table>
  `;
}

coverageForm.addEventListener("submit", async (event) => {
  event.preventDefault();

  coverageResult.textContent = "Running coverage pipeline...";

  let modifications;

  try {
    modifications = JSON.parse(document.getElementById("coverage-modifications").value);
  } catch (error) {
    coverageResult.textContent = "Modifications must be valid JSON.";
    return;
  }

  const payload = {
    parent_smiles: document.getElementById("coverage-parent-smiles").value.trim(),
    modifications: modifications,
    buyable_analogs: linesToSmiles(document.getElementById("coverage-buyable").value),
    makeable_analogs: linesToSmiles(document.getElementById("coverage-makeable").value),
  };

  try {
    const response = await fetch("/api/coverage-pipeline", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(payload),
    });

    const data = await response.json();

    if (!response.ok || data.analog_generation?.error) {
      coverageResult.innerHTML = `
        <div><strong>Coverage failed.</strong></div>
        <pre class="smiles-text">${data.analog_generation?.stderr || data.error || "Unknown error"}</pre>
      `;
      return;
    }

    coverageResult.innerHTML = `
      <div class="status-line">Job ID: ${data.job_id}</div>
      ${summaryTable(data.summary || {})}
      <div class="external-links">
        <a href="${data.downloads.coverage_csv}" download>Download coverage CSV</a>
        <a href="${data.downloads.close_analogs_csv}" download>Download analog CSV</a>
        <a href="${data.downloads.close_analogs_png}" target="_blank">View analog PNG</a>
        <a href="${data.downloads.parent_numbered_png}" target="_blank">View numbered parent</a>
      </div>
    `;
  } catch (error) {
    coverageResult.textContent = `Could not run coverage pipeline: ${error.message}`;
  }
});
