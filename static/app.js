/*
  Created by Benjamin "Ben" Romero (Undergrad Northwestern University '29), Summer 2026.
  Developed in the Irwin Lab at UCSF.
  Released as part of an open-source UCSF research tool; please credit the creator and lab when using or adapting this work.
*/

const form = document.getElementById("job-form");
const result = document.getElementById("result");
const selectedCombos = {};
const loadedBlockLists = {};
const buildingBlockStore = {};
let buildingBlockCounter = 0;
const buildingBlockPageSize = 20;
const enumerationPageSize = 20;
const defaultSearchSettings = {
  source: "db050",
  tc_cutoff: 0.6,
  limit: 50,
  fallback_top_match: true,
};

const appParams = new URLSearchParams(window.location.search);
const prefillSmiles = appParams.get("prefill_smiles");
const synthonReloadTimers = {};

if (prefillSmiles && document.getElementById("user-input")) {
  document.getElementById("user-input").value = prefillSmiles;
}

window.activeSearchSettings = window.activeSearchSettings || defaultSearchSettings;
const synthonSearchSettings = {};

function toggleButton(button) {
  button.classList.toggle("is-active");
}

function svgUrl(path) {
  if (!path) {
    return "";
  }

  const parts = path.split("/");
  const jobIndex = parts.indexOf("jobs");

  if (jobIndex === -1 || parts.length < jobIndex + 4) {
    return "";
  }

  const jobId = parts[jobIndex + 1];
  const filename = parts[parts.length - 1];

  return `/image/${jobId}/${filename}`;
}

function moleculeCell(smiles, svgPath) {
  const url = svgUrl(svgPath);

  if (!smiles && !url) {
    return "";
  }

  const imageHtml = url
    ? `<img class="mol-image-small" src="${url}" alt="Molecule structure" onclick="openImageModal('${url}', 'Molecule structure')">`
    : `<div class="no-structure-image">No stable/valid structure image</div>`;

  const smilesHtml = smiles
    ? `<div class="smiles-text">${smiles}</div>`
    : "";

  return `${imageHtml}${smilesHtml}`;
}

function renderProductRows(products, visibleCount) {
  return products.slice(0, visibleCount).map((product) => `
    <tr>
      <td class="smiles-text">${product.product_smiles || ""}</td>
      <td>${product.tranche || ""}</td>
      <td class="smiles-text">${(product.building_blocks || []).join("<br>")}</td>
    </tr>
  `).join("");
}

function productsToTable(title, products) {
  if (!products || products.length === 0) {
    return `
      <div><strong>${title}</strong></div>
      <div>No products generated.</div>
    `;
  }

  const visibleCount = Math.min(enumerationPageSize, products.length);
  const productListKey = `products-${buildingBlockCounter}`;
  buildingBlockCounter += 1;
  buildingBlockStore[productListKey] = products;

  return `
    <div><strong>${title}</strong></div>

    <div class="synthon-text product-display-status">
      Displaying ${visibleCount}/${products.length}
    </div>

    <table class="synthon-results-table">
      <thead>
        <tr>
          <th>Product SMILES</th>
          <th>Tranche</th>
          <th>Building Blocks</th>
        </tr>
      </thead>
      <tbody>
        ${renderProductRows(products, visibleCount)}
      </tbody>
    </table>

    ${
      visibleCount < products.length
        ? `
          <button
            type="button"
            class="synthon-info-button view-molecule-button"
            data-product-list-key="${productListKey}"
            data-visible-count="${visibleCount}"
            onclick="loadMoreProducts(this)"
          >
            Load more
          </button>
        `
        : ""
    }
  `;
}

function loadMoreProducts(button) {
  const container = button.parentElement;
  const tableBody = container.querySelector("tbody");
  const status = container.querySelector(".product-display-status");

  const key = button.dataset.productListKey;
  const products = buildingBlockStore[key] || [];

  const currentCount = Number(button.dataset.visibleCount || enumerationPageSize);
  const nextCount = Math.min(currentCount + enumerationPageSize, products.length);

  tableBody.innerHTML = renderProductRows(products, nextCount);
  status.textContent = `Displaying ${nextCount}/${products.length}`;
  button.dataset.visibleCount = String(nextCount);

  if (nextCount >= products.length) {
    button.remove();
  }
}

function synthonInfoButton(pairId, reactionId, bbSlot, reagentCount, synthonId, referenceSmiles) {
  if (!synthonId) {
    return "";
  }

  const safePairId = JSON.stringify(pairId);
  const safeReactionId = JSON.stringify(reactionId);
  const safeBbSlot = JSON.stringify(bbSlot);
  const safeReagentCount = JSON.stringify(reagentCount);
  const safeSynthonId = JSON.stringify(synthonId);
  const safeReferenceSmiles = JSON.stringify(referenceSmiles || "");

  return `
    <button
      type="button"
      class="synthon-info-button"
      onclick='findSynthonInfo(${safePairId}, ${safeReactionId}, ${safeBbSlot}, ${safeReagentCount}, ${safeSynthonId}, ${safeReferenceSmiles}, this)'
    >
      Find more info
    </button>
    <div class="synthon-extra-results"></div>
  `;
}

function renderBuildingBlockRows(results, pairId, reactionId, bbSlot, reagentCount, visibleCount) {
  return results.slice(0, visibleCount).map((item) => {
    const blockKey = `block-${buildingBlockCounter}`;
    buildingBlockCounter += 1;
    buildingBlockStore[blockKey] = item;

    const smiles = item.smiles || "";
    const smilesJson = JSON.stringify(smiles);

    const selectedItem =
      selectedCombos[pairId] &&
      selectedCombos[pairId].reagents &&
      selectedCombos[pairId].reagents[bbSlot];

    const isSelected =
      selectedItem &&
      selectedItem.smiles === item.smiles &&
      selectedItem.bb_id === item.bb_id;

    const selectedClass = isSelected ? " is-selected" : "";
    const selectedText = isSelected ? "Selected" : "Select";

    return `
      <tr>
        <td>
          <button
            type="button"
            class="select-building-block-button${selectedClass}"
            onclick="selectBuildingBlock('${pairId}', '${reactionId}', '${bbSlot}', '${reagentCount}', buildingBlockStore['${blockKey}'], this)"
          >
            ${selectedText}
          </button>
        </td>
        <td>${item.tc_score || ""}</td>
        <td>${item.bb_id || ""}</td>
        <td class="smiles-text">${smiles}</td>
        <td>
          <button
            type="button"
            class="view-molecule-button"
            onclick='viewMolecule(${smilesJson}, this)'
          >
            View
          </button>
          <div class="synthon-extra-results"></div>
        </td>
      </tr>
    `;
  }).join("");
}

function loadMoreBuildingBlocks(button) {
  const container = button.closest(".synthon-extra-results");
  const tableBody = container.querySelector("tbody");
  const status = container.querySelector(".building-block-display-status");

  const key = button.dataset.blockListKey;
  const pairId = button.dataset.pairId;
  const reactionId = button.dataset.reactionId;
  const bbSlot = button.dataset.bbSlot;
  const reagentCount = button.dataset.reagentCount;

  const results = buildingBlockStore[key] || [];
  const currentCount = Number(button.dataset.visibleCount || buildingBlockPageSize);
  const nextCount = Math.min(currentCount + buildingBlockPageSize, results.length);

  tableBody.innerHTML = renderBuildingBlockRows(
    results,
    pairId,
    reactionId,
    bbSlot,
    reagentCount,
    nextCount
  );

  status.textContent = `Displaying ${nextCount}/${results.length}`;

  button.dataset.visibleCount = String(nextCount);

  if (nextCount >= results.length) {
    button.remove();
  }
}

function availabilityOptions(selectedAvailability) {
  const options = [
    ["priced", "Priced / In-stock"],
    ["fast", "FAST"],
    ["made", "MADE"],
  ];

  return options.map(([value, label]) => {
    const selected = value === selectedAvailability ? " selected" : "";
    return `<option value="${value}"${selected}>${label}</option>`;
  }).join("");
}

function priceCapOptions(selectedPriceCap) {
  const options = [
    ["db050", "$50 or less"],
    ["db100", "$100 or less"],
    ["db250", "$250 or less"],
    ["db500", "$500 or less"],
    ["db999", "$999 or less"],
  ];

  return options.map(([value, label]) => {
    const selected = value === selectedPriceCap ? " selected" : "";
    return `<option value="${value}"${selected}>${label}</option>`;
  }).join("");
}

function sourceFromSettings(settings) {
  if (settings.availability === "fast") {
    return "bbAFAS";
  }

  if (settings.availability === "made") {
    return "incrMADE";
  }

  return settings.price_cap || "db050";
}

function updateSynthonSetting(pairId, bbSlot, key, value) {
  const settingsKey = `${pairId}-${bbSlot}`;

  if (!synthonSearchSettings[settingsKey]) {
    synthonSearchSettings[settingsKey] = Object.assign({}, window.activeSearchSettings || defaultSearchSettings);
  }

  synthonSearchSettings[settingsKey][key] = value;
}

function updateSynthonSettingAndReload(pairId, reactionId, bbSlot, reagentCount, synthonId, referenceSmiles, key, value, control) {
  updateSynthonSetting(pairId, bbSlot, key, value);

  const settingsKey = `${pairId}-${bbSlot}`;

  if (synthonReloadTimers[settingsKey]) {
    clearTimeout(synthonReloadTimers[settingsKey]);
  }

  synthonReloadTimers[settingsKey] = setTimeout(() => {
    const container = control.closest(".synthon-extra-results");
    const button = container ? container.previousElementSibling : null;

    if (!button) {
      return;
    }

    container.innerHTML = "";
    button.textContent = "Find more info";
    findSynthonInfo(pairId, reactionId, bbSlot, reagentCount, synthonId, referenceSmiles, button);
  }, 350);
}

async function findSynthonInfo(pairId, reactionId, bbSlot, reagentCount, synthonId, referenceSmiles, button) {
  const container = button.nextElementSibling;

  if (container.innerHTML.trim()) {
    container.innerHTML = "";
    button.textContent = "Find more info";
    return;
  }

  button.textContent = "Loading...";
  button.disabled = true;
  container.textContent = "Loading building blocks...";

  try {
    const settingsKey = `${pairId}-${bbSlot}`;
    const activeSettings = synthonSearchSettings[settingsKey] || window.activeSearchSettings || defaultSearchSettings;
    const activeSource = sourceFromSettings(activeSettings);
    const activeTcCutoff = activeSettings.tc_cutoff || 0.6;
    const activeLimit = activeSettings.limit || 50;
    const activeFallbackTopMatch = activeSettings.fallback_top_match ? "true" : "false";

    const url =
      `/synthon/${encodeURIComponent(String(synthonId).trim())}` +
      `?source=${encodeURIComponent(activeSource)}` +
      `&limit=${encodeURIComponent(activeLimit)}` +
      `&reference_smiles=${encodeURIComponent(referenceSmiles || "")}` +
      `&tc_cutoff=${encodeURIComponent(activeTcCutoff)}` +
      `&fallback_top_match=${encodeURIComponent(activeFallbackTopMatch)}`;

    const response = await fetch(url);
    const data = await response.json();
    const results = Array.isArray(data.results) ? data.results : [];

    if (!response.ok) {
      throw new Error(data.error || "Server returned an error");
    }

    if (results.length === 0) {
      container.textContent = `No building blocks returned for synthon ${synthonId}.`;
      button.textContent = "Find more info";
      return;
    }

    if (!loadedBlockLists[pairId]) {
      loadedBlockLists[pairId] = {};
    }

    const sortedResults = results.slice().sort((a, b) => {
      const scoreA = Number(a.tc_score || 0);
      const scoreB = Number(b.tc_score || 0);
      return scoreB - scoreA;
    });

    loadedBlockLists[pairId][bbSlot] = sortedResults;

    const visibleCount = Math.min(buildingBlockPageSize, sortedResults.length);
    const blockListKey = `list-${buildingBlockCounter}`;
    buildingBlockCounter += 1;
    buildingBlockStore[blockListKey] = sortedResults;

    const tableRows = renderBuildingBlockRows(
      sortedResults,
      pairId,
      reactionId,
      bbSlot,
      reagentCount,
      visibleCount
    );

    const inclusionSmarts = Array.isArray(data.in_smarts)
      ? data.in_smarts.join(", ")
      : data.in_smarts || "";

    container.innerHTML = `
      <div class="synthon-settings">
        <label>
          Availability
          <select oninput="updateSynthonSettingAndReload('${pairId}', '${reactionId}', '${bbSlot}', '${reagentCount}', '${synthonId}', ${JSON.stringify(referenceSmiles || "")}, 'availability', this.value, this)">
            ${availabilityOptions(activeSettings.availability || "priced")}
          </select>
        </label>

        <label>
          Price cap
          <select oninput="updateSynthonSettingAndReload('${pairId}', '${reactionId}', '${bbSlot}', '${reagentCount}', '${synthonId}', ${JSON.stringify(referenceSmiles || "")}, 'price_cap', this.value, this)">
            ${priceCapOptions(activeSettings.price_cap || "db050")}
          </select>
        </label>

        <label>
          Tc cutoff
          <input
            type="number"
            min="0"
            max="1"
            step="0.05"
            value="${activeTcCutoff}"
            oninput="updateSynthonSettingAndReload('${pairId}', '${reactionId}', '${bbSlot}', '${reagentCount}', '${synthonId}', ${JSON.stringify(referenceSmiles || "")}, 'tc_cutoff', Number(this.value), this)"
          >
        </label>

        <label>
          Max blocks
          <input
            type="number"
            min="1"
            max="500"
            step="1"
            value="${activeLimit}"
            oninput="updateSynthonSettingAndReload('${pairId}', '${reactionId}', '${bbSlot}', '${reagentCount}', '${synthonId}', ${JSON.stringify(referenceSmiles || "")}, 'limit', Number(this.value), this)"
          >
        </label>

        <label class="checkbox-label">
          <input
            type="checkbox"
            ${activeSettings.fallback_top_match ? "checked" : ""}
            oninput="updateSynthonSettingAndReload('${pairId}', '${reactionId}', '${bbSlot}', '${reagentCount}', '${synthonId}', ${JSON.stringify(referenceSmiles || "")}, 'fallback_top_match', this.checked, this)"
          >
          Use top match if no Tc hits
        </label>
      </div>

      ${data.message ? `<div class="synthon-warning">${data.message}</div>` : ""}

      <div class="synthon-text">
        <strong>Inclusion SMARTS:</strong> ${inclusionSmarts || "Not available"}
      </div>

      <div class="synthon-text building-block-display-status">
        Displaying ${visibleCount}/${sortedResults.length}
      </div>

      <table class="synthon-results-table">
        <thead>
          <tr>
            <th>Enumerate</th>
            <th>Tc Score</th>
            <th>ZINC ID</th>
            <th>SMILES</th>
            <th>Structure</th>
          </tr>
        </thead>
        <tbody>
          ${tableRows}
        </tbody>
      </table>

      ${
        visibleCount < sortedResults.length
          ? `
            <button
              type="button"
              class="synthon-info-button"
              data-block-list-key="${blockListKey}"
              data-pair-id="${pairId}"
              data-reaction-id="${reactionId}"
              data-bb-slot="${bbSlot}"
              data-reagent-count="${reagentCount}"
              data-visible-count="${visibleCount}"
              onclick="loadMoreBuildingBlocks(this)"
            >
              Load more
            </button>
          `
          : ""
      }
    `;

    button.textContent = "Show less info";
  } catch (error) {
    container.textContent = `Could not load building blocks: ${error.message}`;
    button.textContent = "Find more info";
  } finally {
    button.disabled = false;
  }
}

function viewMolecule(smiles, button) {
  const container = button.nextElementSibling;

  if (container.innerHTML.trim()) {
    container.innerHTML = "";
    button.textContent = "View";
    return;
  }

  button.textContent = "Hide";

  container.innerHTML = `
    <img
      class="mol-image-small"
      src="/draw?smiles=${encodeURIComponent(smiles)}"
      alt="Molecule structure"
      onclick="openImageModal('/draw?smiles=${encodeURIComponent(smiles)}', 'Molecule structure')"

    >
  `;
}

function selectBuildingBlock(pairId, reactionId, bbSlot, reagentCount, item, button) {
  const wasSelected = button.classList.contains("is-selected");
  const table = button.closest(".synthon-results-table");

  if (!selectedCombos[pairId]) {
    selectedCombos[pairId] = {
      reactionId: reactionId,
      reagentCount: Number(reagentCount),
      reagents: {},
    };
  }

  if (table) {
    table.querySelectorAll(".select-building-block-button").forEach((otherButton) => {
      otherButton.classList.remove("is-selected");
      otherButton.textContent = "Select";
    });
  }

  if (wasSelected) {
    delete selectedCombos[pairId].reagents[bbSlot];
    updateComboStatus(pairId);
    return;
  }

  selectedCombos[pairId].reagents[bbSlot] = item;
  button.classList.add("is-selected");
  button.textContent = "Selected";
  updateComboStatus(pairId);
}

function updateComboStatus(pairId) {
  const combo = selectedCombos[pairId];

  if (!combo) {
    return;
  }

  const selectedCount = Object.keys(combo.reagents).length;
  const status = document.getElementById(`combo-status-${pairId}`);

  if (!status) {
    return;
  }

  if (selectedCount < combo.reagentCount) {
    status.textContent = `Selected ${selectedCount} of ${combo.reagentCount} building blocks.`;
    return;
  }

  runEnumeration(pairId, status);
}

async function runEnumeration(pairId, status) {
  const combo = selectedCombos[pairId];

  if (!combo) {
    return;
  }

  status.innerHTML = `
    <div><strong>Auto-enumerating selected combo...</strong></div>
  `;

  try {
    const response = await fetch("/enumerate", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        reaction_id: combo.reactionId,
        pair_id: pairId,
        reagents: combo.reagents,
      }),
    });

    const data = await response.json();

    if (!response.ok) {
      throw new Error(data.error || "Enumeration failed");
    }

    status.innerHTML = productsToTable("Manual enumeration products", data.products);
  } catch (error) {
    status.innerHTML = `
      <div><strong>Enumeration failed.</strong></div>
      <div class="smiles-text">${error.message}</div>
    `;
  }
}

async function runEnumerateAllShown(pairId, reactionId, button) {
  const container = button.nextElementSibling;

  if (container.innerHTML.trim()) {
    container.innerHTML = "";
    button.textContent = "Enumerate All Shown Blocks";
    return;
  }

  const blockLists = loadedBlockLists[pairId] || {};
  const loadedSlots = Object.keys(blockLists).filter((slot) => {
    return Array.isArray(blockLists[slot]) && blockLists[slot].length > 0;
  });

  console.log("enumerate-all blockLists", pairId, blockLists);
  console.log("enumerate-all loadedSlots", loadedSlots);

  if (loadedSlots.length < 2) {
    container.innerHTML = `
      <div><strong>Open "Find more info" for at least two synthons first.</strong></div>
      <div class="smiles-text">
        Loaded slots: ${loadedSlots.length ? loadedSlots.join(", ") : "none"}
      </div>
    `;
    return;
  }

  button.textContent = "Hide Enumerated Results";
  container.innerHTML = `
    <div><strong>Enumerating all shown building-block combinations...</strong></div>
  `;

  try {
    const response = await fetch("/enumerate-all", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify({
        reaction_id: reactionId,
        pair_id: pairId,
        block_lists: Object.fromEntries(loadedSlots.map((slot) => [slot, blockLists[slot]])),
        max_combos: 200,
      }),
    });

    const rawText = await response.text();
    let data;

    try {
      data = JSON.parse(rawText);
    } catch (error) {
      throw new Error(rawText || "Server returned non-JSON error");
    }

    if (!response.ok) {
      throw new Error(data.error || "Enumerate all failed");
    }

    container.innerHTML = `
      ${productsToTable("Enumerated products from shown blocks", data.products)}
      <div class="synthon-text">Enumerated ${data.total_combos || 0} possible combinations.</div>
    `;
  } catch (error) {
    container.innerHTML = `
      <div><strong>Enumerate all failed.</strong></div>
      <div class="smiles-text">${error.message}</div>
    `;
  }
}

function outputToTable(output, mode = "full", countData = null) {
  const lines = output
    .trim()
    .split("\n")
    .filter((line) => line && !line.startsWith("#"));

  if (lines.length < 2) {
    return "<p>No matching retro reactions found.</p>";
  }

  const headers = lines[0].split("\t");
  const rows = lines.slice(1).map((line) => line.split("\t"));

  const firstRow = rows[0];
  const originalSmiles = firstRow[headers.indexOf("original_smiles")];
  const originalSvg = firstRow[headers.indexOf("original_svg")];
  const originalUrl = svgUrl(originalSvg);

  const originalHtml = originalUrl
    ? `
      <div class="original-molecule">
        <h3>Submitted Molecule</h3>
        <img class="mol-image-large" src="${originalUrl}" alt="Submitted molecule structure" onclick="openImageModal('${originalUrl}', 'Submitted molecule structure')">
        <div class="smiles-text">${originalSmiles}</div>
      </div>
    `
    : "";

  const visibleHeaders = [
    "Reaction Info",
    "Synthon 1",
    "Synthon 2",
    "Synthon 3",
    "Synthon 4",
  ];

  const resultTables = rows.map((row) => {
    const get = (name) => {
      const index = headers.indexOf(name);
      return index === -1 ? "" : row[index] || "";
    };

    const synthonIds = [
      get("bb1_synthon_id"),
      get("bb2_synthon_id"),
      get("bb3_synthon_id"),
      get("bb4_synthon_id"),
    ].filter(Boolean);

    return `
      <div class="results-table-wrap">
        <table class="results-table">
          <thead>
            <tr>
              ${visibleHeaders.map((header) => `<th>${header}</th>`).join("")}
            </tr>
          </thead>
          <tbody>
            <tr>
              <td class="reaction-info">
                <div><span class="reaction-label">ID:</span> ${get("reaction_id")}</div>
                <div><span class="reaction-label">Name:</span> ${get("reaction_name")}</div>
                <div><span class="reaction-label">Reagents:</span> ${get("reagent_count")}</div>
                ${reactionAvailabilityHtml(countData, synthonIds)}
                ${
                  mode === "overview"
                    ? `
                      <a
                        class="button-link"
                        href="/reaction/${document.body.dataset.jobId}/${encodeURIComponent(get("reaction_id"))}"
                      >
                        Open Reaction
                      </a>

                      <a
                        class="button-link"
                        href="/coverage?parent_smiles=${encodeURIComponent(originalSmiles)}&job_id=${encodeURIComponent(document.body.dataset.jobId || "")}"
                      >
                        Run Analog Coverage
                      </a>
                    `
                    : ""
                }

                <details class="reaction-smarts">
                  <summary>Reaction SMARTS</summary>
                  <div>
                    <span class="reaction-label">Build:</span>
                    <code>${get("build_smarts")}</code>
                  </div>
                  <div>
                    <span class="reaction-label">Breakdown:</span>
                    <code>${get("breakdown_smarts")}</code>
                  </div>
                </details>

                ${
                  mode === "full"
                    ? `
                      <div id="combo-status-${get("pair_id")}" class="synthon-text">
                        Select building blocks to enumerate.
                      </div>

                      <button
                        type="button"
                        class="synthon-info-button"
                        onclick="runEnumerateAllShown('${get("pair_id")}', '${get("reaction_id")}', this)"
                      >
                        Enumerate All Shown Blocks
                      </button>

                      <div class="synthon-extra-results"></div>
                    `
                    : ""
                }
              </td>

              <td>
                ${moleculeCell(get("bb1_smiles"), get("bb1_svg"))}
                <div class="synthon-text">Synthon ID: ${get("bb1_synthon_id")}</div>
                <div class="synthon-text">Name: ${get("bb1_synthon_name")}</div>
                ${renderSingleSynthonCounts(countData, get("bb1_synthon_id"))}
                ${mode === "full" ? synthonInfoButton(get("pair_id"), get("reaction_id"), "bb1", get("reagent_count"), get("bb1_synthon_id"), get("bb1_smiles")) : ""}
              </td>

              <td>
                ${moleculeCell(get("bb2_smiles"), get("bb2_svg"))}
                <div class="synthon-text">Synthon ID: ${get("bb2_synthon_id")}</div>
                <div class="synthon-text">Name: ${get("bb2_synthon_name")}</div>
                ${renderSingleSynthonCounts(countData, get("bb2_synthon_id"))}
                ${mode === "full" ? synthonInfoButton(get("pair_id"), get("reaction_id"), "bb2", get("reagent_count"), get("bb2_synthon_id"), get("bb2_smiles")) : ""}
              </td>

              <td>
                ${moleculeCell(get("bb3_smiles"), get("bb3_svg"))}
                <div class="synthon-text">
                  ${get("bb3_synthon_id") ? `Synthon ID: ${get("bb3_synthon_id")}` : ""}
                </div>
                <div class="synthon-text">
                  ${get("bb3_synthon_name") ? `Name: ${get("bb3_synthon_name")}` : ""}
                </div>
                ${renderSingleSynthonCounts(countData, get("bb3_synthon_id"))}
                ${mode === "full" ? synthonInfoButton(get("pair_id"), get("reaction_id"), "bb3", get("reagent_count"), get("bb3_synthon_id"), get("bb3_smiles")) : ""}
              </td>

              <td>
                ${moleculeCell(get("bb4_smiles"), get("bb4_svg"))}
                <div class="synthon-text">
                  ${get("bb4_synthon_id") ? `Synthon ID: ${get("bb4_synthon_id")}` : ""}
                </div>
                <div class="synthon-text">
                  ${get("bb4_synthon_name") ? `Name: ${get("bb4_synthon_name")}` : ""}
                </div>
                ${renderSingleSynthonCounts(countData, get("bb4_synthon_id"))}
                ${mode === "full" ? synthonInfoButton(get("pair_id"), get("reaction_id"), "bb4", get("reagent_count"), get("bb4_synthon_id"), get("bb4_smiles")) : ""}
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    `;
  }).join("");

  return `
    ${originalHtml}
    ${resultTables}
  `;
}

function openImageModal(src, altText) {
  const existingModal = document.querySelector(".image-modal");

  if (existingModal) {
    existingModal.remove();
  }

  const modal = document.createElement("div");
  modal.className = "image-modal";

  modal.innerHTML = `
    <div class="image-modal-content">
      <img src="${src}" alt="${altText || "Molecule structure"}">
    </div>
  `;

  modal.addEventListener("click", () => {
    modal.remove();
  });

  document.body.appendChild(modal);
}

const pageJobId = document.body.dataset.jobId;
const pageReactionId = document.body.dataset.reactionId;

function formatCountValue(value) {
  if (value === null || value === undefined) {
    return "pending";
  }

  return Number(value).toLocaleString();
}

function renderSingleSynthonCounts(countData, synthonId) {
  const rows = countData && countData.counts
    ? countData.counts[String(synthonId)]
    : null;

  if (!rows) {
    return "";
  }

  return `
    <table class="synthon-count-table compact">
      <thead>
        <tr>
          <th>Tier</th>
          <th>New</th>
          <th>Total</th>
        </tr>
      </thead>
      <tbody>
        ${rows.map((row) => `
          <tr>
            <td>${row.label}</td>
            <td>${row.new_source ? formatCountValue(row.new_count) : "-"}</td>
            <td>${row.total_source ? formatCountValue(row.total_count) : "-"}</td>
          </tr>
        `).join("")}
      </tbody>
    </table>
  `;
}

function reactionAvailabilityHtml(countData, synthonIds) {
  const ids = synthonIds.filter(Boolean);

  if (!countData || ids.length === 0) {
    return `<div class="availability-signal availability-unknown">Availability unknown</div>`;
  }

  const everyHasDb050 = ids.every((id) =>
    (countData.counts[String(id)] || []).some((row) =>
      row.total_source === "db050" && row.total_count > 0
    )
  );

  const everyHasDb999 = ids.every((id) =>
    (countData.counts[String(id)] || []).some((row) =>
      row.total_source === "db999" && row.total_count > 0
    )
  );

  if (everyHasDb050) {
    return `<div class="availability-signal availability-green">Easy purchase</div>`;
  }

  if (everyHasDb999) {
    return `<div class="availability-signal availability-yellow">Purchasable</div>`;
  }

  return `<div class="availability-signal availability-red">On-demand needed</div>`;
}

function renderSynthonCounts(countData) {
  if (!countData || !countData.counts) {
    return "";
  }

  return `
    <div class="synthon-counts-panel">
      <h3>Building Block Availability</h3>
      <div class="synthon-counts-grid">
        ${Object.entries(countData.counts).map(([synthonId, rows]) => `
          <div class="synthon-count-card">
            <h4>Synthon ${synthonId}</h4>
            <table class="synthon-count-table">
              <thead>
                <tr>
                  <th>Availability</th>
                  <th>New</th>
                  <th>Total</th>
                </tr>
              </thead>
              <tbody>
                ${rows.map((row) => `
                  <tr>
                    <td>${row.label}</td>
                    <td>${row.new_source ? formatCountValue(row.new_count) : "-"}</td>
                    <td>${row.total_source ? formatCountValue(row.total_count) : "-"}</td>
                  </tr>
                `).join("")}
              </tbody>
            </table>
          </div>
        `).join("")}
      </div>
    </div>
  `;
}

async function loadSynthonCounts(jobId) {
  try {
    const response = await fetch(`/api/synthon-counts/${jobId}`);
    const data = await response.json();

    if (!response.ok || data.error) {
      return "";
    }

    return data;
  } catch (error) {
    return "";
  }
}

if (pageJobId && !pageReactionId) {
  fetch(`/api/synthons/${pageJobId}`)
    .then((response) => response.json())
    .then((data) => {
      if (data.error) {
        result.textContent = data.error;
        return;
      }

      result.innerHTML = `
        <div class="status-line">Job ID: ${data.job_id}</div>
        ${outputToTable(data.output, "overview", "")}
        <a class="button-link" href="/download/${data.job_id}" download>
          Download output.zip
        </a>
      `;

      loadSynthonCounts(data.job_id).then((countsHtml) => {
        result.innerHTML = `
          <div class="status-line">Job ID: ${data.job_id}</div>
          ${outputToTable(data.output, "overview", countsHtml)}
          <a class="button-link" href="/download/${data.job_id}" download>
            Download output.zip
          </a>
        `;
      });
    })
    .catch((error) => {
      result.textContent = `Could not load synthon results: ${error.message}`;
    });
}

if (pageJobId && pageReactionId) {
  fetch(`/api/reaction/${pageJobId}/${encodeURIComponent(pageReactionId)}`)
    .then((response) => response.json())
    .then((data) => {
      if (data.error) {
        result.textContent = data.error;
        return;
      }

      const headers = Object.keys(data.rows[0]);
      const output = [
        headers.join("\t"),
        ...data.rows.map((row) => headers.map((header) => row[header] || "").join("\t")),
      ].join("\n");

      loadSynthonCounts(pageJobId)
        .then((countsHtml) => {
          result.innerHTML = outputToTable(output, "full", countsHtml);
        })
        .catch(() => {
          result.innerHTML = outputToTable(output, "full", "");
        });
    })
    .catch((error) => {
      result.textContent = `Could not load reaction detail: ${error.message}`;
    });
}

if (form) {
form.addEventListener("submit", async (event) => {
  event.preventDefault();

  const userInput = document.getElementById("user-input").value;
  const email = document.getElementById("email").value;
  const source = document.getElementById("availability-source").value;
  const tcCutoff = Number(document.getElementById("tc-cutoff").value);
  const sortMode = document.getElementById("sort-mode").value;
  const submitButton = document.getElementById("submit-button");

  submitButton.classList.add("is-active");
  submitButton.textContent = "Running...";
  submitButton.disabled = true;

  result.textContent = "Running job...";

  try {
    const response = await fetch("/reaction", {
      method: "POST",
      headers: {
        "Content-Type": "application/json"
      },
      body: JSON.stringify({
        user_input: userInput,
        email: email,
        source: source,
        tc_cutoff: tcCutoff,
        sort_mode: sortMode
      })
    });

    if (!response.ok) {
      throw new Error("Server returned an error");
    }

    const data = await response.json();

    if (data.job_id) {
      window.location.href = `/synthons/${data.job_id}`;
      return;
    }

    result.textContent = "Job finished, but no job ID was returned.";
    submitButton.textContent = "Submit Again";
    submitButton.disabled = false;

    submitButton.textContent = "Submit Again";
    submitButton.disabled = false;
  } catch (error) {
    result.textContent = "Could not connect to the backend server. Make sure it is running.";

    submitButton.classList.remove("is-active");
    submitButton.textContent = "Submit Job";
    submitButton.disabled = false;
  }
});
}


const pastJobForm = document.getElementById("past-job-form");

if (pastJobForm) {
  pastJobForm.addEventListener("submit", (event) => {
    event.preventDefault();

    const input = document.getElementById("past-job-id");
    const jobId = input.value.trim();

    if (!jobId) {
      input.focus();
      return;
    }

    fetch(`/api/past-job/${encodeURIComponent(jobId)}?t=${Date.now()}`)
      .then((response) => response.json())
      .then((data) => {
        if (data.error) {
          alert(data.error);
          return;
        }

        window.open(data.url, "_blank");
      })
      .catch(() => {
        alert("Could not look up that job.");
      });
  });
}

const madeForm = document.getElementById("made-job-form");
const madeResult = document.getElementById("made-result");

function madeWait(milliseconds) {
  return new Promise((resolve) => {
    window.setTimeout(resolve, milliseconds);
  });
}

function madeDownloadLink(files, key, label) {
  const path = files && files[key];

  if (!path) {
    return "";
  }

  return `
    <a
      class="button-link made-download-link"
      href="${path}"
      download
    >
      ${label}
    </a>
  `;
}

function renderMadeResult(state) {
  const workflow = state.result || {};
  const files = state.files || {};
  const routes = workflow.routes || [];

  const routeRows = routes.length
    ? routes.map((route) => `
        <tr>
          <td>${route.reaction_id || ""}</td>
          <td>${route.status || ""}</td>
          <td>${route.reason || ""}</td>
        </tr>
      `).join("")
    : `
        <tr>
          <td colspan="3">No validated routes were reported.</td>
        </tr>
      `;

  madeResult.innerHTML = `
    <div class="made-complete-message">
      MADE workflow complete. Job ID:
      <code>${state.job_id}</code>
    </div>

    <div class="made-summary-grid">
      <div>
        <strong>${workflow.retro_route_count || 0}</strong>
        <span>Retro routes</span>
      </div>
      <div>
        <strong>${workflow.validated_route_count || 0}</strong>
        <span>Validated routes</span>
      </div>
      <div>
        <strong>${workflow.filtered_candidate_count || 0}</strong>
        <span>Filtered candidates</span>
      </div>
      <div>
        <strong>${workflow.product_count || 0}</strong>
        <span>Products</span>
      </div>
      <div>
        <strong>${workflow.match_count || 0}</strong>
        <span>InChIKey matches</span>
      </div>
    </div>

    <div class="made-downloads">
      ${madeDownloadLink(
        files,
        "made_all",
        "Download MADE source manifest"
      )}
      ${madeDownloadLink(
        files,
        "made_filtered",
        "Download filtered MADE candidates"
      )}
      ${madeDownloadLink(
        files,
        "made_products",
        "Download enumerated products"
      )}
      ${madeDownloadLink(
        files,
        "made_matches",
        "Download InChIKey matches"
      )}
    </div>

    <div class="results-table-wrap">
      <table class="results-table made-route-table">
        <thead>
          <tr>
            <th>Reaction</th>
            <th>Status</th>
            <th>Details</th>
          </tr>
        </thead>
        <tbody>
          ${routeRows}
        </tbody>
      </table>
    </div>
  `;
}

async function pollMadeJob(statusUrl) {
  for (let check = 0; check < 360; check += 1) {
    const response = await fetch(statusUrl, {
      cache: "no-store"
    });

    if (!response.ok) {
      throw new Error("Could not read MADE job status");
    }

    const state = await response.json();

    if (state.status === "complete") {
      return state;
    }

    if (state.status === "failed") {
      throw new Error(
        state.error || "The MADE workflow failed"
      );
    }

    madeResult.textContent =
      `MADE workflow ${state.status || "running"}...`;

    await madeWait(5000);
  }

  throw new Error(
    "The MADE job is still running. Save the job ID and check again later."
  );
}

if (madeForm) {
  const originalInput = document.getElementById("user-input");
  const madeInput = document.getElementById("made-smiles");

  if (
    originalInput &&
    madeInput &&
    originalInput.value &&
    !madeInput.value
  ) {
    madeInput.value = originalInput.value;
  }

  madeForm.addEventListener("submit", async (event) => {
    event.preventDefault();

    const smiles = madeInput.value.trim();
    const filterMode = document.getElementById(
      "made-filter-mode"
    ).value;
    const button = document.getElementById(
      "made-submit-button"
    );

    button.disabled = true;
    button.classList.add("is-running");
    button.textContent = "Starting MADE workflow...";
    madeResult.textContent = "Submitting MADE workflow...";

    try {
      const response = await fetch("/api/made/jobs", {
        method: "POST",
        headers: {
          "Content-Type": "application/json"
        },
        body: JSON.stringify({
          smiles: smiles,
          filter_mode: filterMode
        })
      });

      const submission = await response.json();

      if (!response.ok) {
        throw new Error(
          submission.detail ||
          submission.error ||
          "The server rejected the MADE job"
        );
      }

      if (!submission.status_url) {
        throw new Error(
          "The server did not return a MADE status URL"
        );
      }

      madeResult.textContent =
        `MADE job ${submission.job_id} is queued...`;
      button.textContent = "MADE workflow running...";

      const state = await pollMadeJob(
        submission.status_url
      );

      renderMadeResult(state);
      button.textContent = "Run MADE Workflow Again";
    } catch (error) {
      madeResult.textContent =
        `MADE workflow failed: ${error.message}`;
      button.textContent = "Retry MADE Workflow";
    } finally {
      button.disabled = false;
      button.classList.remove("is-running");
      button.classList.remove("is-active");
    }
  });
}

document.addEventListener("click", (event) => {
  const clickedButton = event.target.closest("button, .button-link");

  if (!clickedButton) {
    return;
  }

  if (
    clickedButton.id === "submit-button" ||
    clickedButton.id === "made-submit-button"
  ) {
    return;
  }

  clickedButton.classList.toggle("is-active");
});
