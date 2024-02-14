let index = 0;
function addSequence(sequence) {
  if (index == 0) {
    let tableHead = document.createElement("tr");
    document.getElementById("sequence-table-head").innerHTML = "";
    tableHead.innerHTML = `
        <th>No.</th>
        <th>Sequence</th>
        <th>Start codon index</th>
        <th>Stop codon index</th>
        <th>Status</th>
        <th>Remove</th>
        `;
    document.getElementById("sequence-table-head").appendChild(tableHead);
  }
  index += 1;
  let element = document.createElement("tr");
  element.innerHTML = `
    <tr id="sequence${index}">
    <td>${index}</td>
    <td>
        <input type="text" name="mrna" id="mrna" placeholder="mRNA Sequence" value="${sequence}" >
    </td>
    <td>
        <input type="number" name="start-codon" id="startcodon${index}" placeholder="0" required >
    </td>
    <td>
        <input type="number" name="stop-codon" id="stopcodon${index}" placeholder="0" required >
    </td>
    <td class="text-center">
    </td>
    <td class="text-center" style="cursor: pointer" onclick="removeSequence('sequence${index}')">
    ❌
    </td>
</tr>`;
  document.getElementById("sequences-table").appendChild(element);
}

function removeSequence(id) {
  let removalIndex = Number(id.replace("sequence", ""));
  document
    .getElementById("sequences-table")
    .querySelectorAll("tr")
    .forEach((row) => {
      let cells = row.querySelectorAll("td");
      let id = Number(cells[0].innerText);
      if (id == removalIndex) {
        row.remove();
      } else if (id > removalIndex) {
        cells[0].innerText = id - 1;
        cells[2].id = `startcodon${id - 1}`;
        cells[3].id = `stopcodon${id - 1}`;
        cells[5].setAttribute("onclick", `removeSequence('sequence${id - 1}')`);
      }
    });
  index -= 1;
}

document.getElementById("sequence-adder").addEventListener("click", (_) => {
  addSequence("");
});

document.getElementById("file").addEventListener("change", (_) => {
  let file = document.getElementById("file").files[0];
  let reader = new FileReader();
  reader.readAsText(file, "UTF-8");
  reader.onload = function (e) {
    let data = e.target.result.trim();
    let sequences = [];
    let lines = data.split("\n");
    for (let i = 0; i < lines.length; i++) {
      if (lines[i].startsWith(">")) {
        sequences.push("");
      } else {
        sequences[sequences.length - 1] += lines[i];
      }
    }
    for (let j = 0; j < sequences.length; j++) {
      addSequence(sequences[j]);
    }
  };
});

function is_valid_sequence(gene_sequence, start_codon_index, stop_codon_index) {
  const start_codon = "AUG";
  const stop_codons = ["UAA", "UAG", "UGA"];

  if (
    gene_sequence.slice(start_codon_index - 1, start_codon_index + 2) ===
      start_codon &&
    stop_codons.includes(
      gene_sequence.slice(stop_codon_index - 1, stop_codon_index + 2)
    )
  ) {
    const cds_sequence = gene_sequence.slice(
      start_codon_index + 2,
      stop_codon_index - 1
    );

    if (cds_sequence.length % 3 !== 0) {
      return "Please make sure coding sequence is in triplets";
    }

    return "✅";
  } else {
    return "Invalid start or stop codon index.";
  }
}

document.getElementById("submit-btn").addEventListener("click", (_) => {
  let sequencesData = [
    ...document.getElementById("sequences-table").querySelectorAll("tr"),
  ].map((tr) => {
    let cells = tr.cells;
    if (cells[1].firstElementChild.value == "") {
      removeSequence(`sequence${cells[0].innerText}`);
      return undefined;
    } else {
      cells[4].innerText = is_valid_sequence(
        cells[1].firstElementChild.value,
        Number(cells[2].firstElementChild.value),
        Number(cells[3].firstElementChild.value)
      );
      return [
        cells[1].firstElementChild.value,
        Number(cells[2].firstElementChild.value),
        Number(cells[3].firstElementChild.value),
      ];
    }
  });
  sequencesData = sequencesData.filter((data) => {
    if (data != undefined) {
      return data;
    }
  });
  fetch("/initiation_rate_prediction", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify(sequencesData),
  })
    .then((response) => response.json())
    .then((data) => {
      let outputBody = document.getElementById("output-body");
      outputBody.innerHTML = "";
      for (let i = 0; i < data.length; i++) {
        if (i == 0) {
          let outputHead = document.createElement("tr");
          document.getElementById("output-head").innerHTML = "";
          outputHead.innerHTML = `
          <tr>
            <th>No.</th>
            <th>N1</th>
            <th>N4</th>
            <th>Folding Energy (70 nucleotides)</th>
            <th>Folding Energy (80 nucleotides)</th>
            <th>CDS Length</th>
            <th>In frame AUG</th>
            <th>Kozak Score</th>
            <th>Length of 5' UTR</th>
            <th>Initiation Rate</th>
            </tr>`;
          document.getElementById("output-head").appendChild(outputHead);
        }
        let k = data[i];
        let childCode = `
            <tr>
            <td>${i + 1}</td>
            <td>${k.N1}</td>
            <td>${k.N4}</td>
            <td>${k.folding_energy_70}</td>
            <td>${k.folding_energy_80}</td>
            <td>${k.gene_length}</td>
            <td>${k["in_frame AUG"]}</td>
            <td>${k.kozak_score}</td>
            <td>${k.length_of_5prime_utr}</td>
            <td>${k.initiation_rate}</td>
        </tr>`;
        let child = document.createElement("tr");
        child.innerHTML = childCode;
        outputBody.appendChild(child);
      }
    });
});

document.getElementById("optimize").addEventListener("click", (e) => {
  let targetI = Number(document.getElementById("target-i").value);
  let iterations = Number(document.getElementById("iterations").value);
  let method = Number(document.getElementById("method").value);
  console.log(method);
  // let listOfSequences = [
  //   ...document.getElementById("sequences-table").querySelectorAll("tr"),
  // ].map((tr) => {
  //   let cells = tr.cells;
  //   return {
  //     value: cells[1].firstElementChild.value,
  //     start_codon_index: Number(cells[2].firstElementChild.value),
  //     stop_codon_index: Number(cells[3].firstElementChild.value),
  //     five_prime_utr: Number(cells[2].firstElementChild.value) - 1,
  //   };
  // });
  let listOfSequences = [{
    value: document.getElementById("op-sequence").value,
    start_codon_index: Number(document.getElementById("start_codon_index").value),
    stop_codon_index: Number(document.getElementById("stop_codon_index").value),
    five_prime_utr: Number(document.getElementById("start_codon_index").value)-1,
  }];
  document.getElementById("loading").style.display = "";
  fetch("/optimize", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
    },
    body: JSON.stringify({
      targetI: targetI,
      iterations: iterations,
      method: method,
      sequences: listOfSequences,
      kbt_values: [...document.querySelectorAll(".kbt-value")].map(element => Number(element.value))

    }),
  })
    .then((response) => {
      return response.json();
    })
    .then((data) => {
      for (i = 0; i < data.length; i++) {
        let outputBody = document.getElementById("optimized-output");
        if (i == 0) {
          let geneHead = document.createElement("tr");
          document.getElementById("gene-head").innerHTML = "";
          geneHead.innerHTML = `
            <th>Sr. No</th>
            <th>Target</th>
            <th>I
            </th>
            <th class="gene-td">Gene</th>`;
          document.getElementById("gene-head").appendChild(geneHead);
        }
        outputBody.innerHTML = "";
        for (let i = 0; i < data.length; i++) {
          let k = data[i];
          let childCode = `
                <tr>
                <td>${i + 1}</td>
                <td>${k.tir}</td>
                <td>${k.I}</td>
                <td style="overflow-x: scroll; max-width: 250px">${k.gene}</td>
            </tr>`;
          let child = document.createElement("tr");
          child.innerHTML = childCode;
          outputBody.appendChild(child);
        }
      }
      document.getElementById("loading").style.display = "none";
    });
});

function tableToCSV(id, title) {
  let csv_data = [];
  let rows = document.getElementById(id).getElementsByTagName("tr");
  for (let i = 0; i < rows.length; i++) {
    let cols = rows[i].querySelectorAll("td,th");
    let csvrow = [];
    for (let j = 0; j < cols.length; j++) {
      csvrow.push(cols[j].innerText);
    }
    csv_data.push(csvrow.join(","));
  }
  csv_data = csv_data.join("\n");
  CSVFile = new Blob([csv_data], { type: "text/csv" });
  let temp_link = document.createElement("a");
  temp_link.download = `${title} ${new Date()
    .toISOString()
    .replace(/[-T:.Z]/g, "")
    .replace(
      /(\d{4})(\d{2})(\d{2})(\d{2})(\d{2})(\d{2})/,
      "$1-$3-$2 $4-$5-$6"
    )}.csv`;
  let url = window.URL.createObjectURL(CSVFile);
  temp_link.href = url;
  temp_link.style.display = "none";
  document.body.appendChild(temp_link);
  temp_link.click();
  document.body.removeChild(temp_link);
}
