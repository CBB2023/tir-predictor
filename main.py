import InitiationRate
import OptimizationUTRcodon
import OptimizationUTR
import warnings
from flask import Flask, render_template, request, send_from_directory
warnings.filterwarnings("ignore")

app = Flask(__name__, static_folder='static')


@app.get("/")
def index():
    return render_template("index.html")


@app.get("/tir-predictor")
def tir_predictor():
    return render_template("tir-predictor.html")


@app.post("/initiation_rate_prediction")
def initiation_rate_prediction():
    json = request.get_json()
    response = []
    for data in json:
        sequence = data[0]
        start_codon_index = data[1]
        stop_codon_index = data[2]
        output = InitiationRate.InitiationRate(
            sequence, start_codon_index, stop_codon_index)
        gene_info = output[0]
        gene_info["initiation_rate"] = output[1]
        response.append(gene_info)
    return response


@app.post("/optimize")
def optimize():
    json = request.get_json()
    # print(json)
    iterations = json["iterations"]
    targeti = json["targetI"]
    method = json["method"]
    final = []
    for sequence in json["sequences"]:
        if method == 1:
            df = OptimizationUTR.OptimizationUTR(
                sequence["value"],
                int(sequence["start_codon_index"]),
                int(sequence["stop_codon_index"]),
                targeti,
                iterations,
                kbt_values=json["kbt_values"]
            )
        elif method == 2:
            df = OptimizationUTRcodon.OptimizationUTRcodon(
                sequence["value"],
                sequence["value"][:sequence["five_prime_utr"]],
                int(sequence["start_codon_index"]),
                int(sequence["stop_codon_index"]),
                targeti,
                iterations,
                kbt_values=json["kbt_values"]
            )
        stuff = df.tail(1)
        final.append({
            "tir": stuff.loc[iterations-1, 'tir'],
            "I": stuff.loc[iterations-1, 'I'],
            "gene": stuff.loc[iterations-1, 'gene']
        })
    return final

@app.route('/download/<path:filename>', methods=['GET', 'POST'])
def download(filename):
    # uploads = os.path.join(current_app.root_path, app.config['UPLOAD_FOLDER'])
    return send_from_directory(app.static_folder, filename, as_attachment=True)

if __name__ == "__main__":
    app.run(
        debug=True,
        port=5000,
        host="0.0.0.0"
    )
