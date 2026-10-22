"""
Flask app: takes an NCBI accession or a pasted protein sequence and reports
the AL???LW motifs that fall inside a hydrophobic Kyte-Doolittle window
(window 19, threshold 1.6 unless changed in the settings).
"""

import os
import re
from flask import Flask, redirect, render_template, request, send_file, url_for

import analysis as an
import export
import ncbi

# Vercel only serves static files from public/, so they live in public/static.
# Locally Flask serves the same folder, so /static/... works in both places.
app = Flask(__name__, static_folder="public/static", static_url_path="/static")
app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024  # 2 MB of form data is plenty

# Protein from the example output in the problem statement (481 aa).
# Has ALHAGLW at 297-303, which is what the example plot shows.
EXAMPLE_SEQUENCE = (
    "MPELSPRRRLLVLAICCLSLLIVSLDNTVLNVALPSMQKELHASVAGLQWTIDAYTLVLASLLMLAGSTA"
    "DRIGRRRVFRTGLVLFTIGSVLCSLAPNLESLVAFRMVQAVGGSMLNPVAMSIITNTFTEPRERARAIGV"
    "WGGVVGISMAAGPIVGGLLVDSVGWRSIFWINLPVGLVALLLTWRYVPESRAPKARRPDPVGQLLVIALL"
    "GSLTYAIIEAPSSGAAEILAFGGTALAALIALLWVEPRRAEPLIDLRFFRSVPFSGATVVAISAFAALSG"
    "FLFLSTLYLQNVRGLSALHAGLWMLPMAVMCFVCAPLSGRLVGSRGPRFSLLTAGIAMTASGVLFAAFEA"
    "ETSDVTLVIGYVLFGVGFGFVNAPITNTAVSGMPRAQAGVAAAIASTSRQIGGTLGVAVVGAVLASGVRT"
    "GPYRQVFESAARPGWWIIAGCGLSVLVLGALSSGAWARATARRTAQRLESTEVRDEAGVNA"
)
EXAMPLE_ACCESSION = "Q26473"
# "Try" chips on the home page: (label, accession, tooltip)
EXAMPLE_PROTEINS = [
    ("Semaphorin-7A", "O75326", "O75326, human: ALAASLW at 646-652, in a hydrophobic window"),
    ("Ammonium transporter", "P69681", "P69681, E. coli AmtB: ALVCILW at 75-81, in a hydrophobic window"),
]

def read_settings(values):
    """Read window/threshold/pattern/allowed from the form, using defaults for blanks."""
    def get(name, default):
        value = (values.get(name) or "").strip()
        return value if value else default

    try:
        window = int(get("window", an.DEFAULT_WINDOW))
    except ValueError:
        raise ValueError("Window size must be a whole number.") from None
    try:
        threshold = float(get("threshold", an.DEFAULT_THRESHOLD))
    except ValueError:
        raise ValueError("Threshold must be a number, e.g. 1.6.") from None
    pattern = re.sub(r"\s", "", get("pattern", an.DEFAULT_PATTERN)).upper()
    allowed = an.normalize_allowed(get("allowed", an.DEFAULT_ALLOWED))
    return an.Settings(window=window, threshold=threshold, pattern=pattern, allowed=allowed)


def settings_query(settings):
    """Only the settings that differ from the defaults, for the export link."""
    defaults = an.Settings()
    params = {}
    for name in ("window", "threshold", "pattern", "allowed"):
        value = getattr(settings, name)
        if value != getattr(defaults, name):
            params[name] = f"{value:g}" if isinstance(value, float) else value
    return params


def looks_like_accession(text):
    """Q26473, NP_000509.1 etc. Sequences never contain digits, accessions do."""
    return bool(ncbi.ACCESSION_RE.fullmatch(text)) and any(c.isdigit() for c in text)


def load_input(values):
    """Returns (sequence, record or None, warnings)."""
    mode = values.get("mode")
    accession = (values.get("accession") or "").strip()
    sequence_text = values.get("sequence") or ""

    # same box is used for accessions and sequences
    query = (values.get("query") or "").strip()
    if query:
        if looks_like_accession(query):
            mode, accession = "accession", query
        else:
            mode, sequence_text = "sequence", query

    if mode == "accession" or (mode is None and accession and not sequence_text.strip()):
        if not accession:
            raise ValueError("Enter an NCBI protein accession number, e.g. Q26473.")
        record = ncbi.fetch_record(accession)
        sequence, warnings = an.normalize_sequence(record.sequence)
        return sequence, record, warnings
    if mode is None and accession and sequence_text.strip():
        raise ValueError("Enter either an accession number or a sequence, not both.")
    sequence, warnings = an.normalize_sequence(sequence_text)
    return sequence, None, warnings


def is_fetch():
    """True when the request comes from app.js (it only wants the results HTML)."""
    return request.headers.get("X-Requested-With") == "fetch"


def render_page(values=None, error=None, results=None, status=200):
    """Render the page, or just the results/error part for fetch requests."""
    values = values or {}
    query = values.get("query") or values.get("accession") or values.get("sequence") or ""
    try:
        settings_changed = bool(values) and read_settings(values) != an.Settings()
    except ValueError:
        settings_changed = True     # keep the panel open so the bad value is visible
    if is_fetch():
        if error:
            return render_template("_error.html", error=error), status
        return render_template("_results.html", **results), status
    return render_template(
        "index.html",
        error=error,
        form=values,
        settings_changed=settings_changed,
        query=query,
        defaults=an.Settings(),
        example_accession=EXAMPLE_ACCESSION,
        example_proteins=EXAMPLE_PROTEINS,
        example_sequence=EXAMPLE_SEQUENCE,
        window_limits=(an.MIN_WINDOW, an.MAX_WINDOW),
        results=results,
        **(results or {}),
    ), status


def results_context(result, record, warnings):
    export_url = None
    if record:
        export_url = url_for("export_record", accession=record.accession, **settings_query(result.settings))
    return {
        "a": result,
        "record": record,
        "warnings": warnings + result.warnings,
        "rows": an.sequence_rows(result),
        "chart": an.chart_data(result),
        "aa_names": an.AA_NAMES,
        "export_url": export_url,
        "fasta_header": f"{record.version or record.accession} {record.definition}" if record else "pasted_sequence",
    }


@app.route("/")
def index():
    return render_page()


@app.route("/analyze", methods=["GET", "POST"])
def analyze():
    # Don't show results on a plain GET (reload, old bookmark), only after Analyze.
    if request.method == "GET":
        return redirect(url_for("index"))
    values = request.form
    try:
        settings = read_settings(values)
        sequence, record, warnings = load_input(values)
        result = an.analyze(sequence, settings)
    except ValueError as exc:
        return render_page(values, error=str(exc), status=400)
    except ncbi.NCBIError as exc:
        return render_page(values, error=str(exc), status=502)
    return render_page(values, results=results_context(result, record, warnings))


@app.route("/export/<accession>")
def export_record(accession):
    """Excel file of the NCBI record, plus the analysis sheets."""
    try:
        settings = read_settings(request.args)
        record = ncbi.fetch_record(accession)
        sequence, _ = an.normalize_sequence(record.sequence)
        result = an.analyze(sequence, settings)
    except ValueError as exc:
        return render_page(error=str(exc), status=400)
    except ncbi.NCBIError as exc:
        return render_page(error=str(exc), status=502)

    workbook = export.record_workbook(record.text, result, f"NCBI {record.version or record.accession}")
    return send_file(
        workbook,
        as_attachment=True,
        download_name=f"{safe_name(record.accession)}_NCBI_protein.xlsx",
        mimetype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )


@app.route("/export-analysis", methods=["POST"])
def export_analysis():
    """Excel export for pasted sequences (no NCBI record, so analysis only)."""
    try:
        settings = read_settings(request.form)
        sequence, _ = an.normalize_sequence(request.form.get("sequence", ""))
        result = an.analyze(sequence, settings)
    except ValueError as exc:
        return render_page(error=str(exc), status=400)
    return send_file(
        export.analysis_workbook(result, "Pasted sequence"),
        as_attachment=True,
        download_name="protein_motif_analysis.xlsx",
        mimetype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


@app.route("/health")
def health():
    return {"status": "ok"}


@app.errorhandler(404)
def not_found(_):
    return render_page(error="That page does not exist.", status=404)


@app.errorhandler(413)
def too_large(_):
    return render_page(error="The submitted text is too large (limit 2 MB).", status=413)


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=os.environ.get("FLASK_DEBUG") == "1")
