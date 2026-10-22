# Protein Motif & Hydropathy Analyzer

Small Flask app for the protein sequence analysis assignment. You give it an NCBI protein
accession (or paste a sequence) and it finds every `AL???LW` motif that sits inside a
hydrophobic region of the protein.

"Hydrophobic region" here means a 19-residue window whose average Kyte-Doolittle score is
1.6 or more. Each `?` in the motif can be any of `A V L I M P F Y W S T N C H Q G`.

## What's on the results page

These follow the five points listed under "HTML output page" in the problem statement.

1. **Link to the NCBI search page** for the accession (the "NCBI record" button).
2. **Excel export of the NCBI record.** Same fields as the GenPept page (LOCUS, DEFINITION,
   references, features, ORIGIN...), no HTML.
3. **Hydrophobicity plot** for the whole protein. The problem statement suggests taking this
   from TMpred, but TMpred (ch.embnet.org) has been shut down, so the app draws it itself.
4. **Kyte-Doolittle plot for each motif.** Shows every window that contains the motif and
   whether any of them passes the threshold.
5. **The sequence**, with the motif in blue and the high-scoring residues in red.

A few extra things I added along the way:

- one input box for both accessions and sequences; it figures out which one you typed
- window size, threshold, motif pattern and allowed residues can all be changed under Settings
- motifs that don't pass are listed too, with the best score they reached
- Excel export also works for pasted sequences (analysis sheets only)
- charts show values on hover and can be saved as PNG
- clear messages for things like a wrong accession, DNA pasted instead of protein, or NCBI being down

## How it works

**Getting the sequence.** For an accession the app calls NCBI's `efetch`
(`db=protein&rettype=gp`), which returns the same text as the GenPept page on the NCBI site.
The sequence is read from the ORIGIN section. Pasted text can be plain, FASTA, or copied
straight from GenBank with the numbers still in it.

**Scoring.** Every residue gets its Kyte-Doolittle value (from -4.5 for R up to 4.5 for I).
A 19-residue window slides along the sequence and the average of each window is plotted at
its middle residue.

**Finding the motif.** `AL???LW` becomes the regex `AL[AVLIMPFYWSTNCHQG]{3}LW`. It's run with a
lookahead so overlapping matches aren't missed.

**Deciding if a motif counts.** A motif counts if at least one full window with an average of
1.6 or more contains the whole motif. For a 7-letter motif and a 19-letter window that's 13
possible windows, which is why each motif plot has 13 points.

I checked this against the example in the problem statement. Running the 481-residue example
protein gives the same 13 values as the example plot (0.36, 0.49, 0.72 ... 1.81 around
ALHAGLW at 297-303). The only difference is that the example labels each point one residue to
the right. I've used the middle of the window, which is the usual way to plot it.

## Examples

The "Try" buttons on the home page load these:

- **Semaphorin-7A** (O75326, human): ALAASLW at 646-652
- **Ammonium transporter** (P69681, *E. coli* AmtB): ALVCILW at 75-81
- **Assignment example**: the 481-residue protein from the example output, with ALHAGLW (297-303)
  and ALIALLW (238-244)

All three have a motif that passes. The motif is actually quite rare: out of about 25,000
reviewed human and E. coli proteins in UniProt, only 41 have it inside a hydrophobic window.

Q26473, the accession used in the problem statement, doesn't contain `AL???LW` at all. It
still works fine in the app and is a good one to try the Excel export with.

## Running it locally

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```

Then open http://127.0.0.1:5000. Use `PORT=8000 python app.py` for a different port.

You can optionally set `NCBI_EMAIL` and `NCBI_API_KEY`. NCBI asks for these if you send a lot
of requests, and an API key raises the limit from 3 to 10 requests a second.

## Deploying

It's set up for Vercel. Vercel picks up the Flask `app` in `app.py` on its own, so importing
the repo and clicking Deploy is enough. The CSS/JS are in `public/static/` because Vercel only
serves static files from `public/`. Locally Flask serves the same folder at `/static`.

If NCBI starts rate limiting on Vercel, add `NCBI_API_KEY` under Settings > Environment Variables.

There's also a `Procfile` (gunicorn) if you'd rather use Render or Railway.

## Tests

```bash
pip install -r requirements-dev.txt
pytest
```

The tests check the scores against the example plot, the motif rules, input cleanup, the
GenPept parsing, all the routes and what ends up in the Excel files. They use a saved copy of
the Q26473 record, so they don't need internet.

## Files

```text
app.py           Flask routes
analysis.py      Kyte-Doolittle scoring and motif search
ncbi.py          fetching from NCBI and parsing GenPept records
export.py        Excel export
templates/       HTML templates
public/static/   CSS, JS (charts use Chart.js), logo
tests/           pytest tests
```

## Limitations

- A high hydrophobicity score suggests the region crosses the membrane but doesn't prove it.
  Tools like DeepTMHMM or Phobius are more accurate if you need a proper prediction.
- Looking up an accession needs a working connection to NCBI.

## Author

Madhvi Mishra
