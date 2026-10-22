"""
Fetch protein records from NCBI and parse the GenPept text format.

efetch with rettype=gp gives the same text as the GenPept view on the NCBI
site. The sequence is in the ORIGIN section at the end.
"""

import os
import re
from dataclasses import dataclass
from functools import lru_cache
from urllib.parse import quote

import requests

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
SEARCH_URL = "https://www.ncbi.nlm.nih.gov/protein/?term={}"
RECORD_URL = "https://www.ncbi.nlm.nih.gov/protein/{}"
ACCESSION_RE = re.compile(r"[A-Za-z0-9_.]{2,30}")
TIMEOUT_SECONDS = 20


class NCBIError(Exception):
    """NCBI is down, slow, or rate limiting us."""


@dataclass(frozen=True)
class ProteinRecord:
    accession: str
    text: str           # full GenPept text
    sequence: str
    definition: str
    organism: str
    version: str

    @property
    def search_url(self):
        return search_url(self.accession)

    @property
    def record_url(self):
        return RECORD_URL.format(quote(self.version or self.accession))


def search_url(accession):
    return SEARCH_URL.format(quote(accession))


def clean_accession(accession):
    accession = (accession or "").strip()
    if not ACCESSION_RE.fullmatch(accession):
        raise ValueError(
            "That does not look like an NCBI protein accession "
            "(examples: Q26473, NP_000509.1, AAA29808)."
        )
    return accession


def fetch_record(accession):
    return _fetch_record(clean_accession(accession).upper())


@lru_cache(maxsize=128)
def _fetch_record(accession):
    params = {
        "db": "protein",
        "id": accession,
        "rettype": "gp",
        "retmode": "text",
        "tool": "protein-analysis-tool",
    }
    # optional, NCBI asks for these if you send a lot of requests
    if os.environ.get("NCBI_EMAIL"):
        params["email"] = os.environ["NCBI_EMAIL"]
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]

    try:
        response = requests.get(EFETCH_URL, params=params, timeout=TIMEOUT_SECONDS)
    except requests.Timeout as exc:
        raise NCBIError("NCBI did not respond in time. Please try again in a moment.") from exc
    except requests.RequestException as exc:
        raise NCBIError("Could not connect to NCBI. Check the internet connection and try again.") from exc

    if response.status_code == 429:
        raise NCBIError("NCBI is rate-limiting requests right now. Wait a few seconds and try again.")
    if response.status_code == 400 or (response.ok and "ORIGIN" not in response.text):
        raise ValueError(
            f"NCBI has no protein record for '{accession}'. Check the accession number "
            "(it must be a protein accession, not a nucleotide one)."
        )
    if not response.ok:
        raise NCBIError(f"NCBI returned an error (HTTP {response.status_code}). Please try again later.")

    text = response.text.strip()
    # should only be one record, but just in case keep the first
    text = text.split("\n//", 1)[0] + "\n//"
    sequence = parse_origin(text)
    if not sequence:
        raise ValueError(f"The NCBI record for '{accession}' does not contain a protein sequence.")

    fields = dict(top_level_fields(text))
    return ProteinRecord(
        accession=fields.get("ACCESSION", accession).split()[0],
        text=text,
        sequence=sequence,
        definition=readable_definition(fields.get("DEFINITION", "")),
        organism=fields.get("SOURCE", ""),
        version=fields.get("VERSION", "").split()[0] if fields.get("VERSION") else "",
    )


def readable_definition(definition):
    """
    UniProt-based records look like 'RecName: Full=Semaphorin-1A; AltName: ...'.
    Pull out just the name. Anything else is returned as is.
    """
    definition = definition.strip().rstrip(".")
    match = re.match(r"RecName: Full=([^;]+)", definition)
    return match.group(1).strip() if match else definition


def parse_origin(text):
    if "\nORIGIN" not in text:
        return ""
    origin = text.split("\nORIGIN", 1)[1].split("//", 1)[0]
    return re.sub(r"[^A-Za-z]", "", origin).upper()


def top_level_fields(text):
    """(KEY, value) for the top-level sections like LOCUS, DEFINITION, SOURCE."""
    for entry in parse_genpept(text):
        if entry.level == 0 and entry.key not in ("FEATURES", "ORIGIN"):
            yield entry.key, entry.value


@dataclass
class Entry:
    key: str
    value: str
    level: int          # 0 = LOCUS etc, 1 = AUTHORS etc, 2 = PUBMED
    lines: list


def parse_genpept(text):
    """
    Split a GenPept record into entries. The label is in the first 12
    columns and the value starts at column 13; lines with an empty label
    continue the previous entry. Indentation of the label tells the level
    (REFERENCE > AUTHORS > PUBMED). FEATURES and ORIGIN keep their raw lines.
    """
    entries = []
    current = None
    for line in text.splitlines():
        if line.startswith("//"):
            break
        if not line.strip():
            continue
        label = line[:12]
        key = label.strip()
        in_block = current is not None and current.key in ("FEATURES", "ORIGIN")
        if key and not (in_block and label.startswith(" ")):
            indent = len(label) - len(label.lstrip(" "))
            level = 0 if indent == 0 else 1 if indent <= 2 else 2
            current = Entry(key=key, value="", level=level, lines=[line[12:].rstrip()])
            entries.append(current)
        elif current is not None:
            current.lines.append(line[12:].rstrip() if not in_block else line.rstrip())

    for entry in entries:
        if entry.key == "COMMENT":
            entry.value = "\n".join(l.strip() for l in entry.lines if l.strip())
        else:
            entry.value = " ".join(l.strip() for l in entry.lines if l.strip())
    return entries


@dataclass
class Feature:
    key: str
    location: str
    qualifiers: list    # [(name, value)]


def parse_features(text):
    """FEATURES table -> list of features with their /qualifiers."""
    block = None
    for entry in parse_genpept(text):
        if entry.key == "FEATURES":
            block = entry.lines[1:]
            break
    if not block:
        return []

    features = []
    raw = None      # qualifier currently being read: [name, value]
    for line in block:
        if len(line) > 5 and line[5] != " " and line[:5] == "     ":
            key, _, location = line.strip().partition(" ")
            features.append(Feature(key, location.strip(), []))
            raw = None
            continue
        if not features:
            continue
        content = line.strip()
        feature = features[-1]
        # a line starting with '/' can still be part of an open "..." value
        quote_open = raw is not None and raw[1].count('"') % 2 == 1
        if content.startswith("/") and not quote_open:
            name, _, value = content[1:].partition("=")
            raw = [name, value]
            feature.qualifiers.append(raw)
        elif raw is not None:
            # value continues on this line (translation is joined without spaces)
            raw[1] += ("" if raw[0] == "translation" else " ") + content
        else:
            feature.location += content

    for feature in features:
        feature.qualifiers = [(name, value.strip('"')) for name, value in feature.qualifiers]
    return features


def origin_lines(text):
    """ORIGIN lines as NCBI prints them, e.g. '1 mraalvavaa llwvalhaaa ...'."""
    if "\nORIGIN" not in text:
        return []
    origin = text.split("\nORIGIN", 1)[1].split("//", 1)[0]
    return [line.strip() for line in origin.splitlines() if line.strip()]
