"""
Core sequence analysis: Kyte-Doolittle hydropathy profile and motif search.

Conventions used throughout this module:
  * Positions shown to users are 1-based and inclusive (residue 1 is the first).
  * A window of size w that starts at residue s covers residues s .. s+w-1.
    Its mean Kyte-Doolittle (KD) score is plotted at the window's midpoint,
    s + (w-1)//2, as described in the Kyte & Doolittle (1982) method.
  * A motif "is in a high hydrophobicity window" when at least one complete
    window whose mean KD score is >= threshold contains every motif residue.
"""

import re
from dataclasses import dataclass, field

# Kyte J & Doolittle RF (1982) J Mol Biol 157:105-132.
KD_SCALE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
STANDARD_AA = frozenset(KD_SCALE)

# IUPAC ambiguity / rare residue codes that appear in real NCBI records.
# They are accepted and scored with the closest sensible KD value.
AMBIGUOUS_AA = {
    "B": ("Asn or Asp", (KD_SCALE["N"] + KD_SCALE["D"]) / 2),
    "Z": ("Gln or Glu", (KD_SCALE["Q"] + KD_SCALE["E"]) / 2),
    "J": ("Leu or Ile", (KD_SCALE["L"] + KD_SCALE["I"]) / 2),
    "U": ("selenocysteine, scored as Cys", KD_SCALE["C"]),
    "O": ("pyrrolysine, scored as Lys", KD_SCALE["K"]),
    "X": ("unknown residue, scored as 0", 0.0),
}

AA_NAMES = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartate",
    "C": "Cysteine", "Q": "Glutamine", "E": "Glutamate", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
}

DEFAULT_WINDOW = 19
DEFAULT_THRESHOLD = 1.6
DEFAULT_PATTERN = "AL???LW"
# The 16 residues the assignment allows at each "?" position.
DEFAULT_ALLOWED = "AVLIMPFYWSTNCHQG"

MIN_WINDOW, MAX_WINDOW = 3, 99
MAX_SEQUENCE_LENGTH = 100_000
# Floating-point tolerance so that a window averaging exactly 1.6 passes 1.6.
EPSILON = 1e-9


@dataclass(frozen=True)
class Settings:
    window: int = DEFAULT_WINDOW
    threshold: float = DEFAULT_THRESHOLD
    pattern: str = DEFAULT_PATTERN
    allowed: str = DEFAULT_ALLOWED

    @property
    def is_default(self):
        return self == Settings()


@dataclass
class Window:
    start: int      # first residue, 1-based
    end: int        # last residue, 1-based inclusive
    center: int     # residue the score is plotted at
    score: float

    @property
    def is_high(self):
        return self.score >= self._threshold - EPSILON

    _threshold: float = field(default=0.0, repr=False)


@dataclass
class MotifHit:
    motif: str
    start: int
    end: int
    windows: list           # every complete window containing the whole motif
    qualifying: bool
    best_window: Window | None

    @property
    def high_windows(self):
        return [w for w in self.windows if w.is_high]


@dataclass
class Segment:
    """A run of overlapping high-scoring windows: a candidate hydrophobic region."""
    start: int
    end: int
    peak_score: float
    peak_center: int

    @property
    def length(self):
        return self.end - self.start + 1


@dataclass
class Analysis:
    sequence: str
    settings: Settings
    pattern_regex: str
    windows: list
    residue_scores: list    # index i -> score of the window centred on residue i+1, or None
    high_residues: set      # 1-based positions whose centred window passes the threshold
    hits: list
    segments: list
    warnings: list

    @property
    def qualifying_hits(self):
        return [h for h in self.hits if h.qualifying]

    @property
    def high_window_count(self):
        return sum(1 for w in self.windows if w.is_high)

    @property
    def max_window(self):
        return max(self.windows, key=lambda w: w.score) if self.windows else None


def normalize_sequence(raw):
    """
    Turn pasted text into a clean upper-case sequence.

    Accepts plain sequences, FASTA (header lines starting with '>' or ';'),
    and GenBank ORIGIN-style text with position numbers. A trailing stop
    codon '*' is dropped. Returns (sequence, warnings).
    """
    lines = [line.strip() for line in (raw or "").splitlines()]
    body = [line for line in lines if line and not line.startswith((">", ";"))]
    if sum(1 for line in lines if line.startswith(">")) > 1:
        raise ValueError(
            "Please paste a single protein sequence (the text contains more than one FASTA record)."
        )
    sequence = re.sub(r"[\s\d]", "", "".join(body)).upper().rstrip("*")
    if not sequence:
        raise ValueError("Enter an NCBI accession number or paste a protein sequence.")

    invalid = sorted(set(sequence) - STANDARD_AA - set(AMBIGUOUS_AA))
    if invalid:
        shown = ", ".join(f"'{c}'" for c in invalid[:10])
        raise ValueError(
            f"The sequence contains characters that are not amino-acid codes: {shown}. "
            "Paste a protein (amino-acid) sequence, not DNA/RNA or other text."
        )
    if len(sequence) > MAX_SEQUENCE_LENGTH:
        raise ValueError(
            f"The sequence is {len(sequence):,} residues long; the limit is {MAX_SEQUENCE_LENGTH:,}."
        )

    warnings = []
    for code in sorted(set(sequence) & set(AMBIGUOUS_AA)):
        meaning, value = AMBIGUOUS_AA[code]
        warnings.append(
            f"{sequence.count(code)} '{code}' residue(s) ({meaning}) were scored as {value:g}."
        )
    if set(sequence) <= set("ACGTUN"):
        warnings.append(
            "This sequence contains only A, C, G, T/U and N, so it may be DNA/RNA rather than protein."
        )
    return sequence, warnings


def normalize_allowed(allowed):
    letters = re.sub(r"[\s,;]", "", allowed or "").upper()
    if not letters:
        raise ValueError("The allowed residues for '?' cannot be empty.")
    bad = sorted(set(letters) - STANDARD_AA)
    if bad:
        raise ValueError(f"Allowed residues must be standard amino acids; remove: {', '.join(bad)}.")
    # Keep the user's order, drop duplicates.
    return "".join(dict.fromkeys(letters))


def compile_pattern(pattern, allowed=DEFAULT_ALLOWED):
    """
    Convert a motif such as 'AL???LW' into a regular expression.

    Syntax: a letter matches that residue, '?' matches any residue in
    `allowed`, and '[ILV]' matches any one of the bracketed residues.
    Returns (compiled_regex, motif_length, readable_regex).
    """
    text = re.sub(r"\s", "", pattern or "").upper()
    if not text:
        raise ValueError("The motif pattern cannot be empty.")

    parts, length, i = [], 0, 0
    while i < len(text):
        char = text[i]
        if char == "?":
            parts.append(f"[{allowed}]")
        elif char == "[":
            close = text.find("]", i)
            if close == -1:
                raise ValueError("The motif pattern has a '[' without a matching ']'.")
            group = text[i + 1:close]
            if not group or set(group) - STANDARD_AA:
                raise ValueError("Brackets in the motif must contain amino-acid letters, e.g. [ILV].")
            parts.append(f"[{''.join(dict.fromkeys(group))}]")
            i = close
        elif char in STANDARD_AA:
            parts.append(char)
        else:
            raise ValueError(
                f"'{char}' is not allowed in the motif. Use amino-acid letters, '?' and [..] groups."
            )
        length += 1
        i += 1

    readable = "".join(parts)
    # The lookahead lets overlapping occurrences be reported as well.
    return re.compile(f"(?=({readable}))"), length, readable


def kd_value(residue):
    if residue in KD_SCALE:
        return KD_SCALE[residue]
    return AMBIGUOUS_AA[residue][1]


def sliding_windows(sequence, window, threshold):
    """Mean KD score for every complete window, computed with a running sum."""
    if len(sequence) < window:
        return []
    values = [kd_value(aa) for aa in sequence]
    offset = (window - 1) // 2
    total = sum(values[:window])
    windows = []
    for start0 in range(len(sequence) - window + 1):
        if start0:
            total += values[start0 + window - 1] - values[start0 - 1]
        windows.append(Window(
            start=start0 + 1,
            end=start0 + window,
            center=start0 + 1 + offset,
            score=round(total / window, 4),
            _threshold=threshold,
        ))
    return windows


def find_motifs(sequence, regex):
    return [
        (m.group(1), m.start() + 1, m.start() + len(m.group(1)))
        for m in regex.finditer(sequence)
    ]


def merge_segments(windows):
    segments = []
    for w in (w for w in windows if w.is_high):
        if segments and w.start <= segments[-1].end + 1:
            seg = segments[-1]
            seg.end = max(seg.end, w.end)
            if w.score > seg.peak_score:
                seg.peak_score, seg.peak_center = w.score, w.center
        else:
            segments.append(Segment(w.start, w.end, w.score, w.center))
    return segments


def analyze(sequence, settings=Settings()):
    if not MIN_WINDOW <= settings.window <= MAX_WINDOW:
        raise ValueError(f"Window size must be between {MIN_WINDOW} and {MAX_WINDOW}.")
    if not -4.5 <= settings.threshold <= 4.5:
        raise ValueError("Threshold must be between -4.5 and 4.5 (the range of the Kyte-Doolittle scale).")

    regex, motif_length, readable = compile_pattern(settings.pattern, settings.allowed)
    if motif_length > settings.window:
        raise ValueError(
            f"The motif is {motif_length} residues long but the window is only {settings.window}; "
            "a window must be at least as long as the motif to contain it."
        )

    # A sequence shorter than the window gets no windows; the page explains why.
    windows = sliding_windows(sequence, settings.window, settings.threshold)

    residue_scores = [None] * len(sequence)
    for w in windows:
        residue_scores[w.center - 1] = w.score
    high_residues = {w.center for w in windows if w.is_high}

    hits = []
    last_start = len(sequence) - settings.window + 1
    for motif, start, end in find_motifs(sequence, regex):
        # Windows containing the whole motif start between end-w+1 and start.
        first = max(1, end - settings.window + 1)
        last = min(start, last_start)
        containing = windows[first - 1:last] if last >= first else []
        high = [w for w in containing if w.is_high]
        best = max(containing, key=lambda w: w.score) if containing else None
        hits.append(MotifHit(motif, start, end, containing, bool(high), best))

    return Analysis(
        sequence=sequence,
        settings=settings,
        pattern_regex=readable,
        windows=windows,
        residue_scores=residue_scores,
        high_residues=high_residues,
        hits=hits,
        segments=merge_segments(windows),
        warnings=[],
    )


def sequence_rows(analysis, line_length=60, block=10):
    """
    Lay the sequence out like NCBI's ORIGIN section (60 residues per line,
    blocks of 10) with a highlight class for each residue.
    """
    qualifying = set()
    other_motif = set()
    for hit in analysis.hits:
        target = qualifying if hit.qualifying else other_motif
        target.update(range(hit.start, hit.end + 1))

    rows = []
    seq = analysis.sequence
    for line_start in range(0, len(seq), line_length):
        blocks = []
        for block_start in range(line_start, min(line_start + line_length, len(seq)), block):
            residues = []
            for index in range(block_start, min(block_start + block, len(seq))):
                pos = index + 1
                classes = []
                if pos in qualifying:
                    classes.append("motif")
                elif pos in other_motif:
                    classes.append("motif-weak")
                if pos in analysis.high_residues:
                    classes.append("high")
                residues.append({
                    "aa": seq[index],
                    "pos": pos,
                    "cls": " ".join(classes),
                    "score": analysis.residue_scores[index],
                })
            blocks.append(residues)
        rows.append({"start": line_start + 1, "blocks": blocks})
    return rows


def chart_data(analysis):
    """Plain data for the browser-side charts."""
    seq = analysis.sequence
    return {
        "threshold": analysis.settings.threshold,
        "window": analysis.settings.window,
        "length": len(seq),
        "profile": [{"x": w.center, "y": w.score} for w in analysis.windows],
        "segments": [
            {"start": s.start, "end": s.end, "peak": s.peak_score} for s in analysis.segments
        ],
        "motifs": [
            {
                "motif": h.motif,
                "start": h.start,
                "end": h.end,
                "qualifying": h.qualifying,
                "labels": [f"{seq[w.center - 1]}({w.center})" for w in h.windows],
                "scores": [w.score for w in h.windows],
                "ranges": [f"{w.start}-{w.end}" for w in h.windows],
            }
            for h in analysis.hits
        ],
    }
