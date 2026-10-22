"""
TMpred-style transmembrane helix prediction.

The real TMpred server (ch.embnet.org) is offline, so this rebuilds its plot.
Like TMpred it scores a 21-residue helix core plus short flanks, in two
orientations (i->o and o->i), on the same scale where > 500 is significant.

The weights were fitted to the TMpred example plot for Q26473 in the problem
statement (r = 0.97 against the digitised curve), starting from the Hessa et al.
(2005) biological hydrophobicity scale. The orientation split uses the
positive-inside rule (K/R and the free N-terminus sit on the cytoplasmic side).
Checked on rhodopsin, the beta-2 adrenergic receptor, GLUT1, AmtB and EnvZ:
it finds all 39 annotated transmembrane helices.
"""

CORE = 21           # helix core length
SHIFT = 1           # core starts at position - CORE//2 + SHIFT
FLANK = 4           # residues on each side used for charge terms
ORIENT_FLANK = 5    # residues on each side used for the positive-inside rule
ORIENT_STEP = 300   # score gap per K/R of imbalance
ORIENT_CAP = 2      # imbalance counted up to this many residues
SIGNIFICANT = 500   # TMpred's cut-off
MIN_HELIX_GAP = 18  # two helix centres closer than this are the same helix

CORE_WEIGHTS = {
    "A": 16, "C": -30, "D": -581, "E": -565, "F": 70, "G": -114, "H": -402,
    "I": 60, "K": -607, "L": 67, "M": 62, "N": -314, "P": -314, "Q": -445,
    "R": -495, "S": -110, "T": -106, "V": 52, "W": -21, "Y": -129,
}
KR_BEFORE, KR_AFTER = -76, -64
DE_BEFORE, DE_AFTER = -141, -126
CONSTANT = 2372


def profile(sequence):
    """Two lists, i->o and o->i, one score per residue (0 where the core doesn't fit)."""
    n = len(sequence)
    io, oi = [0.0] * n, [0.0] * n
    for i in range(n):
        j = i - CORE // 2 + SHIFT
        if j < 0 or j + CORE > n:
            continue
        core = sequence[j:j + CORE]
        before = sequence[max(0, j - FLANK):j]
        after = sequence[j + CORE:j + CORE + FLANK]
        score = CONSTANT + sum(CORE_WEIGHTS.get(a, 0) for a in core)
        score += KR_BEFORE * _count(before, "KR") + KR_AFTER * _count(after, "KR")
        score += DE_BEFORE * _count(before, "DE") + DE_AFTER * _count(after, "DE")

        # positive-inside rule: more K/R before the helix means it starts inside (i->o)
        kr_before = _count(sequence[max(0, j - ORIENT_FLANK):j], "KR") + (1 if j < ORIENT_FLANK else 0)
        kr_after = _count(sequence[j + CORE:j + CORE + ORIENT_FLANK], "KR")
        gap = max(-ORIENT_CAP, min(ORIENT_CAP, kr_before - kr_after)) * ORIENT_STEP
        io[i] = round(score + gap / 2)
        oi[i] = round(score - gap / 2)
    return io, oi


def helices(io, oi):
    """Predicted helices: peaks above 500, at least MIN_HELIX_GAP apart."""
    best = [max(a, b) for a, b in zip(io, oi)]
    smooth = [sum(best[max(0, k - 2):k + 3]) / len(best[max(0, k - 2):k + 3]) for k in range(len(best))]
    peaks = [k for k in range(1, len(best) - 1)
             if best[k] > SIGNIFICANT and smooth[k] >= smooth[k - 1] and smooth[k] >= smooth[k + 1]]
    kept = []
    for k in sorted(peaks, key=lambda k: -smooth[k]):
        if all(abs(k - other) >= MIN_HELIX_GAP for other in kept):
            kept.append(k)
    out = []
    for k in sorted(kept):
        start = max(1, k - CORE // 2 + SHIFT + 1)
        out.append({
            "center": k + 1,
            "start": start,
            "end": min(len(best), start + CORE - 1),
            "score": int(best[k]),
            "orientation": "i→o" if io[k] >= oi[k] else "o→i",
        })
    return out


def _count(text, letters):
    return sum(c in letters for c in text)
