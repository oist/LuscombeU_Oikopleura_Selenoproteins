#!/usr/bin/env python3
"""Spliced leader (SL) at the 5' end of the transcripts drawn in the machinery figures.

For each gene, the transcripts drawn below the gold-standard models (rows "transcript ... exons" of
figures/inputs/reannotation_<gene>.tsv, make machinery-reannotation-inputs) are read from the transcript file of
their assembly (stage1/samplesheet.csv, column transcripts: OKI2018 and OSKA2016 transcript assemblies, OdB3 cDNA
collection). A transcript starts with the SL when its sequence carries the 3' end of the 40-nt O. dioica SL
ACTCATCCCATTTTTGAGTCCGATTTCGATTGTCTAACAG (Ganot et al. 2004) within the first --window nucleotides (at most one
mismatch), or the reverse complement within the last --window nucleotides (transcript assemblies are not oriented):
the criterion of make operons (stage2/operons.py). Masked bases (lower case or N) are counted at both ends, because
the OdB3 cDNA collection masks the 5' end of its reads; for those, the GenBank record of the read decides
(curation/spliced_leader_transcripts.tsv).

Output (--out): one row per drawn transcript: gene, assembly, transcript, evidence (empty when no SL), sl (sense,
antisense or none), sl_end (length of the transcript, read as the mRNA, up to the end of the SL), length, masked_head
and masked_tail (masked bases among the first and last --window bases of the record), head and tail (the first and the
last 50 bases of the record, as stored). figures/reannotation_page.py marks the transcripts whose sl is not none.
"""
import argparse, collections, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gold_standard import read_fasta, rc   # noqa: E402

SL_3P = "GATTGTCTAACAG"


def sl_match(seq, sl, window):
    """Position just after the SL at the 5' end of seq (sense), or None; at most one mismatch (as stage2/operons.py)."""
    n = len(sl)
    for i in range(0, max(0, min(window, len(seq)) - n + 1)):
        if sum(a != b for a, b in zip(seq[i:i + n], sl)) <= 1: return i + n
    return None


def masked(s): return sum(1 for c in s if c.islower() or c in "Nn")


def drawn(inputs, gene):
    """(assembly, transcript) drawn in the figure of a gene, in order."""
    out = []
    path = f"{inputs}/reannotation_{gene}.tsv"
    if not os.path.exists(path): return out
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["kind"] == "transcript" and r["field"] == "exons": out.append((r["assembly"], r["name"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--inputs", default="figures/inputs")
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--sl", default=SL_3P, help="3' part of the SL to search (default: as make operons)")
    ap.add_argument("--window", type=int, default=80)
    ap.add_argument("--out", default="results/tables/sl_transcripts.tsv")
    a = ap.parse_args()

    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    todo = [(g, asm, t) for g in a.genes for asm, t in drawn(a.inputs, g)]
    names = collections.defaultdict(set)
    for g, asm, t in todo: names[asm].add(t)
    seqs = {}
    for asm, keep in names.items():
        path = sheet.get(asm, {}).get("transcripts", "")
        if not path or not os.path.exists(path):
            sys.exit(f"no transcript file for {asm} ({path or 'samplesheet column empty'})")
        for k, s in read_fasta(path, keep): seqs[(asm, k)] = s
    rows = []
    for g, asm, t in todo:
        s = seqs.get((asm, t))
        if s is None:
            rows.append({"gene": g, "assembly": asm, "transcript": t, "evidence": "", "sl": "not found in the file",
                         "sl_end": "", "length": "", "masked_head": "", "masked_tail": "", "head": "", "tail": ""})
            continue
        u = s.upper()
        p, how = sl_match(u, a.sl, a.window), "sense"
        if p is None: p, how = sl_match(rc(u), a.sl, a.window), "antisense"
        if p is None: how = "none"
        ev = "" if p is None else (f"the SL (3' part, {a.sl}, at most one mismatch) at the 5' end of the transcript "
                                   f"({'reverse complement of the record' if how == 'antisense' else 'record'}), "
                                   f"make sl-transcripts")
        rows.append({"gene": g, "assembly": asm, "transcript": t, "evidence": ev, "sl": how,
                     "sl_end": "" if p is None else p, "length": len(s), "masked_head": masked(s[:a.window]),
                     "masked_tail": masked(s[-a.window:]), "head": s[:50], "tail": s[-50:]})
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["gene"], delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    n = collections.Counter((r["gene"], r["sl"] not in ("none", "not found in the file")) for r in rows)
    for g in a.genes:
        print(f"{g}: {n[(g, True)]} of {n[(g, True)] + n[(g, False)]} drawn transcripts start with the SL", file=sys.stderr)


if __name__ == "__main__":
    main()
