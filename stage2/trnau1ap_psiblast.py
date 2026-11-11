#!/usr/bin/env python3
"""Do the appendicularian proteins found by a PSI-BLAST search for SECp43 (TRNAU1AP) belong to SECp43?

Input: the sequences of a PSI-BLAST search of the NCBI proteins of tunicates seeded with SECp43, downloaded as FASTA
from the NCBI web site (curation/trnau1ap_blast/psiblast.fa, search by Charles Plessy, 2026-11-11; parameters in
curation/trnau1ap_blast/searches.tsv). Such a search returns SECp43 first, then, as the iterations broaden the
profile, other proteins with RNA-recognition motifs (poly(A)-binding, ELAV-like, CELF, Musashi, TIA1-like and others),
so that the presence of an appendicularian protein in the result says nothing by itself.

Each sequence is labelled SECp43 when its NCBI name says so ("tRNA selenocysteine 1-associated protein"); an unnamed
or uncharacterised sequence ("unnamed protein product", "uncharacterized protein", "hypothetical protein") takes the
label of its best blastp hit among the named sequences of other species. Each sequence of --target (Oikopleura) is
then compared by blastp (legacy BLAST, as in the container) with the sequences of the other species: best bit score to
a SECp43 and to any other protein, and the margin (other minus SECp43; positive: closer to another protein). The
SECp43 of the other species, compared in the same way with the sequences of the other species (their own species
left out), are the positive controls.
Output (--out, TSV): species, accession, name, label, best SECp43 hit and score, best other hit and score, margin.
"""
import argparse, os, re, subprocess, tempfile

UNNAMED = re.compile(r"unnamed protein product|uncharacterized protein|hypothetical protein|LOW QUALITY PROTEIN: uncharacterized")
SECP43 = re.compile(r"tRNA selenocysteine 1-associated protein", re.I)


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        l = l.rstrip("\n")
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:], []
        elif l.strip(): buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def species(name):
    m = re.search(r"\[([^\]]+)\]\s*$", name)
    return m.group(1) if m else ""


def blast(queries, db, d):
    """Best hits: (query, subject) -> bit score, for all pairs reported by blastp."""
    with open(f"{d}/db.faa", "w") as o:
        for acc, s in db.items(): o.write(f">{acc}\n{s}\n")
    with open(f"{d}/q.faa", "w") as o:
        for acc, s in queries.items(): o.write(f">{acc}\n{s}\n")
    subprocess.run(["formatdb", "-i", "db.faa", "-p", "T"], cwd=d, check=True, stdout=subprocess.DEVNULL)
    out = subprocess.run(["blastall", "-p", "blastp", "-d", "db.faa", "-i", "q.faa", "-e", "10", "-m", "8",
                          "-v", "1000", "-b", "1000"], cwd=d, check=True, capture_output=True, text=True).stdout
    best = {}
    for l in out.splitlines():
        f = l.split("\t")
        k = (f[0], f[1]); b = float(f[11])
        if b > best.get(k, 0): best[k] = b
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", default="curation/trnau1ap_blast/psiblast.fa")
    ap.add_argument("--target", default="Oikopleura", help="genus whose proteins are tested")
    ap.add_argument("--out", default="results/tables/trnau1ap_psiblast.tsv")
    a = ap.parse_args()
    recs = [(n.split()[0], n.split(" ", 1)[1] if " " in n else "", species(n), s) for n, s in read_fasta(a.fasta)]
    name = {acc: nm for acc, nm, sp, s in recs}
    sp = {acc: x for acc, nm, x, s in recs}
    seq = {acc: s for acc, nm, x, s in recs}
    label = {acc: ("SECp43" if SECP43.search(nm) else "other") for acc, nm, x, s in recs if not UNNAMED.search(nm)}
    with tempfile.TemporaryDirectory() as d:
        unnamed = {acc: seq[acc] for acc, nm, x, s in recs if UNNAMED.search(nm)}
        if unnamed:   # label of the best hit among the named sequences of other species
            hits = blast(unnamed, {acc: seq[acc] for acc in label}, d)
            for q in unnamed:
                cand = [(b, s) for (qq, s), b in hits.items() if qq == q and sp[s] != sp[q]]
                label[q] = label[max(cand)[1]] if cand else "other"
        need = {acc: seq[acc] for acc in seq if a.target in sp[acc] or label[acc] == "SECp43"}
        hits = blast(need, seq, d)
    rows = []
    for acc, nm, x, s in recs:
        test = a.target in x
        if not test and label[acc] != "SECp43": continue
        best = {}
        for (q, t), b in hits.items():
            if q != acc or sp[t] == x or a.target in sp[t]: continue   # other species, outside the target genus
            k = label[t]
            if b > best.get(k, (0, ""))[0]: best[k] = (b, t)
        bs, bo = best.get("SECp43", (0, "")), best.get("other", (0, ""))
        rows.append(dict(role="tested" if test else "control", species=x, accession=acc, name=nm.rsplit(" [", 1)[0],
                         label=label[acc], secp43_hit=bs[1], secp43_bits=f"{bs[0]:.1f}", other_hit=bo[1],
                         other_name=name.get(bo[1], "").rsplit(" [", 1)[0], other_bits=f"{bo[0]:.1f}",
                         margin=f"{bo[0] - bs[0]:.1f}"))
    cols = list(rows[0])
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(r[c] for c in cols) + "\n")
    t = [r for r in rows if r["role"] == "tested"]; c = [r for r in rows if r["role"] == "control"]
    print(f"{len(recs)} sequences, {sum(1 for v in label.values() if v == 'SECp43')} labelled SECp43 "
          f"({len({sp[k] for k, v in label.items() if v == 'SECp43'})} species)")
    print(f"tested ({a.target}): {len(t)}; labelled SECp43: {sum(r['label'] == 'SECp43' for r in t)}; margins "
          f"{min(float(r['margin']) for r in t):.0f} to {max(float(r['margin']) for r in t):.0f} bits; best SECp43 score "
          f"{max(float(r['secp43_bits']) for r in t):.0f} bits")
    print(f"controls: {len(c)}; margins {min(float(r['margin']) for r in c):.0f} to "
          f"{max(float(r['margin']) for r in c):.0f} bits; SECp43 scores {min(float(r['secp43_bits']) for r in c):.0f} "
          f"to {max(float(r['secp43_bits']) for r in c):.0f}")


if __name__ == "__main__":
    main()
