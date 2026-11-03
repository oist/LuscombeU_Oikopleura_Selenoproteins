#!/usr/bin/env python3
"""Intron positions of curation/intron_positions.py placed on a multiple alignment.

The proteins of the models of --positions (RefSeq proteins from --proteins, matched by accession; gold-standard
proteins from <gold-dir>/<genome>.proteins.faa) and supporting sequences of the same family (--extra: names of
--proteins starting with these prefixes) are aligned with MAFFT L-INS-i (--localpair --maxiterate 1000). Each intron
gets the alignment column of the residue that follows it (phase 0) or that contains it (phase 1 or 2), and the
residue of the reference protein at that column. Introns of different species at the same column and phase are
candidate shared introns. A translation of the model that differs from the protein of --proteins is reported
(column note) and left out.
Outputs: --out (TSV: model, id, intron, phase, residue, column, reference_residue, note) and --alignment (FASTA).
"""
import argparse, csv, os, subprocess, tempfile


def read_fasta(path):
    out, name = {}, None
    for l in open(path):
        if l.startswith(">"): name = l[1:].split()[0]; out[name] = []
        elif name: out[name].append(l.strip())
    return {k: "".join(v).replace("-", "").replace("*", "") for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positions", default="results/tables/pstk_intron_positions.tsv")
    ap.add_argument("--proteins", default="results/trees/pstk.faa")
    ap.add_argument("--gene", default="PSTK")
    ap.add_argument("--gold-dir", default="curation/gold_standard_machinery")
    ap.add_argument("--extra", nargs="*", default=["Acanthaster", "Apostichopus", "Lingula", "Crassostrea", "Capitella",
                                                   "Orbicella", "Exaiptasia", "Xenopus_tropicalis", "Alligator"])
    ap.add_argument("--out", default="results/tables/pstk_intron_columns.tsv")
    ap.add_argument("--alignment", default="results/trees/pstk_introns.aln")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.positions), delimiter="\t"))
    prot = read_fasta(a.proteins)
    seqs, ref = {}, None
    for r in rows:
        key = f"{r['model']}_{r['id']}"
        if key in seqs: continue
        if r["source"] == "gold":
            name, s = None, None
            for l in open(f"{a.gold_dir}/{r['model']}.proteins.faa"):
                if l.startswith(">"): name = l[1:].split("|")[0]
                elif name == r["id"]: s = l.strip()
        else:
            s = next((v for k, v in prot.items() if k.endswith(r["id"])), None)
        if s is None or len(s) != int(r["protein_length"]):
            seqs[key] = None; continue
        seqs[key] = s
        if r["source"] == "reference": ref = key
    extra = {k: v for k, v in prot.items() if any(k.startswith(x) for x in a.extra)}
    with tempfile.TemporaryDirectory() as d:
        with open(f"{d}/in.fa", "w") as o:
            for k, s in list(seqs.items()) + list(extra.items()):
                if s: o.write(f">{k}\n{s}\n")
        aln = subprocess.run(["mafft", "--quiet", "--localpair", "--maxiterate", "1000", f"{d}/in.fa"],
                             capture_output=True, text=True, check=True).stdout
    os.makedirs(os.path.dirname(a.alignment) or ".", exist_ok=True)
    open(a.alignment, "w").write(aln)
    msa = {}
    for block in aln.split(">")[1:]:
        lines = block.strip().split("\n"); msa[lines[0].split()[0]] = "".join(lines[1:])

    def column(s, idx):
        k = -1
        for c, ch in enumerate(s):
            if ch != "-":
                k += 1
                if k == idx: return c

    out = []
    for r in rows:
        if not r["intron"]: continue
        key = f"{r['model']}_{r['id']}"
        if not seqs.get(key):
            out.append(dict(model=r["model"], id=r["id"], intron=r["intron"], phase=r["phase"], residue=r["residue"],
                            column="", reference_residue="", note="protein differs from the translation of the model")); continue
        x = float(r["residue"]); n = int(x)
        c = column(msa[key], n)
        refres = sum(1 for ch in msa[ref][:c + 1] if ch != "-") if ref else ""
        out.append(dict(model=r["model"], id=r["id"], intron=r["intron"], phase=r["phase"], residue=r["residue"],
                        column=c + 1, reference_residue=refres, note=""))
    cols = ["model", "id", "intron", "phase", "residue", "column", "reference_residue", "note"]
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in sorted(out, key=lambda r: (int(r["column"] or 0), r["model"])): o.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"{len(out)} introns written to {a.out}; alignment {a.alignment}")


if __name__ == "__main__":
    main()
