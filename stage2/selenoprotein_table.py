#!/usr/bin/env python3
"""Table SELENO-ORTHOLOGY of the manuscript: the evidence for the orthology of each gold-standard gene of O. dioica
(round 2 of step 8: trees of the stage 2 final run with the gold-standard proteins).

Inputs (all in the repository):
  results/tables/stage2_selenoprotein_orthology.tsv  appendicularian clades of the ML and PMSF trees
                                                     (selenoprotein_orthology.py)
  results/tables/stage2_selenoprotein_rbh.tsv        best human hits of the final-run predictions (machinery_rbh.py)
  results/trees/selenoproteins/<family>.removed.tsv  final-run predictions replaced by each gold-standard protein
For each gene (gold-standard name, one row; the three assemblies pooled):
  best human hit   best human gene of the final-run predictions that the gold-standard proteins replaced, with the
                   number of assemblies, how many of them are reciprocal, and the margins to the next human gene
                   (stage2/machinery_rbh.py); for a best hit outside the human genes of the family (MIEN1 for the
                   chromosome 2 SELENOW), how far it scores above them; "none" when the human proteome gave no hit;
                   assemblies without such a prediction are not counted
  clade            the appendicularian clade of the ML tree that holds the gene: sequences, genomes, ultrafast
                   bootstrap support under the ML model and under LG+C60+F+G (PMSF; same leaves, or the PMSF clade
                   holding the gene when the leaves differ, marked *)
  sister group     taxonomic groups of the sister group (unrooted reading of selenoprotein_orthology.py), ML; PMSF
  vertebrate genes human genes of the nearest side of the ML tree that holds vertebrate sequences
Outputs: --markdown and --table.
"""
import argparse, collections, csv, glob, os, re

GENOMES = {"OKI": "OKI2018", "OSKA": "OSKA2016", "OdB3": "OdB3"}
ORDER = ["DIO_PAR", "DIO_chr1", "GPX_PAR", "GPX_chr1a", "GPX_chr1b", "GPX_chr2", "SELENOF", "SELENOK", "SELENOO",
         "SELENOS", "SELENOT", "SELENOW_chr1a", "SELENOW_pair1", "SELENOW_pair2", "SELENOW_PAR", "SELENOW_chr2",
         "TXNRD"]


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def replaced(trees):
    """(genome, gene) -> final-run predictions that the gold-standard protein replaced."""
    out = collections.defaultdict(list)
    for f in glob.glob(f"{trees}/*.removed.tsv"):
        for r in rows(f):
            m = re.fullmatch(r"gold-standard gene (\S+) used", r["reason"])
            if m:
                g, _, p = r["sequence"].partition("_")
                out[(g, m.group(1))].append(p)
    return out


def hit_text(hits):
    """hits: [(gene, status, margin)] of the assemblies -> text."""
    if not hits: return "not in the best-hit table"
    genes = collections.Counter(h[0] for h in hits)
    out = []
    for gene, n in genes.most_common():
        sub = [h for h in hits if h[0] == gene]
        if gene in ("-", ""):
            out.append(f"none ({n})"); continue
        rec = sum(h[1].startswith("RBH") for h in sub)
        margins = sorted(abs(float(h[2])) for h in sub if h[2] not in ("-", ""))
        span = (f"{margins[0]:g}" + (f" to {margins[-1]:g}" if margins[-1] != margins[0] else "")) if margins else ""
        if all(h[1].startswith("other gene") for h in sub):   # best hit outside the genes of the family
            out.append(f"{gene} ({n}; outside the family" + (f", {span} bits above its genes" if span else "") + ")")
            continue
        out.append(f"{gene} ({n}; {rec} reciprocal" + (f"; margin {span} bits" if span else "") + ")")
    return "; ".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orthology", default="results/tables/stage2_selenoprotein_orthology.tsv")
    ap.add_argument("--rbh", default="results/tables/stage2_selenoprotein_rbh.tsv")
    ap.add_argument("--trees", default="results/trees/selenoproteins")
    ap.add_argument("--markdown", default="results/tables/stage2_selenoprotein_table.md")
    ap.add_argument("--table", default="results/tables/stage2_selenoprotein_table.tsv")
    a = ap.parse_args()
    rbh = {(r["genome"], r["prediction"]): r for r in rows(a.rbh)}
    rep = replaced(a.trees)
    orth = rows(a.orthology)
    out = []
    genes = sorted({g.split(":", 1)[1] for r in orth for g in r["curated O. dioica genes"].split() if ":" in g},
                   key=lambda g: (ORDER.index(g) if g in ORDER else 99, g))
    for gene in genes:
        def clade(tree):
            for r in orth:
                if r["tree"] == tree and any(x.split(":", 1)[1] == gene for x in r["curated O. dioica genes"].split()
                                             if ":" in x):
                    return r
        ml, pm = clade("ML"), clade("PMSF")
        hits = []
        for g in GENOMES:
            for p in rep.get((g, gene), []):
                r = rbh.get((g, p))
                if r: hits.append((r["best_human_gene"], r["status"], r["margin"]))
        same = pm is not None and ml is not None and pm["members"] == ml["members"]
        support = (f"{ml['clade support']} / {pm['clade support']}{'' if same else '*'}" if ml and pm else
                   ml["clade support"] if ml else "-")
        out.append({
            "gene": gene,
            "assemblies": ", ".join(GENOMES[g] for g in GENOMES
                                    if any(f"{g}:{gene}" == x for r in orth if r["tree"] == "ML"
                                           for x in r["curated O. dioica genes"].split())),
            "best human hit": hit_text(hits),
            "clade sequences": ml["sequences"] if ml else "-",
            "clade genomes": len(ml["genomes"].split(",")) if ml else "-",
            "clade support (ML / PMSF)": support,
            "sister group (ML)": f"{ml['sister group']} ({ml['sister support']})" if ml else "-",
            "sister group (PMSF)": f"{pm['sister group']} ({pm['sister support']})" if pm else "-",
            "nearest vertebrate genes (ML)": ml["nearest vertebrate genes"] if ml else "-",
        })
    cols = list(out[0])
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in out: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(a.markdown, "w") as o:
        o.write("| Gene | Best human hit | Appendicularian clade (sequences, genomes) | Support ML / PMSF | "
                "Sister group, ML | Sister group, PMSF | Nearest vertebrate genes |\n|---|---|---|---|---|---|---|\n")
        for r in out:
            o.write(f"| {r['gene']} | {r['best human hit']} | {r['clade sequences']}, {r['clade genomes']} | "
                    f"{r['clade support (ML / PMSF)']} | {r['sister group (ML)']} | {r['sister group (PMSF)']} | "
                    f"{r['nearest vertebrate genes (ML)']} |\n")
    print(open(a.markdown).read())


if __name__ == "__main__":
    main()
