#!/usr/bin/env python3
"""Table MACHINERY of the manuscript: the census of the selenocysteine machinery per genome with the evidence for
the orthology of each gene.

Inputs (all in the repository):
  results/tables/stage2_machinery_census.tsv   presence per gene and genome (stage2/machinery_census.py)
  results/tables/stage2_machinery_rbh.tsv      reciprocal best hits against the human proteome (machinery_rbh.py)
  results/tables/stage2_machinery_lba.tsv      supports and AU tests of the trees (machinery_lba.py)
  results/trees/SEPHS2.{treefile,groups.tsv}, results/trees/C20/SEPHS2.treefile, results/trees/lba/SEPHS2.*.treefile
Rows added to the census:
  "SEPHS, second copy": appendicularian SEPHS predictions outside the selenocysteine clade of the maximum-likelihood
  tree and without selenocysteine (the copy that groups with the tunicate glycine forms); a genome has it when such a
  prediction exists.
The orthology column is written from these results: reciprocal best hits for SEPSECS, EEFSEC and TRNAU1AP; the
trees and AU tests for PSTK, SEPHS and the SECIS-binding protein. The SECIS-binding protein row is named
SECISBP2/SECISBP2L (co-ortholog of both vertebrate paralogs).
Outputs: --markdown and --table.
"""
import argparse, collections, csv, glob, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import machinery_lba as lba   # taxon assignment, tree reading, appendicularian clade

GENOMES = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]
NAMES = {"OKI": "*O. dioica* OKI2018", "OSKA": "*O. dioica* OSKA2016", "OdB3": "*O. dioica* OdB3",
         "Olon": "*O. longicauda*", "Oalb": "*O. albicans*", "Ovan": "*O. vanhoeffeni*", "Bsty": "*B. stygius*",
         "Mery": "*M. erythrocephalus*", "Fbor": "*F. borealis*"}
RBH_FAMILY = {"SEPSECS": "SecS", "EEFSEC": "eEFsec", "TRNAU1AP": "secp43"}
ROW_NAMES = {"SECISBP2": "SECISBP2/SECISBP2L", "TRNAU1AP": "TRNAU1AP (associated factor)"}


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def p_au(text):
    m = re.search(r"=\s*([\d.eE-]+)", text)
    return float(m.group(1)) if m else None


def fmt_p(p):
    return f"{p:.2g}"   # two significant digits


def rbh_statement(rbh, family):
    by_genome = collections.defaultdict(set)
    for r in rbh:
        if r["family"] == family: by_genome[r["genome"]].add(r["status"].split(" (")[0])
    with_rbh = [g for g in GENOMES if "RBH" in by_genome.get(g, ())]
    forward = [g for g in GENOMES if g not in with_rbh and "forward" in by_genome.get(g, ())]
    if not with_rbh: return "no reciprocal best hit"
    if len(with_rbh) == 1 and not forward: return f"reciprocal best hit ({NAMES[with_rbh[0]]})"
    text = "reciprocal best hit in all nine genomes" if len(with_rbh) == len(GENOMES) else "reciprocal best hit"
    if forward: text += "; best hit only in " + " and ".join(NAMES[g] for g in forward)
    return text


def sephs_second_copy(trees):
    """Genomes with a second SEPHS copy, and the range of supports of the glycine-form clade it falls in."""
    groups = dict(l.rstrip("\n").split("\t") for l in open(f"{trees}/SEPHS2.groups.tsv"))
    taxa = lba.load_taxa("config/taxon_groups.tsv")
    ml = f"{trees}/SEPHS2.treefile"
    names = {x.name for x in lba.read_tree(ml).get_terminals()}
    tax = {n: lba.taxon(n, groups, taxa) for n in names}
    sec_clade = lba.appendicularian_clade(ml, tax)
    second = {n for n in names if tax[n] == "Appendicularians" and n not in sec_clade
              and groups.get(n) == "appendicularian" and not n.endswith(".selenocysteine")}
    genomes = {n.split("_")[0] for n in second}
    glycine = {n for n in names if tax[n] == "Other tunicates" and re.search(lba.ORTHOLOG_EXCLUDE["SEPHS2"], n)}
    clade = glycine | second
    supports = []
    files = [ml, f"{trees}/C20/SEPHS2.treefile"] + sorted(glob.glob(f"{trees}/lba/SEPHS2.pmsf.treefile")) + \
        sorted(glob.glob(f"{trees}/lba/SEPHS2.slowfast*.treefile"))
    for f in files:
        if not os.path.exists(f): continue
        t = lba.read_tree(f)
        leaves = {x.name for x in t.get_terminals()}
        side = clade & leaves
        for c in t.get_nonterminals():
            s = {x.name for x in c.get_terminals()}
            if (s == side or s == leaves - side) and c.confidence is not None:
                supports.append(c.confidence); break
    return genomes, supports


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--census", default="results/tables/stage2_machinery_census.tsv")
    ap.add_argument("--rbh", default="results/tables/stage2_machinery_rbh.tsv")
    ap.add_argument("--lba", default="results/tables/stage2_machinery_lba.tsv")
    ap.add_argument("--trees", default="results/trees")
    ap.add_argument("--markdown", default="results/tables/stage2_machinery_table.md")
    ap.add_argument("--table", default="results/tables/stage2_machinery_table.tsv")
    a = ap.parse_args()
    census = rows(a.census)
    rbh = rows(a.rbh)
    lba_rows = rows(a.lba)
    lt = collections.defaultdict(dict)   # family -> analysis -> row (the first row of a name; AU tests below)
    for r in lba_rows: lt[r["family"]].setdefault(r["analysis"], r)
    def au_values(fam):   # all AU tests of a family, under every model
        return [p_au(r["appendicularian clade"]) for r in lba_rows if r["family"] == fam and r["analysis"].startswith("AU test")]
    second_genomes, gly = sephs_second_copy(a.trees)

    pstk_out = min(int(r["family vs outgroup"]) for k, r in lt["pstk"].items()
                   if not k.startswith("AU") and r["family vs outgroup"].isdigit())
    orthology = {
        "tRNA-Sec": "not tested (RNA gene)",
        "SEPHS2": f"selenocysteine at the catalytic site; outside the tunicate SEPHS in the trees "
                  f"(AU p at most {fmt_p(max(au_values('SEPHS2')))})",
        "SEPHS, second copy": f"ortholog of the tunicate glycine form ({min(gly):.0f} to {max(gly):.0f}%)" if gly else
                              "not placed",
        "PSTK": f"PSTK, not KTI12 ({pstk_out}%); position within PSTK unresolved "
                f"(AU p = {fmt_p(au_values('pstk')[0])})",
        "SECISBP2": f"co-ortholog of SECISBP2 and SECISBP2L; with the other tunicates "
                    f"({lt['SBP2']['reference']['appendicularians+tunicates']}%, "
                    f"C60 {lt['SBP2']['PMSF LG+C60+F+G']['appendicularians+tunicates']}%)",
    }
    for gene, fam in RBH_FAMILY.items(): orthology[gene] = rbh_statement(rbh, fam)

    out = []
    for r in census:
        out.append([ROW_NAMES.get(r["gene"], r["gene"]), orthology.get(r["gene"], "")] + [r[g] for g in GENOMES])
        if r["gene"] == "SEPHS2":
            out.append(["SEPHS, second copy", orthology["SEPHS, second copy"]] +
                       ["found" if g in second_genomes else "not found" for g in GENOMES])
    header = ["Gene", "Orthology"] + [NAMES[g] for g in GENOMES]
    with open(a.table, "w") as o:
        o.write("\t".join(["gene", "orthology"] + GENOMES) + "\n")
        for r in out: o.write("\t".join(r) + "\n")
    with open(a.markdown, "w") as o:
        o.write("| " + " | ".join(header) + " |\n|" + "---|" * len(header) + "\n")
        for r in out: o.write("| " + " | ".join(r) + " |\n")
    print(open(a.markdown).read())


if __name__ == "__main__":
    main()
