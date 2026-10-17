#!/usr/bin/env python3
"""Census of the selenocysteine machinery per genome, from:
  - results/tables/stage1_trnasec.tsv (tRNA-Sec loci supported by two tools and G73),
  - the predictions of the appendicularian machinery profiles (work/stage2/machinery_test/<genome>.gtf),
  - the fragment analysis (work/stage2/machinery_fragments/fragments.tsv).
Status per gene and genome:
  selenocysteine     SEPHS2 predicted with a selenocysteine codon
  found              predicted (homologue label, or for SEPHS2 another residue at the selenocysteine position)
  partial            predicted only as "pseudo" or "unaligned" (in-frame stop, frameshift or selenocysteine region
                     not aligned)
  fragment           dropped by Selenoprofiles but supported by the fragment analysis
  not found          none of the above
Outputs: a TSV (--table) and a Markdown table (--markdown).
"""
import argparse, collections, csv, os, re

GENOMES = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]
NAMES = {"OKI": "*O. dioica* OKI2018", "OSKA": "*O. dioica* OSKA2016", "OdB3": "*O. dioica* OdB3",
         "Olon": "*O. longicauda*", "Oalb": "*O. albicans*", "Ovan": "*O. vanhoeffeni*", "Bsty": "*B. stygius*",
         "Mery": "*M. erythrocephalus*", "Fbor": "*F. borealis*"}
GENES = [("tRNA-Sec", None), ("SEPHS2", "SEPHS2"), ("PSTK", "pstk"), ("SEPSECS", "SecS"), ("EEFSEC", "eEFsec"),
         ("SECISBP2", "SBP2"), ("TRNAU1AP", "secp43")]
ASSOCIATED = {"TRNAU1AP"}  # SECp43: associated factor, not part of the machinery (not required, not specific)


def labels(gtf):
    c = collections.defaultdict(collections.Counter)
    seen = set()
    for l in open(gtf):
        m = re.search(r'gene_id "([^":]+)"', l)
        if not m or m.group(1) in seen: continue
        seen.add(m.group(1))
        c[m.group(1).split(".")[0].removesuffix("_app")][m.group(1).rsplit(".", 1)[1]] += 1
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trnasec", default="results/tables/stage1_trnasec.tsv")
    ap.add_argument("--test-dir", default="work/stage2/machinery_test")
    ap.add_argument("--fragments", default="work/stage2/machinery_fragments/fragments.tsv")
    ap.add_argument("--table", default="results/tables/stage2_machinery_census.tsv")
    ap.add_argument("--markdown", default="results/tables/stage2_machinery_census.md")
    a = ap.parse_args()
    trna = {r["genome"]: int(r["loci_supported"] or 0) for r in csv.DictReader(open(a.trnasec), delimiter="\t")}
    frag = collections.defaultdict(set)
    if os.path.exists(a.fragments):
        for r in csv.DictReader(open(a.fragments), delimiter="\t"):
            if r["status"] == "fragment of the gene": frag[r["genome"]].add(r["family"])
    status = {}
    for g in GENOMES:
        lab = labels(f"{a.test_dir}/{g}.gtf") if os.path.exists(f"{a.test_dir}/{g}.gtf") else {}
        for gene, fam in GENES:
            if fam is None:
                n = trna.get(g, 0); status[(gene, g)] = f"found ({n} loci)" if n > 1 else ("found" if n else "not found")
                continue
            c = lab.get(fam, collections.Counter())
            if fam == "SEPHS2" and c.get("selenocysteine"): s = "selenocysteine"
            elif any(k not in ("pseudo", "unaligned") for k in c): s = "found"
            elif c: s = "partial"
            elif fam in frag[g]: s = "fragment"
            else: s = "not found"
            status[(gene, g)] = s
    with open(a.table, "w") as o:
        o.write("gene\t" + "\t".join(GENOMES) + "\n")
        for gene, _ in GENES: o.write(gene + "\t" + "\t".join(status[(gene, g)] for g in GENOMES) + "\n")
    with open(a.markdown, "w") as o:
        o.write("| Gene | " + " | ".join(NAMES[g] for g in GENOMES) + " |\n")
        o.write("|---|" + "---|" * len(GENOMES) + "\n")
        for gene, _ in GENES:
            label = f"{gene} (associated factor)" if gene in ASSOCIATED else gene
            o.write(f"| {label} | " + " | ".join(status[(gene, g)] for g in GENOMES) + " |\n")
    print(open(a.markdown).read())


if __name__ == "__main__":
    main()
