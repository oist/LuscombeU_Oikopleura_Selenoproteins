#!/usr/bin/env python3
"""Machinery genes of the gold-standard assemblies, as input of curation/gold_standard.py --anchor sense.

For each machinery prediction of the final stage 2 run (results/tables/stage2_final_predictions.tsv) in the
gold-standard assemblies, one row: the gene name, the prediction, and an anchor codon, a sense codon in the middle
of the coding sequence of the prediction (<gtf-dir>/<genome>.selenoprofiles.gtf), read on the genome. The gold-standard
script builds the model of the gene around this codon, as it does around the selenocysteine codon of the
selenoprotein genes.
  - Families and gene names: pstk PSTK, SecS SEPSECS, eEFsec EEFSEC, SBP2 SBP2 (the SECIS-binding protein, co-ortholog
    of SECISBP2 and SECISBP2L), secp43 TRNAU1AP, SEPHS2 SEPHS_b (the SEPHS copy without selenocysteine: the
    selenocysteine copy, SEPHS2, is a gene of the selenoprotein gold standard).
  - Left out: predictions labelled pseudo (by Selenoprofiles: frameshift or in-frame stop codon), and predictions
    whose coding span overlaps a gene of the selenoprotein gold standard on the same strand (--gold: the SEPHS2
    predictions of the selenocysteine copy). In genomes without a selenoprotein gold standard (O. longicauda,
    O. albicans, O. vanhoeffeni), the SEPHS2 predictions kept as SEPHS_b are those labelled with the residue found at
    the selenocysteine position (alanine, serine, cysteine...): predictions labelled selenocysteine or unaligned may be
    the selenocysteine copy, which stays a prediction in the trees.
  - Several predictions of one gene in an assembly are named <gene>, <gene>_2, ... by decreasing protein length;
    predictions shorter than --min-fraction (default 0.5) of the longest one of the gene in the assembly are left out
    (fragments on short contigs, for example six short eEFsec predictions in O. longicauda).
  - Frame: the coding rows of the prediction are joined in transcript order and translated; the frame with the fewest
    stop codons (a terminal stop codon not counted) is used (the first one when tied), and the anchor is the sense
    codon nearest to the middle of that frame whose three bases are contiguous in the genome.
Output (--out): genome, gene, family, contig, anchor_codon (1-based lowest coordinate of the codon), strand, status
("prediction"), context (empty: no context check), prediction, label, protein_length, notes.
With --human, also the reference proteins of each gene (--references): the human RefSeq proteins of its genes
(REFERENCE_GENES; longest NP_ isoform, else the longest isoform; selenocysteine written X), header
<gene>|Human_<symbol>|<protein id>; they rank the candidate models of curation/gold_standard.py --anchor sense and
are projected on the genomes by curation/machinery_projections.py.
"""
import argparse, csv, gzip, os, re

FAMILY_GENE = {"pstk": "PSTK", "SecS": "SEPSECS", "eEFsec": "EEFSEC", "SBP2": "SBP2", "secp43": "TRNAU1AP",
               "SEPHS2": "SEPHS_b"}
REFERENCE_GENES = {"PSTK": ["PSTK"], "SEPSECS": ["SEPSECS"], "EEFSEC": ["EEFSEC"], "SBP2": ["SECISBP2", "SECISBP2L"],
                   "TRNAU1AP": ["TRNAU1AP"], "SEPHS_b": ["SEPHS1", "SEPHS2"]}
RESIDUE_LABELS = {"alanine", "serine", "cysteine", "glycine", "threonine", "arginine", "leucine", "valine",
                  "glutamine", "asparagine", "aspartate", "glutamate", "lysine", "histidine", "isoleucine",
                  "methionine", "phenylalanine", "proline", "tryptophan", "tyrosine"}
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s): return s.translate(COMP)[::-1]


def read_fasta(path, keep):
    out, name, buf = {}, None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name in keep: out[name] = "".join(buf).upper()
                name, buf = l[1:].split()[0], []
            elif name in keep: buf.append(l.strip())
    if name in keep: out[name] = "".join(buf).upper()
    return out


def gtf_cds(path):
    """prediction id -> [(start0, end0 inclusive)] of its CDS rows, and its contig and strand."""
    out = {}
    if not os.path.exists(path): return out
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        m = re.search(r'transcript_id "([^"]+)"', f[8]) or re.search(r'gene_id "([^"]+)"', f[8])
        r = out.setdefault(m.group(1), {"contig": f[0], "strand": f[6], "rows": []})
        r["rows"].append((int(f[3]) - 1, int(f[4]) - 1))
    return out


def gold_spans(path):
    """(contig, strand, lo0, hi0) of the coding spans of the genes of a gold-standard GFF3."""
    spans = {}
    if not os.path.exists(path): return []
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        k = re.search(r"Parent=([^;]+)", f[8]).group(1)
        s = spans.setdefault(k, [f[0], f[6], int(f[3]) - 1, int(f[4]) - 1])
        s[2], s[3] = min(s[2], int(f[3]) - 1), max(s[3], int(f[4]) - 1)
    return [tuple(v) for v in spans.values()]


def anchor(genome, rows, strand):
    """Lowest 0-based coordinate of the anchor codon, and a note on the frame, or (None, reason)."""
    pos = [p for a, b in sorted(rows) for p in range(a, b + 1)]
    if strand == "-": pos = pos[::-1]
    seq = "".join(genome[p] for p in pos)
    if strand == "-": seq = seq.translate(COMP)
    stops = [sum(seq[i:i + 3] in STOPS for i in range(fr, len(seq) - 2, 3)) for fr in range(3)]
    if len(seq) >= 3 and seq[-3:] in STOPS and len(seq) % 3 == 0: stops[0] -= 1   # terminal stop codon of frame 0
    fr = min(range(3), key=lambda k: (stops[k], k))
    codons = list(range(fr, len(seq) - 2, 3))
    mid = len(codons) // 2
    for i in sorted(codons, key=lambda c: abs(codons.index(c) - mid)):
        p3 = pos[i:i + 3]
        if seq[i:i + 3] in STOPS or seq[i:i + 3].count("N"): continue
        if p3 != list(range(p3[0], p3[0] + 3)) and p3 != list(range(p3[0], p3[0] - 3, -1)): continue
        note = "" if fr == 0 and stops[0] == 0 else f"frame {fr} of the joined coding rows, {stops[fr]} stop codon(s)"
        return min(p3), note
    return None, "no sense codon with contiguous bases"


def write_references(human, out):
    gene_of = {}
    with gzip.open(f"{human}/annotation.gff3.gz", "rt") as fh:
        for l in fh:
            if "\tCDS\t" not in l: continue
            m, g = re.search(r"protein_id=([^;\n]+)", l), re.search(r"gene=([^;\n]+)", l)
            if m and g: gene_of[m.group(1)] = g.group(1)
    prot, name, buf = {}, None, []
    with gzip.open(f"{human}/protein.faa.gz", "rt") as fh:
        for l in fh:
            if l.startswith(">"):
                if name: prot[name] = "".join(buf)
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
    if name: prot[name] = "".join(buf)
    with open(out, "w") as o:
        for gene, symbols in REFERENCE_GENES.items():
            for sym in symbols:
                ids = [p for p, gg in gene_of.items() if gg == sym and p in prot]
                if not ids:
                    print(f"WARNING: no human RefSeq protein for {sym}"); continue
                p = max(ids, key=lambda x: (x.startswith("NP_"), len(prot[x])))
                o.write(f">{gene}|Human_{sym}|{p}\n{prot[p].replace('U', 'X').rstrip('*')}\n")
    print(f"references written to {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictions", default="results/tables/stage2_final_predictions.tsv")
    ap.add_argument("--gtf-dir", default="work/stage2/pipeline/final/selenoprofiles")
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3"])
    ap.add_argument("--out", default="curation/machinery_genes.tsv")
    ap.add_argument("--min-fraction", type=float, default=0.5)
    ap.add_argument("--human", default="", help="directory with the human RefSeq protein.faa.gz and annotation.gff3.gz")
    ap.add_argument("--references", default="curation/machinery_references.faa")
    a = ap.parse_args()
    if a.human: write_references(a.human, a.references)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    preds = [r for r in csv.DictReader(open(a.predictions), delimiter="\t")
             if r["genome"] in a.genomes and r["family"] in FAMILY_GENE]
    cols = ["genome", "gene", "family", "contig", "anchor_codon", "strand", "status", "context", "prediction", "label",
            "protein_length", "notes"]
    rows, left_out = [], []
    for g in a.genomes:
        mine = [r for r in preds if r["genome"] == g]
        if not mine: continue
        cds = gtf_cds(f"{a.gtf_dir}/{g}.selenoprofiles.gtf")
        gold = gold_spans(f"{a.gold}/{g}.gff3")
        genome = read_fasta(sheet[g]["fasta"], {r["contig"] for r in mine})
        named, longest = {}, {}
        for r in mine:
            if r["label"] != "pseudo":
                longest[r["family"]] = max(longest.get(r["family"], 0), int(r["protein_length"]))
        for r in sorted(mine, key=lambda r: (r["family"], -int(r["protein_length"]))):
            lo, hi = int(r["start"]) - 1, int(r["end"]) - 1
            if r["label"] == "pseudo":
                left_out.append(f"{g} {r['prediction']}: labelled pseudo"); continue
            hit = [s for s in gold if s[0] == r["contig"] and s[1] == r["strand"] and lo <= s[3] and hi >= s[2]]
            if hit:
                left_out.append(f"{g} {r['prediction']}: overlaps a gold-standard selenoprotein gene"); continue
            if r["family"] == "SEPHS2" and not gold and r["label"] not in RESIDUE_LABELS:
                left_out.append(f"{g} {r['prediction']}: SEPHS2 labelled {r['label']} in a genome without "
                                f"selenoprotein gold standard (selenocysteine copy possible)"); continue
            if int(r["protein_length"]) < a.min_fraction * longest[r["family"]]:
                left_out.append(f"{g} {r['prediction']}: {r['protein_length']} aa, shorter than {a.min_fraction} of the "
                                f"longest prediction of the family ({longest[r['family']]} aa)"); continue
            c = cds.get(r["prediction"])
            if not c:
                left_out.append(f"{g} {r['prediction']}: no CDS rows in {a.gtf_dir}/{g}.selenoprofiles.gtf"); continue
            p0, note = anchor(genome[c["contig"]], c["rows"], c["strand"])
            if p0 is None:
                left_out.append(f"{g} {r['prediction']}: {note}"); continue
            base = FAMILY_GENE[r["family"]]
            k = named[base] = named.get(base, 0) + 1
            rows.append({"genome": g, "gene": base if k == 1 else f"{base}_{k}", "family": base, "contig": c["contig"],
                         "anchor_codon": p0 + 1, "strand": c["strand"], "status": "prediction", "context": "",
                         "prediction": r["prediction"], "label": r["label"], "protein_length": r["protein_length"],
                         "notes": note or "-"})
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    for x in left_out: print("left out:", x)
    print(f"{len(rows)} genes written to {a.out}")


if __name__ == "__main__":
    main()
