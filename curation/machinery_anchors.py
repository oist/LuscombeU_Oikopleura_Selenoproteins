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
    predictions of the selenocysteine copy).
  - Several predictions of one gene in an assembly are named <gene>, <gene>_2, ... by decreasing protein length.
  - Frame: the coding rows of the prediction are joined in transcript order and translated; the frame with the fewest
    stop codons is used (the first one when tied), and the anchor is the sense codon nearest to the middle of that
    frame whose three bases are contiguous in the genome.
Output (--out): genome, gene, family, contig, anchor_codon (1-based lowest coordinate of the codon), strand, status
("prediction"), context (empty: no context check), prediction, label, protein_length, notes.
"""
import argparse, csv, gzip, os, re

FAMILY_GENE = {"pstk": "PSTK", "SecS": "SEPSECS", "eEFsec": "EEFSEC", "SBP2": "SBP2", "secp43": "TRNAU1AP",
               "SEPHS2": "SEPHS_b"}
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictions", default="results/tables/stage2_final_predictions.tsv")
    ap.add_argument("--gtf-dir", default="work/stage2/pipeline/final/selenoprofiles")
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3"])
    ap.add_argument("--out", default="curation/machinery_genes.tsv")
    a = ap.parse_args()
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
        named = {}
        for r in sorted(mine, key=lambda r: (r["family"], -int(r["protein_length"]))):
            lo, hi = int(r["start"]) - 1, int(r["end"]) - 1
            if r["label"] == "pseudo":
                left_out.append(f"{g} {r['prediction']}: labelled pseudo"); continue
            hit = [s for s in gold if s[0] == r["contig"] and s[1] == r["strand"] and lo <= s[3] and hi >= s[2]]
            if hit:
                left_out.append(f"{g} {r['prediction']}: overlaps a gold-standard selenoprotein gene"); continue
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
