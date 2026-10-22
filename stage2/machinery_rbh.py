#!/usr/bin/env python3
"""Reciprocal best hits of the machinery predictions, or of the selenoprotein predictions, against the human RefSeq
proteome.

For each Selenoprofiles prediction of a family of --families in each genome (by default the machinery genes and the
associated factor SECp43; --families selenoproteins for the selenoprotein families of SELENOPROTEIN_GENES):
  forward  blastp of the predicted protein against the human RefSeq proteins; the gene of each protein comes from the
           RefSeq annotation (CDS protein_id -> gene). The best hit gives the human gene; the margin (bits) is the
           score of the best expected gene minus that of the best other gene (closest paralog or related family).
           A selenoprotein family can have several expected human genes (paralogs, for example GPX1 to GPX8), or
           none when the family is absent from human (SELENOL, SELENOU).
  reverse  tblastn of the best human protein against the genome; its best locus must overlap the prediction.
Status:
  RBH         the best human gene is an expected one, and the reverse search returns to the prediction
  forward     the best human gene is an expected one, but the human protein hits another locus of the genome better
              (a paralog or a better copy elsewhere)
  other gene  the best human gene is not an expected one (for example SECISBP2L for SECISBP2)
  best hit only  the family has no expected human gene; the best hit and the reverse search are reported
  no hit      no forward hit
Small margins (below --ambiguous bits) are flagged for phylogeny.
Inputs: Selenoprofiles predictions (<dir>/<genome>.selenoprofiles.faa, headers with chromosome:, strand: and
positions:), data/genomes/Hsap/protein.faa.gz and annotation.gff3.gz, data/genomes/<genome>/genome.fna.gz.
Output: --table (one row per prediction).
"""
import argparse, collections, glob, gzip, os, re, shutil, subprocess, tempfile

MACHINERY_GENES = {"SEPHS2": ["SEPHS2"], "pstk": ["PSTK"], "SecS": ["SEPSECS"], "eEFsec": ["EEFSEC"],
                   "SBP2": ["SECISBP2"], "secp43": ["TRNAU1AP"]}
# Selenoprofiles selenoprotein families -> human genes of the family (HGNC symbols, as in the RefSeq annotation);
# empty when the family has no human member
SELENOPROTEIN_GENES = {
    "DIO": ["DIO1", "DIO2", "DIO3"],
    "GPX": ["GPX1", "GPX2", "GPX3", "GPX4", "GPX5", "GPX6", "GPX7", "GPX8"],
    "TXNRD": ["TXNRD1", "TXNRD2", "TXNRD3"],
    "SELENOW": ["SELENOW", "SELENOV"],
    "SELENOF": ["SELENOF"], "SELENOM": ["SELENOM"], "SELENOK": ["SELENOK"], "SELENOO": ["SELENOO"],
    "SELENOS": ["SELENOS"], "SELENOT": ["SELENOT"], "SELENOH": ["SELENOH"], "SELENOI": ["SELENOI"],
    "SELENON": ["SELENON"], "SELENOP": ["SELENOP"], "MSRB": ["MSRB1", "MSRB2", "MSRB3"], "MsrA": ["MSRA"],
    "SELENOL": [], "SELENOU": [], "SELENOJ": [], "SELENOE": [], "FrnE": [], "AhpC": [], "SelKi": []}
FAMILY_SETS = {"machinery": MACHINERY_GENES, "selenoproteins": SELENOPROTEIN_GENES}


def read_fasta(path):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = l[1:].rstrip("\n"), []
            else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def run(cmd, cwd): return subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=cwd).stdout


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--predictions", nargs="+", default=sorted(glob.glob("work/stage1/pipeline/selenoprofiles/*.selenoprofiles.faa")),
                    help="Selenoprofiles protein FASTA files, one per genome (<genome>.<anything>.faa or <genome>.faa)")
    ap.add_argument("--families", nargs="+", default=["machinery"],
                    help="families, or the sets 'machinery' and 'selenoproteins'")
    ap.add_argument("--human", default="data/genomes/Hsap")
    ap.add_argument("--genomes", default="data/genomes")
    ap.add_argument("--ambiguous", type=float, default=10.0)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--table", default="results/tables/stage2_machinery_rbh.tsv")
    a = ap.parse_args()
    EXPECTED = {}
    for f in a.families:
        if f in FAMILY_SETS: EXPECTED.update(FAMILY_SETS[f])
        else: EXPECTED[f] = {**MACHINERY_GENES, **SELENOPROTEIN_GENES}.get(f, [])
    d = tempfile.mkdtemp()
    # human proteins and their genes
    gene_of = {}
    with gzip.open(f"{a.human}/annotation.gff3.gz", "rt") as fh:
        for l in fh:
            if "\tCDS\t" not in l: continue
            m, g = re.search(r"protein_id=([^;\n]+)", l), re.search(r"gene=([^;\n]+)", l)
            if m and g: gene_of[m.group(1)] = g.group(1)
    missing = sorted({x for v in EXPECTED.values() for x in v} - set(gene_of.values()))
    if missing: print("WARNING: expected genes absent from the human annotation:", " ".join(missing), flush=True)
    with open(f"{d}/human.faa", "w") as o:
        for n, s in read_fasta(f"{a.human}/protein.faa.gz"): o.write(f">{n.split()[0]}\n{s}\n")
    human = {n.split()[0]: s for n, s in read_fasta(f"{d}/human.faa")}
    run(["formatdb", "-i", "human.faa", "-p", "T"], d)
    rows = []
    for faa in a.predictions:
        g = os.path.basename(faa).split(".")[0]
        preds = []
        for n, s in read_fasta(faa):
            pid = n.split()[0]; fam = pid.split(".")[0].removesuffix("_app")
            if fam not in EXPECTED: continue
            loc = dict(re.findall(r"(chromosome|strand|positions):(\S+)", n))
            pos = [int(x) for x in re.findall(r"\d+", loc.get("positions", ""))]
            preds.append((pid, fam, s.replace("*", "").replace("U", "X"), loc.get("chromosome", ""),
                          min(pos) if pos else 0, max(pos) if pos else 0))
        if not preds: continue
        with open(f"{d}/q.faa", "w") as o:
            for p in preds: o.write(f">{p[0]}\n{p[2]}\n")
        hits = collections.defaultdict(dict)   # query -> gene -> (bits, protein)
        for l in run(["blastall", "-p", "blastp", "-d", "human.faa", "-i", "q.faa", "-e", "1e-3", "-m", "8",
                      "-b", "200", "-v", "200", "-a", str(a.cpus)], d).splitlines():
            f = l.split("\t"); gene = gene_of.get(f[1], f[1]); bits = float(f[11])
            if bits > hits[f[0]].get(gene, (0, ""))[0]: hits[f[0]][gene] = (bits, f[1])
        # reverse: best human protein of each prediction against the genome
        genome = f"{a.genomes}/{g}/genome.fna.gz"
        reverse = {}
        best_prot = {p[0]: max(hits[p[0]].values())[1] for p in preds if hits.get(p[0])}
        if best_prot and os.path.exists(genome):
            with gzip.open(genome, "rb") as i, open(f"{d}/g.fa", "wb") as o: shutil.copyfileobj(i, o)
            run(["formatdb", "-i", "g.fa", "-p", "F"], d)
            with open(f"{d}/h.faa", "w") as o:
                for prot in sorted(set(best_prot.values())): o.write(f">{prot}\n{human[prot]}\n")
            for l in run(["blastall", "-p", "tblastn", "-d", "g.fa", "-i", "h.faa", "-e", "1e-3", "-m", "8",
                          "-a", str(a.cpus)], d).splitlines():
                f = l.split("\t"); bits = float(f[11])
                if bits > reverse.get(f[0], (0,))[0]:
                    reverse[f[0]] = (bits, f[1], min(int(f[8]), int(f[9])), max(int(f[8]), int(f[9])))
        for pid, fam, s, chrom, start, end in preds:
            h = hits.get(pid, {})
            expected = EXPECTED[fam]
            exp_text = ",".join(expected) or "-"
            if not h:
                rows.append((g, pid, fam, exp_text, "-", "-", "-", "-", "-", "-", "no hit")); continue
            ranked = sorted(h.items(), key=lambda x: -x[1][0])
            best_gene, (best_bits, best_p) = ranked[0]
            # reference gene for the margin: the best expected gene, or the best gene when none is expected
            ref = max(expected, key=lambda x: h.get(x, (0, ""))[0]) if expected else best_gene
            exp_bits = h.get(ref, (0, ""))[0]
            other = max((v[0] for k, v in h.items() if k != ref), default=0)
            other_gene = next((k for k, v in ranked if k != ref), "-")
            margin = exp_bits - other
            r = reverse.get(best_p)
            back = bool(r and r[1] == chrom and r[2] <= end + 1000 and r[3] >= start - 1000)
            reverse_run = os.path.exists(genome)
            if not expected: status = "best hit only" + ("" if back or not reverse_run else ", reverse elsewhere")
            elif best_gene not in expected: status = "other gene"
            elif back: status = "RBH"
            elif reverse_run: status = "forward"
            else: status = "forward (no reverse search: genome not available)"
            if status != "other gene" and margin < a.ambiguous: status += " (ambiguous)"
            rows.append((g, pid, fam, exp_text, best_gene, f"{best_bits:.0f}", other_gene, f"{margin:.1f}",
                         f"{r[1]}:{r[2]}-{r[3]}" if r else "-", ("yes" if back else "no") if reverse_run else "-", status))
        print(g, len(preds), "predictions", flush=True)
    os.makedirs(os.path.dirname(a.table) or ".", exist_ok=True)
    with open(a.table, "w") as o:
        o.write("genome\tprediction\tfamily\texpected_gene\tbest_human_gene\tbest_bits\tnext_gene\tmargin\t"
                "reverse_best_locus\treverse_returns\tstatus\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    shutil.rmtree(d)


if __name__ == "__main__":
    main()
