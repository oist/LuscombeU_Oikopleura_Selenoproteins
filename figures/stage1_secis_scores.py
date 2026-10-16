#!/usr/bin/env python3
"""Stage 1 SECIS scores: are the selenoprotein 3' regions below the RF00031 gathering threshold, or inside the
background?

For each genome, scores with Rfam RF00031 (Infernal cmsearch, default filters, all hits of at least 0 bits):
  - the 3' sequences written by Selenoprofiles after each prediction, separately for predictions labelled
    selenocysteine and for the other predictions (cysteine homologs, unaligned, machinery...);
  - background: random genomic windows of the same length (6,000 nt), on a random strand.
The best score per sequence goes to results/tables/stage1_secis_scores.tsv; the figure shows its distribution per
genome, with the gathering threshold. Sequences without any reported hit are drawn below 0.
"""
import argparse, gzip, os, random, subprocess, tempfile
from collections import defaultdict

GENOMES = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]
NAMES = {"OKI": "O. dioica OKI2018", "OSKA": "O. dioica OSKA2016", "OdB3": "O. dioica OdB3", "Olon": "O. longicauda",
         "Oalb": "O. albicans", "Ovan": "O. vanhoeffeni", "Bsty": "B. stygius", "Mery": "M. erythrocephalus",
         "Fbor": "F. borealis"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = l[1:].rstrip("\n"), []
            else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def label(header):
    """Selenoprofiles label: last dotted field of the prediction name (region_of_N_bp_downstream_of:FAM.k.label)."""
    return header.split()[0].split(":", 1)[1].rsplit(".", 1)[1]


def background(genome, n, length, seed):
    rng = random.Random(seed)
    contigs = [(name.split()[0], seq) for name, seq in read_fasta(genome) if len(seq) >= length]
    total = sum(len(s) for _, s in contigs)
    out = []
    for i in range(n):
        x = rng.randrange(total)
        for name, seq in contigs:
            if x < len(seq): break
            x -= len(seq)
        start = min(x, len(seq) - length)
        w = seq[start:start + length]
        if rng.random() < 0.5: w = w.translate(COMP)[::-1]
        out.append((f"bg{i}_{name}_{start}", w))
    return out


def curated(table, genome):
    """6,000 nt downstream of each curated Sec codon (strand-aware), from the genome."""
    import csv
    rows = [r for r in csv.DictReader(open(table), delimiter="\t") if r["genome"] == genome and r["status"] == "sec"]
    if not rows: return []
    want = {r["contig"] for r in rows}
    seqs = {n.split()[0]: s for n, s in read_fasta(f"data/genomes/{genome}/genome.fna.gz") if n.split()[0] in want}
    out = []
    for r in rows:
        q, p = seqs[r["contig"]], int(r["sec_codon"])
        w = q[p - 1:p - 1 + 6000] if r["strand"] == "+" else q[max(0, p + 2 - 6000):p + 2].translate(COMP)[::-1]
        out.append((f"curated:{r['gene']}", w))
    return out


def best_scores(cm, records, cpus):
    with tempfile.TemporaryDirectory() as d:
        fa, tbl = f"{d}/in.fa", f"{d}/out.tbl"
        with open(fa, "w") as o:
            for i, (_, s) in enumerate(records): o.write(f">r{i}\n{s}\n")
        subprocess.run(["cmsearch", "--toponly", "-T", "0", "--cpu", str(cpus), "--tblout", tbl, cm, fa],
                       check=True, stdout=subprocess.DEVNULL)
        best = {}
        for l in open(tbl):
            if l.startswith("#"): continue
            f = l.split(); i = int(f[0][1:]); sc = float(f[14])
            if sc > best.get(i, (-1,))[0]: best[i] = (sc, int(f[7]), int(f[8]))
    return [best.get(i, (None, None, None)) for i in range(len(records))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline", default="work/stage1/pipeline")
    ap.add_argument("--cm", default="data/references/RF00031.cm")
    ap.add_argument("--background", type=int, default=200, help="random windows per genome")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--table", default="results/tables/stage1_secis_scores.tsv")
    ap.add_argument("--figure", default="figures/stage1_secis_scores")
    ap.add_argument("--model-name", default="RF00031")
    ap.add_argument("--baseline", default="", help="curation table: adds the 3' sequences of the curated Sec genes")
    ap.add_argument("--seed", type=int, default=0, help="offset for the random windows")
    ap.add_argument("--replot", action="store_true", help="only redraw the figure from --table")
    a = ap.parse_args()
    if a.replot:
        ga = next((float(l.split()[1]) for l in open(a.cm) if l.startswith("GA ")), 0.0)
        rows = [(r["genome"], r["group"], r["sequence"], int(r["length"]), float(r["best_bits"]) if r["best_bits"] else None,
                 None, None) for r in __import__("csv").DictReader(open(a.table), delimiter="\t")]
        return plot(rows, ga, a.figure, a.model_name)
    ga = next((float(l.split()[1]) for l in open(a.cm) if l.startswith("GA ")), 0.0)
    rows = []
    for g in GENOMES:
        tp = list(read_fasta(f"{a.pipeline}/selenoprofiles/{g}.three_prime.fa"))
        bg = background(f"data/genomes/{g}/genome.fna.gz", a.background, 6000, seed=a.seed + GENOMES.index(g))
        recs = [(h, s, "selenocysteine" if label(h) == "selenocysteine" else "other prediction") for h, s in tp] + \
               [(h, s, "random window") for h, s in bg]
        if a.baseline:
            recs += [(h, s, "curated Sec gene") for h, s in curated(a.baseline, g)]
        for (h, s, grp), (sc, b, e) in zip(recs, best_scores(a.cm, [(h, s) for h, s, _ in recs], a.cpus)):
            rows.append((g, grp, h.split()[0], len(s), sc, b, e))
        print(g, len(tp), "3' sequences,", len(bg), "random windows", flush=True)
    with open(a.table, "w") as o:
        o.write("genome\tgroup\tsequence\tlength\tbest_bits\tstart\tend\n")
        for r in rows: o.write("\t".join("" if x is None else str(x) for x in r) + "\n")
    plot(rows, ga, a.figure, a.model_name)


def plot(rows, ga, prefix, model_name="RF00031"):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    groups = [g for g in ["curated Sec gene", "selenocysteine", "other prediction", "random window"]
              if any(r[1] == g for r in rows)]
    colours = {"curated Sec gene": "#e67e22", "selenocysteine": "#c0392b", "other prediction": "#7f8c8d",
               "random window": "#2c3e50"}
    ticks = {"curated Sec gene": "curated\nSec genes", "selenocysteine": "Sec\npredictions",
             "other prediction": "other\npredictions", "random window": "random\nwindows"}
    fig, axes = plt.subplots(3, 3, figsize=(10, 9), sharey=True)
    top = max(50, max((r[4] or 0) for r in rows) + 10)
    rng = random.Random(0)
    data = defaultdict(list)
    for g, grp, _, length, sc, _, _ in rows: data[(g, grp)].append(-3 + rng.uniform(-1.5, 1.5) if sc is None else sc)
    for ax, g in zip(axes.flat, GENOMES):
        for k, grp in enumerate(groups):
            ys = data[(g, grp)]
            ax.scatter([k + rng.uniform(-0.25, 0.25) for _ in ys], ys, s=6, alpha=0.6, color=colours[grp], linewidths=0)
            above = sum(1 for y in ys if y >= 20)
            if ys: ax.text(k, top - 4, f"{above}/{len(ys)}", ha="center", fontsize=7)
        if ga: ax.axhline(ga, color="black", lw=0.8, ls="--")
        ax.axhline(20, color="black", lw=0.5, ls=":")
        ax.axhspan(-5, -0.5, color="#ecf0f1", zorder=0)
        ax.set_title(NAMES[g], fontsize=9, style="italic"); ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([ticks[x] for x in groups], fontsize=7)
        ax.set_ylim(-5, top)
    for ax in axes[:, 0]: ax.set_ylabel(f"best {model_name} score (bits)", fontsize=8)
    thr = f"dashed: gathering threshold, {ga:g} bits; " if ga else ""
    fig.suptitle(f"{model_name} scores of the Selenoprofiles 3′ sequences and of random genomic windows\n"
                 f"({thr}dotted: 20 bits; grey band: no hit; numbers: sequences ≥ 20 bits)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for ext in ("png", "pdf"): fig.savefig(f"{prefix}.{ext}", dpi=200)


if __name__ == "__main__":
    main()
