#!/usr/bin/env python3
"""Build an appendicularian SECIS covariance model from the stage 1 outputs, by iterative bootstrap.

Training sequences: the 3' sequences that Selenoprofiles wrote after its predictions labelled selenocysteine, in the
training genomes. Background: random genomic windows of 6,000 nt from the same genomes. Held-out genomes are not used
for training and serve as controls.

Round 0 takes as seeds the best Rfam RF00031 hit of each training sequence when it scores at least --seed-bits.
Each round then aligns the seeds to the current model (cmalign), builds a new model (cmbuild), calibrates it
(cmcalibrate), searches the training sequences and the background, and takes as next seeds the best hit of each
training sequence that scores above the best background score of the round. Rounds stop when the seed set no longer
changes, or after --rounds rounds. The last model is calibrated fully.

Outputs in --out: round_<k>.sto and .cm, rounds.tsv (seeds and background maximum per round), final model
appendicularia_secis.cm and its seed alignment appendicularia_secis.sto.
"""
import argparse, gzip, os, random, shutil, subprocess

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def background(genome, n, length, seed):
    rng = random.Random(seed)
    contigs = [(name, seq) for name, seq in read_fasta(genome) if len(seq) >= length]
    total = sum(len(s) for _, s in contigs)
    for i in range(n):
        x = rng.randrange(total)
        for name, seq in contigs:
            if x < len(seq): break
            x -= len(seq)
        start = min(x, len(seq) - length)
        w = seq[start:start + length]
        yield f"bg|{name}|{start}", (w.translate(COMP)[::-1] if rng.random() < 0.5 else w)


def run(cmd, **kw):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, **kw)


def search(cm, fasta, out, cpus):
    run(["cmsearch", "--toponly", "-T", "0", "--cpu", str(cpus), "--tblout", out, cm, fasta])
    best = {}
    for l in open(out):
        if l.startswith("#"): continue
        f = l.split(); sc, a, b = float(f[14]), int(f[7]), int(f[8])
        if sc > best.get(f[0], (-1,))[0]: best[f[0]] = (sc, a, b)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline", default="work/stage1/pipeline")
    ap.add_argument("--genomes", default="data/genomes")
    ap.add_argument("--rfam", default="data/references/RF00031.cm")
    ap.add_argument("--train", nargs="+", default=["OKI", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery"])
    ap.add_argument("--seed-bits", type=float, default=15.0)
    ap.add_argument("--background", type=int, default=200, help="random windows per training genome")
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--quick-calibration", type=float, default=0.2, help="cmcalibrate -L for intermediate rounds (Mb)")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", default="work/stage2/secis_model")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    o = lambda f: os.path.join(a.out, f)

    # training sequences and background
    train = {}
    for g in a.train:
        for name, seq in read_fasta(f"{a.pipeline}/selenoprofiles/{g}.three_prime.fa"):
            if name.rsplit(".", 1)[1] == "selenocysteine":
                train[f"{g}|{name.split(':', 1)[1]}"] = seq
    with open(o("train.fa"), "w") as fh:
        for k, s in train.items(): fh.write(f">{k}\n{s}\n")
    if not os.path.exists(o("background.fa")):
        with open(o("background.fa"), "w") as fh:
            for i, g in enumerate(a.train):
                for k, s in background(f"{a.genomes}/{g}/genome.fna.gz", a.background, 6000, seed=100 + i):
                    fh.write(f">{g}|{k}\n{s}\n")

    # round 0: RF00031 hits
    best = search(a.rfam, o("train.fa"), o("round_0.train.tbl"), a.cpus)
    seeds = {k: v for k, v in best.items() if v[0] >= a.seed_bits}
    log = [("round", "seeds", "seeds_by_genome", "background_max_bits", "training_sequences_above")]
    log.append((0, len(seeds), ",".join(f"{g}:{sum(k.startswith(g + '|') for k in seeds)}" for g in a.train), "", ""))
    model, previous = a.rfam, None
    for r in range(1, a.rounds + 1):
        with open(o(f"round_{r}.seeds.fa"), "w") as fh:
            for k, (sc, s, e) in sorted(seeds.items()):
                fh.write(f">{k}/{s}-{e}\n{train[k][s - 1:e]}\n")
        run(["cmalign", "--noprob", "-o", o(f"round_{r}.sto"), model, o(f"round_{r}.seeds.fa")])
        run(["cmbuild", "-F", "-n", "appendicularia_SECIS", o(f"round_{r}.cm"), o(f"round_{r}.sto")])
        run(["cmcalibrate", "-L", str(a.quick_calibration), "--cpu", str(a.cpus), o(f"round_{r}.cm")])
        model = o(f"round_{r}.cm")
        bg = search(model, o("background.fa"), o(f"round_{r}.background.tbl"), a.cpus)
        bgmax = max([v[0] for v in bg.values()] or [0])
        best = search(model, o("train.fa"), o(f"round_{r}.train.tbl"), a.cpus)
        new = {k: v for k, v in best.items() if v[0] > bgmax}
        log.append((r, len(seeds), ",".join(f"{g}:{sum(k.startswith(g + '|') for k in seeds)}" for g in a.train),
                    f"{bgmax:.1f}", len(new)))
        print(*log[-1], flush=True)
        stable = set(new) == set(seeds)
        seeds = new
        if stable: break
    # final model from the last seed set, fully calibrated
    with open(o("final.seeds.fa"), "w") as fh:
        for k, (sc, s, e) in sorted(seeds.items()):
            fh.write(f">{k}/{s}-{e}\n{train[k][s - 1:e]}\n")
    run(["cmalign", "--noprob", "-o", o("appendicularia_secis.sto"), model, o("final.seeds.fa")])
    run(["cmbuild", "-F", "-n", "appendicularia_SECIS", o("appendicularia_secis.cm"), o("appendicularia_secis.sto")])
    run(["cmcalibrate", "--cpu", str(a.cpus), o("appendicularia_secis.cm")])
    with open(o("rounds.tsv"), "w") as fh:
        for row in log: fh.write("\t".join(map(str, row)) + "\n")
    print("final seeds:", len(seeds))


if __name__ == "__main__":
    main()
