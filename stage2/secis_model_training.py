#!/usr/bin/env python3
"""Step 3 of the plan (docs/plan.md, 2026-10-24): a SECIS covariance model trained on the six Oikopleura assemblies,
evaluated by leaving out one family or one species at a time, and compared with Rfam RF00031, RF00031 merged with
appendicularian sequences, and the current appendicularian model.

Positive windows: the --window nucleotides that follow the stop codon of each selenoprotein gene, on its strand,
shortened at the next selenoprotein coding sequence downstream on the same strand (so that the SELENOW pair of chr1
does not give one SECIS element twice), and left out when shorter than --min-window (contig ends):
  - OKI2018, OSKA2016 and OdB3: the gold-standard models of step 2 (curation/gold_standard/<genome>.gff3);
  - O. longicauda, O. albicans and O. vanhoeffeni: the stage 1 predictions labelled selenocysteine and kept by the
    first round of hit validation (results/tables/stage2_secis_hits.tsv), whose frame is extended on the genome to
    the next stop codon (at most --max-extension nucleotides; otherwise the window starts at the end of the
    prediction, flagged).
Negative windows: random genomic windows (--random per genome, uniform over the contigs at least --window long, on a
random strand; a window that overlaps a positive window is drawn again, so that the SECIS elements of known
selenoproteins do not count as false positives) of the six assemblies and of the non-target tunicates (--nontarget),
drawn with one random seed for training and another for evaluation; each positive window shuffled with its dinucleotide composition preserved
(Altschul and Erickson 1985); and the windows after the stop codons of the stage 1 cysteine homologs kept by the
first round, except those overlapping a positive window.

Training (each fold, and the final model on all positive windows): bootstrap from RF00031. Round 0 seeds are the
best RF00031 hit of each training window with at least --seed-bits. Each round aligns the seeds to the current model
(cmalign), builds a model (cmbuild) and calibrates it quickly (cmcalibrate -L --quick-calibration), searches the
training windows and the training background (the training random windows of the training genomes), and takes as
next seeds the best hit of each training window scoring above the best background hit; then, in decreasing order of
score, a seed is dropped when it overlaps an earlier seed in the genome (one SECIS element per locus) or is at least
--identity identical to an earlier seed (identical bases over the length of the shorter sequence, in a global
alignment with free end gaps). Rounds stop when the seed set no longer changes, or after --rounds rounds. The merged
model of a fold is built from the sequences of the RF00031 seed alignment and the final seeds of the fold, aligned
together to RF00031 (cmalign), so that both models of a fold see the same appendicularian sequences.

Evaluation: one fold per family (its positive windows in the six assemblies held out) and one per species
(O. dioica, its three assemblies held out together because they share near-identical SECIS elements; O. longicauda;
O. albicans; O. vanhoeffeni). In each fold, all windows are scored with the fold's new and merged models, RF00031
and the current model (results/models/appendicularia_secis.cm, trained on sequences of most held-out sets: its
results are optimistic). For each model, the training threshold is the best score of the training random windows of
the fold's training genomes; the held-out positive windows are counted above it and above the best score of the
evaluation negatives of the fold (evaluation random windows of the held-out genomes, or of all six genomes for a
family fold, and of the non-target tunicates; shuffled held-out windows; held-out cysteine windows), and the
evaluation negatives above the training threshold are counted per megabase.

Outputs: in --out, the windows (windows.fa, windows.tsv), one directory per fold (rounds, seeds, models, scores) and
final/; in results/tables, stage2_secis_v2_windows.tsv, stage2_secis_v2_rounds.tsv, stage2_secis_v2_summary.tsv
(fold, model and pooled rows), stage2_secis_v2_heldout.tsv (one row per positive window: its held-out scores) and
stage2_secis_v2_negatives.tsv (the evaluation negatives above the training threshold, per fold and model, with their
coordinates: candidate SECIS elements of selenoproteins not yet known);
the final model is copied to results/models/oikopleura_secis.{cm,sto} with its seeds.
"""
import argparse, collections, concurrent.futures, csv, gzip, os, random, shutil, subprocess

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = {"TAA", "TAG", "TGA"}
SPECIES = {"OKI": "O. dioica", "OSKA": "O. dioica", "OdB3": "O. dioica", "Olon": "O. longicauda",
           "Oalb": "O. albicans", "Ovan": "O. vanhoeffeni"}
species = lambda g: SPECIES.get(g, g)


def rc(s): return s.translate(COMP)[::-1]


def read_fasta(path):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def attributes(field):
    out = {}
    for x in field.strip().strip(";").split(";"):
        x = x.strip()
        if not x: continue
        if "=" in x: k, v = x.split("=", 1)
        else: k, _, v = x.partition(" ")
        out[k] = v.strip('"')
    return out


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# ---------------------------------------------------------------------------------------------------------------
# windows

def window(seq, contig, strand, stop_end, length, limit=None):
    """Window of up to length nt after a stop codon ending at stop_end (0-based exclusive on +, 0-based start of the
    stop codon on -), shortened at limit (the next coding sequence downstream: its first base on +, its last base
    + 1 on -). Returns (start, end, sequence on the gene strand), genomic 0-based half-open."""
    if strand == "+":
        a, b = stop_end, min(len(seq), stop_end + length, limit if limit is not None else len(seq))
        return a, b, seq[a:b].upper()
    a, b = max(0, stop_end - length, limit if limit is not None else 0), stop_end
    return a, b, rc(seq[a:b].upper())


def gold_standard_cds(path):
    """Gold-standard GFF3: gene -> (contig, strand, [(start0, end0 exclusive)] of its CDS rows)."""
    out = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        gene = attributes(f[8])["Parent"].split("_", 1)[1]
        out.setdefault(gene, (f[0], f[6], []))[2].append((int(f[3]) - 1, int(f[4])))
    return out


def read_gtf(path):
    """Selenoprofiles GTF: prediction -> (contig, strand, [(start0, end0 exclusive, phase)])."""
    out = {}
    if not os.path.exists(path): return out
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        at = attributes(f[8])
        pid = at.get("transcript_id") or at.get("gene_id")
        out.setdefault(pid, (f[0], f[6], []))[2].append((int(f[3]) - 1, int(f[4]), 0 if f[7] in (".", "") else int(f[7])))
    return out


def extend_to_stop(seq, strand, cds, max_extension, need_tga):
    """End of the stop codon of a predicted coding sequence (CDS rows, possibly partial): the frame is the phase of
    the first row, or the frame with fewest stop codons that holds a TGA when need_tga; it is extended on the genome
    beyond the last row to the next stop codon. Returns (stop_end as in window(), found)."""
    rows = sorted(cds, key=lambda x: x[0], reverse=(strand == "-"))
    pos = [p for a, b, _ in rows for p in (range(a, b) if strand == "+" else range(b - 1, a - 1, -1))]
    nt = "".join(seq[p] for p in pos).upper()
    if strand == "-": nt = nt.translate(COMP)
    codons = lambda f: [nt[i:i + 3] for i in range(f, len(nt) - 2, 3)]
    frame = rows[0][2] % 3
    if need_tga:
        ok = [f for f in (0, 1, 2) if "TGA" in codons(f)]
        if ok and frame not in ok:
            frame = min(ok, key=lambda f: sum(c in STOPS for c in codons(f)))
    last = codons(frame)
    if last and last[-1] in STOPS:   # the prediction ends with its stop codon
        end_pos = pos[len(nt) - (len(nt) - frame) % 3 - 1]
        return (end_pos + 1 if strand == "+" else end_pos), True
    tail = (len(nt) - frame) % 3   # bases of an incomplete codon at the end of the prediction
    edge = pos[-1]
    for k in range(-tail, max_extension, 3):
        if strand == "+":
            a = edge + 1 + k; c = seq[a:a + 3].upper()
            if len(c) < 3: break
            if c in STOPS: return a + 3, True
        else:
            a = edge - 1 - k; c = rc(seq[a - 2:a + 1].upper()) if a - 2 >= 0 else ""
            if len(c) < 3: break
            if c in STOPS: return a - 2, True
    return (edge + 1 if strand == "+" else edge), False


def random_windows(contigs, n, length, seed, exclude=None):
    """n random windows (uniform over the contigs at least length long, random strand); windows that overlap a span
    of exclude (contig -> [(start, end)], the positive windows) are drawn again."""
    rng = random.Random(seed)
    usable = [(name, s) for name, s in contigs.items() if len(s) >= length]
    total = sum(len(s) for _, s in usable)
    i = tries = 0
    while total and i < n and tries < 100 * n:
        tries += 1
        x = rng.randrange(total)
        for name, s in usable:
            if x < len(s): break
            x -= len(s)
        a = min(x, len(s) - length)
        strand = "-" if rng.random() < 0.5 else "+"
        if exclude and any(p < a + length and a < q for p, q in exclude.get(name, [])): continue
        w = s[a:a + length].upper()
        yield i, name, strand, a, a + length, (rc(w) if strand == "-" else w)
        i += 1


def dinucleotide_shuffle(seq, rng):
    """Random sequence with the dinucleotide counts and the first and last letters of seq (Altschul and Erickson
    1985: Eulerian walk with a random spanning arborescence of last exits)."""
    s = seq.upper()
    if len(s) < 3: return s
    edges = collections.defaultdict(list)
    for x, y in zip(s, s[1:]): edges[x].append(y)
    last = s[-1]
    while True:
        exit_ = {v: rng.choice(e) for v, e in edges.items() if v != last}
        ok = True
        for v in edges:
            seen, u = set(), v
            while u != last:
                if u in seen or u not in exit_: ok = False; break
                seen.add(u); u = exit_[u]
            if not ok: break
        if ok: break
    lists = {}
    for v, e in edges.items():
        e = list(e)
        if v in exit_:
            e.remove(exit_[v]); rng.shuffle(e); e.append(exit_[v])
        else: rng.shuffle(e)
        lists[v] = e
    out, u, idx = [s[0]], s[0], collections.Counter()
    for _ in range(len(s) - 1):
        v = lists[u][idx[u]]; idx[u] += 1; out.append(v); u = v
    return "".join(out)


def build_windows(a):
    """All windows, as a list of dicts (id, genome, set, family, species, contig, strand, start, end, seq, note)."""
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    fam = {(r["genome"], r["gene"]): r["family"] for r in csv.DictReader(open(a.baseline), delimiter="\t")}
    kept = {}
    for r in csv.DictReader(open(a.hits), delimiter="\t"):
        if r["decision"] == "kept": kept[(r["genome"], r["prediction"])] = r["label"]
    rng = random.Random(a.shuffle_seed)
    wins = []
    genomes = a.genomes
    for i, g in enumerate(genomes + a.nontarget):
        contigs = {k: v for k, v in read_fasta(sheet[g]["fasta"])}
        print(g, len(contigs), "contigs", flush=True)
        pos_spans = collections.defaultdict(list)
        if g in genomes:
            # selenoprotein genes: (name, family, contig, strand, stop_end, all CDS rows, note)
            genes = []
            if g in a.gold_genomes:
                for gene, (c, st, cds) in gold_standard_cds(f"{a.gold}/{g}.gff3").items():
                    stop_end = max(b for _, b in cds) if st == "+" else min(x for x, _ in cds)
                    genes.append((gene, fam.get((g, gene), gene.split("_")[0]), c, st, stop_end, cds, "gold standard"))
            gtf = read_gtf(f"{a.stage1}/{g}.selenoprofiles.gtf")
            cys = []
            for pid, (c, st, cds) in gtf.items():
                label = kept.get((g, pid))
                if label is None or c not in contigs: continue
                stop_end, found = extend_to_stop(contigs[c], st, cds, a.max_extension, label == "selenocysteine")
                rows = [(x, y) for x, y, _ in cds]
                note = "stage 1 prediction" + ("" if found else ", no stop codon found: window after the prediction")
                if label == "selenocysteine" and g not in a.gold_genomes:
                    genes.append((pid, pid.split(".")[0], c, st, stop_end, rows, note))
                elif label == "cysteine":
                    cys.append((pid, pid.split(".")[0], c, st, stop_end, rows, note))
            # next selenoprotein coding sequence downstream on the same strand
            def limit(c, st, stop_end, own):
                spans = [(min(x for x, _ in rows), max(y for _, y in rows)) for _, _, c2, st2, _, rows, _ in genes
                         if c2 == c and st2 == st and rows is not own]
                if st == "+":
                    nxt = [x for x, _ in spans if x >= stop_end]
                    return min(nxt) if nxt else None
                nxt = [y for _, y in spans if y <= stop_end]
                return max(nxt) if nxt else None
            for name, family, c, st, stop_end, rows, note in genes:
                s, e, w = window(contigs[c], c, st, stop_end, a.window, limit(c, st, stop_end, rows))
                if len(w) < a.min_window:
                    print(f"  {g} {name}: window of {len(w)} nt, left out", flush=True); continue
                wins.append({"id": f"{g}|{name}", "genome": g, "set": "positive", "family": family,
                             "species": species(g), "contig": c, "strand": st, "start": s, "end": e, "seq": w,
                             "note": note})
                pos_spans[c].append((s, e))
                wins.append({"id": f"{g}|{name}|shuffled", "genome": g, "set": "shuffled", "family": family,
                             "species": species(g), "contig": c, "strand": st, "start": s, "end": e,
                             "seq": dinucleotide_shuffle(w, rng), "note": "dinucleotide shuffle"})
            for name, family, c, st, stop_end, rows, note in cys:
                s, e, w = window(contigs[c], c, st, stop_end, a.window)
                if len(w) < a.min_window or any(x < e and s < y for x, y in pos_spans[c]): continue
                wins.append({"id": f"{g}|{name}", "genome": g, "set": "cysteine", "family": family,
                             "species": species(g), "contig": c, "strand": st, "start": s, "end": e, "seq": w,
                             "note": note})
        for kind, seed in (("random_train", a.train_seed + i), ("random_eval", a.eval_seed + i)):
            if kind == "random_train" and g not in genomes: continue
            for k, c, st, s, e, w in random_windows(contigs, a.random, a.window, seed, pos_spans if g in genomes else None):
                wins.append({"id": f"{g}|{kind}|{k}", "genome": g, "set": kind, "family": "-",
                             "species": species(g) if g in genomes else "non-target", "contig": c, "strand": st, "start": s, "end": e,
                             "seq": w, "note": ""})
        del contigs
    return wins


# ---------------------------------------------------------------------------------------------------------------
# models

def search(cm, fasta, tbl, cpus):
    """Best hit of each sequence: id -> (bits, from, to), positions on the sequence (top strand)."""
    run(["cmsearch", "--toponly", "-T", "0", "--cpu", str(cpus), "--tblout", tbl, cm, fasta])
    best = {}
    for l in open(tbl):
        if l.startswith("#"): continue
        f = l.split(); sc, x, y = float(f[14]), int(f[7]), int(f[8])
        if sc > best.get(f[0], (-1,))[0]: best[f[0]] = (sc, x, y)
    return best


def genomic(w, x, y):
    """Genomic 0-based half-open interval of positions x..y (1-based, inclusive) of a window."""
    return (w["start"] + x - 1, w["start"] + y) if w["strand"] == "+" else (w["end"] - y, w["end"] - x + 1)


def identity(s, t, aligner):
    """Identical bases over the length of the shorter sequence, global alignment with free end gaps."""
    bs, bt = aligner.align(s, t)[0].aligned
    same = sum(s[i] == t[j] for (a0, a1), (b0, b1) in zip(bs, bt) for i, j in zip(range(a0, a1), range(b0, b1)))
    return same / min(len(s), len(t))


def select_seeds(hits, wins, threshold, max_identity, aligner):
    """Best hit of each window above threshold, then one per locus and below max_identity to higher-scoring seeds."""
    chosen = []
    for wid, (sc, x, y) in sorted(hits.items(), key=lambda kv: (-kv[1][0], kv[0])):
        if sc <= threshold: continue
        w = wins[wid]
        iv = genomic(w, x, y)
        seq = w["seq"][x - 1:y]
        if any(wins[o]["genome"] == w["genome"] and wins[o]["contig"] == w["contig"] and ov[0] < iv[1] and iv[0] < ov[1]
               for o, ov, _ in chosen): continue
        if any(identity(seq, os_, aligner) >= max_identity for _, _, os_ in chosen): continue
        chosen.append((wid, iv, seq))
    return {wid: hits[wid] for wid, _, _ in chosen}


def write_seeds(path, seeds, wins):
    with open(path, "w") as fh:
        for wid, (sc, x, y) in sorted(seeds.items()):
            fh.write(f">{wid}/{x}-{y}\n{wins[wid]['seq'][x - 1:y]}\n")


def train(name, train_ids, background_ids, wins, a, cpus):
    """Bootstrap from RF00031 on the training windows; returns the paths of the new and merged models, the rounds
    and the seeds."""
    from Bio import Align
    aligner = Align.PairwiseAligner(mode="global", match_score=1, mismatch_score=-1, open_gap_score=-2,
                                    extend_gap_score=-1)
    aligner.end_insertion_score = aligner.end_deletion_score = 0
    d = os.path.join(a.out, name); os.makedirs(d, exist_ok=True)
    o = lambda f: os.path.join(d, f)
    with open(o("train.fa"), "w") as fh:
        for k in train_ids: fh.write(f">{k}\n{wins[k]['seq']}\n")
    with open(o("background.fa"), "w") as fh:
        for k in background_ids: fh.write(f">{k}\n{wins[k]['seq']}\n")
    best = search(a.rfam, o("train.fa"), o("round_0.train.tbl"), cpus)
    seeds = select_seeds(best, wins, a.seed_bits - 1e-9, a.identity, aligner)
    rounds = [(name, 0, len(seeds), "", "")]
    model = a.rfam
    for r in range(1, a.rounds + 1):
        write_seeds(o(f"round_{r}.seeds.fa"), seeds, wins)
        run(["cmalign", "--noprob", "-o", o(f"round_{r}.sto"), model, o(f"round_{r}.seeds.fa")])
        run(["cmbuild", "-F", "-n", "oikopleura_SECIS", o(f"round_{r}.cm"), o(f"round_{r}.sto")])
        run(["cmcalibrate", "-L", str(a.quick_calibration), "--cpu", str(cpus), o(f"round_{r}.cm")])
        model = o(f"round_{r}.cm")
        bg = search(model, o("background.fa"), o(f"round_{r}.background.tbl"), cpus)
        bgmax = max([v[0] for v in bg.values()] or [0])
        hits = search(model, o("train.fa"), o(f"round_{r}.train.tbl"), cpus)
        new = select_seeds(hits, wins, bgmax, a.identity, aligner)
        rounds.append((name, r, len(new), f"{bgmax:.1f}", sum(v[0] > bgmax for v in hits.values())))
        stable = set(new) == set(seeds)
        seeds = new
        if stable: break
    write_seeds(o("final.seeds.fa"), seeds, wins)
    run(["cmalign", "--noprob", "-o", o("new.sto"), model, o("final.seeds.fa")])
    run(["cmbuild", "-F", "-n", "oikopleura_SECIS", o("new.cm"), o("new.sto")])
    with open(o("merged.fa"), "w") as fh:
        for k, v in stockholm_sequences(a.rfam_seed).items(): fh.write(f">rfam|{k}\n{v}\n")
        fh.write(open(o("final.seeds.fa")).read())
    run(["cmalign", "--noprob", "-o", o("merged.sto"), a.rfam, o("merged.fa")])
    run(["cmbuild", "-F", "-n", "RF00031_oikopleura", o("merged.cm"), o("merged.sto")])
    mb = (a.final_calibration if name == "final" else a.quick_calibration)
    cal = ["-L", str(mb)] if mb else []
    for m in ("new", "merged"):
        run(["cmcalibrate"] + cal + ["--cpu", str(cpus), o(f"{m}.cm")])
    return {"new": o("new.cm"), "merged": o("merged.cm")}, rounds, seeds


def fold_job(args):
    """Train one fold and score all windows with its models (run in a worker process)."""
    name, train_ids, background_ids, a, cpus, wins = args
    models, rounds, seeds = train(name, train_ids, background_ids, wins, a, cpus)
    scores = {m: search(p, os.path.join(a.out, "windows.fa"), os.path.join(a.out, name, f"{m}.windows.tbl"), cpus)
              for m, p in models.items()}
    return name, rounds, {wid: (sc, x, y) for wid, (sc, x, y) in seeds.items()}, \
        {m: {k: v[0] for k, v in s.items()} for m, s in scores.items()}


# ---------------------------------------------------------------------------------------------------------------

def evaluation_negatives(wlist, kind, label, held, eval_genomes, nontarget):
    """Evaluation negatives of a fold: evaluation random windows of its evaluation genomes and of the non-target
    genomes, shuffled held-out positive windows, and cysteine windows of the held-out family or species."""
    out = []
    for w in wlist:
        if w["set"] == "random_eval" and (w["genome"] in eval_genomes or w["genome"] in nontarget): out.append(w["id"])
        elif w["set"] == "shuffled" and w["id"][:-len("|shuffled")] in held: out.append(w["id"])
        elif w["set"] == "cysteine" and w[kind] == label: out.append(w["id"])
    return out


def stockholm_sequences(path):
    """Ungapped sequences of a Stockholm alignment (U written as T)."""
    seqs = collections.OrderedDict()
    for l in open(path):
        if not l.strip() or l.startswith("#") or l.startswith("//"): continue
        name, part = l.split()[:2]
        seqs[name] = seqs.get(name, "") + part
    return {k: "".join(c for c in v.upper() if c.isalpha()).replace("U", "T") for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--baseline", default="curation/baseline_sec_genes.tsv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--stage1", default="work/stage1/pipeline/selenoprofiles")
    ap.add_argument("--hits", default="results/tables/stage2_secis_hits.tsv")
    ap.add_argument("--rfam", default="data/references/RF00031.cm")
    ap.add_argument("--rfam-seed", default="data/references/RF00031.sto")
    ap.add_argument("--current", default="results/models/appendicularia_secis.cm")
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan"])
    ap.add_argument("--gold-genomes", nargs="*", default=["OKI", "OSKA", "OdB3"],
                    help="genomes whose positive windows come from the gold standard (the others: stage 1 predictions)")
    ap.add_argument("--nontarget", nargs="*", default=["Cint", "Csav", "Scla", "Tdem"])
    ap.add_argument("--window", type=int, default=2000)
    ap.add_argument("--min-window", type=int, default=100, help="shortest positive window (a SECIS element is 60-80 nt)")
    ap.add_argument("--max-extension", type=int, default=600)
    ap.add_argument("--random", type=int, default=200, help="random windows per genome, for training and for evaluation")
    ap.add_argument("--train-seed", type=int, default=300)
    ap.add_argument("--eval-seed", type=int, default=400)
    ap.add_argument("--shuffle-seed", type=int, default=500)
    ap.add_argument("--seed-bits", type=float, default=15.0)
    ap.add_argument("--identity", type=float, default=0.9)
    ap.add_argument("--rounds", type=int, default=8)
    ap.add_argument("--quick-calibration", type=float, default=0.2, help="cmcalibrate -L for the fold models (Mb)")
    ap.add_argument("--final-calibration", type=float, default=0, help="cmcalibrate -L for the final models (0: full)")
    ap.add_argument("--folds", nargs="*", help="only these folds (names as in the summary), for tests")
    ap.add_argument("--no-final", action="store_true", help="skip the final model, for tests")
    ap.add_argument("--jobs", type=int, default=4, help="folds trained in parallel")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", default="work/stage2/secis_v2")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--models", default="results/models")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); os.makedirs(a.tables, exist_ok=True)
    o = lambda f: os.path.join(a.out, f)

    wlist = build_windows(a)
    wins = {w["id"]: w for w in wlist}
    with open(o("windows.fa"), "w") as fh:
        for w in wlist: fh.write(f">{w['id']}\n{w['seq']}\n")
    cols = ["id", "genome", "set", "family", "species", "contig", "strand", "start", "end", "length", "note"]
    with open(o("windows.tsv"), "w") as fh, open(f"{a.tables}/stage2_secis_v2_windows.tsv", "w") as pub:
        fh.write("\t".join(cols) + "\n"); pub.write("\t".join(cols) + "\n")
        for w in wlist:
            row = "\t".join(str(len(w["seq"])) if c == "length" else (str(w[c] + 1) if c == "start" else str(w[c]))
                            for c in cols) + "\n"
            fh.write(row)
            if w["set"] in ("positive", "cysteine"): pub.write(row)
    pos = [w for w in wlist if w["set"] == "positive"]
    print(len(pos), "positive windows;", collections.Counter(w["genome"] for w in pos), flush=True)

    # folds: training positives, training background, held-out positives, evaluation negatives
    genomes = a.genomes
    folds = []
    for fam in sorted({w["family"] for w in pos}):
        held = [w["id"] for w in pos if w["family"] == fam]
        folds.append((f"family:{fam}", [w["id"] for w in pos if w["family"] != fam], genomes, held, genomes))
    for sp in sorted({species(g) for g in genomes}):
        gs = [g for g in genomes if species(g) == sp]
        folds.append((f"species:{sp}", [w["id"] for w in pos if w["species"] != sp], [g for g in genomes if g not in gs],
                      [w["id"] for w in pos if w["species"] == sp], gs))
    if a.folds: folds = [f for f in folds if f[0] in a.folds]
    jobs = []
    for name, tr, train_genomes, held, eval_genomes in folds:
        bg = [w["id"] for w in wlist if w["set"] == "random_train" and w["genome"] in train_genomes]
        jobs.append((name.replace(":", "_").replace(" ", "_").replace(".", ""), tr, bg, a, max(1, a.cpus // a.jobs), wins))
    if not a.no_final:
        bg = [w["id"] for w in wlist if w["set"] == "random_train"]
        jobs.append(("final", [w["id"] for w in pos], bg, a, max(1, a.cpus // a.jobs), wins))

    # fixed models, scored once
    fixed = {m: {k: v[0] for k, v in search(p, o("windows.fa"), o(f"{m}.windows.tbl"), a.cpus).items()}
             for m, p in (("RF00031", a.rfam), ("current", a.current))}
    results = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for name, rounds, seeds, scores in ex.map(fold_job, jobs):
            results[name] = (rounds, seeds, scores)
            print(name, "done:", len(seeds), "seeds", flush=True)

    # rounds
    with open(f"{a.tables}/stage2_secis_v2_rounds.tsv", "w") as fh:
        fh.write("fold\tround\tseeds\tbackground_max_bits\ttraining_windows_above\n")
        for name in results:
            for r in results[name][0]: fh.write("\t".join(map(str, r)) + "\n")
    # summary per fold and model, and pooled
    summary = [("fold", "model", "held-out positive windows", "training threshold (bits)",
                "detected above the training threshold", "best evaluation negative (bits)",
                "detected above the best evaluation negative", "evaluation negatives above the training threshold",
                "per Mb", "negative windows", "negative Mb")]
    pooled = collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0.0])
    heldout = collections.defaultdict(dict)
    negrows = [("fold", "model", "window", "set", "genome", "contig", "strand", "start", "end", "bits", "training threshold (bits)")]
    for (fname, tr, train_genomes, held, eval_genomes), job in zip(folds, jobs):
        key = job[0]
        if key not in results: continue
        kind, label = fname.split(":", 1)
        neg = evaluation_negatives(wlist, kind, label, set(held), eval_genomes, a.nontarget)
        train_bg = job[2]
        models = dict(results[key][2]); models.update(fixed)
        for m, sc in models.items():
            thr = max([sc.get(k, 0.0) for k in train_bg] or [0.0])
            negmax = max([sc.get(k, 0.0) for k in neg] or [0.0])
            det = sum(sc.get(k, 0.0) > thr for k in held)
            det_strict = sum(sc.get(k, 0.0) > negmax for k in held)
            fp = sum(sc.get(k, 0.0) > thr for k in neg)
            for k in neg:
                if sc.get(k, 0.0) > thr:
                    w = wins[k]
                    negrows.append((fname, m, k, w["set"], w["genome"], w["contig"], w["strand"], w["start"] + 1, w["end"],
                                    f"{sc[k]:.1f}", f"{thr:.1f}"))
            mb = sum(len(wins[k]["seq"]) for k in neg) / 1e6
            summary.append((fname, m, len(held), f"{thr:.1f}", det, f"{negmax:.1f}", det_strict, fp,
                            f"{fp / mb:.2f}" if mb else "", len(neg), f"{mb:.2f}"))
            p = pooled[(kind, m)]
            p[0] += len(held); p[1] += det; p[2] += det_strict; p[3] += fp; p[4] += len(neg); p[5] += mb
            for k in held: heldout[k][f"{kind} {m}"] = f"{sc.get(k, 0.0):.1f}"
    for (kind, m), p in sorted(pooled.items()):
        summary.append((f"all {kind} folds", m, p[0], "", p[1], "", p[2], p[3], f"{p[3] / p[5]:.2f}" if p[5] else "",
                        p[4], f"{p[5]:.2f}"))
    with open(f"{a.tables}/stage2_secis_v2_summary.tsv", "w") as fh:
        for r in summary: fh.write("\t".join(map(str, r)) + "\n")
    with open(f"{a.tables}/stage2_secis_v2_negatives.tsv", "w") as fh:
        for r in negrows: fh.write("\t".join(map(str, r)) + "\n")
    hcols = sorted({c for v in heldout.values() for c in v})
    with open(f"{a.tables}/stage2_secis_v2_heldout.tsv", "w") as fh:
        fh.write("\t".join(["window", "genome", "family"] + hcols) + "\n")
        for k in sorted(heldout):
            fh.write("\t".join([k, wins[k]["genome"], wins[k]["family"]] + [heldout[k].get(c, "") for c in hcols]) + "\n")
    print(open(f"{a.tables}/stage2_secis_v2_summary.tsv").read())
    # final model
    if "final" in results:
        os.makedirs(a.models, exist_ok=True)
        for src, dst in (("new.cm", "oikopleura_secis.cm"), ("new.sto", "oikopleura_secis.sto"),
                         ("final.seeds.fa", "oikopleura_secis.seeds.fa")):
            shutil.copy(os.path.join(a.out, "final", src), os.path.join(a.models, dst))
        print("final model:", len(results["final"][1]), "seeds; copied to", a.models)


if __name__ == "__main__":
    main()
