#!/usr/bin/env python3
"""Description of the Oikopleura SECIS model of step 3 (results/models/oikopleura_secis.cm, seeds in .seeds.fa).

The seeds are realigned to the model with cmalign, and so are two other sets, with the same command:
- decoys: in dinucleotide-shuffled copies of the windows after the stop codons of the gold-standard models
  (curation/gold_standard/<genome>.after_stop.fa, shortened as in training; --shuffles copies of each), the best hit
  of the model on the top strand of each copy (cmsearch --max); the --decoys best of these (by default as many as the
  seeds). They show how much pairing a structural alignment gives to sequences that are not SECIS elements.
- gold-standard hits: the best hit of the model in each of these windows (cmsearch --max); among them the genes that
  joined the gold standard after training (--extra: TXNL1.b, KKSU), and the orthologous hits of the three O. dioica
  assemblies, compared by sequence (see ortholog_events).
The cmsearch tables are kept in --work and reused when the model and the sequences are unchanged.

Outputs (results/tables/, prefix stage2_secis_v2_model):
- _seeds.tsv: per seed and per hit of the new genes: genome, gene or prediction, family, position of the hit after
  the stop codon, length, score (new genes, cmsearch --max), identity to the nearest seed (consensus columns where
  both have a residue) and that seed, cluster (single linkage at --cluster-identity), core motifs (UGA on the 5'
  strand and GA on the 3' strand of the internal loop between helices 1 and 2; AAR in the loop between helices 2 and
  3, or in the apical loop), canonical pairs per helix.
- _identities.tsv: identity between each two seeds, and whether they are orthologous (same gene, two O. dioica
  assemblies).
- _positions.tsv: per consensus position: consensus residue, structure, residue counts in the seeds, fraction of the
  most frequent residue, information content of the seeds (2 minus the Shannon entropy, bits), and relative entropy
  of the match emissions of the model (model_weights).
- _pairs.tsv: per base pair of the consensus structure, in the seeds and in the decoys: canonical pairs (Watson-Crick
  and G-U), pair types, mutual information (bits) between the two columns and its p-value by permutation
  (--permutations shuffles of one column among the cluster representatives, so that the near-identical seeds of
  orthologous genes count once), Benjamini-Hochberg q-value over the pairs of the set, and the cluster
  representatives whose pair differs from the consensus pair at both positions and is canonical.
- _orthologs.tsv: per gene and base pair, the pairs of the three O. dioica assemblies and the change between them.
- _summary.tsv: identities, clusters, motifs and pairs per helix (seeds and decoys), and the changes between orthologs
  with their expected fractions.
"""
import argparse, collections, csv, hashlib, itertools, math, os, random, re, statistics, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uga_secis_screen import dinucleotide_shuffle   # noqa: E402

CANON = {"AU", "UA", "GC", "CG", "GU", "UG"}
HELICES = {1: "helix 1", 2: "helix 2", 3: "helix 3"}


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        l = l.rstrip()
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].split()[0], []
        elif l: buf.append(l)
    if name is not None: yield name, "".join(buf)


def read_stockholm(text):
    seqs, gc = {}, {}
    for l in text.splitlines():
        if not l.strip() or l.startswith("//") or l.startswith("# STOCKHOLM") or l.startswith("#=GF") or l.startswith("#=GS"):
            continue
        if l.startswith("#=GC"):
            _, tag, s = l.split(None, 2); gc[tag] = gc.get(tag, "") + s.strip()
        elif l.startswith("#=GR") or l.startswith("#"): continue
        else:
            n, s = l.split(); seqs[n] = seqs.get(n, "") + s
    return seqs, gc


def run(cmd, **kw):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kw)
    if p.returncode: sys.exit(f"command failed: {' '.join(cmd)}\n{p.stderr[-2000:]}")
    return p.stdout


def cmalign(cm, seqs, tmp, tag, maps=None):
    """Aligned residues at the consensus positions (gap '-'), one string per sequence, and the consensus line and
    structure at those positions; with maps (a dict), also the index in the sequence of the residue at each consensus
    position (None for a gap)."""
    fa = os.path.join(tmp, f"{tag}.fa")
    with open(fa, "w") as o:
        for k, s in seqs.items(): o.write(f">{k}\n{s.upper().replace('T', 'U')}\n")
    sto, gc = read_stockholm(run(["cmalign", "--noprob", "--outformat", "Pfam", cm, fa]))
    rf, ss = gc["RF"], gc["SS_cons"]
    cols = [i for i, c in enumerate(rf) if c not in ".~"]
    aln = {k: "".join(v[i] if v[i] not in "." else "-" for i in cols).upper() for k, v in sto.items()}
    if maps is not None:
        cset = set(cols)
        for k, v in sto.items():
            idx, m = 0, []
            for i, c in enumerate(v):
                res = c.isalpha()
                if i in cset: m.append(idx if res else None)
                if res: idx += 1
            maps[k] = m
    return aln, "".join(rf[i] for i in cols), "".join(ss[i] for i in cols)


def pairs_of(ss):
    st, out = [], []
    for i, c in enumerate(ss):
        if c in "<([{": st.append(i)
        elif c in ">)]}": out.append((st.pop(), i))
    return sorted(out)


def helices(pairs):
    """Helix number of each pair: consecutive stacked pairs (i+1, j-1) form a helix, numbered 5' to 3'."""
    h, k, prev = {}, 0, None
    for i, j in pairs:
        if prev is None or not (i == prev[0] + 1 and j == prev[1] - 1): k += 1
        h[(i, j)] = k; prev = (i, j)
    return h


def identity(a, b):
    both = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
    return sum(x == y for x, y in both) / len(both) if both else 0.0


def clusters(names, pid, t):
    parent = {n: n for n in names}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for (a, b), x in pid.items():
        if x >= t: parent[find(a)] = find(b)
    groups = collections.defaultdict(list)
    for n in names: groups[find(n)].append(n)
    return list(groups.values())


def mutual_information(xs, ys):
    n = len(xs)
    if not n: return 0.0
    cx, cy, cxy = collections.Counter(xs), collections.Counter(ys), collections.Counter(zip(xs, ys))
    return sum(c / n * math.log2(c * n / (cx[x] * cy[y])) for (x, y), c in cxy.items())


def bh(ps):
    """Benjamini-Hochberg q-values, in the order of ps."""
    m = len(ps)
    order = sorted(range(m), key=lambda k: ps[k])
    q, best = [0.0] * m, 1.0
    for r in range(m - 1, -1, -1):
        k = order[r]; best = min(best, ps[k] * m / (r + 1)); q[k] = best
    return q


def pair_stats(aln, reps, cons, pairs, hel, perms, rng):
    out = []
    for (i, j) in pairs:
        both = [(aln[n][i], aln[n][j]) for n in aln if aln[n][i] != "-" and aln[n][j] != "-"]
        rb = [(aln[n][i], aln[n][j]) for n in reps if aln[n][i] != "-" and aln[n][j] != "-"]
        canon = sum(x + y in CANON for x, y in both)
        types = collections.Counter(x + y for x, y in rb if x + y in CANON)
        xs, ys = [x for x, _ in rb], [y for _, y in rb]
        mi = mutual_information(xs, ys)
        ge = sum(mutual_information(xs, rng.sample(ys, len(ys))) >= mi - 1e-12 for _ in range(perms))
        ci, cj = cons[i].upper(), cons[j].upper()
        comp = sum(x != ci and y != cj and x + y in CANON for x, y in rb)
        out.append({"helix": hel[(i, j)], "i": i + 1, "j": j + 1, "consensus": f"{cons[i]}-{cons[j]}", "sequences": len(both),
                    "canonical": canon, "canonical_fraction": canon / len(both) if both else 0.0,
                    "representatives": len(rb), "pair_types": ",".join(f"{k}:{v}" for k, v in types.most_common()),
                    "mi_bits": mi, "p": (ge + 1) / (perms + 1), "compensatory": comp})
    for r, q in zip(out, bh([r["p"] for r in out])): r["q"] = q
    return out


def motifs(s, pos):
    """Core motifs and helix 3 of one aligned sequence (consensus positions given 1-based in pos)."""
    at = lambda k: s[k - 1]
    loop23 = "".join(at(k) for k in pos["loop23"])
    apical = "".join(at(k) for k in pos["apical"])
    return {"core_5p": "".join(at(k) for k in pos["core5"]), "core_3p": "".join(at(k) for k in pos["core3"]),
            "AAR_loop_2_3": bool(re.search("AA[AG]", loop23.replace("-", ""))),
            "AAR_apical": bool(re.search("AA[AG]", apical.replace("-", ""))), "loop_2_3": loop23, "apical_loop": apical}


def model_weights(cm, n):
    """Relative entropy (bits) of the match emissions of the model at each consensus position, against its null model
    (uniform in this model): ML and MR states for unpaired positions; for a base pair, the relative entropy of the 16
    pair emissions of the MP state, half to each position. Assumes a model without bifurcation (consensus positions
    then follow the nodes: left positions from the 5' end, right positions from the 3' end)."""
    text = open(cm).read()
    null = [float(x) for x in re.search(r"^NULL\s+(.*)$", text, re.M).group(1).split()]
    q = [0.25 * 2 ** x for x in null]
    q = [x / sum(q) for x in q]
    body = text.split("\nCM\n", 1)[1].split("\n//", 1)[0]
    if "[ BIF" in body: sys.exit("model_weights: bifurcations are not handled")
    w, left, right, node = [0.0] * n, 0, n + 1, None
    for line in body.splitlines():
        m = re.match(r"\s+\[ (\w+)\s+\d+ \]", line)
        if m:
            node = m.group(1)
            if node in ("MATP", "MATL"): left += 1
            if node in ("MATP", "MATR"): right -= 1
            continue
        f = line.split()
        if not f or f[0] not in ("MP", "ML", "MR"): continue
        if f[0] == "MP" and node == "MATP":
            bg = [qa * qb for qa in q for qb in q]
            p = [b * 2 ** float(x) for b, x in zip(bg, f[-16:])]
            t = sum(p); p = [x / t for x in p]
            re_ = sum(x * math.log2(x / b) for x, b in zip(p, bg) if x > 0)
            w[left - 1] += re_ / 2; w[right - 1] += re_ / 2
        elif (f[0] == "ML" and node == "MATL") or (f[0] == "MR" and node == "MATR"):
            p = [b * 2 ** float(x) for b, x in zip(q, f[-4:])]
            t = sum(p); p = [x / t for x in p]
            w[(left if f[0] == "ML" else right) - 1] += sum(x * math.log2(x / b) for x, b in zip(p, q) if x > 0)
    return w


def ortholog_events(real, rmaps, windows, pairs, hel, npos, core, genomes=("OKI", "OSKA", "OdB3"), flank=50):
    """Substitutions between the orthologous SECIS hits of the three O. dioica assemblies. The hit of the first genome
    is placed on the consensus positions by cmalign; the windows of the other genomes are aligned to it by sequence
    (local alignment of the hit and flank nt on each side, match 2, mismatch -3, gaps -5/-2), so that the structure
    does not drive the alignment of the orthologs. Three-taxon parsimony: at a position or pair where two genomes agree
    and the third differs, one change from the shared state (that of the internal node) to the third. A change at a
    pair whose shared state is canonical either keeps a canonical pair (single change, or compensatory when both bases
    differ) or not; expected fractions of kept pairs if the new base(s) were random: for a single change, the fraction
    of the three other bases that pair with the unchanged partner; for a double change, the fraction of the nine
    combinations that are canonical, each new base weighted by the fraction of transitions among the changes at the
    other unpaired positions (transitions: A-G, C-U; each of the two transversions gets half of the rest). Changes from
    a non-canonical shared pair are counted apart (to canonical: restoring).
    Sites (positions of the hits): changes per site at paired positions, at the core motifs (core) and at the other
    unpaired positions."""
    from Bio import Align
    al = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-2)
    by_gene = collections.defaultdict(dict)
    for k in real:
        name, span = k.rsplit("/", 1)
        genome, gene = name.split("|")
        by_gene[gene][genome] = (k, int(span.split("-")[0]))
    paired = {p for ij in pairs for p in ij}
    rows, events, sites = [], collections.Counter(), collections.Counter()
    pair_changes, unpaired_changes = [], []
    for gene, hits in sorted(by_gene.items()):
        if not all(g in hits for g in genomes): continue
        ref = genomes[0]
        kA, frmA = hits[ref]
        wA = windows[f"{gene}|{ref}"]
        absA = [None if x is None else frmA - 1 + x for x in rmaps[kA]]
        lo = max(0, min(x for x in absA if x is not None) - flank)
        hi = max(x for x in absA if x is not None) + flank + 1
        cols = {ref: [wA[x] if x is not None else "-" for x in absA]}
        ok = True
        for g in genomes[1:]:
            wB = windows[f"{gene}|{g}"]
            aln = al.align(wA[lo:hi], wB)[0]
            m = {}
            for (a0, a1), (b0, b1) in zip(*aln.aligned):
                for d in range(a1 - a0): m[lo + a0 + d] = b0 + d
            mapped = [x for x in absA if x is not None and x in m]
            ident = sum(wA[x] == wB[m[x]] for x in mapped) / max(1, len([x for x in absA if x is not None]))
            if ident < 0.5: ok = False
            cols[g] = [wB[m[x]] if x is not None and x in m else "-" for x in absA]
        if not ok:
            rows.append({"gene": gene, "position": "", "helix": "", "states": "", "event": "orthologs not aligned (identity < 0.5)"})
            continue
        cols = {g: "".join(v).replace("T", "U") for g, v in cols.items()}
        for (i, j) in pairs:
            st = [cols[g][i] + cols[g][j] for g in genomes]
            if any("-" in x for x in st): ev = "gap"
            else:
                c = collections.Counter(st)
                if len(c) == 1: ev = "none"
                elif len(c) == 3: ev = "ambiguous (three states)"
                else:
                    anc, new = c.most_common()[0][0], c.most_common()[1][0]
                    both = anc[0] != new[0] and anc[1] != new[1]
                    if anc in CANON:
                        if both: ev = "compensatory" if new in CANON else "double change, pair lost"
                        else: ev = "single change, pair kept" if new in CANON else "single change, pair lost"
                        pair_changes.append((anc, new, both))
                    else:
                        ev = "restoring" if new in CANON else "non-canonical to non-canonical"
            events[ev] += 1
            rows.append({"gene": gene, "position": f"{i + 1}-{j + 1}", "helix": hel[(i, j)], "states": ",".join(st), "event": ev})
        for k in range(npos):
            st = [cols[g][k] for g in genomes]
            if "-" in st: continue
            cat = "paired" if k in paired else ("core" if k + 1 in core else "other unpaired")
            sites[f"{cat} sites"] += 1
            if len(set(st)) > 1: sites[f"{cat} sites with a change"] += 1
            c = collections.Counter(st)
            if cat == "other unpaired" and len(c) == 2: unpaired_changes.append((c.most_common()[0][0], c.most_common()[1][0]))
    for cat in ("paired", "core", "other unpaired"): sites[f"{cat} sites with a change"] += 0
    ts = {("A", "G"), ("G", "A"), ("C", "U"), ("U", "C")}
    t = sum(x in ts for x in unpaired_changes) / len(unpaired_changes) if unpaired_changes else 1 / 3
    sites["transition fraction at the other unpaired positions"] = f"{t:.3f} ({len(unpaired_changes)} changes)"
    pnew = lambda a, b: t if (a, b) in ts else (1 - t) / 2
    exp = collections.Counter()
    for anc, new, both in pair_changes:
        if both:
            exp["double"] += sum(pnew(anc[0], x) * pnew(anc[1], y) for x in "ACGU" for y in "ACGU"
                                 if x != anc[0] and y != anc[1] and x + y in CANON)
        elif anc[0] != new[0]:
            exp["single"] += sum(pnew(anc[0], b) for b in "ACGU" if b != anc[0] and b + anc[1] in CANON)
        else:
            exp["single"] += sum(pnew(anc[1], b) for b in "ACGU" if b != anc[1] and anc[0] + b in CANON)
    return rows, events, sites, exp


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cm", default="results/models/oikopleura_secis.cm")
    ap.add_argument("--seeds", default="results/models/oikopleura_secis.seeds.fa")
    ap.add_argument("--windows", nargs="+", default=[f"curation/gold_standard/{g}.after_stop.fa" for g in ("OKI", "OSKA", "OdB3")])
    ap.add_argument("--extra", nargs="+", default=["TXNL1.b", "KKSU"], help="genes added to the gold standard after training")
    ap.add_argument("--baseline", default="curation/baseline_sec_genes.tsv", help="gene families")
    ap.add_argument("--training-windows", default="results/tables/stage2_secis_v2_windows.tsv",
                    help="lengths of the training windows (some are shortened at the next gene)")
    ap.add_argument("--shuffles", type=int, default=20, help="shuffled copies of each window")
    ap.add_argument("--decoys", type=int, default=0, help="decoys kept (0: as many as the seeds)")
    ap.add_argument("--cluster-identity", type=float, default=0.7)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--prefix", default="stage2_secis_v2_model")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 2)))
    ap.add_argument("--work", default="work/stage2/secis_model_description", help="cmsearch tables, reused")
    ap.add_argument("--reference-genome", default="OKI", choices=["OKI", "OSKA", "OdB3"],
                    help="assembly whose hit is placed on the model in the comparison of orthologs")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    family = {}
    for r in csv.DictReader(open(a.baseline), delimiter="\t"):
        family[r["gene"]] = r.get("family", r["gene"])
    fam_of = lambda gene: family.get(gene, re.split(r"[._]", gene)[0])

    with tempfile.TemporaryDirectory() as tmp:
        seeds = dict(read_fasta(a.seeds))
        aln, cons, ss = cmalign(a.cm, seeds, tmp, "seeds")
        pairs = pairs_of(ss)
        hel = helices(pairs)
        # positions of the motifs, from the consensus structure: the internal loop between helices 1 and 2, and the
        # loop between helices 2 and 3 (5' strand) and the apical loop
        h = collections.defaultdict(list)
        for p, k in hel.items(): h[k].append(p)
        h1, h2, h3 = (sorted(h[k]) for k in (1, 2, 3))
        loop12_5 = list(range(h1[-1][0] + 2, h2[0][0] + 1))          # 1-based, between the last 5' base of h1 and h2
        loop12_3 = list(range(h2[0][1] + 2, h1[-1][1] + 1))
        pos = {"core5": loop12_5[-4:-1], "core3": loop12_3[1:3],
               "loop23": list(range(h2[-1][0] + 2, h3[0][0] + 1)), "apical": list(range(h3[-1][0] + 2, h3[-1][1] + 1))}
        # windows of the gold-standard genes, for the decoys and the new genes
        windows = {}
        for f in a.windows:
            for k, s in read_fasta(f): windows[k] = s.upper()
        # shortened as in training (at the next selenoprotein coding sequence: SELENOW_pair2 before SELENOW_pair1)
        for r in csv.DictReader(open(a.training_windows), delimiter="\t"):
            if r["set"] != "positive": continue
            g, gene = r["id"].split("|", 1)
            for k in windows:
                if k.split("|")[:2] == [gene, g]: windows[k] = windows[k][:int(r["length"])]
        decoy_fa = os.path.join(tmp, "decoys.fa")
        with open(decoy_fa, "w") as o:
            for k, s in windows.items():
                for c in range(a.shuffles): o.write(f">{k.split('|')[1]}|{k.split('|')[0]}|shuffled{c}\n{dinucleotide_shuffle(s, rng)}\n")

        def best_hits(fa, tag):
            # --max: no filter, so that every window has its best hit; the search is kept in --work and reused when
            # the model and the sequences are the same
            os.makedirs(a.work, exist_ok=True)
            key = hashlib.md5(open(a.cm, "rb").read() + open(fa, "rb").read()).hexdigest()
            tbl = os.path.join(a.work, f"{tag}.{key}.tbl")
            if not os.path.exists(tbl):
                run(["cmsearch", "--max", "--toponly", "-T", "0", "--cpu", str(a.cpus), "--tblout", tbl + ".tmp", "-o",
                     os.devnull, a.cm, fa])
                os.replace(tbl + ".tmp", tbl)
            best = {}
            for l in open(tbl):
                if l.startswith("#"): continue
                f = l.split()
                t, frm, to, bits = f[0], int(f[7]), int(f[8]), float(f[14])
                if t not in best or bits > best[t][2]: best[t] = (frm, to, bits)
            seqs = dict(read_fasta(fa))
            return {f"{t}/{frm}-{to}": (seqs[t][frm - 1:to], bits) for t, (frm, to, bits) in best.items()}

        dec = best_hits(decoy_fa, "decoys")
        n_dec = a.decoys or len(seeds)
        dec = dict(sorted(dec.items(), key=lambda kv: -kv[1][1])[:n_dec])
        daln, _, _ = cmalign(a.cm, {k: v[0] for k, v in dec.items()}, tmp, "decoys")
        dmot = [motifs(x, pos) for x in daln.values()]
        # best hit in each gold-standard window: the new genes, and the orthologous SECIS of the O. dioica assemblies
        real_fa = os.path.join(tmp, "real.fa")
        with open(real_fa, "w") as o:
            for k, s in windows.items(): o.write(f">{k.split('|')[1]}|{k.split('|')[0]}\n{s}\n")
        real = best_hits(real_fa, "real")
        rmaps = {}
        raln, _, _ = cmalign(a.cm, {k: v[0] for k, v in real.items()}, tmp, "real", rmaps)
        ext = {k: v for k, v in real.items() if k.split("|")[1].split("/")[0] in a.extra}
        ealn = {k: raln[k] for k in ext}

    names = list(aln)
    pid = {(x, y): identity(aln[x], aln[y]) for x, y in itertools.combinations(names, 2)}
    get = lambda x, y: pid.get((x, y), pid.get((y, x)))
    cl = sorted(clusters(names, pid, a.cluster_identity), key=lambda c: (-len(c), names.index(c[0])))
    cluster_of = {n: k + 1 for k, c in enumerate(cl) for n in c}
    reps = [c[0] for c in cl]
    dnames = list(daln)
    dpid = {(x, y): identity(daln[x], daln[y]) for x, y in itertools.combinations(dnames, 2)}
    dreps = [c[0] for c in clusters(dnames, dpid, a.cluster_identity)]

    def helix_pairs(s):
        return {k: sum(s[i] + s[j] in CANON for (i, j) in pairs if hel[(i, j)] == k) for k in (1, 2, 3)}

    # seeds and new genes
    rows = []
    for n in names + list(ealn):
        s = aln.get(n) or ealn[n]
        genome, rest = n.split("|", 1)
        gene, rng_ = rest.rsplit("/", 1)
        start = int(rng_.split("-")[0])
        others = [(identity(s, aln[m]), m) for m in names if m != n]
        best = max(others)
        hp = helix_pairs(s)
        row = {"sequence": n, "set": "seed" if n in aln else "new gene", "genome": genome, "gene": gene,
               "family": fam_of(gene.split(".")[0]) if n in aln else gene, "after_stop": start - 1,
               "length": len(s.replace("-", "")), "bits": f"{ext[n][1]:.1f}" if n in ext else "",
               "nearest_seed_identity": f"{best[0]:.3f}", "nearest_seed": best[1], "cluster": cluster_of.get(n, ""),
               "helix1_canonical": f"{hp[1]}/5", "helix2_canonical": f"{hp[2]}/{len(h2)}", "helix3_canonical": f"{hp[3]}/{len(h3)}"}
        row.update(motifs(s, pos))
        rows.append(row)
    os.makedirs(a.tables, exist_ok=True)
    cols = list(rows[0])
    with open(f"{a.tables}/{a.prefix}_seeds.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r[c]) for c in cols) + "\n")

    # identities between seeds; orthologous: the same gene in two O. dioica assemblies
    dioica = {"OKI", "OSKA", "OdB3"}
    gene_of = lambda n: n.split("|")[1].rsplit("/", 1)[0]
    with open(f"{a.tables}/{a.prefix}_identities.tsv", "w") as o:
        o.write("seed_1\tseed_2\tidentity\torthologous\n")
        for (x, y), v in pid.items():
            orth = x.split("|")[0] in dioica and y.split("|")[0] in dioica and gene_of(x) == gene_of(y)
            o.write(f"{x}\t{y}\t{v:.3f}\t{'yes' if orth else 'no'}\n")

    # positions
    weights = model_weights(a.cm, len(cons))
    with open(f"{a.tables}/{a.prefix}_positions.tsv", "w") as o:
        o.write("position\tconsensus\tstructure\tA\tC\tG\tU\tgap\tmost_frequent_fraction\tinformation_bits\t"
                "model_relative_entropy_bits\n")
        for k in range(len(cons)):
            c = collections.Counter(aln[n][k] for n in names)
            res = sum(c[x] for x in "ACGU")
            top = max("ACGU", key=lambda x: c[x])
            info = 2 + sum(c[x] / res * math.log2(c[x] / res) for x in "ACGU" if c[x])
            o.write(f"{k + 1}\t{cons[k]}\t{ss[k]}\t{c['A']}\t{c['C']}\t{c['G']}\t{c['U']}\t{c['-']}\t{c[top] / res:.3f}\t"
                    f"{info:.3f}\t{weights[k]:.3f}\n")

    # pairs
    ps = pair_stats(aln, reps, cons, pairs, hel, a.permutations, rng)
    pd = pair_stats(daln, dreps, cons, pairs, hel, a.permutations, rng)
    pcols = ["helix", "i", "j", "consensus", "sequences", "canonical", "canonical_fraction", "representatives", "pair_types",
             "mi_bits", "p", "q", "compensatory"]
    with open(f"{a.tables}/{a.prefix}_pairs.tsv", "w") as o:
        o.write("set\t" + "\t".join(pcols) + "\n")
        for tag, rs in (("seeds", ps), ("decoys", pd)):
            for r in rs:
                o.write(tag + "\t" + "\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in pcols) + "\n")

    # summary
    v = sorted(pid.values())
    nn = [max(get(n, m) for m in names if m != n) for n in names]
    sd = [r for r in rows if r["set"] == "seed"]
    summ = [("seeds", len(names)), ("decoys (best hits in shuffled windows)", len(dnames)),
            ("decoy scores (bits)", f"{min(x[1] for x in dec.values()):.1f}-{max(x[1] for x in dec.values()):.1f}"),
            ("seed pairs", len(v)), ("pairwise identity, median", f"{statistics.median(v):.3f}"),
            ("pairwise identity, quartiles", "{:.3f}-{:.3f}".format(*statistics.quantiles(v, n=4)[::2])),
            ("pairwise identity, maximum", f"{v[-1]:.3f}")]
    summ += [(f"seed pairs with identity >= {t}", sum(x >= t for x in v)) for t in (0.6, 0.7, 0.8, 0.9)]
    summ += [("nearest-seed identity, median", f"{statistics.median(nn):.3f}"),
             (f"clusters at {a.cluster_identity}", len(cl)),
             (f"largest cluster at {a.cluster_identity}", len(cl[0])),
             ("decoy pairwise identity, median", f"{statistics.median(dpid.values()):.3f}"),
             ("seeds with core UGA (5' strand)", sum(r["core_5p"] == "UGA" for r in sd)),
             ("seeds with core GA (3' strand)", sum(r["core_3p"] == "GA" for r in sd)),
             ("seeds with AAR between helices 2 and 3", sum(r["AAR_loop_2_3"] for r in sd)),
             ("seeds with AAR in the apical loop only", sum(r["AAR_apical"] and not r["AAR_loop_2_3"] for r in sd)),
             ("seeds with at least 2 of 3 canonical pairs in helix 3", sum(int(r["helix3_canonical"][0]) >= 2 for r in sd))]
    for tag, rs, al in (("seeds", ps, aln), ("decoys", pd, daln)):
        for k in (1, 2, 3):
            hr = [r for r in rs if r["helix"] == k]
            summ += [(f"{tag}: {HELICES[k]} pairs", len(hr)),
                     (f"{tag}: {HELICES[k]} canonical fraction (mean over pairs)", f"{statistics.mean(r['canonical_fraction'] for r in hr):.3f}"),
                     (f"{tag}: {HELICES[k]} pairs with MI q < 0.05", sum(r["q"] < 0.05 for r in hr)),
                     (f"{tag}: {HELICES[k]} representatives with a canonical pair differing at both positions from the consensus pair", sum(r["compensatory"] for r in hr))]
        summ += [(f"{tag}: representatives", len(reps if tag == "seeds" else dreps))]
    core = set(pos["core5"]) | set(pos["core3"]) | set(pos["loop23"][1:4])
    orows, oev, osites, oexp = ortholog_events(real, rmaps, {k.split("|")[0] + "|" + k.split("|")[1]: v for k, v in windows.items()},
                                               pairs, hel, len(cons), core,
                                               genomes=tuple([a.reference_genome] + [g for g in ("OKI", "OSKA", "OdB3") if g != a.reference_genome]))
    with open(f"{a.tables}/{a.prefix}_orthologs.tsv", "w") as o:
        o.write(f"gene\tposition\thelix\tstates ({','.join([a.reference_genome] + [g for g in ('OKI', 'OSKA', 'OdB3') if g != a.reference_genome])})\tevent\n")
        for r in orows: o.write(f"{r['gene']}\t{r['position']}\t{r['helix']}\t{r['states']}\t{r['event']}\n")
    from scipy.stats import binomtest
    summ += [("orthologs: genes compared (three assemblies)", len({r["gene"] for r in orows if r["position"]}))]
    summ += [(f"orthologs: pair events, {k}", v) for k, v in sorted(oev.items())]
    for kind, kept, lost in (("single", "single change, pair kept", "single change, pair lost"),
                             ("double", "compensatory", "double change, pair lost")):
        n = oev[kept] + oev[lost]
        if n:
            e = oexp[kind] / n
            summ += [(f"orthologs: {kind} changes from a canonical pair, kept canonical", f"{oev[kept]}/{n}"),
                     (f"orthologs: {kind} changes, expected fraction kept if random", f"{e:.3f}"),
                     (f"orthologs: {kind} changes, binomial p (kept more than expected)",
                      f"{binomtest(oev[kept], n, e, alternative='greater').pvalue:.2g}")]
    summ += [(f"orthologs: {k}", v) for k, v in sorted(osites.items())]
    summ += [("decoys with core UGA (5' strand)", sum(m["core_5p"] == "UGA" for m in dmot)),
             ("decoys with core GA (3' strand)", sum(m["core_3p"] == "GA" for m in dmot)),
             ("decoys with AAR between helices 2 and 3", sum(m["AAR_loop_2_3"] for m in dmot))]
    with open(f"{a.tables}/{a.prefix}_summary.tsv", "w") as o:
        o.write("measure\tvalue\n")
        for k, x in summ: o.write(f"{k}\t{x}\n")
    for k, x in summ: print(f"{k}\t{x}")


if __name__ == "__main__":
    main()
