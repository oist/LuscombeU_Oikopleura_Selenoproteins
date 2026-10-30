#!/usr/bin/env python3
"""Stem-loop screen for SECIS elements that the covariance models miss (2026-10-31; first use: SELENOI).

The screen does not start from a SECIS model. Step 1 lists the stable unbranched stem-loops of each sequence; step 2
describes each with the features of SECIS elements; step 3 scores them with the covariance models without their
filters and down to very low scores.

1. Stem-loops. Local minimum free energy structures (ViennaRNA, windows of --span nt, as RNALfold) are split into
   unbranched stem-loop domains: the stretch closed by a base pair that encloses one hairpin loop and no multiloop.
   Domains of --min-length to --max-length nt with at least --min-pairs base pairs are refolded alone, and their
   minimum free energy is compared with that of --shuffles dinucleotide-preserving shuffles of their sequence
   (z-score). Overlapping domains are all kept (a SECIS element can be part of a larger domain).
2. SECIS features of each domain, read on its structure:
   - quartet: a helix of the domain whose closing base pair (i, j) is preceded on the 5' side by UGA (ending at
     i - 1 or i - 2) and followed on the 3' side by GA (starting at j + 1 or j + 2), the arrangement of the
     non-Watson-Crick quartet below helix 2 (5'-AUGA ... GA-3', with the G·A pairs not formed in the predicted
     structure); also GGGA instead of UGA [Novoselov et al. 2007];
   - helix_2: base pairs from that closing pair to the first unpaired base on both strands (bulges of 1 nt allowed);
   - AAR: two or three unpaired A (AA or AAA, then R) on the 5' strand or in the loop, 9 to 16 nt after the 5' UGA
     (the AAR of form 1 and form 2 elements); CC: two unpaired C in the same place [Korotkov et al. 2002].
3. Best scores of the Oikopleura SECIS model and of RF00031 (and of any other --models) on the domain with --flank nt
   on each side, from cmsearch --max (no filters) --toponly -T --min-bits.
Outputs: --out (one row per domain). Coordinates are 1-based on the input sequences, which are read as given (sense
strand).
"""
import argparse, collections, concurrent.futures, os, random, re, subprocess, tempfile

import RNA


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].split()[0], []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def dinuc_shuffle(s, rng):
    """Dinucleotide-preserving shuffle (Altschul and Erickson 1985; random Eulerian path)."""
    if len(s) < 3: return s
    edges = collections.defaultdict(list)
    for a, b in zip(s, s[1:]): edges[a].append(b)
    last = s[-1]
    while True:   # last edges forming an arborescence towards the last letter
        lastedge = {v: rng.choice(es) for v, es in edges.items() if v != last}
        ok = True
        for v in lastedge:
            seen, u = set(), v
            while u != last:
                if u in seen or u not in lastedge: ok = False; break
                seen.add(u); u = lastedge[u]
            if not ok: break
        if ok: break
    lists = {}
    for v, es in edges.items():
        es = list(es)
        if v in lastedge: es.remove(lastedge[v])
        rng.shuffle(es)
        if v in lastedge: es.append(lastedge[v])
        lists[v] = es
    out, u = [s[0]], s[0]
    idx = collections.Counter()
    for _ in range(len(s) - 1):
        u = lists[u][idx[u]]; idx[out[-1]] += 1; out.append(u)
    return "".join(out)


def pairs_of(db):
    st, p = [], [-1] * len(db)
    for i, c in enumerate(db):
        if c == "(": st.append(i)
        elif c == ")": j = st.pop(); p[i], p[j] = j, i
    return p


def domains(db, offset, max_length):
    """Unbranched stem-loop domains of a dot-bracket structure, one per hairpin loop: from the hairpin outwards
    through stacks, bulges and interior loops, up to the outermost pair that keeps the domain within max_length nt;
    (start, end) 0-based inclusive, shifted by offset."""
    p = pairs_of(db)
    n, out = len(db), []
    for k in range(n):
        l = p[k]
        if l <= k or any(p[x] != -1 for x in range(k + 1, l)): continue   # not a hairpin-closing pair
        i, j = k, l
        while True:
            a, b = i - 1, j + 1
            while a >= 0 and p[a] == -1: a -= 1
            while b < n and p[b] == -1: b += 1
            if a < 0 or b >= n or p[a] != b or b - a + 1 > max_length: break
            i, j = a, b
        out.append((i + offset, j + offset))
    return out


def features(seq, db):
    """SECIS features read on the structure of a domain (sequence in RNA letters)."""
    p = pairs_of(db)
    n = len(seq)
    best = None
    for i in range(n):
        j = p[i]
        if j <= i: continue
        if i > 0 and p[i - 1] == j + 1: continue     # not the closing pair of a helix
        for s5 in (1, 2):
            # 5' UGA ending at i - s5 (s5 = 1: UGA immediately before i; 2: one spacer)
            u = seq[i - s5 - 2:i - s5 + 1] if i - s5 - 2 >= 0 else ""
            g = seq[i - s5 - 3:i - s5 + 1] if i - s5 - 3 >= 0 else ""
            kind = "UGA" if u == "UGA" else "GGGA" if g == "GGGA" else None
            if kind is None: continue
            for s3 in (1, 2):
                if seq[j + s3:j + s3 + 2] != "GA": continue
                # helix 2: pairs from (i, j) inwards, 1-nt bulges allowed
                k, l, bp = i, j, 1
                while True:
                    if k + 1 < l and p[k + 1] == l - 1: k, l = k + 1, l - 1
                    elif k + 2 < l and p[k + 2] == l - 1 and p[k + 1] == -1: k, l = k + 2, l - 1
                    elif k + 1 < l - 1 and p[k + 1] == l - 2 and p[l - 1] == -1: k, l = k + 1, l - 2
                    else: break
                    bp += 1
                q = i - s5   # last base of the 5' UGA
                window = seq[q + 9:q + 17]
                unp = "".join(c if p[q + 9 + x] == -1 else "." for x, c in enumerate(window)) if q + 17 <= n else ""
                aar = bool(re.search(r"AA[AG]", unp)) or bool(re.search(r"AA", unp))
                cc = "CC" in unp
                cand = (bp, aar or cc, kind, i, j, s5, s3, aar, cc)
                if best is None or (cand[1], cand[0]) > (best[1], best[0]): best = cand
    if best is None: return {"quartet": "", "helix_2_pairs": "", "AAR": "", "CC": ""}
    bp, _, kind, i, j, s5, s3, aar, cc = best
    return {"quartet": f"{kind}@{i - s5 - (2 if kind == 'UGA' else 3) + 1}/GA@{j + s3 + 1}", "helix_2_pairs": bp,
            "AAR": "yes" if aar else "", "CC": "yes" if cc else ""}


def cm_scores(cands, models, flank, seqs, min_bits, cpus, work):
    """Best score of each model on each domain (with flanks): cmsearch --max --toponly."""
    out = collections.defaultdict(dict)
    fa = os.path.join(work, "domains.fa")
    with open(fa, "w") as o:
        for k, c in enumerate(cands):
            s = seqs[c["sequence"]]
            o.write(f">d{k}\n{s[max(0, c['start'] - 1 - flank):c['end'] + flank]}\n")
    for m, cm in models.items():
        tbl = os.path.join(work, f"{m}.tbl")
        subprocess.run(["cmsearch", "--max", "--toponly", "-T", str(min_bits), "--cpu", str(cpus), "--tblout", tbl,
                        "-o", os.devnull, cm, fa], check=True)
        for l in open(tbl):
            if l.startswith("#"): continue
            f = l.split()
            k = int(f[0][1:])
            out[k][m] = max(out[k].get(m, -1e9), float(f[14]))
    return out


def stem_loops(job):
    """Stem-loop domains of one sequence, with z-scores and features."""
    name, s, a = job
    rng = random.Random(f"{a.seed}|{name}")
    md = RNA.md(); md.window_size = a.span; md.max_bp_span = a.span
    found, out, cache = set(), [], {}
    fc = RNA.fold_compound(s, md, RNA.OPTION_WINDOW)
    fc.mfe_window_cb(lambda i, j, db, e, *rest: found.update(domains(db, i - 1, a.max_length)), None)
    for d0, d1 in sorted(found):
        L = d1 - d0 + 1
        if not a.min_length <= L <= a.max_length: continue
        sub = s[d0:d1 + 1]
        db, e = RNA.fold(sub)
        if db.count("(") < a.min_pairs: continue
        if sub not in cache:
            sh = [RNA.fold(dinuc_shuffle(sub, rng))[1] for _ in range(a.shuffles)]
            m = sum(sh) / len(sh); sd = (sum((x - m) ** 2 for x in sh) / len(sh)) ** 0.5
            cache[sub] = (e - m) / sd if sd > 0 else 0.0
        c = {"sequence": name, "start": d0 + 1, "end": d1 + 1, "length": L, "pairs": db.count("("),
             "mfe": f"{e:.1f}", "z": f"{cache[sub]:.2f}", "structure": db, "rna": sub}
        c.update(features(sub, db))
        out.append(c)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", nargs="+")
    ap.add_argument("--span", type=int, default=150)
    ap.add_argument("--min-length", type=int, default=40)
    ap.add_argument("--max-length", type=int, default=150)
    ap.add_argument("--min-pairs", type=int, default=12)
    ap.add_argument("--shuffles", type=int, default=100)
    ap.add_argument("--flank", type=int, default=10)
    ap.add_argument("--models", nargs="*", default=["oikopleura=results/models/oikopleura_secis.cm",
                                                     "RF00031=data/references/RF00031.cm"])
    ap.add_argument("--min-bits", type=float, default=-50.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    seqs = {}
    for f in a.fasta:
        for k, s in read_fasta(f): seqs[k] = s.upper().replace("T", "U")
    cands = []
    with concurrent.futures.ProcessPoolExecutor(a.cpus) as ex:
        for name, cs in zip(seqs, ex.map(stem_loops, [(k, v, a) for k, v in seqs.items()])):
            cands += cs
            print(f"{name}: {len(cs)} stem-loops", flush=True)
    models = dict(m.split("=", 1) for m in a.models)
    models = {k: v for k, v in models.items() if os.path.exists(v)}
    with tempfile.TemporaryDirectory() as work:
        sc = cm_scores(cands, models, a.flank, seqs, a.min_bits, a.cpus, work) if cands and models else {}
    cols = ["sequence", "start", "end", "length", "pairs", "mfe", "z", "quartet", "helix_2_pairs", "AAR", "CC"] + \
           [f"bits_{m}" for m in models] + ["structure", "rna"]
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for k, c in enumerate(cands):
            for m in models:
                v = sc.get(k, {}).get(m)
                c[f"bits_{m}"] = f"{v:.1f}" if v is not None else ""
            o.write("\t".join(str(c.get(x, "")) for x in cols) + "\n")
    print(f"{len(cands)} stem-loops -> {a.out}")


if __name__ == "__main__":
    main()
