#!/usr/bin/env python3
"""A wider search for the SECIS element of SELENOI (2026-10-31), starting from stem-loops rather than from a SECIS
model. Reads repository files only: the regions of make selenoi-regions (results/sequences/selenoi_regions.fa), the
candidate protein (curation/candidates.faa), the stage 1 selenocysteine sites, the SECIS models and, as controls,
the gold-standard 3' sequences (curation/gold_standard/<genome>.after_stop.fa).

1. Candidate codon in each appendicularian region (one per genome; for O. dioica, the region centred on the codon
   placed by the transcripts): tblastn (legacy BLAST, no filter, E <= 1) of the --peptide residues before the first
   U of the candidate protein, followed by X; in the best HSP that reaches the last of these residues, on the sense
   strand, the codon after it, with the 8 codons before and the codons to the next stop codon (TAA, TAG or TGA).
   Vertebrates and C. intestinalis (positive controls of the searches): the stage 1 selenocysteine codon.
2. Transcripts of the three O. dioica assemblies (the SELENOI transcripts of the regions file) aligned to the region
   of their assembly (minimap2 -x splice:hq -un): whether each holds the candidate codon in an exon, or splices out a
   segment that contains it (length, end dinucleotides, position of the codon in it), and where the transcript ends
   after the codon.
3. 3' sequences: from the candidate codon to --after nt downstream (appendicularians, C. intestinalis), or
   --after-vertebrates nt (vertebrates, whose SECIS can lie far in a long 3' UTR).
4. Searches on each 3' sequence:
   a. cmsearch --max (no filters) --toponly -T 0 with each model (Oikopleura, first appendicularian model, RF00031):
      best hit in the whole sequence and in its first --near nt;
   b. background: --shuffles dinucleotide-preserving shuffles of each appendicularian 3' sequence, best hit of each
      model in each (median, 95th percentile and maximum, pooled);
   c. stem-loop screen (stage2/secis_stemloop_screen.py) on the first --near nt: stem-loops, and those with the SECIS
      features (UGA ... GA arrangement around a helix, unpaired AA or CC above it); all written to --stemloops;
   d. core patterns without folding (secis_stemloop_screen.core_patterns) in the first --near nt, leaving out
      patterns whose UGA is the candidate codon itself; background rate from the shuffles;
   e. comparative: the first --compare nt of the appendicularian 3' sequences with a TGA aligned with MAFFT
      (L-INS-i), consensus structure and energy with RNAalifold (ViennaRNA), pairwise identities; alignment written to
      --alignment (Stockholm, with the consensus structure).
5. Calibration of 4c and 4d on the gold-standard 3' sequences: known SECIS = best Oikopleura-model hit of at least
   20 bits starting in the first 700 nt; a known SECIS is recovered by a stem-loop that overlaps at least half of
   it (with or without the features), or by a core pattern that overlaps it; background = nt 1001-2000 of the
   sequences without such a hit there.
Outputs: --table (one row per 3' sequence), --transcripts-table, --stemloops, --calibration, --alignment.
"""
import argparse, collections, concurrent.futures, csv, glob, os, random, re, shutil, statistics, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import RNA                                                                          # noqa: E402
from secis_stemloop_screen import core_patterns, dinuc_shuffle, stem_loops         # noqa: E402

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
APPENDICULARIANS = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]
VERTEBRATES = ["Hsap", "Mmus", "Drer"]


def tr(s): return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].strip(), []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def run(cmd, **kw): return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def cm_best(seqs, models, cpus, work, tag):
    """{name: {model: [(bits, from, to)]}} from cmsearch --max --toponly -T 0."""
    fa = f"{work}/{tag}.fa"
    with open(fa, "w") as o:
        for k, s in seqs.items(): o.write(f">{k}\n{s}\n")
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for m, cm in models.items():
        tbl = f"{work}/{tag}.{m}.tbl"
        run(["cmsearch", "--max", "--toponly", "-T", "0", "--cpu", str(cpus), "--tblout", tbl, "-o", os.devnull, cm, fa])
        for l in open(tbl):
            if l.startswith("#"): continue
            f = l.split()
            out[f[0]][m].append((float(f[14]), int(f[7]), int(f[8])))
    return out


def best_in(hits, lo, hi):
    hs = [h for h in hits if lo <= h[1] <= hi]
    return max(hs) if hs else None


def fmt(h): return f"{h[0]:.1f} ({h[1]}-{h[2]})" if h else "no hit"


def screen_args(a):
    return argparse.Namespace(seed=a.seed, span=150, min_length=40, max_length=150, min_pairs=12,
                              shuffles=a.stemloop_shuffles)


def with_features(c): return bool(c["quartet"]) and bool(c["AAR"] or c["CC"])


def calibration(a, models, work):
    seqs = {}
    for f in sorted(glob.glob(f"{a.gold}/*.after_stop.fa")):
        for k, s in read_fasta(f): seqs[k.split()[0]] = s.upper()
    hits = cm_best(seqs, {"oikopleura": models["oikopleura"]}, a.cpus, work, "gold")
    known = {k: max(hs) for k, d in hits.items() for hs in [[h for h in d["oikopleura"] if h[0] >= 20 and h[1] <= 700]]
             if hs}
    late = {k for k, d in hits.items() if any(h[0] >= 20 and h[1] > 1000 for h in d["oikopleura"])}
    pos = {k: s[:700] for k, s in seqs.items() if k in known}
    neg = {k: s[1000:2000] for k, s in seqs.items() if k not in late and len(s) >= 2000}
    ns = screen_args(a)
    res = {}
    with concurrent.futures.ProcessPoolExecutor(a.cpus) as ex:
        for tag, d in (("pos", pos), ("neg", neg)):
            for k, cs in zip(d, ex.map(stem_loops, [(k, v.replace("T", "U"), ns) for k, v in d.items()])):
                res[(tag, k)] = cs
    rec = collections.Counter()
    for k, (bits, s, e) in known.items():
        cs = [c for c in res[("pos", k)] if min(c["end"], e) - max(c["start"], s) + 1 >= 0.5 * (e - s + 1)]
        rec["stem-loop"] += bool(cs)
        rec["stem-loop with the features"] += any(with_features(c) for c in cs)
        rec["core pattern"] += any(s - 20 <= p + 1 <= e for p, q, h, mm in core_patterns(pos[k]))
    kb = sum(len(v) for v in neg.values()) / 1000
    bg_sl = sum(len(res[("neg", k)]) for k in neg) / kb
    bg_f = sum(sum(with_features(c) for c in res[("neg", k)]) for k in neg) / kb
    bg_core = sum(len({p for p, *_ in core_patterns(v)}) for v in neg.values()) / kb
    rows = [("known SECIS (Oikopleura model >= 20 bits, first 700 nt)", len(known), ""),
            ("recovered by a stem-loop", rec["stem-loop"], f"{bg_sl:.1f} stem-loops per kb in background"),
            ("recovered by a stem-loop with the features", rec["stem-loop with the features"],
             f"{bg_f:.2f} per kb in background"),
            ("recovered by a core pattern", rec["core pattern"], f"{bg_core:.2f} per kb in background"),
            ("background (nt 1001-2000 of the gold-standard 3' sequences)", len(neg), f"{kb:.0f} kb")]
    with open(a.calibration, "w") as o:
        o.write("measure\tcount\tbackground\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", default="results/sequences/selenoi_regions.fa")
    ap.add_argument("--proteins", default="curation/candidates.faa")
    ap.add_argument("--gene", default="SELENOI")
    ap.add_argument("--stage1", default="results/tables/stage1_sec_sites.tsv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--models", nargs="*", default=["oikopleura=results/models/oikopleura_secis.cm",
                                                     "appendicularia=results/models/appendicularia_secis.cm",
                                                     "RF00031=data/references/RF00031.cm"])
    ap.add_argument("--peptide", type=int, default=45)
    ap.add_argument("--after", type=int, default=3000)
    ap.add_argument("--after-vertebrates", type=int, default=8000)
    ap.add_argument("--near", type=int, default=600)
    ap.add_argument("--compare", type=int, default=450)
    ap.add_argument("--shuffles", type=int, default=20)
    ap.add_argument("--stemloop-shuffles", type=int, default=50)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-calibration", action="store_true")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--table", default="results/tables/stage2_selenoi_secis_search.tsv")
    ap.add_argument("--transcripts-table", default="results/tables/stage2_selenoi_transcripts.tsv")
    ap.add_argument("--stemloops", default="results/tables/stage2_selenoi_stemloops.tsv")
    ap.add_argument("--calibration", default="results/tables/stage2_secis_stemloop_calibration.tsv")
    ap.add_argument("--alignment", default="results/sequences/selenoi_3prime.sto")
    a = ap.parse_args()
    models = {k: v for k, v in (m.split("=", 1) for m in a.models) if os.path.exists(v)}
    missing = [m.split("=", 1)[0] for m in a.models if m.split("=", 1)[0] not in models]
    if missing: print(f"models not found, left out: {', '.join(missing)}")
    regions = {k: s.upper() for k, s in read_fasta(a.regions)}
    prot = {k.split()[0]: s for k, s in read_fasta(a.proteins)}[a.gene]
    u = prot.find("U")
    work = tempfile.mkdtemp(prefix="selenoi_secis_")
    # 1. candidate codons
    genomic = {}
    for k in regions:
        g, kind = k.split("|")[:2]
        if g in APPENDICULARIANS and kind in ("codon", "tblastn_locus"):
            if g not in genomic or kind == "codon": genomic[g] = k
    with open(f"{work}/regions.fa", "w") as o:
        for g, k in genomic.items(): o.write(f">{g}\n{regions[k]}\n")
    with open(f"{work}/q.faa", "w") as o: o.write(f">q\n{prot[u - a.peptide:u].replace('U', 'X')}X\n")
    run(["formatdb", "-i", "regions.fa", "-p", "F"], cwd=work)
    best = {}
    for l in run(["blastall", "-p", "tblastn", "-d", "regions.fa", "-i", "q.faa", "-F", "F", "-e", "1", "-m", "8"],
                 cwd=work).splitlines():
        f = l.split("\t")
        qe, s0, s1, bits = int(f[7]), int(f[8]), int(f[9]), float(f[11])
        if qe >= a.peptide and s1 > s0 and (f[1] not in best or bits > best[f[1]][0]): best[f[1]] = (bits, s1, f[2])
    rows, three = [], {}
    for g, k in genomic.items():
        r = {"genome": g, "region": k.split("|")[2], "source": k.split("|")[1]}
        if g not in best:
            rows.append(dict(r, codon="not found")); continue
        bits, end, ident = best[g]
        s, c = regions[k], end            # 0-based start of the codon after the hit
        n = c + 3
        while n + 3 <= len(s) and s[n:n + 3] not in ("TAA", "TAG", "TGA"): n += 3
        r.update({"peptide_hit": f"{bits:.1f} bits, {ident}% identity", "codon": s[c:c + 3],
                  "context": tr(s[c - 24:c]) + "[" + s[c:c + 3] + "]" + tr(s[c + 3:n + 3]),
                  "codons_to_stop": (n - c - 3) // 3})
        rows.append(r)
        if s[c:c + 3] == "TGA": three[g] = (s[c:c + a.after], c)
    # vertebrates and C. intestinalis: stage 1 codon
    st1 = [x for x in csv.DictReader(open(a.stage1), delimiter="\t")
           if x["family"] == a.gene and x["genome"] not in APPENDICULARIANS]
    st1.sort(key=lambda x: x["selenoprofiles_label"] != "selenocysteine")
    for x in st1:
        g = x["genome"]
        if g in three: continue
        for k in regions:
            kg, kind, loc, strand = k.split("|")[:4]
            if kg != g or kind != "stage1_site": continue
            contig, span = loc.rsplit(":", 1); s0, e0 = map(int, span.split("-"))
            if contig != x["contig"] or strand != x["strand"]: continue
            p = int(x["sec_codon_start"])
            if not s0 <= p <= e0: continue
            off = p - s0 if strand == "+" else e0 - p - 2
            s = regions[k]
            if s[off:off + 3] != "TGA": continue
            L = a.after_vertebrates if g in VERTEBRATES else a.after
            three[g] = (s[off:off + L], off)
            rows.append({"genome": g, "region": loc, "source": f"stage 1 site ({x['sources']})", "codon": "TGA",
                         "context": tr(s[off - 24:off]) + "[TGA]" + tr(s[off + 3:off + 30])})
            break
    # 2. transcripts
    trows = []
    for g in ("OKI", "OSKA", "OdB3"):
        if g not in genomic: continue
        ref = f"{work}/{g}.ref.fa"
        open(ref, "w").write(f">{g}\n{regions[genomic[g]]}\n")
        txs = {k.split("|")[2]: s for k, s in regions.items() if k.startswith(f"{g}|transcript|")}
        with open(f"{work}/{g}.tx.fa", "w") as o:
            for k, s in txs.items(): o.write(f">{k}\n{s}\n")
        c = three[g][1] if g in three else None
        for l in run(["minimap2", "-x", "splice:hq", "-un", "-c", "--secondary=no", ref, f"{work}/{g}.tx.fa"]).splitlines():
            f = l.split("\t")
            if f[4] != "+" or c is None: continue
            cg = next(x[5:] for x in f[12:] if x.startswith("cg:Z:"))
            t = int(f[7])
            state = "transcript ends before the codon" if int(f[8]) <= c else "outside the alignment"
            for n, op in re.findall(r"(\d+)([MIDN])", cg):
                n = int(n)
                if op in "MD":
                    if t <= c < t + n: state = "in an exon"
                    t += n
                elif op == "N":
                    if t <= c < t + n:
                        s = regions[genomic[g]]
                        state = (f"spliced out: segment of {n} nt ({s[t:t + 2]}...{s[t + n - 2:t + n]}), codon at "
                                 f"{c - t + 1}")
                    t += n
            trows.append({"genome": g, "transcript": f[0], "length": f[1], "aligned": f"{int(f[2]) + 1}-{f[3]}",
                          "codon": state, "transcript_end_after_codon": int(f[8]) - c - 3})
    with open(a.transcripts_table, "w") as o:
        cols = ["genome", "transcript", "length", "aligned", "codon", "transcript_end_after_codon"]
        o.write("\t".join(cols) + "\n")
        for r in trows: o.write("\t".join(str(r[x]) for x in cols) + "\n")
    # 4a. covariance models
    seqs = {g: v[0] for g, v in three.items()}
    hits = cm_best(seqs, models, a.cpus, work, "three")
    # 4b. background
    rng = random.Random(a.seed)
    shuf = {f"{g}_{i}": dinuc_shuffle(s, rng) for g, s in seqs.items() if g in APPENDICULARIANS
            for i in range(a.shuffles)}
    sh = cm_best(shuf, models, a.cpus, work, "shuffled") if shuf else {}
    bg = {}
    for m in models:
        shb = sorted(max((h[0] for h in sh.get(k, {}).get(m, [])), default=0.0) for k in shuf)
        if shb: bg[m] = (f"median {statistics.median(shb):.1f}, 95th percentile {shb[int(0.95 * (len(shb) - 1))]:.1f}, "
                         f"maximum {shb[-1]:.1f} ({len(shb)} shuffles)")
    core_bg = (sum(len({p for p, *_ in core_patterns(v[:a.near]) if p > 0}) for v in shuf.values()) /
               (len(shuf) * a.near / 1000)) if shuf else 0
    # 4c, 4d
    ns = screen_args(a)
    with concurrent.futures.ProcessPoolExecutor(a.cpus) as ex:
        sls = dict(zip(seqs, ex.map(stem_loops, [(g, s[:a.near].replace("T", "U"), ns) for g, s in seqs.items()])))
    with open(a.stemloops, "w") as o:
        cols = ["sequence", "start", "end", "length", "pairs", "mfe", "z", "quartet", "helix_2_pairs", "AAR", "CC",
                "structure", "rna"]
        o.write("\t".join(cols) + "\n")
        for g, cs in sls.items():
            for c in cs: o.write("\t".join(str(c.get(x, "")) for x in cols) + "\n")
    for r in rows:
        g = r["genome"]
        if g not in seqs: continue
        r["searched_nt"] = len(seqs[g])
        for m in models:
            r[f"best_{m}"] = fmt(best_in(hits.get(g, {}).get(m, []), 1, len(seqs[g])))
            r[f"best_{m}_first_{a.near}"] = fmt(best_in(hits.get(g, {}).get(m, []), 1, a.near))
        cs = sls[g]
        r[f"stemloops_first_{a.near}"] = len(cs)
        r["stemloops_with_features"] = "; ".join(f"{c['start']}-{c['end']} (z {c['z']})" for c in cs
                                                 if with_features(c)) or "none"
        core = sorted({p for p, *_ in core_patterns(seqs[g][:a.near]) if p > 0})
        r[f"core_patterns_first_{a.near}"] = len(core)
        if g in APPENDICULARIANS:
            for m in models: r[f"background_{m}"] = bg.get(m, "")
            r["background_core_patterns_per_kb"] = f"{core_bg:.2f}"
    # 4e. comparative
    comp = {g: seqs[g][:a.compare] for g in APPENDICULARIANS if g in seqs}
    if len(comp) >= 2:
        with open(f"{work}/comp.fa", "w") as o:
            for g, s in comp.items(): o.write(f">{g}\n{s}\n")
        aln = dict(read_fasta_text(run(["mafft", "--quiet", "--localpair", "--maxiterate", "1000", f"{work}/comp.fa"])))
        names = list(aln)
        ss, e = RNA.alifold([aln[n].upper().replace("T", "U") for n in names])
        os.makedirs(os.path.dirname(a.alignment) or ".", exist_ok=True)
        with open(a.alignment, "w") as o:
            o.write("# STOCKHOLM 1.0\n")
            o.write(f"#=GF CC first {a.compare} nt from the candidate TGA of SELENOI; MAFFT L-INS-i; RNAalifold "
                    f"consensus, {e:.2f} kcal/mol\n")
            for n in names: o.write(f"{n:12s} {aln[n].upper()}\n")
            o.write(f"{'#=GC SS_cons':12s} {ss}\n//\n")
        idn = []
        for i, x in enumerate(names):
            for y in names[i + 1:]:
                m = sum(1 for p, q in zip(aln[x], aln[y]) if p == q and p != "-")
                n = sum(1 for p, q in zip(aln[x], aln[y]) if p != "-" and q != "-")
                idn.append(f"{x}-{y} {100 * m / n:.0f}%")
        pairs = ss.count("(")
        for r in rows:
            if r["genome"] in comp:
                r["alifold_consensus"] = f"{pairs} base pairs, {e:.2f} kcal/mol ({len(comp)} species)"
        print("pairwise identity:", ", ".join(idn))
    cols = ["genome", "source", "region", "peptide_hit", "codon", "context", "codons_to_stop", "searched_nt"] + \
           [f"best_{m}" for m in models] + [f"best_{m}_first_{a.near}" for m in models] + \
           [f"stemloops_first_{a.near}", "stemloops_with_features", f"core_patterns_first_{a.near}",
            ] + [f"background_{m}" for m in models] + ["background_core_patterns_per_kb", "alifold_consensus"]
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r.get(x, "")) for x in cols) + "\n")
    print(open(a.table).read())
    # 5. calibration
    if not a.no_calibration:
        for r in calibration(a, models, work): print("\t".join(map(str, r)))
    shutil.rmtree(work, ignore_errors=True)


def read_fasta_text(text):
    name, buf = None, []
    for l in text.splitlines():
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].strip(), []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


if __name__ == "__main__":
    main()
