#!/usr/bin/env python3
"""Content of introns: secondary structure, pairing of their two ends, and copies elsewhere in the genome.

Introns: those of every gold-standard model (coding exons; --gold-dirs, genomes --genomes), read on the genome
sequences of the samplesheet (gene strand). For each intron:
  - length, GC fraction, splice sites as named by the model, and the bases repeated at both ends (number of
    placements that leave the spliced sequence unchanged, minus one);
  - mfe (introns up to --max-fold nt): minimum free energy of the intron (RNAfold, ViennaRNA, default parameters), and its z-score against
    --shuffles dinucleotide-preserving shuffles of the intron (Altschul and Erickson 1985; random generator seeded
    with --seed, the genome, the gene and the intron number, so that the z-scores of an intron do not change when other
    models are added): z < -2 means more
    structure than expected from the composition; introns that hold an assembly gap (a run of 10 N or more; single N
    are ambiguous bases, frequent in the O. albicans assembly) are not folded, and their ends, when outside the
    gap, are compared with shuffles of the intron without its N;
  - ends: energy of the duplex between the first and the last --end nt of the intron (RNAduplex), a pairing that
    would bring the two splice sites together, and its z-score against the same shuffles;
  - copies: BLASTN (legacy BLAST 2.2.26, -e --evalue, no filter) of the intron against its genome: the number of hits of
    at least --min-hit nt outside the intron itself, the number of loci (hits closer than 1 kb merged), the best such
    hit (identity, length) and the fraction of the intron covered by such hits; a repeated element (transposon,
    satellite) inside an intron gives many loci.
Output (--out, TSV): genome, gene, intron (number in gene order), contig, start, end (1-based, genome), strand, length,
gc, sites, repeated_bases, mfe, mfe_shuffled_mean, mfe_z, ends_energy, ends_z, hits, loci, best_hit, covered, focus
("yes" for the genes of --focus).
"""
import argparse, collections, csv, gzip, os, random, re, subprocess, tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s): return s.translate(COMP)[::-1]


def read_contigs(path, names):
    out, name, buf = {}, None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name in names: out[name] = "".join(buf).upper()
                name, buf = l[1:].split()[0], []
            elif name in names: buf.append(l.strip())
    if name in names: out[name] = "".join(buf).upper()
    return out


def models(path):
    m = collections.OrderedDict()
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) == 9 and f[2] == "CDS":
            p = re.search(r"Parent=rna-([^;]+)", f[8]).group(1)
            m.setdefault(p, []).append((f[0], int(f[3]), int(f[4]), f[6]))
    return m


def dinuc_shuffle(s, rng):
    """Dinucleotide-preserving shuffle (Altschul and Erickson 1985): random Eulerian path with the same last edges."""
    if len(s) < 3: return s
    edges = collections.defaultdict(list)
    for a, b in zip(s, s[1:]): edges[a].append(b)
    last = s[-1]
    while True:
        # a random last edge out of every vertex except the last one, forming a tree towards the last vertex
        lastedge = {}
        for v in edges:
            if v != last: lastedge[v] = rng.choice(edges[v])
        ok = True
        for v in lastedge:            # every vertex must reach the last one through the last edges
            seen, u = set(), v
            while u != last:
                if u in seen or u not in lastedge: ok = False; break
                seen.add(u); u = lastedge[u]
            if not ok: break
        if ok: break
    lists = {}
    for v, outs in edges.items():
        rest = list(outs)
        if v in lastedge: rest.remove(lastedge[v])
        rng.shuffle(rest)
        lists[v] = rest + ([lastedge[v]] if v in lastedge else [])
    out, u = [s[0]], s[0]
    for _ in range(len(s) - 1):
        u = lists[u].pop(0); out.append(u)
    return "".join(out)


def rnafold(seqs):
    p = subprocess.run(["RNAfold", "--noPS"], input="".join(f">s{i}\n{x}\n" for i, x in enumerate(seqs)),
                       capture_output=True, text=True, check=True).stdout
    return [float(m) for m in re.findall(r"\(\s*(-?\d+\.\d+)\)\s*$", p, flags=re.M)]


def rnaduplex(pairs):
    p = subprocess.run(["RNAduplex"], input="".join(f"{a}\n{b}\n" for a, b in pairs), capture_output=True, text=True,
                       check=True).stdout
    return [float(m) for m in re.findall(r"\(\s*(-?\d+\.\d+)\)\s*$", p, flags=re.M)]


def zscore(x, ctrl):
    if not ctrl: return "", ""
    m = sum(ctrl) / len(ctrl)
    sd = (sum((c - m) ** 2 for c in ctrl) / max(len(ctrl) - 1, 1)) ** 0.5
    return round(m, 2), (round((x - m) / sd, 2) if sd > 0 else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gold-dirs", nargs="+", default=["curation/gold_standard", "curation/gold_standard_machinery"])
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3", "Oalb"])
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--focus", nargs="*", default=["PSTK"])
    ap.add_argument("--shuffles", type=int, default=100)
    ap.add_argument("--end", type=int, default=25)
    ap.add_argument("--max-fold", type=int, default=2500, help="introns longer than this are not folded")
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--min-hit", type=int, default=40)
    ap.add_argument("--work", default="work/intron_content")
    ap.add_argument("--out", default="results/tables/intron_content.tsv")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    os.makedirs(a.work, exist_ok=True)
    rows = []
    for g in a.genomes:
        ms = collections.OrderedDict()
        for d in a.gold_dirs:
            if os.path.exists(f"{d}/{g}.gff3"):
                for k, v in models(f"{d}/{g}.gff3").items(): ms.setdefault(k, v)
        if not ms: continue
        seq = read_contigs(sheet[g]["fasta"], {v[0][0] for v in ms.values()})
        mine = []
        for name, cds in ms.items():
            contig, strand = cds[0][0], cds[0][3]
            gene = name[len(g) + 1:] if name.startswith(g + "_") else name
            ex = sorted((x, y) for _, x, y, _ in cds)
            introns = [(b1 + 1, a2 - 1) for (a1, b1), (a2, b2) in zip(ex, ex[1:]) if a2 > b1 + 1]
            if strand == "-": introns = introns[::-1]
            for i, (lo, hi) in enumerate(introns):
                w = 30
                s = seq[contig][lo - 1 - w:hi + w]
                if strand == "-": s = rc(s)
                A, B = w, len(s) - w
                kmax = 0
                while kmax < w and s[A + kmax] == s[B + kmax]: kmax += 1
                kmin = 0
                while kmin < w and s[A - 1 - kmin] == s[B - 1 - kmin]: kmin += 1
                intr = s[A:B]
                mine.append(dict(genome=g, gene=gene, intron=i + 1, contig=contig, start=lo, end=hi, strand=strand,
                                 length=len(intr), gc=round((intr.count("G") + intr.count("C")) / len(intr), 3),
                                 sites=intr[:2] + "-" + intr[-2:], repeated_bases=kmax + kmin, seq=intr,
                                 focus="yes" if gene in a.focus else ""))
        # structure
        for r in mine:
            if r["length"] > a.max_fold: continue        # folding of long introns: too slow for the shuffles
            s = r["seq"].replace("T", "U")
            n = min(a.end, len(s) // 2)
            # an assembly gap (a run of 10 N or more): the fold of the intron is unknown; the pairing of its ends is
            # scored when they lie outside the gap, against shuffles of the intron without its N (shuffled N would
            # fall into the ends of the controls and weaken their pairing)
            gap = "N" * 10 in s
            if gap and ("N" in s[:n] or "N" in s[-n:]): continue
            # shuffles seeded per intron, so that its z-scores do not depend on the other models scored
            rng = random.Random(f"{a.seed}:{g}:{r['gene']}:{r['intron']}")
            sh = [dinuc_shuffle(s.replace("N", "") if gap else s, rng) for _ in range(a.shuffles)]
            if not gap:
                e = rnafold([s] + sh)
                r["mfe"] = e[0]; r["mfe_shuffled_mean"], r["mfe_z"] = zscore(e[0], e[1:])
            d = rnaduplex([(s[:n], s[-n:])] + [(x[:n], x[-n:]) for x in sh])
            r["ends_energy"] = d[0]; _, r["ends_z"] = zscore(d[0], d[1:])
        # copies in the genome
        fa = f"{a.work}/{g}.genome.fa"
        if not os.path.exists(fa + ".nsq") and not os.path.exists(fa + ".00.nsq"):
            with (gzip.open(sheet[g]["fasta"], "rt") if sheet[g]["fasta"].endswith(".gz") else open(sheet[g]["fasta"])) as fh, open(fa, "w") as o:
                for l in fh: o.write(l.split()[0] + "\n" if l.startswith(">") else l)
            subprocess.run(["formatdb", "-i", os.path.basename(fa), "-p", "F", "-o", "F"], check=True, cwd=a.work)
        q = f"{a.work}/{g}.introns.fa"
        with open(q, "w") as o:
            for k, r in enumerate(mine): o.write(f">{k}\n{r['seq']}\n")
        out = subprocess.run(["blastall", "-p", "blastn", "-d", fa, "-i", q, "-e", a.evalue, "-F", "F", "-m", "8",
                              "-b", "5000", "-v", "5000"], capture_output=True, text=True, check=True).stdout
        hits = collections.defaultdict(list)
        for l in out.splitlines():
            f = l.split("\t")
            k, t, pid, alen, qs, qe, ts, te = int(f[0]), f[1], float(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8]), int(f[9])
            r = mine[k]
            lo, hi = min(ts, te), max(ts, te)
            if t == r["contig"] and lo <= r["end"] + 50 and hi >= r["start"] - 50: continue      # the intron itself
            if alen < a.min_hit: continue
            hits[k].append((t, lo, hi, pid, alen, qs, qe))
        for k, r in enumerate(mine):
            h = sorted(hits[k], key=lambda x: (x[0], x[1]))
            loci, prev = 0, None
            for t, lo, hi, *_ in h:
                if prev is None or t != prev[0] or lo > prev[1] + 1000: loci += 1; prev = [t, hi]
                else: prev[1] = max(prev[1], hi)
            cov = set()
            for *_, qs, qe in h: cov.update(range(min(qs, qe), max(qs, qe) + 1))
            best = max(h, key=lambda x: (x[4] * x[3]), default=None)
            r.update(hits=len(h), loci=loci, best_hit=f"{best[3]:.0f}% over {best[4]} nt ({best[0]}:{best[1]}-{best[2]})" if best else "",
                     covered=round(len(cov) / r["length"], 2))
        rows += mine
    cols = ["genome", "gene", "intron", "contig", "start", "end", "strand", "length", "gc", "sites", "repeated_bases", "mfe",
            "mfe_shuffled_mean", "mfe_z", "ends_energy", "ends_z", "hits", "loci", "best_hit", "covered", "focus"]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"{len(rows)} introns written to {a.out}")


if __name__ == "__main__":
    main()
