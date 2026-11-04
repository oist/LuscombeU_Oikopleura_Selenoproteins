#!/usr/bin/env python3
"""An intron of a gold-standard model that is an insertion: its ends, terminal repeats, reading frames and copies.

The intron (--intron contig:start-end:strand, 1-based, gene strand) is read with --flank nucleotides on each side from
the genome of --genome (--samplesheet files: the gold-standard assemblies and those of the other individuals). The same locus without the insertion is looked for in --empty (FASTA files of
the region of the gene in other assemblies, as written by curation/gene_regions.py): the record where both flanks
align (legacy BLASTN, same strand, the two alignments at most --max-gap nt apart). From the last base of the left flank
and the first base of the right flank, exact matches are extended into the intron on the genome and into the empty
site; the bases of the empty site matched from both sides are the target-site duplication, and the bases matched by
neither side are the insertion.
The insertion is then described by:
  - the exact terminal inverted repeat (bases of its start equal to the reverse complement of its end) and the
    longest alignment of its first part to its last part on the same strand starting and ending at its ends (long
    terminal repeats, LTR) or on the other strand (terminal inverted repeats);
  - its open reading frames of at least --min-orf codons (ATG to stop, both strands) and two motifs in them: YXDD
    (the catalytic motif of reverse transcriptases) and H-X3-7-H-X20-40-C-X2-C (an HHCC zinc-binding motif, as at the
    N-terminus of integrases);
  - its copies (BLASTN, E at most --evalue, no filter) in the genomes of --copies-in (databases shared with
    curation/intron_content.py in --work) and in the sequences of --copies-fasta, the insertion itself left out: hits
    of the same sequence and strand closer than --merge nt form one copy, with the fraction of the insertion covered,
    the identity (weighted by alignment length), the parts of the insertion aligned (with LTRs, a part at one end can
    be either LTR), and whether the copy reaches the end of its sequence (partial copies at contig ends).
Outputs: --out (TSV, field and value) and --copies (TSV, one row per copy).
"""
import argparse, collections, csv, gzip, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from intron_content import read_contigs, rc   # noqa: E402

CODE = {a + b + c: aa for (a, b, c), aa in zip(
    [(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
MOTIFS = {"YXDD": r"Y[A-Z]DD", "HHCC": r"H[A-Z]{3,7}H[A-Z]{20,40}C[A-Z]{2}C"}


def read_fasta(path):
    out, name, buf = {}, None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name: out[name] = "".join(buf).upper()
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
    if name: out[name] = "".join(buf).upper()
    return out


def blast(query_seqs, db_fasta, work, evalue, max_hits=5000):
    """Legacy BLASTN of {name: sequence} against a FASTA file (formatted in place when needed); tabular rows."""
    if not os.path.exists(db_fasta + ".nsq") and not os.path.exists(db_fasta + ".00.nsq"):
        subprocess.run(["formatdb", "-i", os.path.basename(db_fasta), "-p", "F", "-o", "F"], check=True,
                       cwd=os.path.dirname(db_fasta) or ".", stdout=subprocess.DEVNULL)
    q = f"{work}/insertion_query.fa"
    with open(q, "w") as o:
        for k, s in query_seqs.items(): o.write(f">{k}\n{s}\n")
    out = subprocess.run(["blastall", "-p", "blastn", "-d", db_fasta, "-i", q, "-e", str(evalue), "-F", "F", "-m", "8",
                          "-b", str(max_hits), "-v", str(max_hits)], capture_output=True, text=True, check=True).stdout
    rows = []
    for l in out.splitlines():
        f = l.split("\t")
        rows.append(dict(q=f[0], s=f[1], pid=float(f[2]), alen=int(f[3]), qs=int(f[6]), qe=int(f[7]), ss=int(f[8]),
                         se=int(f[9]), ev=float(f[10])))
    return rows


def plain_fasta(seqs, path):
    with open(path, "w") as o:
        for k, s in seqs.items(): o.write(f">{k}\n{s}\n")
    return path


def empty_site(locus, F, L, empties, a):
    """Record, orientation and junction of the empty site, with the insertion and TSD as locus indices."""
    seqs = {}
    for p in empties:
        for k, s in read_fasta(p).items(): seqs[f"{os.path.basename(p)}:{k}"] = s
    if not seqs: return None
    names = {f"e{i}": k for i, k in enumerate(seqs)}
    db = plain_fasta({f"e{i}": seqs[k] for i, k in enumerate(seqs)}, f"{a.work}/insertion_empty.fa")
    for ext in (".nhr", ".nin", ".nsq"):
        if os.path.exists(db + ext): os.remove(db + ext)
    hits = blast({"left": locus[:F], "right": locus[F + L:]}, db, a.work, "1e-10")
    best = None
    for e in names:
        for hl in (h for h in hits if h["s"] == e and h["q"] == "left"):
            for hr in (h for h in hits if h["s"] == e and h["q"] == "right"):
                plus_l, plus_r = hl["ss"] < hl["se"], hr["ss"] < hr["se"]
                if plus_l != plus_r: continue
                E = seqs[names[e]] if plus_l else rc(seqs[names[e]])
                n = len(E)
                tr = (lambda x: x - 1) if plus_l else (lambda x: n - x)      # 1-based subject -> 0-based index of E
                eL = tr(hl["se"]) + (F - hl["qe"])                         # index in E of the last left-flank base
                eR = tr(hr["ss"]) - (hr["qs"] - 1)                         # index in E of the first right-flank base
                if not (0 <= eL < n and 0 <= eR < n) or abs(eR - eL - 1) > a.max_gap: continue
                score = hl["alen"] + hr["alen"]
                if best is None or score > best[0]: best = (score, names[e], plus_l, E, eL, eR)
    if best is None: return None
    _, name, plus, E, eL, eR = best
    i, j = F, eL + 1                       # extend the left match into the intron
    while i < len(locus) and j < len(E) and locus[i] == E[j]: i, j = i + 1, j + 1
    left_end, left_end_E = i, j
    i, j = F + L - 1, eR - 1               # extend the right match backwards into the intron
    while i >= 0 and j >= 0 and locus[i] == E[j]: i, j = i - 1, j - 1
    right_start, right_start_E = i + 1, j + 1
    tsd = E[right_start_E:left_end_E] if right_start_E < left_end_E else ""
    return dict(record=name, orientation="+" if plus else "-", ins_start=left_end, ins_end=right_start, tsd=tsd)


def terminal_repeats(ins, a):
    """Exact terminal inverted repeat; best end-to-end alignment of the first and last parts (LTR or TIR)."""
    k = 0
    r = rc(ins)
    while k < len(ins) // 2 and ins[k] == r[k]: k += 1
    part = min(3000, len(ins) // 2)
    db = plain_fasta({"last": ins[-part:]}, f"{a.work}/insertion_last.fa")
    for ext in (".nhr", ".nin", ".nsq"):
        if os.path.exists(db + ext): os.remove(db + ext)
    best = None
    for h in blast({"first": ins[:part]}, db, a.work, "1e-5"):
        if h["qs"] > 3: continue
        if h["ss"] < h["se"] and h["se"] >= part - 3: kind = "LTR"          # same strand, to the end
        elif h["ss"] > h["se"] and h["ss"] >= part - 3: kind = "TIR"        # other strand, from the end
        else: continue
        if best is None or h["alen"] > best[1]: best = (kind, h["alen"], h["pid"], h["qe"])
    return k, best


def orfs(ins, min_codons):
    out = []
    for strand, s in (("+", ins), ("-", rc(ins))):
        for f in range(3):
            p = "".join(CODE.get(s[x:x + 3], "X") for x in range(f, len(s) - 2, 3))
            pos = 0
            for seg in p.split("*"):
                m = seg.find("M")
                if m >= 0 and len(seg) - m >= min_codons:
                    st = f + 3 * (pos + m)
                    prot = seg[m:]
                    a_, b_ = (st + 1, st + 3 * len(prot)) if strand == "+" else (len(s) - st - 3 * len(prot) + 1, len(s) - st)
                    out.append(dict(codons=len(prot), strand=strand, start=a_, end=b_,
                                    motifs=";".join(f"{k}@{x.start() + 1}" for k, rx in MOTIFS.items()
                                                    for x in re.finditer(rx, prot))))
                pos += len(seg) + 1
    return sorted(out, key=lambda o: -o["codons"])


def copies(hits, Lins, lengths, self_loc, merge):
    """Copies from BLAST rows: merged per subject and strand; self_loc (subject, lo, hi) left out."""
    by = collections.defaultdict(list)
    for h in hits:
        lo, hi = min(h["ss"], h["se"]), max(h["ss"], h["se"])
        if self_loc and h["s"] == self_loc[0] and lo <= self_loc[2] and hi >= self_loc[1]: continue
        by[(h["s"], h["ss"] < h["se"])].append((lo, hi, h))
    out = []
    for (s, plus), hs in by.items():
        hs.sort(key=lambda x: x[0])
        groups = [[hs[0]]]
        for x in hs[1:]:
            if x[0] <= max(y[1] for y in groups[-1]) + merge: groups[-1].append(x)
            else: groups.append([x])
        for g in groups:
            merged = []
            for q1, q2 in sorted((x[2]["qs"], x[2]["qe"]) for x in g):
                if merged and q1 <= merged[-1][1] + 1: merged[-1][1] = max(merged[-1][1], q2)
                else: merged.append([q1, q2])
            cov = sum(q2 - q1 + 1 for q1, q2 in merged)
            lo, hi = min(x[0] for x in g), max(x[1] for x in g)
            n = lengths.get(s)
            out.append(dict(sequence=s, start=lo, end=hi, strand="+" if plus else "-", covered=round(cov / Lins, 3),
                            identity=round(sum(x[2]["pid"] * x[2]["alen"] for x in g) / sum(x[2]["alen"] for x in g), 2),
                            parts=",".join(f"{q1}-{q2}" for q1, q2 in merged),
                            sequence_end="yes" if n and (lo <= 100 or hi >= n - 100) else "no"))
    return sorted(out, key=lambda c: -c["covered"])


def lengths_of(path):
    out, name = {}, None
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"): name = l[1:].split()[0]; out[name] = 0
            else: out[name] += len(l.strip())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--intron", required=True, help="contig:start-end:strand (1-based, gene strand)")
    ap.add_argument("--empty", nargs="+", default=[], help="FASTA files of the locus in other assemblies")
    ap.add_argument("--copies-in", nargs="*", default=[], help="genomes (samplesheet names) searched for copies")
    ap.add_argument("--copies-fasta", nargs="*", default=[], help="other sequences searched for copies")
    ap.add_argument("--samplesheet", nargs="+", default=["stage1/samplesheet.csv", "stage2/run/samplesheet_final.csv"],
                    help="samplesheets (sample and fasta columns; missing files are skipped)")
    ap.add_argument("--flank", type=int, default=300)
    ap.add_argument("--max-gap", type=int, default=20)
    ap.add_argument("--min-orf", type=int, default=150)
    ap.add_argument("--evalue", default="1e-20")
    ap.add_argument("--merge", type=int, default=1000)
    ap.add_argument("--work", default="work/intron_content")
    ap.add_argument("--out", required=True)
    ap.add_argument("--copies", required=True)
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    sheet = {}
    for p in a.samplesheet:
        if os.path.exists(p):
            for r in csv.DictReader(open(p)): sheet.setdefault(r["sample"], r)
    m = re.fullmatch(r"(.+):(\d+)-(\d+):([+-])", a.intron)
    contig, start, end, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
    seq = read_contigs(sheet[a.genome]["fasta"], {contig})[contig]
    F, L = a.flank, end - start + 1
    locus = seq[start - 1 - F:end + F]
    if strand == "-": locus = rc(locus)
    to_genome = (lambda i: start - F + i) if strand == "+" else (lambda i: end + F - i)   # locus index -> 1-based
    rows = [("genome", a.genome), ("intron", a.intron), ("intron_length", L)]
    site = empty_site(locus, F, L, a.empty, a)
    if site:
        s_, e_ = site["ins_start"], site["ins_end"]
        lo, hi = sorted((to_genome(s_), to_genome(e_ - 1)))
        rows += [("empty_site", f"{site['record']} ({site['orientation']})"), ("insertion", f"{contig}:{lo}-{hi}"),
                 ("insertion_length", e_ - s_), ("target_site_duplication", site["tsd"] or "none"),
                 ("tsd_length", len(site["tsd"]))]
        ins = locus[s_:e_]
    else:
        rows += [("empty_site", "not found: the whole intron is described")]
        ins = locus[F:F + L]
    k, rep = terminal_repeats(ins, a)
    rows += [("insertion_start", ins[:20]), ("insertion_end", ins[-20:]),
             ("exact_terminal_inverted_repeat", f"{k} ({ins[:k]})" if k else "0")]
    if rep:
        kind, alen, pid, qe = rep
        rows += [("terminal_repeat", kind), ("terminal_repeat_length", qe), ("terminal_repeat_identity", pid)]
    else:
        rows += [("terminal_repeat", "none")]
    of = orfs(ins, a.min_orf)
    rows += [("orfs", len(of))] + [(f"orf{i + 1}", f"{o['codons']} codons, strand {o['strand']}, {o['start']}-{o['end']}"
                                    + (f", {o['motifs']}" if o["motifs"] else "")) for i, o in enumerate(of)]
    crow = []
    for g in a.copies_in:
        fa = f"{a.work}/{g}.genome.fa"
        if not os.path.exists(fa):
            with (gzip.open(sheet[g]["fasta"], "rt") if sheet[g]["fasta"].endswith(".gz") else open(sheet[g]["fasta"])) as fh, open(fa, "w") as o:
                for l in fh: o.write(l.split()[0] + "\n" if l.startswith(">") else l)
        hits = blast({"insertion": ins}, fa, a.work, a.evalue)
        self_loc = (contig, start - 50, end + 50) if g == a.genome else None
        cs = copies(hits, len(ins), lengths_of(fa), self_loc, a.merge)
        crow += [dict(source=g, **c) for c in cs]
        rows.append((f"copies_{g}", f"{len(cs)} ({sum(c['covered'] >= 0.9 for c in cs)} covering at least 90%)"))
    for p in a.copies_fasta:
        fa = plain_fasta(read_fasta(p), f"{a.work}/insertion_{os.path.basename(p)}")
        for ext in (".nhr", ".nin", ".nsq"):
            if os.path.exists(fa + ext): os.remove(fa + ext)
        cs = copies(blast({"insertion": ins}, fa, a.work, a.evalue), len(ins), lengths_of(fa), None, a.merge)
        crow += [dict(source=p, **c) for c in cs]
        rows.append((f"copies_{p}", f"{len(cs)} ({sum(c['covered'] >= 0.9 for c in cs)} covering at least 90%)"))
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as o:
        o.write("field\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in rows))
    cols = ["source", "sequence", "start", "end", "strand", "covered", "identity", "parts", "sequence_end"]
    with open(a.copies, "w") as o:
        o.write("\t".join(cols) + "\n" + "".join("\t".join(str(c[x]) for x in cols) + "\n" for c in crow))
    print(open(a.out).read(), end="")


if __name__ == "__main__":
    main()
