#!/usr/bin/env python3
"""Equivalent placements of introns of a gold-standard gene, for the side panel of figures/reannotation_page.py.

An intron whose ends repeat the same bases (gene strand: the first bases of the intron equal the first bases of the
next exon, or the last bases of the previous exon equal the last bases of the intron) can be moved by those bases
without changing the spliced sequence: the transcripts cannot tell these placements apart, and each placement reads
other splice-site dinucleotides. curation/gold_standard.py names one of them (an intron class of O. dioica when one
exists, otherwise the smallest move from the alignment). For each intron given, this script lists all equivalent
placements (shift k in nucleotides, donor-acceptor dinucleotides, phase of the intron in the coding sequence).

Assemblies (--assembly NAME=gold or NAME=lift:SOURCE): gold, the exons of the gold-standard model
(<gold-dir>/<NAME>.gff3); lift, the exons of SOURCE carried base by base onto NAME through a minimap2 alignment
(-x asm20, cs tag) of the two regions. Sequences: results/gene_regions/<gene>/<NAME>.fa (make pstk-regions), the first
record, in gene orientation as its header says.
Output (--out, TSV): assembly, source, intron, length, window (gene strand: --flank bases of exon, the repeated bases
in brackets, --flank bases of intron, " ... ", --flank bases of intron, the repeated bases in brackets, --flank bases
of exon), placements (k:donor-acceptor:phase, ";"-separated), named_k (the placement of the gold-standard model,
0), protein (translation of the exons, for the check of lifted models) and note.
"""
import argparse, os, re, subprocess, tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODE = {a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[16 * i + 4 * j + k]
        for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}


def rc(s): return s.translate(COMP)[::-1]


def translate(s): return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def region(path):
    head, seq = None, []
    for l in open(path):
        if l.startswith(">"):
            if head: break
            head = l[1:].strip()
        else: seq.append(l.strip())
    m = re.match(r"[^|]+\|([^:]+):(\d+)-(\d+)\|([+\-])", head)
    seq = "".join(seq).upper()
    r = dict(contig=m.group(1), start=int(m.group(2)), end=int(m.group(3)), strand=m.group(4))
    r["gs"] = seq if r["strand"] == "+" else rc(seq)
    return r


def gold_exons(path, gene, r):
    """Coding exons of the gene as gene-strand intervals [a, b) of the region, in gene order."""
    rows = []
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) == 9 and f[2] == "CDS" and re.search(rf"Parent=rna-[^;]*_{re.escape(gene)}(;|$)", f[8]):
            rows.append((int(f[3]), int(f[4])))
    to = (lambda p: p - r["start"]) if r["strand"] == "+" else (lambda p: r["end"] - p)
    return sorted((min(to(a), to(b)), max(to(a), to(b)) + 1) for a, b in rows)


def lift_exons(src, tgt, exons):
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/q.fa", "w").write(f">q\n{src}\n"); open(f"{d}/t.fa", "w").write(f">t\n{tgt}\n")
        out = subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs=long", f"{d}/t.fa", f"{d}/q.fa"],
                             capture_output=True, text=True, check=True).stdout
    best = max((l.split("\t") for l in out.splitlines()), key=lambda f: int(f[9]))
    if best[4] != "+": raise SystemExit("lifted region not in gene orientation")
    qi, ti, m = int(best[2]), int(best[7]), {}
    cs = next(x[5:] for x in best[12:] if x.startswith("cs:Z:"))
    for op, s in re.findall(r"([=*+\-~])([A-Za-z0-9]+)", cs):
        if op == "=":
            for k in range(len(s)): m[qi + k] = ti + k
            qi += len(s); ti += len(s)
        elif op == "*": m[qi] = ti; qi += 1; ti += 1
        elif op == "+": qi += len(s)
        elif op == "-": ti += len(s)
    return [(m[a], m[b - 1] + 1) for a, b in exons]


def placements(gs, a, b, maxshift=20):
    out = []
    for k in range(-maxshift, maxshift + 1):
        if k > 0 and gs[a:a + k] != gs[b:b + k]: continue
        if k < 0 and gs[a + k:a] != gs[b + k:b]: continue
        out.append((k, gs[a + k:a + k + 2] + "-" + gs[b + k - 2:b + k]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--assembly", nargs="+", required=True)
    ap.add_argument("--introns", nargs="+", type=int, required=True, help="intron numbers in gene order (1-based)")
    ap.add_argument("--gold-dir", default="curation/gold_standard_machinery")
    ap.add_argument("--regions", default="results/gene_regions")
    ap.add_argument("--flank", type=int, default=7)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    exons, seqs, rows = {}, {}, []
    for spec in a.assembly:
        name, how = spec.split("=")
        r = region(f"{a.regions}/{a.gene}/{name}.fa")
        seqs[name] = r["gs"]
        if how == "gold":
            exons[name], source = gold_exons(f"{a.gold_dir}/{name}.gff3", a.gene, r), "gold standard"
        else:
            src = how.split(":")[1]
            exons[name], source = lift_exons(seqs[src], r["gs"], exons[src]), f"exons carried from {src}"
        gs, ex = seqs[name], exons[name]
        prot = translate("".join(gs[x:y] for x, y in ex))
        cum = 0
        for i in range(len(ex) - 1):
            cum += ex[i][1] - ex[i][0]
            if i + 1 not in a.introns: continue
            s, e = ex[i][1], ex[i + 1][0]
            pl = placements(gs, s, e)
            k0, k1 = min(k for k, _ in pl), max(k for k, _ in pl)
            f = a.flank
            window = (gs[s + k0 - f:s + k0].lower() + "[" + gs[s + k0:s + k1] + "]" + gs[s + k1:s + k1 + f] + " ... " +
                      gs[e + k0 - f:e + k0] + "[" + gs[e + k0:e + k1] + "]" + gs[e + k1:e + k1 + f].lower())
            rows.append(dict(assembly=name, source=source, intron=i + 1, length=e - s, window=window,
                             placements=";".join(f"{k}:{d}:{(cum + k) % 3}" for k, d in pl), named_k=0,
                             protein=prot.rstrip("*"), note=""))
    cols = ["assembly", "source", "intron", "length", "window", "placements", "named_k", "protein", "note"]
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"{len(rows)} rows written to {a.out}")


if __name__ == "__main__":
    main()
