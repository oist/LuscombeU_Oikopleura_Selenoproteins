#!/usr/bin/env python3
"""TXNL1 homologs in the three O. dioica assemblies: reciprocal best hits with human TXNL1, the next hits of O. dioica
TXNL1 in its own assembly, and the hits of the two parts of human TXNL1 (thioredoxin domain, PITH region).

The O. dioica TXNL1 (gold standard, curation/gold_standard/<genome>.proteins.faa; 168 residues, selenocysteine read as
X) aligns with the PITH region of human TXNL1 only. For each assembly of --genomes:
1. forward: blastp of O. dioica TXNL1 against the human RefSeq proteome; hits grouped by gene (RefSeq annotation:
   CDS protein_id -> gene): the --top best genes, each with its best protein, and the margin (bits) of each gene to
   the next.
2. reverse: human TXNL1 (--human-txnl1), its thioredoxin part (residues 1 to --split - 1; Pfam PF00085 at 12-104 in
   UniProt O43396) and its PITH part (--split to the end; PF06201 at 118-279, then the C-terminal tail) against the
   assembly:
   - blastp against its proteins: the annotated coding sequences, translated from the genome (phase of the first
     segment; internal stop codons read as X), and the gold-standard proteins;
   - tblastn against its genome; a locus groups colinear HSPs on one strand within --locus-span nt (function loci),
     and loci are ranked by the sum of the bits of their HSPs.
   For each of the --top best subjects: whether it overlaps the TXNL1 gold-standard locus (same strand), the annotated
   transcripts it overlaps, and its best human gene, by blastp of the protein (proteins) or blastx of the locus with
   --flank nt on each side (genome); reciprocal when that gene is TXNL1.
3. self: O. dioica TXNL1 against the same proteins and genome. Subjects at the TXNL1 locus are marked; the others are
   the next hits (paralogs, other haplotypes of a haplotype-split assembly such as OdB3, or unannotated copies), with
   identity, query coverage and best human gene.
Legacy BLAST (blastall, formatdb), no low-complexity filter, E <= --evalue. Databases are built in --work and removed
at the end unless --keep.
Inputs: stage1/samplesheet.csv (genome FASTA and annotation), curation/gold_standard/<genome>.gff3 and .proteins.faa,
data/genomes/Hsap/protein.faa.gz and annotation.gff3.gz.
Outputs: results/tables/stage2_txnl1_homologs.tsv (every row) and stage2_txnl1_homologs_summary.tsv (answers per
assembly).
"""
import argparse, collections, csv, os, re, shutil, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uga_secis_screen import read_fasta, rc, opener, attributes   # noqa: E402

CODE = {a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[16 * i + 4 * j + k]
        for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}


def run(cmd, cwd=None):
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode: sys.exit(f"command failed ({p.returncode}) in {cwd or '.'}: {' '.join(cmd)}\n{p.stderr[-2000:]}")
    return p.stdout


def clean(seq):
    """Protein for BLAST: selenocysteine and internal stops as X, no terminal stop."""
    return seq.upper().rstrip("*").replace("U", "X").replace("*", "X")


def translate(s): return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


class Db:
    """A legacy BLAST database with numbered subjects (s1, s2, ...), to avoid the parsing of identifiers."""

    def __init__(self, work, name, records, protein):
        self.dir = os.path.join(work, name)
        os.makedirs(self.dir, exist_ok=True)
        self.path, self.protein = os.path.join(self.dir, "db.fa"), protein
        self.name, self.length = {}, {}
        with open(self.path, "w") as o:
            for n, (k, s) in enumerate(records, 1):
                sid = f"s{n}"
                self.name[sid], self.length[sid] = k, len(s)
                o.write(f">{sid}\n{s}\n")
        run(["formatdb", "-i", "db.fa", "-p", "T" if protein else "F", "-o", "F"], self.dir)


def blast(program, db, queries, work, tag, evalue, cpus, max_hits=50):
    """-m 8 rows of a search of queries (name -> sequence) against a Db; subject names translated back."""
    fa = os.path.join(work, f"{tag}.query.fa")
    names = {}
    with open(fa, "w") as o:
        for n, (k, s) in enumerate(queries.items(), 1):
            names[f"q{n}"] = k
            o.write(f">q{n}\n{s}\n")
    out = run(["blastall", "-p", program, "-d", db.path, "-i", fa, "-e", str(evalue), "-F", "F", "-m", "8",
               "-v", str(max_hits), "-b", str(max_hits), "-a", str(cpus)])
    rows = []
    for l in out.splitlines():
        f = l.split("\t")
        if len(f) < 12: continue
        rows.append({"q": names[f[0]], "s": db.name[f[1]], "sid": f[1], "ident": float(f[2]), "len": int(f[3]),
                     "qs": int(f[6]), "qe": int(f[7]), "ss": int(f[8]), "se": int(f[9]), "evalue": float(f[10]),
                     "bits": float(f[11])})
    return rows


def coverage(ranges, qlen):
    """Fraction of the query covered by the union of (start, end) ranges."""
    cov = set()
    for a, b in ranges: cov.update(range(min(a, b), max(a, b) + 1))
    return len(cov) / qlen if qlen else 0


def annotated_proteins(fasta, gff):
    """Annotated coding sequences translated from the genome: id -> {"contig", "strand", "lo", "hi" (0-based, half-open,
    span of the coding sequence), "protein", "gene"}; id = protein_id when present, else the Parent of the CDS rows."""
    cds = collections.defaultdict(list)
    meta = {}
    for l in opener(gff):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        at = attributes(f[8])
        k = at.get("protein_id") or at.get("Parent", at.get("ID", "?")).split(",")[0]
        cds[k].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] in "012" else 0))
        label = " ".join(x for x in (at.get("gene") or at.get("locus_tag") or at.get("Name", ""),
                                      at.get("product", "").replace("%2C", ",")) if x)
        meta.setdefault(k, {"contig": f[0], "strand": f[6], "gene": label})
    by_contig = collections.defaultdict(list)
    for k, m in meta.items(): by_contig[m["contig"]].append(k)
    out = {}
    for contig, seq in read_fasta(fasta):
        if contig not in by_contig: continue
        seq = seq.upper()
        for k in by_contig[contig]:
            m, segs = meta[k], sorted(cds[k])
            if m["strand"] == "+":
                nt, phase = "".join(seq[a:b] for a, b, _ in segs), segs[0][2]
            else:
                nt, phase = "".join(rc(seq[a:b]) for a, b, _ in segs[::-1]), segs[-1][2]
            out[k] = dict(m, lo=segs[0][0], hi=segs[-1][1], protein=clean(translate(nt[phase:])))
    return out


def gold_loci(gff):
    """Gold-standard models: name -> (contig, strand, lo, hi) of the mRNA (0-based, half-open)."""
    out = {}
    for l in open(gff):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "mRNA": continue
        at = attributes(f[8])
        out[at.get("Name", at.get("ID"))] = (f[0], f[6], int(f[3]) - 1, int(f[4]))
    return out


def overlaps(loc, contig, strand, lo, hi):
    return loc is not None and loc[0] == contig and loc[1] == strand and loc[2] < hi and lo < loc[3]


def loci(rows, span, max_overlap=0.5):
    """tblastn HSPs grouped into loci per query: an HSP joins a locus on the same sequence and strand when it lies
    within span nt of it, overlaps none of its HSPs on the query by more than max_overlap of the shorter one (HSPs
    extended across the short introns of O. dioica overlap on the query), and keeps their order along the query
    (midpoints; colinear); loci ranked by the sum of the bits of their HSPs (overlaps counted twice)."""
    out = collections.defaultdict(list)
    for r in sorted(rows, key=lambda r: -r["bits"]):
        strand = "+" if r["se"] >= r["ss"] else "-"
        lo, hi = min(r["ss"], r["se"]) - 1, max(r["ss"], r["se"])
        for L in out[r["q"]]:
            if L["contig"] != r["s"] or L["strand"] != strand or not (lo < L["hi"] + span and L["lo"] - span < hi):
                continue
            if any(min(r["qe"], x["qe"]) - max(r["qs"], x["qs"]) + 1 >
                   max_overlap * min(r["qe"] - r["qs"] + 1, x["qe"] - x["qs"] + 1) for x in L["hsps"]): continue
            qm, gm = r["qs"] + r["qe"], lo + hi
            if all((qm < x["qs"] + x["qe"]) == ((gm < x["gm"]) if strand == "+" else (gm > x["gm"])) for x in L["hsps"]):
                L["hsps"].append(dict(r, gm=gm))
                L["lo"], L["hi"] = min(L["lo"], lo), max(L["hi"], hi)
                L["bits"] += r["bits"]
                break
        else:
            out[r["q"]].append({"contig": r["s"], "strand": strand, "lo": lo, "hi": hi, "bits": r["bits"],
                                "hsps": [dict(r, gm=lo + hi)]})
    for q in out: out[q].sort(key=lambda L: -L["bits"])
    return out


def human_genes(rows):
    """Per query: [(gene, protein, bits, evalue)] by decreasing bits, one entry per gene."""
    out = collections.defaultdict(list)
    for r in sorted(rows, key=lambda r: -r["bits"]):
        if any(g == r["gene"] for g, *_ in out[r["q"]]): continue
        out[r["q"]].append((r["gene"], r["s"], r["bits"], r["evalue"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3"])
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--gene", default="TXNL1.b", help="gold-standard gene name")
    ap.add_argument("--human", default="data/genomes/Hsap")
    ap.add_argument("--human-gene", default="TXNL1")
    ap.add_argument("--human-txnl1", default="NP_004777.1")
    ap.add_argument("--split", type=int, default=118, help="first residue of the PITH part of the human protein")
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--locus-span", type=int, default=5000)
    ap.add_argument("--flank", type=int, default=100)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--work", default="work/stage2/txnl1_homologs")
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--out", default="results/tables/stage2_txnl1_homologs.tsv")
    ap.add_argument("--summary", default="results/tables/stage2_txnl1_homologs_summary.tsv")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}

    # human proteome, genes and descriptions
    gene_of = {}
    for l in opener(os.path.join(a.human, "annotation.gff3.gz")):
        if "\tCDS\t" not in l: continue
        m, g = re.search(r"protein_id=([^;\n]+)", l), re.search(r"gene=([^;\n]+)", l)
        if m and g: gene_of[m.group(1)] = g.group(1)
    human, desc = {}, {}
    with opener(os.path.join(a.human, "protein.faa.gz")) as fh:
        name = None
        for l in fh:
            if l.startswith(">"):
                name, _, d = l[1:].strip().partition(" ")
                desc[name], human[name] = re.sub(r" \[[^]]+\]$", "", d), []
            elif name: human[name].append(l.strip())
    human = {k: clean("".join(v)) for k, v in human.items()}
    if a.human_txnl1 not in human: sys.exit(f"{a.human_txnl1} not in the human proteome")
    hdb = Db(a.work, "human", human.items(), True)
    h = human[a.human_txnl1]
    hq = {f"human {a.human_gene} (1-{len(h)})": h,
          f"human {a.human_gene} thioredoxin part (1-{a.split - 1})": h[:a.split - 1],
          f"human {a.human_gene} PITH part ({a.split}-{len(h)})": h[a.split - 1:]}

    def to_human(seqs, tag, program="blastp"):
        """Best human genes of sequences (name -> sequence)."""
        if not seqs: return {}
        rows = blast(program, hdb, seqs, a.work, tag, a.evalue, a.cpus, 100)
        for r in rows: r["gene"] = gene_of.get(r["s"], r["s"])
        return human_genes(rows)

    def gene_cells(genes):
        b = genes[0] if genes else None
        n = genes[1] if len(genes) > 1 else None
        return [b[0] if b else "", f"{b[2]:.1f}" if b else "", n[0] if n else "", f"{n[2]:.1f}" if n else ""]

    cols = ["genome", "search", "query", "database", "rank", "subject", "location", f"at_{a.gene}_locus",
            "annotation", "evalue", "bits", "identity", "query_range", "query_coverage", "subject_length",
            "best_human_gene", "best_human_bits", "next_human_gene", "next_human_bits", "reciprocal"]
    table, summary = [], []
    for g in a.genomes:
        if g not in sheet: sys.exit(f"{g} not in {a.samplesheet}")
        gold = dict(read_fasta(os.path.join(a.gold, f"{g}.proteins.faa")))
        gname = next((k for k in gold if k.split("|")[0] == a.gene), None)
        if gname is None:
            print(f"{g}: no gold-standard {a.gene}", file=sys.stderr)
            continue
        oq = {f"{g} {a.gene}": clean(gold[gname])}
        gloc = gold_loci(os.path.join(a.gold, f"{g}.gff3"))
        tloc = gloc.get(a.gene)
        fasta, gff = sheet[g]["fasta"], sheet[g].get("annotation", "")
        ann = annotated_proteins(fasta, gff) if gff and os.path.exists(gff) else {}
        records, where = [], {}
        for k, r in ann.items():
            if r["protein"]:
                records.append((f"annotation:{k}", r["protein"]))
                where[f"annotation:{k}"] = (r["contig"], r["strand"], r["lo"], r["hi"], r["gene"])
        for k, s in gold.items():
            n = k.split("|")[0]
            records.append((f"gold:{n}", clean(s)))
            if n in gloc: where[f"gold:{n}"] = gloc[n] + (n,)
        pdb = Db(a.work, f"{g}.proteins", records, True)
        seq_of = dict(records)
        gseq = {k: s.upper() for k, s in read_fasta(fasta)}
        ndb = Db(a.work, f"{g}.genome", gseq.items(), False)

        def at_locus(contig, strand, lo, hi):
            return "yes" if overlaps(tloc, contig, strand, lo, hi) else "no"

        def ann_at(contig, strand, lo, hi):
            names = [k.split(":", 1)[1] + (f" ({w[4]})" if w[4] else "") for k, w in where.items()
                     if w[0] == contig and w[1] == strand and w[2] < hi and lo < w[3]]
            return ", ".join(sorted(names))

        # 1. forward
        fw = blast("blastp", hdb, oq, a.work, f"{g}.forward", a.evalue, a.cpus, 100)
        for r in fw: r["gene"] = gene_of.get(r["s"], r["s"])
        fgenes = human_genes(fw)[f"{g} {a.gene}"]
        best_fw = {}
        for r in fw:
            if r["gene"] not in best_fw: best_fw[r["gene"]] = r
        for i, (gene, prot, bits, ev) in enumerate(fgenes[:a.top], 1):
            r = best_fw[gene]
            nxt = fgenes[i] if i < len(fgenes) else None
            table.append([g, "forward", f"{g} {a.gene}", "human proteins", i, prot, "", "", f"{gene}: {desc.get(prot, '')}",
                          f"{ev:.1e}", f"{bits:.1f}", f"{r['ident']:.0f}", f"{r['qs']}-{r['qe']}",
                          f"{coverage([(x['qs'], x['qe']) for x in fw if x['s'] == prot], len(oq[f'{g} {a.gene}'])):.2f}",
                          len(human[prot]), gene, f"{bits:.1f}", nxt[0] if nxt else "", f"{nxt[2]:.1f}" if nxt else "", ""])
        fbest = fgenes[0][0] if fgenes else ""
        margin = f"{fgenes[0][2] - fgenes[1][2]:.1f}" if len(fgenes) > 1 else "no other gene"

        # 2-3. reverse and self, against proteins and genome
        queries = dict(hq, **oq)
        prow = blast("blastp", pdb, queries, a.work, f"{g}.proteins", a.evalue, a.cpus)
        nrow = blast("tblastn", ndb, queries, a.work, f"{g}.genome", a.evalue, a.cpus)
        # proteins: best HSP per subject, then the human genes of the top subjects
        tops = {}
        for q in queries:
            best = {}
            for r in sorted((r for r in prow if r["q"] == q), key=lambda r: -r["bits"]):
                if r["s"] not in best: best[r["s"]] = r
            tops[q] = list(best.values())[:a.top]
        need = {r["s"]: seq_of[r["s"]] for q in tops for r in tops[q]}
        hg = to_human(need, f"{g}.proteins_to_human")
        # genome: loci, then blastx of the top loci
        L = loci(nrow, a.locus_span)
        regions = {}
        for q in queries:
            for i, x in enumerate(L[q][:a.top], 1):
                lo, hi = max(0, x["lo"] - a.flank), min(len(gseq[x["contig"]]), x["hi"] + a.flank)
                s = gseq[x["contig"]][lo:hi]
                x["key"] = f"{x['contig']}:{lo + 1}-{hi}{x['strand']}"
                regions[x["key"]] = s if x["strand"] == "+" else rc(s)
        hx = to_human(regions, f"{g}.loci_to_human", "blastx")
        answers = {}
        for q in queries:
            kind = "self" if q in oq else "reverse"
            qlen = len(queries[q])
            for i, r in enumerate(tops[q], 1):
                w = where.get(r["s"])
                loc = f"{w[0]}:{w[2] + 1}-{w[3]}{w[1]}" if w else ""
                at = at_locus(*w[:4]) if w else ""
                genes = hg.get(r["s"], [])
                recip = ("yes" if genes and genes[0][0] == a.human_gene else "no") if kind == "reverse" else ""
                cov = coverage([(x["qs"], x["qe"]) for x in prow if x["q"] == q and x["s"] == r["s"]], qlen)
                table.append([g, kind, q, "proteins", i, r["s"], loc, at, w[4] if w else "", f"{r['evalue']:.1e}",
                              f"{r['bits']:.1f}", f"{r['ident']:.0f}", f"{r['qs']}-{r['qe']}", f"{cov:.2f}",
                              len(seq_of[r["s"]])] + gene_cells(genes) + [recip])
                answers.setdefault((q, "proteins"), []).append((r, at, genes, recip, w, cov))
            for i, x in enumerate(L[q][:a.top], 1):
                b = x["hsps"][0]
                genes = hx.get(x["key"], [])
                at = at_locus(x["contig"], x["strand"], x["lo"], x["hi"])
                recip = ("yes" if genes and genes[0][0] == a.human_gene else "no") if kind == "reverse" else ""
                cov = coverage([(h["qs"], h["qe"]) for h in x["hsps"]], qlen)
                table.append([g, kind, q, "genome", i, f"locus {i} ({len(x['hsps'])} HSP, {x['bits']:.1f} bits in all)",
                              f"{x['contig']}:{x['lo'] + 1}-{x['hi']}{x['strand']}", at,
                              ann_at(x["contig"], x["strand"], x["lo"], x["hi"]), f"{b['evalue']:.1e}", f"{b['bits']:.1f}",
                              f"{b['ident']:.0f}", f"{b['qs']}-{b['qe']}", f"{cov:.2f}", ""] + gene_cells(genes) + [recip])
                answers.setdefault((q, "genome"), []).append(
                    ({"s": f"{x['contig']}:{x['lo'] + 1}-{x['hi']}{x['strand']} ({len(x['hsps'])} HSP, {x['bits']:.1f} bits in all)",
                      "evalue": b["evalue"], "bits": b["bits"],
                      "ident": b["ident"]}, at, genes, recip,
                     (x["contig"], x["strand"], x["lo"], x["hi"], ann_at(x["contig"], x["strand"], x["lo"], x["hi"])), cov))

        # summary
        def describe(item):
            r, at, genes, recip, w, cov = item
            dist = ""
            if w and tloc and at != "yes":
                if w[0] == tloc[0]:
                    d = max(w[2] - tloc[3], tloc[2] - w[3], 0)
                    dist = f"; {d} nt from the {a.gene} locus ({'same' if w[1] == tloc[1] else 'other'} strand)"
                else: dist = f"; on another sequence than {a.gene} ({w[0]})"
            gtxt = f"best human gene {genes[0][0]} ({genes[0][2]:.1f} bits" + (
                f"; next {genes[1][0]}, {genes[1][2]:.1f})" if len(genes) > 1 else ")") if genes else "no human hit"
            name = r["s"] + (f" [{w[4]}]" if w and w[4] else "")
            return (f"{name}; E = {r['evalue']:.1e}, {r['bits']:.1f} bits, {r['ident']:.0f}% identity, query coverage "
                    f"{cov:.2f}; at the {a.gene} locus: {at}; {gtxt}" + (f"; reciprocal: {recip}" if recip else "") + dist)

        summary.append([g, f"forward: best human gene of {g} {a.gene}",
                        f"{fbest} ({fgenes[0][2]:.1f} bits; margin to the next gene {margin}"
                        + (f", {fgenes[1][0]}" if len(fgenes) > 1 else "") + ")" if fgenes else "no hit"])
        for db in ("proteins", "genome"):
            q = f"human {a.human_gene} (1-{len(h)})"
            items = answers.get((q, db), [])
            summary.append([g, f"reverse ({db}): best hit of human {a.human_gene}",
                            describe(items[0]) if items else "no hit"])
            if items:
                rbh = fbest == a.human_gene and items[0][1] == "yes"
                summary.append([g, f"reciprocal best hits ({db})", "yes" if rbh else "no"])
            for part in ("thioredoxin part", "PITH part"):
                q = next(k for k in hq if part in k)
                items = answers.get((q, db), [])
                summary.append([g, f"reverse ({db}): best hit of the {part} of human {a.human_gene}",
                                describe(items[0]) if items else "no hit"])
                other = next((x for x in items if x[1] != "yes"), None)
                if part == "thioredoxin part" and other and other is not items[0]:
                    summary.append([g, f"reverse ({db}): best hit of the {part} away from the {a.gene} locus", describe(other)])
            q = f"{g} {a.gene}"
            items = [x for x in answers.get((q, db), []) if x[1] != "yes"]
            summary.append([g, f"self ({db}): best hit of {g} {a.gene} away from its locus",
                            describe(items[0]) if items else "none"])
            if len(items) > 1:
                summary.append([g, f"self ({db}): second hit away from its locus", describe(items[1])])
        if not a.keep:
            for d in (pdb.dir, ndb.dir): shutil.rmtree(d, ignore_errors=True)
    if not a.keep: shutil.rmtree(hdb.dir, ignore_errors=True)
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in table: o.write("\t".join(str(x) for x in r) + "\n")
    with open(a.summary, "w") as o:
        o.write("genome\tquestion\tanswer\n")
        for r in summary: o.write("\t".join(str(x) for x in r) + "\n")
    print(f"{len(table)} rows -> {a.out}; {len(summary)} answers -> {a.summary}")


if __name__ == "__main__":
    main()
