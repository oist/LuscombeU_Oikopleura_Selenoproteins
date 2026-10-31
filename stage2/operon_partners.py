#!/usr/bin/env python3
"""Functions of the genes that share an operon with a selenoprotein gene, in the three O. dioica assemblies.

Question (2026-10-31): do the genes in the same operon as the selenoprotein genes have anything to do with redox,
selenium or selenoproteins?

Loci and operons: the same rules as stage2/operons.py (annotated transcripts with a coding sequence grouped into loci
by overlapping coding spans; loci overlapping a gold-standard selenoprotein gene replaced by it; consecutive loci on
the same strand with coding spans at most --max-gap nt apart form an operon).
Proteins: each annotated locus is represented by the protein of its longest coding sequence, translated from the
genome (internal stop codons read as X); the gold-standard genes are not searched. --proteins-only writes them
(<work>/<genome>.proteins.faa) and stops.
Search, either of:
  - DIAMOND against the Swiss-Prot database of the OIST Bioinformatics User Group (module DB/diamondDB/sprot/2022-07, or another release,
    https://github.com/oist/BioinfoUgrp/blob/master/DataBases.md), human entries only (--taxonlist 9606), outside the
    container (make operon-search), giving <work>/<genome>.diamond.tsv (qseqid sseqid pident evalue bitscore stitle);
  - or, when that file is absent, legacy blastp against --human-fasta (formatdb), giving <work>/<genome>.blastp.tsv.
The best hit (bit score, E <= 1e-5) of each protein gives the human protein. Flags:
  - from the hit's protein name (Swiss-Prot title or --human-table), always: "oxidoreductase (name)" for names with
    reductase, oxidase, dehydrogenase, oxygenase, peroxidase, peroxiredoxin, thioredoxin, glutaredoxin, catalase,
    superoxide dismutase, desaturase or hydroxylase; "selenium (name)" for names with seleno or selenium;
  - from UniProt (--human-table, make human-uniprot; current release, matched by accession), when available:
    keywords Oxidoreductase, Redox-active center, Antioxidant, Selenium, Selenocysteine; directly annotated GO terms
    GO:0045454 cell redox homeostasis, GO:0006979 response to oxidative stress, GO:0098869 cellular oxidant
    detoxification, GO:0016209 antioxidant activity, GO:0004601 peroxidase activity, GO:0001887 selenium compound
    metabolic process, GO:0016259 selenocysteine metabolic process (not propagated);
  - redox_or_selenium: any of them.
Sets, in each assembly: partners (the other loci of the operons that hold a selenoprotein gene), immediate neighbours
(the loci just before and just after a selenoprotein gene in its operon), other loci in operons (operons without a
selenoprotein gene), loci outside operons. Partners and neighbours are compared with the other loci in operons by
Fisher's exact test. The three assemblies are reported separately (their genes are orthologous, so pooling them
would count the same genes three times).
Outputs: --table (one row per partner), --summary (counts and tests per set), and per assembly in --work the proteins,
the search table and the best hit of every locus (<genome>.loci.tsv, not committed: several MB).
"""
import argparse, collections, csv, gzip, math, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from operons import DIOICA, Everything, gold_genes, loci_of, tr   # noqa: E402
from operons import read_annotation, read_fasta, rc               # noqa: E402  (from curation/gold_standard.py)

KEYWORDS = ["Oxidoreductase", "Redox-active center", "Antioxidant", "Selenium", "Selenocysteine"]
GO_TERMS = {"GO:0045454": "cell redox homeostasis", "GO:0006979": "response to oxidative stress",
            "GO:0098869": "cellular oxidant detoxification", "GO:0016209": "antioxidant activity",
            "GO:0004601": "peroxidase activity", "GO:0001887": "selenium compound metabolic process",
            "GO:0016259": "selenocysteine metabolic process"}
NAME_FLAGS = {"oxidoreductase (name)": re.compile(
    r"reductase|oxidase|dehydrogenase|oxygenase|peroxidase|peroxiredoxin|thioredoxin|glutaredoxin|catalase|"
    r"superoxide dismutase|desaturase|hydroxylase", re.I),
    "selenium (name)": re.compile(r"seleno|selenium", re.I)}
FLAGS = list(NAME_FLAGS) + KEYWORDS + list(GO_TERMS.values())


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def build_operons(genes, ann, max_gap):
    """Loci (contig -> list sorted by span) with gold genes in place, and operon -> members (stage2/operons.py)."""
    loci = loci_of(ann)
    over = {name: [l for l in loci[x["contig"]] if l["strand"] == x["strand"] and l["span"][0] <= x["span"][1]
                   and x["span"][0] <= l["span"][1]] for name, x in genes.items()}
    for name, x in genes.items():
        L = loci[x["contig"]]
        for l in over[name]:
            if l in L: L.remove(l)
        L.append({"strand": x["strand"], "span": x["span"], "ids": [name], "rep": None, "gold": name})
    for c in loci: loci[c].sort(key=lambda l: l["span"])
    members = collections.defaultdict(list)
    for c, L in loci.items():
        op = 0
        for i, l in enumerate(L):
            if i and L[i - 1]["strand"] == l["strand"] and 0 <= l["span"][0] - L[i - 1]["span"][1] - 1 <= max_gap:
                l["operon"] = L[i - 1]["operon"]
            else:
                op += 1; l["operon"] = f"{c}:{op}"
            members[l["operon"]].append(l)
    return loci, members


def order_of(m):
    """Members of an operon in transcription order."""
    return m if m[0]["strand"] == "+" else m[::-1]


def accession(sid):
    """UniProt accession from a sequence identifier (sp|P12345|NAME_HUMAN, sp|P12345.3|..., P12345.3 or P12345)."""
    acc = sid.split("|")[1] if sid.count("|") >= 2 else sid
    return acc.split(".")[0]


def read_human_table(path):
    """UniProt accession -> {gene, protein, keywords, go} (TSV of make human-uniprot); {} if absent."""
    info = {}
    if not path or not os.path.exists(path): return info
    with opener(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            info[r["Entry"]] = {"gene": r.get("Gene Names (primary)", ""), "protein": r.get("Protein names", ""),
                                "keywords": {k.strip() for k in (r.get("Keywords") or "").split(";")},
                                "go": {g.strip() for g in (r.get("Gene Ontology IDs") or "").split(";")}}
    return info


def title_names(stitle):
    """(gene, protein name) from a UniProt FASTA title: 'sp|P10599|THIO_HUMAN Thioredoxin OS=... GN=TXN PE=1 SV=3'."""
    m = re.search(r"\sGN=(\S+)", stitle)
    rest = stitle.split(None, 1)[1] if " " in stitle else ""
    return (m.group(1) if m else ""), rest.split(" OS=")[0]


def flags_of(acc, protein, info):
    fl = [k for k, rx in NAME_FLAGS.items() if rx.search(protein or "")]
    h = info.get(acc)
    if h:
        fl += [k for k in KEYWORDS if k in h["keywords"]] + [v for k, v in GO_TERMS.items() if k in h["go"]]
    return fl


def best_hits(path, fmt):
    """query -> {acc, pident, evalue, bits, gene, protein}: best hit by bit score."""
    best = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if fmt == "diamond":
            if len(f) < 5: continue
            q, s, pid, ev, bits = f[0], f[1], float(f[2]), float(f[3]), float(f[4])
            gene, prot = title_names(f[5]) if len(f) > 5 else ("", "")
        else:
            if len(f) < 12: continue
            q, s, pid, ev, bits = f[0], f[1], float(f[2]), float(f[10]), float(f[11])
            gene, prot = "", ""
        if q not in best or bits > best[q]["bits"]:
            best[q] = {"acc": accession(s), "pident": pid, "evalue": ev, "bits": bits, "gene": gene, "protein": prot}
    return best


def fisher(a, n, b, m):
    """Two-sided Fisher's exact test (sum of tables no more probable than the observed one), a of n against b of m."""
    if n == 0 or m == 0: return ""
    k, N = a + b, n + m
    def p(x): return math.comb(n, x) * math.comb(m, k - x) / math.comb(N, k)
    p0 = p(a)
    tot = sum(p(x) for x in range(max(0, k - m), min(k, n) + 1) if p(x) <= p0 * (1 + 1e-7))
    return f"{min(1.0, tot):.3g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--genomes", nargs="+", default=DIOICA)
    ap.add_argument("--max-gap", type=int, default=500)
    ap.add_argument("--proteins-only", action="store_true", help="write the proteins of the loci and stop")
    ap.add_argument("--human-fasta", default="", help="fallback search: human proteins for legacy blastp")
    ap.add_argument("--human-table", default="data/references/uniprot_human_sp.tsv.gz")
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--work", default="work/stage2/operon_partners")
    ap.add_argument("--table", default="results/tables/stage2_operon_partners.tsv")
    ap.add_argument("--summary", default="results/tables/stage2_operon_partners_summary.tsv")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    info = read_human_table(a.human_table)
    print(f"UniProt annotation table: {len(info)} entries" if info else
          f"no UniProt annotation table ({a.human_table}): flags from protein names only", flush=True)
    part_rows, summary = [], []
    for g in a.genomes:
        genes = gold_genes(f"{a.gold}/{g}.gff3")
        ann = read_annotation(sheet[g].get("annotation", ""), Everything())
        if not ann:
            print(f"{g}: no annotation", flush=True); continue
        loci, members = build_operons(genes, ann, a.max_gap)
        faa = os.path.join(a.work, f"{g}.proteins.faa")
        prot = {}
        contigs = {k: v.upper() for k, v in read_fasta(sheet[g]["fasta"], set(loci))}
        for c, L in loci.items():
            for l in L:
                if l["gold"] or l["rep"] is None or c not in contigs: continue
                r = l["rep"]
                cds = "".join(contigs[c][s:e + 1] for s, e in sorted(r["cds"]))
                p = tr(cds if r["strand"] == "+" else rc(cds)).rstrip("*").replace("*", "X")
                if len(p) >= 10: prot[l["ids"][0]] = p
        del contigs
        if a.proteins_only or not os.path.exists(faa):
            with open(faa, "w") as o:
                for k, p in prot.items(): o.write(f">{k}\n{p}\n")
        if a.proteins_only:
            print(f"{g}: {len(prot)} proteins written to {faa}", flush=True); continue
        dmd, bl = os.path.join(a.work, f"{g}.diamond.tsv"), os.path.join(a.work, f"{g}.blastp.tsv")
        if os.path.exists(dmd):
            best = best_hits(dmd, "diamond")
        elif a.human_fasta:
            db = os.path.join(a.work, "human.faa")
            if not os.path.exists(db + ".psq"):
                with open(db, "w") as o:
                    for name, s in read_fasta(a.human_fasta): o.write(f">{accession(name)}\n{s}\n")
                subprocess.run(["formatdb", "-i", db, "-p", "T", "-o", "F", "-l", os.path.join(a.work, "formatdb.log")],
                               check=True)
            if not os.path.exists(bl):
                subprocess.run(["blastall", "-p", "blastp", "-d", db, "-i", faa, "-e", a.evalue, "-m", "8", "-b", "5",
                                "-a", str(a.cpus), "-o", bl + ".part"], check=True)
                os.replace(bl + ".part", bl)
            best = best_hits(bl, "blast")
        else:
            sys.exit(f"{g}: no search results ({dmd}): run make operon-search first (or give --human-fasta)")
        for q, b in best.items():      # names from the UniProt table when available, else from the search title
            h = info.get(b["acc"])
            if h: b["gene"], b["protein"] = h["gene"] or b["gene"], h["protein"] or b["protein"]
            b["flags"] = flags_of(b["acc"], b["protein"], info)
        print(f"{g}: {sum(len(L) for L in loci.values())} loci, {len(prot)} proteins, {len(best)} with a human hit",
              flush=True)

        def hit(l):
            b = best.get(l["ids"][0])
            return (b, b["flags"]) if b else (None, [])
        sel_ops = {l["operon"] for L in loci.values() for l in L if l["gold"]}
        sets = collections.defaultdict(list)
        with open(os.path.join(a.work, f"{g}.loci.tsv"), "w") as o:
            o.write("locus\tcontig\tstrand\tspan\toperon_size\tset\thuman_hit\thuman_gene\tbits\tflags\n")
            for c, L in loci.items():
                for l in L:
                    if l["gold"]: continue
                    m = members[l["operon"]]
                    s = "partners" if l["operon"] in sel_ops else "other loci in operons" if len(m) > 1 \
                        else "loci outside operons"
                    sets[s].append(l)
                    b, fl = hit(l)
                    o.write(f"{l['ids'][0]}\t{c}\t{l['strand']}\t{l['span'][0] + 1}-{l['span'][1] + 1}\t{len(m)}\t{s}\t"
                            f"{b['acc'] if b else ''}\t{b['gene'] if b else ''}\t{b['bits'] if b else ''}\t"
                            f"{','.join(fl)}\n")
        for op in sorted(sel_ops):
            order = order_of(members[op])
            if len(order) < 2: continue
            gold_here = [(k, l["gold"]) for k, l in enumerate(order) if l["gold"]]
            for k, l in enumerate(order):
                if l["gold"]: continue
                rel = ";".join(f"{'upstream' if k < gk else 'downstream'} of {gn}" for gk, gn in gold_here
                               if abs(k - gk) == 1) or "other"
                if rel != "other": sets["immediate neighbours"].append(l)
                b, fl = hit(l)
                part_rows.append({"genome": g, "operon": op, "operon_size": len(order),
                                  "selenoprotein_genes": ",".join(f"{gn} ({gk + 1})" for gk, gn in gold_here),
                                  "partner": l["ids"][0], "partner_place": k + 1, "relation": rel,
                                  "human_hit": b["acc"] if b else "", "human_gene": b["gene"] if b else "",
                                  "human_protein": b["protein"] if b else "",
                                  "bits": f"{b['bits']:.0f}" if b else "", "evalue": f"{b['evalue']:.2g}" if b else "",
                                  "identity": f"{b['pident']:.1f}" if b else "", "flags": ",".join(fl)})
        ref_n = [l for l in sets["other loci in operons"] if l["ids"][0] in prot]
        for sname in ("partners", "immediate neighbours", "other loci in operons", "loci outside operons"):
            with_p = [l for l in sets[sname] if l["ids"][0] in prot]
            row = {"genome": g, "set": sname, "loci": len(sets[sname]), "with_protein": len(with_p),
                   "with_human_hit": sum(1 for l in with_p if best.get(l["ids"][0]))}
            any_ = sum(1 for l in with_p if hit(l)[1])
            row["redox_or_selenium"] = any_
            for f in FLAGS: row[f] = sum(1 for l in with_p if f in hit(l)[1])
            if sname in ("partners", "immediate neighbours"):
                row["fisher_p_redox_or_selenium"] = fisher(any_, len(with_p), sum(1 for l in ref_n if hit(l)[1]),
                                                           len(ref_n))
                ox = "Oxidoreductase" if info else "oxidoreductase (name)"
                row["fisher_p_oxidoreductase"] = fisher(row[ox], len(with_p), sum(1 for l in ref_n if ox in hit(l)[1]),
                                                        len(ref_n))
            summary.append(row)
    if a.proteins_only: return
    cols = ["genome", "operon", "operon_size", "selenoprotein_genes", "partner", "partner_place", "relation",
            "human_hit", "human_gene", "human_protein", "bits", "evalue", "identity", "flags"]
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in part_rows: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    scols = ["genome", "set", "loci", "with_protein", "with_human_hit", "redox_or_selenium"] + FLAGS + \
            ["fisher_p_redox_or_selenium", "fisher_p_oxidoreductase"]
    with open(a.summary, "w") as o:
        o.write("\t".join(scols) + "\n")
        for r in summary: o.write("\t".join(str(r.get(c, "")) for c in scols) + "\n")
    print(f"wrote {a.table} ({len(part_rows)} partners) and {a.summary}", flush=True)


if __name__ == "__main__":
    main()
