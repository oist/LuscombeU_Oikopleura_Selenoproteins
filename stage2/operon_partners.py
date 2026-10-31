#!/usr/bin/env python3
"""Functions of the genes that share an operon with a selenoprotein gene, in the three O. dioica assemblies.

Question (2026-10-31): do the genes in the same operon as the selenoprotein genes have anything to do with redox,
selenium or selenoproteins?

Loci and operons: the same rules as stage2/operons.py (annotated transcripts with a coding sequence grouped into loci
by overlapping coding spans; loci overlapping a gold-standard selenoprotein gene replaced by it; consecutive loci on
the same strand with coding spans at most --max-gap nt apart form an operon).
Proteins: each annotated locus is represented by the protein of its longest coding sequence, translated from the
genome (internal stop codons read as X); the gold-standard genes are not searched.
Function: best hit of each protein in the reviewed human proteome of UniProt (Swiss-Prot, proteome UP000005640;
make human-uniprot), legacy blastp (E <= --evalue), by bit score; the hit's UniProt keywords and Gene Ontology terms
(directly annotated, not propagated) give the flags:
  keywords Oxidoreductase, Redox-active center, Antioxidant, Selenium, Selenocysteine;
  GO:0045454 cell redox homeostasis, GO:0006979 response to oxidative stress, GO:0098869 cellular oxidant
  detoxification, GO:0016209 antioxidant activity, GO:0004601 peroxidase activity, GO:0001887 selenium compound
  metabolic process, GO:0016259 selenocysteine metabolic process;
and redox_or_selenium, any of them.
Sets, in each assembly: partners (the other loci of the operons that hold a selenoprotein gene), immediate neighbours
(the loci just before and just after a selenoprotein gene in its operon), other loci in operons (operons without a
selenoprotein gene), all loci outside selenoprotein operons. Partners and neighbours are compared with the other loci
in operons by Fisher's exact test. The three assemblies are reported separately (their genes are orthologous, so
pooling them would count the same genes three times).
Outputs: --table (one row per partner), --summary (counts and tests per set), and per assembly in --work the proteins,
the blastp table and the best hit of every locus (<genome>.loci.tsv, not committed: several MB).
"""
import argparse, collections, csv, gzip, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from operons import DIOICA, Everything, gold_genes, loci_of, tr   # noqa: E402
from operons import read_annotation, read_fasta, rc               # noqa: E402  (from curation/gold_standard.py)

KEYWORDS = ["Oxidoreductase", "Redox-active center", "Antioxidant", "Selenium", "Selenocysteine"]
GO_TERMS = {"GO:0045454": "cell redox homeostasis", "GO:0006979": "response to oxidative stress",
            "GO:0098869": "cellular oxidant detoxification", "GO:0016209": "antioxidant activity",
            "GO:0004601": "peroxidase activity", "GO:0001887": "selenium compound metabolic process",
            "GO:0016259": "selenocysteine metabolic process"}
FLAGS = KEYWORDS + list(GO_TERMS.values())


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


def read_human(fasta, table):
    """UniProt accession -> {gene, protein, flags}; and the plain FASTA (accessions as names) for formatdb."""
    info = {}
    with opener(table) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kws = {k.strip() for k in (r.get("Keywords") or "").split(";")}
            gos = {g.strip() for g in (r.get("Gene Ontology IDs") or "").split(";")}
            flags = [k for k in KEYWORDS if k in kws] + [v for k, v in GO_TERMS.items() if k in gos]
            info[r["Entry"]] = {"gene": r.get("Gene Names (primary)", ""), "protein": r.get("Protein names", ""),
                                "flags": flags}
    seqs = []
    for name, s in read_fasta(fasta):
        acc = name.split("|")[1] if name.count("|") >= 2 else name
        seqs.append((acc, s))
    return info, seqs


def best_hits(path):
    """query -> (subject, pident, evalue, bits) of the best hit by bit score (blastall -m 8)."""
    best = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 12: continue
        q, s, pid, ev, bits = f[0], f[1], float(f[2]), float(f[10]), float(f[11])
        if q not in best or bits > best[q][3]: best[q] = (s, pid, ev, bits)
    return best


def fisher(a, n, b, m):
    """Fisher's exact test, a of n against b of m."""
    from scipy.stats import fisher_exact
    if n == 0 or m == 0: return ""
    return f"{fisher_exact([[a, n - a], [b, m - b]])[1]:.3g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--genomes", nargs="+", default=DIOICA)
    ap.add_argument("--max-gap", type=int, default=500)
    ap.add_argument("--human-fasta", default="data/references/uniprot_human_sp.fasta.gz")
    ap.add_argument("--human-table", default="data/references/uniprot_human_sp.tsv.gz")
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--work", default="work/stage2/operon_partners")
    ap.add_argument("--table", default="results/tables/stage2_operon_partners.tsv")
    ap.add_argument("--summary", default="results/tables/stage2_operon_partners_summary.tsv")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    info, hseqs = read_human(a.human_fasta, a.human_table)
    db = os.path.join(a.work, "human_sp.faa")
    if not os.path.exists(db + ".psq"):
        with open(db, "w") as o:
            for acc, s in hseqs: o.write(f">{acc}\n{s}\n")
        subprocess.run(["formatdb", "-i", db, "-p", "T", "-o", "F", "-l", os.path.join(a.work, "formatdb.log")], check=True)
    print(f"human proteome: {len(hseqs)} proteins, {len(info)} annotated", flush=True)
    part_rows, summary = [], []
    for g in a.genomes:
        genes = gold_genes(f"{a.gold}/{g}.gff3")
        ann = read_annotation(sheet[g].get("annotation", ""), Everything())
        if not ann:
            print(f"{g}: no annotation", flush=True); continue
        loci, members = build_operons(genes, ann, a.max_gap)
        contigs = {k: v.upper() for k, v in read_fasta(sheet[g]["fasta"], set(loci))}
        faa = os.path.join(a.work, f"{g}.proteins.faa")
        prot = {}
        for c, L in loci.items():
            for l in L:
                if l["gold"] or l["rep"] is None or c not in contigs: continue
                r = l["rep"]
                cds = "".join(contigs[c][s:e + 1] for s, e in sorted(r["cds"]))
                p = tr(cds if r["strand"] == "+" else rc(cds)).rstrip("*").replace("*", "X")
                if len(p) >= 10: prot[l["ids"][0]] = p
        with open(faa, "w") as o:
            for k, p in prot.items(): o.write(f">{k}\n{p}\n")
        out = os.path.join(a.work, f"{g}.blastp.tsv")
        if not os.path.exists(out) or os.path.getmtime(out) < os.path.getmtime(faa):
            tmp = out + ".part"
            subprocess.run(["blastall", "-p", "blastp", "-d", db, "-i", faa, "-e", a.evalue, "-m", "8", "-b", "5",
                            "-a", str(a.cpus), "-o", tmp], check=True)
            os.replace(tmp, out)
        best = best_hits(out)
        print(f"{g}: {sum(len(L) for L in loci.values())} loci, {len(prot)} proteins, {len(best)} with a human hit",
              flush=True)

        def hit(l):
            b = best.get(l["ids"][0])
            if b is None: return None, []
            return b, info.get(b[0], {}).get("flags", [])
        # sets
        sel_ops = {l["operon"] for L in loci.values() for l in L if l["gold"]}
        sets = collections.defaultdict(list)
        with open(os.path.join(a.work, f"{g}.loci.tsv"), "w") as o:
            o.write("locus\tcontig\tstrand\tspan\toperon_size\tset\thuman_hit\tbits\tflags\n")
            for c, L in loci.items():
                for l in L:
                    if l["gold"]: continue
                    m = members[l["operon"]]
                    if l["operon"] in sel_ops: s = "partners"
                    elif len(m) > 1: s = "other loci in operons"
                    else: s = "loci outside operons"
                    sets[s].append(l)
                    b, fl = hit(l)
                    o.write(f"{l['ids'][0]}\t{c}\t{l['strand']}\t{l['span'][0] + 1}-{l['span'][1] + 1}\t{len(m)}\t{s}\t"
                            f"{b[0] if b else ''}\t{b[3] if b else ''}\t{','.join(fl)}\n")
        # partners table and immediate neighbours
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
                h = info.get(b[0], {}) if b else {}
                part_rows.append({"genome": g, "operon": op, "operon_size": len(order),
                                  "selenoprotein_genes": ",".join(f"{gn} ({gk + 1})" for gk, gn in gold_here),
                                  "partner": l["ids"][0], "partner_place": k + 1, "relation": rel,
                                  "human_hit": b[0] if b else "", "human_gene": h.get("gene", ""),
                                  "human_protein": h.get("protein", ""),
                                  "bits": f"{b[3]:.0f}" if b else "", "evalue": f"{b[2]:.2g}" if b else "",
                                  "identity": f"{b[1]:.1f}" if b else "", "flags": ",".join(fl)})
        ref = sets["other loci in operons"]
        ref_n = [l for l in ref if l["ids"][0] in prot]
        for sname in ("partners", "immediate neighbours", "other loci in operons", "loci outside operons"):
            L = sets[sname]
            with_p = [l for l in L if l["ids"][0] in prot]
            hits = [l for l in with_p if best.get(l["ids"][0])]
            row = {"genome": g, "set": sname, "loci": len(L), "with_protein": len(with_p), "with_human_hit": len(hits)}
            counts = {f: sum(1 for l in with_p if f in hit(l)[1]) for f in FLAGS}
            any_ = sum(1 for l in with_p if hit(l)[1])
            row["redox_or_selenium"] = any_
            for f in FLAGS: row[f] = counts[f]
            if sname in ("partners", "immediate neighbours"):
                ref_any = sum(1 for l in ref_n if hit(l)[1])
                ref_ox = sum(1 for l in ref_n if "Oxidoreductase" in hit(l)[1])
                row["fisher_p_redox_or_selenium"] = fisher(any_, len(with_p), ref_any, len(ref_n))
                row["fisher_p_oxidoreductase"] = fisher(counts["Oxidoreductase"], len(with_p), ref_ox, len(ref_n))
            summary.append(row)
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
