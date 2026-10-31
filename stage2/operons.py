#!/usr/bin/env python3
"""Step 10 of the plan: selenoprotein genes, operons and SECIS elements in the three O. dioica assemblies.

O. dioica transcribes many genes as operons, whose polycistronic transcripts are resolved into mRNAs by spliced-leader
(SL) trans-splicing at the 5' end of each downstream mRNA (Ganot et al. 2004). A selenoprotein mRNA needs its SECIS
element in its own 3' UTR; for the first gene of an operon, that 3' UTR ends before the next gene. The questions:
which selenoprotein genes lie in operons, at which position, where their SECIS lies relative to the next gene,
whether their mRNAs carry the SL, whether their neighbours are conserved between the cryptic species, and how this
compares with the other genes of each assembly.

Genes. The public annotation of each assembly (stage1/samplesheet.csv: RefSeq for OKI2018, ANISEED for OSKA2016,
GenBank for OdB3): transcripts with a coding sequence, grouped into loci when their coding spans overlap on the same
strand. Annotated loci whose coding span overlaps that of a gold-standard selenoprotein gene on its strand are
replaced by the gold-standard gene (curation/gold_standard/<genome>.gff3); an annotated mRNA that only holds the gene
in its 3' UTR stays a locus of its own (its own coding sequence). Coding spans are used because the three annotations
do not all have UTRs.
Operons: consecutive loci (in contig order) on the same strand whose coding spans are at most --max-gap nt apart
(the 500-nt intergenic distance of Plessy et al. 2024, applied here to coding spans, which makes it more inclusive);
position in transcription order: first, internal, last, or single (not in an operon).
For each gold-standard gene:
  - its neighbours on each side (any strand), their strand and distance (coding span to coding span), and the
    operon it belongs to (members, position);
  - the SECIS element: hits of the Oikopleura model (results/models/oikopleura_secis.cm, cmsearch, sense strand, at
    least --secis-bits) in the 2,000 nt after the stop codon (curation/gold_standard/<genome>.after_stop.fa); the
    SECIS of the gene is the best hit that ends before the coding sequence of the next gene on the same strand (the
    best hit of the window can be the SECIS of that gene), or else the best hit (noted); its distance from the stop
    codon, and the space left before the coding sequence of the next gene downstream on the same strand;
  - 3' UTR length and polyadenylation signal of the gold-standard model (results/tables/gold_standard.tsv);
  - SL: transcripts of the assembly (alignments of make gold-standard, work/stage2/gold_standard/<genome>.transcripts.paf;
    sequences from the samplesheet) that overlap the gene and start with the 3' end (--sl) of the 40-nt O. dioica SL
    ACTCATCCCATTTTTGAGTCCGATTTCGATTGTCTAACAG (Ganot et al. 2004, Fig. 1A; transcript assemblies can lack its 5'
    part, for example OSKA2016 transcripts starting with TTTCGATTGTCTAACAG), at most one mismatch, within the first
    --sl-window nt, or its reverse complement at the end; the genome position where each such transcript starts, its
    distance upstream of the start codon, whether the two genome bases before it are AG (the trans-splice acceptor;
    Ganot et al. 2004: consensus TTT(C/T/A)AG), and whether it lies between the coding sequences of the upstream
    neighbour and the gene; and the transcripts that run from the coding sequence of the upstream neighbour into the
    gene (polycistronic). Ganot et al. 2004 found intercistronic regions of 23 to 30 nt between the polyadenylation
    site of a gene and the trans-splice site of the next, and defined candidate operons by less than 300 bp between
    coding sequences; the gaps are reported so that this stricter definition can be applied;
  - neighbours conserved between assemblies: the proteins of the neighbours (translated from the representative
    coding sequence of each annotated locus) compared by blastp (E <= 1e-10) with the neighbours of the same gene in
    the other two assemblies;
  - for the next gene on the same strand: the trans-splice sites of its SL transcripts that lie between the two coding
    sequences (positions after the stop codon of the selenoprotein gene; * when the acceptor is not AG), and the
    space between the end of the SECIS element and the first of them: the polyadenylation site of the selenoprotein
    mRNA, and hence the end of its 3' UTR, lies before that site (23 to 30 nt before it in Ganot et al. 2004).
Background: the same operon definition for all loci of each assembly; selenoprotein genes against the others,
Fisher's exact test (operon or not; first position or not).
3' UTRs and intercistronic space (--utr-table, --utr-bins): for the loci followed by another gene of their operon
(first or internal), the distance from the stop codon to the next coding sequence and to the first trans-splice site
of the next gene (SL transcripts of the whole transcript set, placed as above), and, for all loci, the 3' UTR length
(annotated exons beyond the coding sequence; for the gold-standard genes, the 3' UTR of the gold-standard model,
from transcripts): selenoprotein genes against the others (Mann-Whitney U test; medians and quartiles), with the
distributions of the others in bins for the figure. The OdB3 annotation has no UTRs and its ESTs have masked 5' ends:
only the distance to the next coding sequence is available there.
Outputs: --table (one row per gene and assembly), --summary (one row per assembly), --utr-table, --utr-bins.
"""
import argparse, collections, csv, os, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "curation"))
from gold_standard import attributes, read_annotation, read_fasta, rc   # noqa: E402

DIOICA = ["OKI", "OSKA", "OdB3"]
CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


class Everything:
    def __contains__(self, x): return True


def tr(s): return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def run(cmd, **kw): return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def gold_genes(path):
    """name -> {contig, strand, cds: [(s0, e0 inclusive)], mrna: (s0, e0 inclusive)}."""
    out = {}
    for l in open(path):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9: continue
        at = attributes(f[8])
        if f[2] == "mRNA":
            out.setdefault(at["Name"], {"contig": f[0], "strand": f[6], "cds": []})["mrna"] = (int(f[3]) - 1, int(f[4]) - 1)
        elif f[2] == "CDS":
            name = at["Parent"].split("-", 1)[1].split("_", 1)[1]
            out.setdefault(name, {"contig": f[0], "strand": f[6], "cds": []})["cds"].append((int(f[3]) - 1, int(f[4]) - 1))
    for g in out.values():
        g["span"] = (min(s for s, _ in g["cds"]), max(e for _, e in g["cds"]))
    return out


def loci_of(ann):
    """Annotated loci: contig -> [{strand, span, ids, rep}] (transcripts with overlapping coding spans merged)."""
    by = collections.defaultdict(list)
    for tid, r in ann.items():
        if not r["cds"]: continue
        s, e = min(x for x, _ in r["cds"]), max(y for _, y in r["cds"])
        by[(r["contig"], r["strand"])].append((s, e, tid, r))
    out = collections.defaultdict(list)
    for (contig, strand), xs in by.items():
        xs.sort()
        cur = None
        for s, e, tid, r in xs:
            if cur and s <= cur["span"][1]:
                cur["span"] = (cur["span"][0], max(cur["span"][1], e)); cur["ids"].append(r["protein"])
                if sum(y - x + 1 for x, y in r["cds"]) > sum(y - x + 1 for x, y in cur["rep"]["cds"]): cur["rep"] = r
            else:
                cur = {"strand": strand, "span": (s, e), "ids": [r["protein"]], "rep": r, "gold": None}
                out[contig].append(cur)
    return out


def sl_match(seq, sl, window):
    """Position just after the SL at the 5' end of seq (sense), or None; at most one mismatch."""
    n = len(sl)
    for i in range(0, max(0, min(window, len(seq)) - n + 1)):
        if sum(a != b for a, b in zip(seq[i:i + n], sl)) <= 1: return i + n
    return None


def utr3_of(r):
    """Annotated 3' UTR length of a transcript (exons beyond the coding sequence), or None without exon rows."""
    if not r or not r["exons"] or not r["cds"]: return None
    if r["strand"] == "+":
        end = max(e for _, e in r["cds"])
        return sum(max(0, e - max(s, end + 1)) for s, e in r["exons"])
    start = min(s for s, _ in r["cds"])
    return sum(max(0, min(e, start) - s) for s, e in r["exons"])


def sl_junctions_all(alns, path, sl, window):
    """Trans-splice sites of all SL transcripts of an assembly: (contig, strand) -> sorted genome positions (0-based
    first base of the mRNA body); the strand is that of the mRNA."""
    names = {al["q"] for v in alns.values() for al in v}
    by_q = collections.defaultdict(list)
    for c, v in alns.items():
        for al in v: by_q[al["q"]].append((c, al))
    out = collections.defaultdict(list)
    for k, t in read_fasta(path, names):
        t = t.upper()
        p, same = sl_match(t, sl, window), True
        if p is None:
            p, same = sl_match(rc(t), sl, window), False
        if p is None: continue
        for c, al in by_q[k]:
            strand = al["strand"] if same else ("-" if al["strand"] == "+" else "+")
            qs_sense = al["qs"] if same else al["qlen"] - al["qe"]
            j = al["ts"] + (p - qs_sense) if strand == "+" else al["te"] - (p - qs_sense)
            out[(c, strand)].append(j)
    for k in out: out[k].sort()
    return out


def quartiles(v):
    import statistics
    if not v: return "", "", ""
    if len(v) < 4: return statistics.median(v), min(v), max(v)
    q = statistics.quantiles(v, n=4)
    return statistics.median(v), q[0], q[2]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--gold-table", default="results/tables/gold_standard.tsv")
    ap.add_argument("--paf-dir", default="work/stage2/gold_standard")
    ap.add_argument("--model", default="results/models/oikopleura_secis.cm")
    ap.add_argument("--secis-bits", type=float, default=20.0)
    ap.add_argument("--max-gap", type=int, default=500)
    ap.add_argument("--sl", default="GATTGTCTAACAG")
    ap.add_argument("--sl-window", type=int, default=80)
    ap.add_argument("--genomes", nargs="+", default=DIOICA)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--table", default="results/tables/stage2_operons.tsv")
    ap.add_argument("--summary", default="results/tables/stage2_operons_summary.tsv")
    ap.add_argument("--utr-table", default="results/tables/stage2_operons_utr.tsv")
    ap.add_argument("--utr-bins", default="results/tables/stage2_operons_utr_bins.tsv")
    ap.add_argument("--bin", type=int, default=25)
    ap.add_argument("--bin-max", type=int, default=3000)
    a = ap.parse_args()
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    gtab = {(r["genome"], r["gene"]): r for r in csv.DictReader(open(a.gold_table), delimiter="\t")}
    work = tempfile.mkdtemp(prefix="operons_")
    rows, summary, neighbour_prot, utr_rows, utr_bins = [], [], {}, [], []
    for g in a.genomes:
        genes = gold_genes(f"{a.gold}/{g}.gff3")
        ann = read_annotation(sheet[g].get("annotation", ""), Everything())
        if not ann:
            print(f"{g}: no annotation"); continue
        loci = loci_of(ann)
        # gold genes replace overlapping annotated loci (an annotated model can join two gold genes)
        over = {name: [l for l in loci[x["contig"]] if l["strand"] == x["strand"] and l["span"][0] <= x["span"][1]
                       and x["span"][0] <= l["span"][1]] for name, x in genes.items()}
        for name, x in genes.items():
            L = loci[x["contig"]]
            for l in over[name]:
                if l in L: L.remove(l)
            L.append({"strand": x["strand"], "span": x["span"], "ids": [name], "rep": None, "gold": name,
                      "annotated_as": ",".join(i for l in over[name] for i in l["ids"])})
        for c in loci: loci[c].sort(key=lambda l: l["span"])
        # operons
        for c, L in loci.items():
            op = 0
            for i, l in enumerate(L):
                if i and L[i - 1]["strand"] == l["strand"] and 0 <= l["span"][0] - L[i - 1]["span"][1] - 1 <= a.max_gap:
                    l["operon"] = L[i - 1]["operon"]
                else:
                    op += 1; l["operon"] = f"{c}:{op}"
        members = collections.defaultdict(list)
        for c, L in loci.items():
            for l in L: members[l["operon"]].append(l)
        def position(l):
            m = members[l["operon"]]
            if len(m) == 1: return "single"
            order = m if l["strand"] == "+" else m[::-1]
            k = order.index(l)
            return "first" if k == 0 else "last" if k == len(order) - 1 else "internal"
        # SECIS after the stop codon
        hits = collections.defaultdict(list)
        tbl = f"{work}/{g}.tbl"
        run(["cmsearch", "--toponly", "-T", str(a.secis_bits), "--cpu", str(a.cpus), "--tblout", tbl, "-o", os.devnull,
             a.model, f"{a.gold}/{g}.after_stop.fa"])
        for l in open(tbl):
            if l.startswith("#"): continue
            f = l.split()
            name = f[0].split("|")[0]
            hits[name].append((float(f[14]), int(f[7]), int(f[8])))
        # genome sequence for neighbour proteins
        contigs = {k: v.upper() for k, v in read_fasta(sheet[g]["fasta"], {x["contig"] for x in genes.values()})}
        # transcripts for the SL
        paf = f"{a.paf_dir}/{g}.transcripts.paf"
        alns = collections.defaultdict(list)
        if os.path.exists(paf):
            for l in open(paf):
                f = l.split("\t")
                if "tp:A:P" not in l: continue
                alns[f[5]].append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4],
                                   "ts": int(f[7]), "te": int(f[8]) - 1})
        want = set()
        for x in genes.values():
            lo, hi = x["mrna"]
            for al in alns.get(x["contig"], []):
                if al["ts"] <= hi + 2000 and al["te"] >= lo - 2000: want.add(al["q"])
        txseq = {k: v.upper() for k, v in read_fasta(sheet[g]["transcripts"], want)} if want and sheet[g].get("transcripts") else {}
        for name, x in sorted(genes.items()):
            L = loci[x["contig"]]
            me = next(l for l in L if l["gold"] == name)
            i = L.index(me)
            plus = x["strand"] == "+"
            up = L[i - 1] if plus and i > 0 else (L[i + 1] if not plus and i + 1 < len(L) else None)
            down = L[i + 1] if plus and i + 1 < len(L) else (L[i - 1] if not plus and i > 0 else None)
            def gap(n):
                if n is None: return ""
                return (n["span"][0] - x["span"][1] - 1) if n["span"][0] > x["span"][1] else (x["span"][0] - n["span"][1] - 1)
            def desc(n):
                if n is None: return "none (contig end)"
                return f"{n['ids'][0]} ({'same' if n['strand'] == x['strand'] else 'opposite'} strand, {gap(n)} nt)"
            for side, n in (("up", up), ("down", down)):
                if n is not None and n["rep"] is not None and n["rep"]["contig"] in contigs:
                    r = n["rep"]
                    cds = "".join(contigs[r["contig"]][s:e + 1] for s, e in sorted(r["cds"]))
                    neighbour_prot[f"{g}|{name}|{side}"] = tr(cds if r["strand"] == "+" else rc(cds)).rstrip("*")
            m = members[me["operon"]]
            down_same = down if down is not None and down["strand"] == x["strand"] else None
            # SECIS: the best hit that ends before the coding sequence of the next gene on the same strand (the best
            # hit in the 2,000 nt can be the SECIS of that gene: SELENOW_pair2 and SELENOW_pair1), else the best hit
            hs = hits.get(name, [])
            limit = gap(down_same) if down_same is not None else None
            before = [h for h in hs if limit is None or h[2] <= limit]
            sec = max(before) if before else (max(hs) if hs else None)
            sec_note = "" if before or not hs else "only beyond the next gene"
            gt = gtab.get((g, name), {})
            # SL transcripts
            lo, hi = x["mrna"]
            cs, ce = x["span"]
            start = cs if plus else ce          # first base of the coding sequence (gene strand)
            up_end = (up["span"][1] if plus else up["span"][0]) if up is not None and up["strand"] == x["strand"] else None
            def junctions(lo_, hi_):
                """SL transcripts overlapping [lo_, hi_] on the gene strand: (transcript, trans-splice site: first
                genome base after the SL, 0-based, acceptor is AG)."""
                out = []
                for al in alns.get(x["contig"], []):
                    if not (al["ts"] <= hi_ and al["te"] >= lo_) or al["q"] not in txseq: continue
                    t = txseq[al["q"]]
                    sense = (al["strand"] == "+") == plus     # stored orientation is the sense of the gene
                    pos = sl_match(t if sense else rc(t), a.sl, a.sl_window)
                    if pos is None: continue
                    qs_sense = al["qs"] if sense else al["qlen"] - al["qe"]
                    g5 = al["ts"] if plus else al["te"]
                    j = g5 + (pos - qs_sense) if plus else g5 - (pos - qs_sense)
                    cs_ = contigs.get(x["contig"], "")
                    acc = (cs_[j - 2:j] if plus else rc(cs_[j + 1:j + 3])) if cs_ else ""
                    out.append((al["q"], j, acc == "AG"))
                return out
            sl_sites = []
            for q, j, ag in junctions(cs, ce):
                dist = (start - j) if plus else (j - start)
                between = up_end is None or (j > up_end if plus else j < up_end)
                sl_sites.append((q, dist, between, ag))
            through = [al["q"] for al in alns.get(x["contig"], []) if al["ts"] <= ce and al["te"] >= cs and up_end is not None
                       and ((plus and al["ts"] <= up_end) or (not plus and al["te"] >= up_end))]
            # trans-splice sites of the next gene on the same strand, between the two coding sequences
            next_sl = []
            if down_same is not None:
                dlo, dhi = down_same["span"]
                for q, j, ag in junctions(dlo, dhi):
                    off = (j - ce) if plus else (cs - j)        # 1-based position after the stop codon
                    if 0 < off <= gap(down_same) + 1: next_sl.append((off, ag))
            first_next = min((o for o, _ in next_sl), default=None)
            rows.append({
                "genome": g, "gene": name, "contig": x["contig"], "strand": x["strand"],
                "coding_span": f"{cs + 1}-{ce + 1}", "annotated_as": me.get("annotated_as", "") or "not annotated",
                "upstream_neighbour": desc(up), "downstream_neighbour": desc(down),
                "operon_size": len(m), "operon_position": position(me),
                "operon_members": ",".join(l["ids"][0] for l in (m if plus else m[::-1])),
                "secis_bits": f"{sec[0]:.1f}" if sec else "none",
                "secis_after_stop": f"{sec[1]}-{sec[2]}" if sec else "",
                "secis_note": sec_note,
                "secis_hits_after_stop": ";".join(f"{h[1]}-{h[2]}:{h[0]:.1f}" for h in sorted(hs, key=lambda h: h[1])),
                "next_gene_same_strand_after_stop": gap(down_same) if down_same is not None else "",
                "next_gene_sl_sites_after_stop": ",".join(f"{o}{'' if ag else '*'}" for o, ag in sorted(next_sl)),
                "space_secis_to_next_sl_site": (first_next - sec[2] - 1) if (sec and first_next is not None) else "",
                "space_after_secis": (gap(down_same) - sec[2]) if (sec and down_same is not None) else "",
                "utr3_length": gt.get("utr3_length", ""), "polyA_signal": gt.get("polyA_signal", ""),
                "sl_transcripts": len(sl_sites),
                "sl_start_upstream_of_cds": ",".join(str(d) for _, d, *_ in sorted(sl_sites, key=lambda s: s[1]))[:200],
                "sl_acceptor_AG": sum(1 for *_, ag in sl_sites if ag),
                "sl_between_neighbour_and_gene": sum(1 for _, _, b, _ in sl_sites if b),
                "transcripts_from_upstream_neighbour": len(through),
            })
        n_all = sum(len(L) for L in loci.values())
        op_all = sum(1 for L in loci.values() for l in L if position(l) != "single")
        sel = [r for r in rows if r["genome"] == g]
        op_sel = sum(1 for r in sel if r["operon_position"] != "single")
        first_sel = sum(1 for r in sel if r["operon_position"] == "first")
        first_all = sum(1 for L in loci.values() for l in L if position(l) == "first")
        from scipy.stats import fisher_exact
        p_op = fisher_exact([[op_sel, len(sel) - op_sel], [op_all - op_sel, n_all - len(sel) - (op_all - op_sel)]])[1]
        p_first = fisher_exact([[first_sel, op_sel - first_sel],
                                [first_all - first_sel, (op_all - op_sel) - (first_all - first_sel)]])[1] if op_sel else ""
        summary.append({"genome": g, "loci": n_all, "loci_in_operons": op_all,
                        "fraction_in_operons": f"{op_all / n_all:.3f}", "first_in_operon": first_all,
                        "selenoprotein_genes": len(sel), "selenoprotein_in_operons": op_sel,
                        "selenoprotein_first": first_sel,
                        "fisher_p_operon": f"{p_op:.3g}", "fisher_p_first_among_operonic": f"{p_first:.3g}" if p_first != "" else ""})
        print(f"{g}: {n_all} loci, {op_all} in operons; selenoprotein genes {op_sel} of {len(sel)} in operons", flush=True)
        # 3' UTR and intercistronic space of the genes followed by another gene of their operon (first or internal),
        # and 3' UTR of the last genes and of the genes outside operons: selenoprotein genes against the others
        junc = sl_junctions_all(alns, sheet[g]["transcripts"], a.sl, a.sl_window) if alns and sheet[g].get("transcripts") else {}
        vals = collections.defaultdict(list)     # (measure, context, group) -> values
        for c, L in loci.items():
            for l in L:
                pos = position(l)
                grp = "selenoprotein" if l["gold"] else "other"
                u = int(gtab[(g, l["gold"])]["utr3_length"]) if l["gold"] and gtab.get((g, l["gold"]), {}).get("utr3_length", "") not in ("", "0") \
                    else (None if l["gold"] else utr3_of(l["rep"]))
                ctx = "followed in operon" if pos in ("first", "internal") else "last or single"
                if u is not None and u > 0: vals[("3' UTR", ctx, grp)].append(u)
                if pos not in ("first", "internal"): continue
                m = members[l["operon"]]
                order = m if l["strand"] == "+" else m[::-1]
                n = order[order.index(l) + 1]
                if l["strand"] == "+":
                    stop, nxt = l["span"][1], n["span"][0]
                    gap_ = nxt - stop - 1
                    js = [j for j in junc.get((c, "+"), []) if stop < j <= nxt]
                    first_sl = (min(js) - stop) if js else None
                else:
                    stop, nxt = l["span"][0], n["span"][1]
                    gap_ = stop - nxt - 1
                    js = [j for j in junc.get((c, "-"), []) if nxt <= j < stop]
                    first_sl = (stop - max(js)) if js else None
                vals[("stop to next coding sequence", ctx, grp)].append(gap_)
                if first_sl is not None: vals[("stop to next trans-splice site", ctx, grp)].append(first_sl)
        from scipy.stats import mannwhitneyu
        for (meas, ctx) in sorted({(k[0], k[1]) for k in vals}):
            o, sp = vals.get((meas, ctx, "other"), []), vals.get((meas, ctx, "selenoprotein"), [])
            pval = mannwhitneyu(sp, o).pvalue if sp and o else ""
            for grp, v in (("other", o), ("selenoprotein", sp)):
                med, q1, q3 = quartiles(v)
                utr_rows.append({"genome": g, "measure": meas, "context": ctx, "group": grp, "n": len(v),
                                 "median": med, "q1": q1, "q3": q3,
                                 "mann_whitney_p": (f"{pval:.3g}" if pval != "" else "") if grp == "selenoprotein" else ""})
            nb = a.bin_max // a.bin
            counts = [0] * (nb + 1)
            for x in o: counts[min(max(x, 0) // a.bin, nb)] += 1
            for b, n_ in enumerate(counts):
                utr_bins.append({"genome": g, "measure": meas, "context": ctx, "bin_start": b * a.bin, "count": n_})
    # neighbours conserved between assemblies
    if neighbour_prot:
        fa = f"{work}/nb.faa"
        with open(fa, "w") as o:
            for k, s in neighbour_prot.items():
                if len(s) >= 30: o.write(f">{k}\n{s}\n")
        run(["formatdb", "-i", fa, "-p", "T"], cwd=work)
        hit = collections.defaultdict(set)
        for l in run(["blastall", "-p", "blastp", "-d", fa, "-i", fa, "-e", "1e-10", "-m", "8", "-a", str(a.cpus)],
                     cwd=work).splitlines():
            q, s = l.split("\t")[:2]
            qg, qn, qside = q.split("|"); sg, sn, sside = s.split("|")
            if qn == sn and qg != sg: hit[(qg, qn, qside)].add(f"{sg}:{sside}")
        for r in rows:
            for side in ("up", "down"):
                r[f"{side}stream_neighbour_conserved_in"] = ",".join(sorted(hit.get((r["genome"], r["gene"], side), []))) or "-"
    cols = list(rows[0])
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    with open(a.summary, "w") as o:
        cols = list(summary[0])
        o.write("\t".join(cols) + "\n")
        for r in summary: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    if utr_rows:
        cols = list(utr_rows[0])
        with open(a.utr_table, "w") as o:
            o.write("\t".join(cols) + "\n")
            for r in utr_rows: o.write("\t".join(str(r[c]) for c in cols) + "\n")
        cols = list(utr_bins[0])
        with open(a.utr_bins, "w") as o:
            o.write("# counts of the other (non-selenoprotein) loci in bins of --bin nt; the last bin holds the values from --bin-max\n")
            o.write("\t".join(cols) + "\n")
            for r in utr_bins: o.write("\t".join(str(r[c]) for c in cols) + "\n")
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
