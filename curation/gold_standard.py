#!/usr/bin/env python3
"""Step 2 of the plan (docs/plan.md, 2026-10-24): gold-standard gene models of the curated O. dioica selenoprotein
genes in OKI2018, OSKA2016 and OdB3.

For each curated selenocysteine codon (curation/baseline_sec_genes.tsv):
  1. Transcripts. The transcripts of the assembly (stage1/samplesheet.csv: the OKI2018 Trinity transcriptome, the
     OSKA2016 transcript assembly, the OdB3 cDNA collection) are aligned to it with minimap2 (-x splice:hq -un, one
     primary alignment each; cached in --work). -un turns off the preference for GT-AG splice sites, which otherwise
     moves the junctions of non-canonical introns (frequent in Oikopleura, for example GA-AG) and creates artificial
     indels; junctions between identical bases stay ambiguous, so the splice sites reported are indicative. A transcript supports the gene when an aligned exon holds the three bases of
     the codon and they read TGA on the gene strand.
  2. Open reading frame. In each supporting transcript (in the sense of the gene), the frame of the TGA is extended
     upstream to the previous stop codon (TAA, TAG or TGA) and the coding sequence starts at the first ATG after it;
     downstream it ends at the next stop codon. The TGA is read as selenocysteine (U). A frame that reaches the end
     of the transcript is flagged as open at that end.
  3. Choice. The model is the supporting transcript with a complete frame (start and stop) and no insertion or
     deletion against the assembly inside the coding sequence, then the longest coding sequence, then the longest
     3' UTR. The furthest 3' end of the transcripts that share the stop codon is also given, as a genomic extent
     after the stop codon (utr3_max_extent, introns included).
  4. Without a supporting transcript, the Selenoprofiles prediction of stage 1 that holds the codon
     (<stage1>/<genome>.selenoprofiles.gtf) gives the coding sequence, without UTRs.
Checks and comparisons, reported per gene: the protein context of the curated codon; introns and their splice
sites (non-canonical ones listed, for example GA-AG); insertions, deletions and mismatches between the chosen
transcript and the assembly inside the coding sequence (an indel means that the assembly carries another allele or
an error: the protein given is the transcript's, flagged); and how the public annotation of the assembly treats the
gene (overlapping annotated coding sequences: selenocysteine declared, codon inside the coding sequence without
declaration, coding sequence ending at the codon or starting after it, no model).
Outputs in --out: <genome>.gff3 (gene, mRNA, exon and CDS features; transl_except=(pos:..,aa:Sec) on the CDS rows),
<genome>.proteins.faa (U for selenocysteine), <genome>.utr3.fa (3' UTR of the chosen transcript),
<genome>.after_stop.fa (--after-stop nucleotides of genome after the stop codon, on the gene strand), and --table.
"""
import argparse, collections, csv, gzip, os, re, statistics, subprocess

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}


def rc(s): return s.translate(COMP)[::-1]


def translate(s, sec=None):
    """Translation; the codon starting at index sec (if given) is read as U."""
    return "".join("U" if i == sec else CODE.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))


def opener(path): return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path, keep=None):
    name, buf = None, []
    with opener(path) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None and (keep is None or name in keep): yield name, "".join(buf)
                name, buf = l[1:].split()[0], []
            elif keep is None or name in keep: buf.append(l.strip())
    if name is not None and (keep is None or name in keep): yield name, "".join(buf)


def attributes(field):
    out = {}
    for x in field.strip().strip(";").split(";"):
        x = x.strip()
        if not x: continue
        if "=" in x: k, v = x.split("=", 1)
        else: k, _, v = x.partition(" ")
        out[k] = v.strip('"')
    return out


# ---------------------------------------------------------------------------------------------------------------
# alignments

def paf_alignments(paf, wanted):
    """Primary alignments of the PAF that overlap a codon of wanted ({contig: [(start0, end0)]})."""
    with open(paf) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            t, ts, te = f[5], int(f[7]), int(f[8])
            if t not in wanted or not any(ts <= a and b <= te for a, b in wanted[t]): continue
            tags = {x.split(":", 2)[0]: x.split(":", 2)[2] for x in f[12:] if x.count(":") >= 2}
            if tags.get("tp", "P") != "P" or "cg" not in tags: continue
            yield {"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "t": t,
                   "ts": ts, "te": te, "cigar": tags["cg"]}


def walk(aln):
    """Map of an alignment, with the transcript in the orientation of the genome forward strand:
    gpos -> qpos for aligned bases, blocks (exons: target spans without intron), introns, inserted query positions
    and deleted target positions."""
    q = aln["qs"] if aln["strand"] == "+" else aln["qlen"] - aln["qe"]
    t = aln["ts"]
    g2q, q2g, ins, dels, introns, blocks = {}, {}, [], [], [], []
    bstart = t
    for n, op in re.findall(r"(\d+)([MIDN=X])", aln["cigar"]):
        n = int(n)
        if op in "M=X":
            for i in range(n): g2q[t + i] = q + i; q2g[q + i] = t + i
            t += n; q += n
        elif op == "I":
            ins += [(q + i, t) for i in range(n)]; q += n
        elif op == "D":
            dels += [t + i for i in range(n)]; t += n
        elif op == "N":
            blocks.append((bstart, t)); introns.append((t, t + n)); t += n; bstart = t
    blocks.append((bstart, t))
    return {"g2q": g2q, "q2g": q2g, "ins": ins, "dels": dels, "introns": introns, "blocks": blocks}


# ---------------------------------------------------------------------------------------------------------------
# models

def orf(mrna, s):
    """Frame through the TGA at s: (start, stop) indices (stop = first base of the stop codon), and open ends."""
    i = s - 3
    while i >= 0 and mrna[i:i + 3].upper() not in STOPS: i -= 3
    up_stop = i if i >= 0 else None
    j = (up_stop + 3) if up_stop is not None else s % 3
    start = next((k for k in range(j, s, 3) if mrna[k:k + 3].upper() == "ATG"), None)
    k = s + 3
    while k + 3 <= len(mrna) and mrna[k:k + 3].upper() not in STOPS: k += 3
    stop = k if k + 3 <= len(mrna) else None
    return start, stop, up_stop is None, stop is None


def splice_sites(genome, introns, strand):
    out = []
    for a, b in introns:
        site = f"{genome[a:a + 2]}-{genome[b - 2:b]}".upper() if strand == "+" else \
            f"{rc(genome[b - 2:b])}-{rc(genome[a:a + 2])}".upper()
        out.append(site)
    return out


def transcript_model(aln, tx, genome, codon0, strand):
    """Model from one transcript, or None when it does not hold the codon as TGA in an aligned exon."""
    w = walk(aln)
    fwd = tx if aln["strand"] == "+" else rc(tx)
    L = len(fwd)
    qs = [w["g2q"].get(codon0 + i) for i in range(3)]
    if None in qs or qs != [qs[0], qs[0] + 1, qs[0] + 2]: return None
    if strand == "+":
        mrna, s = fwd, qs[0]
        m2g = lambda i: w["q2g"].get(i)
    else:
        mrna, s = rc(fwd), L - 1 - qs[2]
        m2g = lambda i: w["q2g"].get(L - 1 - i)
    if mrna[s:s + 3].upper() != "TGA": return None
    start, stop, open5, open3 = orf(mrna, s)
    m = {"tx": aln["q"], "start": start, "stop": stop, "open5": open5, "open3": open3, "sec_idx": s,
         "blocks": w["blocks"], "introns": w["introns"], "utr3": "", "cds": "", "protein": ""}
    if start is None: return m
    end = (stop + 3) if stop is not None else len(mrna) - (len(mrna) - s) % 3
    m["cds"] = mrna[start:end]
    m["protein"] = translate(m["cds"], s - start).rstrip("*")
    m["utr3"] = mrna[end:] if stop is not None else ""
    m["utr5_len"] = start
    gpos = [m2g(i) for i in range(start, end)]
    m["cds_genome"] = [p for p in gpos if p is not None]
    lo, hi = min(m["cds_genome"]), max(m["cds_genome"])
    m["cds_insertions"] = sum(p is None for p in gpos)
    m["cds_deletions"] = sum(lo <= d <= hi for d in w["dels"])
    m["cds_mismatches"] = sum(genome[p].upper() != (mrna[start + i].upper() if strand == "+" else rc(mrna[start + i]).upper())
                              for i, p in enumerate(gpos) if p is not None)
    m["stop_genome"] = (hi if strand == "+" else lo) if stop is not None else None
    tx_end_q = len(mrna) - 1
    m["tx_end_genome"] = m2g(tx_end_q) if m2g(tx_end_q) is not None else (w["blocks"][-1][1] - 1 if strand == "+" else w["blocks"][0][0])
    m["sites"] = splice_sites(genome, [x for x in w["introns"]], strand)
    return m


def segments(blocks, lo, hi):
    """Exon blocks clipped to [lo, hi] (0-based, inclusive hi)."""
    return [(max(a, lo), min(b - 1, hi)) for a, b in blocks if b - 1 >= lo and a <= hi]


def prediction_model(gtf_rows, genome, codon0, strand):
    """Model from the Selenoprofiles prediction whose CDS holds the codon."""
    for pid, cds in gtf_rows.items():
        if not any(a <= codon0 and codon0 + 2 <= b for a, b, _ in cds): continue
        if cds[0][2] != strand: continue
        segs = sorted((a, b) for a, b, _ in cds)
        seq = "".join(genome[a:b + 1] for a, b in segs)
        if strand == "-": seq = rc(seq)
        # index of the codon in the CDS
        off, idx = 0, None
        order = segs if strand == "+" else segs[::-1]
        for a, b in order:
            if a <= codon0 <= b:
                idx = off + (codon0 - a if strand == "+" else b - (codon0 + 2)); break
            off += b - a + 1
        if idx is None or seq[idx:idx + 3].upper() != "TGA": continue
        return {"tx": f"prediction {pid}", "cds": seq, "protein": translate(seq, idx).rstrip("*"),
                "start": 0, "stop": None, "open5": not seq[:3].upper() == "ATG",
                "open3": seq[-3:].upper() not in STOPS, "blocks": [(a, b + 1) for a, b in segs], "introns": [],
                "cds_genome": [p for a, b in segs for p in range(a, b + 1)], "utr3": "", "utr5_len": 0,
                "cds_insertions": 0, "cds_deletions": 0, "cds_mismatches": 0, "sites": [], "from_prediction": True,
                "stop_genome": None, "tx_end_genome": None}
    return None


def read_gtf(path):
    """Selenoprofiles GTF: prediction -> [(start0, end0 inclusive, strand)] of its CDS rows."""
    out = collections.defaultdict(list)
    if not path or not os.path.exists(path): return {}
    with open(path) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS": continue
            at = attributes(f[8])
            pid = at.get("transcript_id") or at.get("gene_id")
            out[(f[0], pid)].append((int(f[3]) - 1, int(f[4]) - 1, f[6]))
    by_contig = collections.defaultdict(dict)
    for (c, pid), rows in out.items(): by_contig[c][pid] = rows
    return by_contig


def read_annotation(path, contigs):
    """Annotated mRNAs on the contigs: id -> {"contig", "strand", "cds": [(s0, e0)], "protein", "sec"}."""
    out = {}
    if not path or not os.path.exists(path): return out
    with opener(path) as fh:
        for l in fh:
            if l.startswith("#"): continue
            f = l.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in contigs or f[2] != "CDS": continue
            at = attributes(f[8])
            parent = at.get("Parent", at.get("ID", "?"))
            r = out.setdefault(parent, {"contig": f[0], "strand": f[6], "cds": [], "protein": at.get("protein_id", parent),
                                        "sec": False})
            r["cds"].append((int(f[3]) - 1, int(f[4]) - 1))
            if "Sec" in at.get("transl_except", ""): r["sec"] = True
    return out


def annotation_status(ann, contig, strand, codon0, lo, hi):
    """How the annotation treats the gene: status and the overlapping annotated proteins."""
    over = []
    for mid, r in ann.items():
        if r["contig"] != contig or r["strand"] != strand: continue
        a, b = min(x for x, _ in r["cds"]), max(y for _, y in r["cds"])
        if b < lo or a > hi: continue
        over.append((mid, r, a, b))
    if not over: return "no annotated coding sequence", "-"
    ids = ",".join(sorted({r["protein"] for _, r, _, _ in over}))
    statuses = []
    for _, r, a, b in over:
        inside = any(x <= codon0 and codon0 + 2 <= y for x, y in r["cds"])
        if inside: statuses.append("selenocysteine declared" if r["sec"] else "codon inside the coding sequence, not declared")
        elif (strand == "+" and b < codon0) or (strand == "-" and a > codon0 + 2): statuses.append("coding sequence ends before the codon")
        else: statuses.append("coding sequence starts after the codon")
    order = ["selenocysteine declared", "codon inside the coding sequence, not declared",
             "coding sequence ends before the codon", "coding sequence starts after the codon"]
    return next(s for s in order if s in statuses), ids


# ---------------------------------------------------------------------------------------------------------------

def gff_lines(g, gene, contig, strand, m, sec_genome):
    """GFF3 rows of one model (0-based inclusive inputs, 1-based output)."""
    lo_cds, hi_cds = min(m["cds_genome"]), max(m["cds_genome"])
    exons = [(a, b - 1) for a, b in m["blocks"]]
    if m.get("from_prediction"): exons = [(a, b - 1) for a, b in m["blocks"]]
    glo, ghi = min(a for a, _ in exons), max(b for _, b in exons)
    gid, rid = f"gene-{g}_{gene}", f"rna-{g}_{gene}"
    note = f"evidence={m['tx']}"
    rows = [(contig, "gold_standard", "gene", glo + 1, ghi + 1, ".", strand, ".", f"ID={gid};Name={gene}"),
            (contig, "gold_standard", "mRNA", glo + 1, ghi + 1, ".", strand, ".", f"ID={rid};Parent={gid};Name={gene};{note}")]
    for k, (a, b) in enumerate(sorted(exons), 1):
        rows.append((contig, "gold_standard", "exon", a + 1, b + 1, ".", strand, ".", f"ID=exon-{g}_{gene}-{k};Parent={rid}"))
    cds = segments([(a, b + 1) for a, b in exons], lo_cds, hi_cds)
    order = cds if strand == "+" else cds[::-1]
    done = 0
    te = f"transl_except=(pos:{sec_genome[0] + 1}..{sec_genome[1] + 1},aa:Sec)" if strand == "+" else \
        f"transl_except=(pos:complement({sec_genome[0] + 1}..{sec_genome[1] + 1}),aa:Sec)"
    phased = []
    for a, b in order:
        phased.append((a, b, (3 - done % 3) % 3)); done += b - a + 1
    for a, b, ph in sorted(phased):
        rows.append((contig, "gold_standard", "CDS", a + 1, b + 1, ".", strand, str(ph), f"ID=cds-{g}_{gene};Parent={rid};{te}"))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3"])
    ap.add_argument("--baseline", default="curation/baseline_sec_genes.tsv")
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--stage1", default="work/stage1/pipeline/selenoprofiles")
    ap.add_argument("--work", default="work/stage2/gold_standard")
    ap.add_argument("--out", default="curation/gold_standard")
    ap.add_argument("--table", default="results/tables/gold_standard.tsv")
    ap.add_argument("--after-stop", type=int, default=2000)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True); os.makedirs(a.out, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    base = [r for r in csv.DictReader(open(a.baseline), delimiter="\t") if r["genome"] in a.genomes]
    report = []
    for g in a.genomes:
        genes = [r for r in base if r["genome"] == g]
        if not genes: continue
        contigs = {r["contig"] for r in genes}
        genome = dict(read_fasta(sheet[g]["fasta"], contigs))
        wanted = collections.defaultdict(list)
        for r in genes: wanted[r["contig"]].append((int(r["sec_codon"]) - 1, int(r["sec_codon"]) + 2))
        alns = []
        tx_path = sheet[g].get("transcripts", "")
        if tx_path and os.path.exists(tx_path):
            paf = f"{a.work}/{g}.transcripts.paf"
            if not os.path.exists(paf):
                with open(paf + ".tmp", "w") as o:
                    subprocess.run(["minimap2", "-x", "splice:hq", "-un", "--secondary=no", "-c", "-t", str(a.cpus),
                                    sheet[g]["fasta"], tx_path], stdout=o, check=True)
                os.replace(paf + ".tmp", paf)
            alns = list(paf_alignments(paf, wanted))
        txseq = dict(read_fasta(tx_path, {x["q"] for x in alns})) if alns else {}
        gtf = read_gtf(f"{a.stage1}/{g}.selenoprofiles.gtf")
        ann = read_annotation(sheet[g].get("annotation", ""), contigs)
        gff, prot, utr, after = ["##gff-version 3"], [], [], []
        for r in genes:
            contig, strand, gene = r["contig"], r["strand"], r["gene"]
            codon0 = int(r["sec_codon"]) - 1
            gs = genome[contig]
            models = []
            for x in alns:
                if x["t"] != contig or not (x["ts"] <= codon0 and codon0 + 3 <= x["te"]): continue
                m = transcript_model(x, txseq[x["q"]], gs, codon0, strand)
                if m: models.append(m)
            usable = [m for m in models if m["cds"]]
            rank = lambda m: (not m["open5"] and not m["open3"], -(m["cds_insertions"] + m["cds_deletions"]),
                              len(m["cds"]), len(m["utr3"]))
            best = max(usable, key=rank) if usable else prediction_model(gtf.get(contig, {}), gs, codon0, strand)
            row = {"genome": g, "gene": gene, "family": r["family"], "contig": contig, "sec_codon": r["sec_codon"],
                   "strand": strand, "baseline_status": r["status"], "supporting_transcripts": len(models)}
            if best is None:
                row.update({"evidence": "none", "notes": "no supporting transcript and no prediction holding the codon"})
                report.append(row); continue
            same_stop = [m for m in usable if best.get("stop_genome") is not None and m.get("stop_genome") == best["stop_genome"]]
            ends = [m["tx_end_genome"] for m in same_stop if m.get("tx_end_genome") is not None]
            utr3_max = max((abs(e - best["stop_genome"]) for e in ends), default="")
            ctx = r.get("context", "").split("*")[0]
            notes = []
            if best["open5"]: notes.append("5' end open (start may lie upstream)")
            if best["open3"]: notes.append("3' end open (no stop codon in the transcript)")
            if best["cds_insertions"] or best["cds_deletions"]:
                notes.append("indel between transcript and assembly in the coding sequence: protein from the transcript")
            noncanon = [s for s in best["sites"] if s not in ("GT-AG", "GC-AG", "AT-AC")]
            lo_cds, hi_cds = min(best["cds_genome"]), max(best["cds_genome"])
            status, ids = annotation_status(ann, contig, strand, codon0, lo_cds, hi_cds)
            row.update({"evidence": best["tx"], "protein_length": len(best["protein"]),
                        "complete": "yes" if not (best["open5"] or best["open3"]) else "no",
                        "cds_start": (lo_cds + 1) if strand == "+" else (hi_cds + 1),
                        "cds_end": (hi_cds + 1) if strand == "+" else (lo_cds + 1),
                        "exons": len(best["blocks"]), "introns": len(best["introns"]),
                        "non_canonical_introns": ",".join(noncanon) or "-",
                        "utr5_length": best.get("utr5_len", ""), "utr3_length": len(best["utr3"]),
                        "utr3_max_extent": utr3_max,
                        "polyA_signal": "yes" if re.search("AATAAA|ATTAAA", best["utr3"][-60:].upper()) else "no",
                        "cds_insertions": best["cds_insertions"], "cds_deletions": best["cds_deletions"],
                        "cds_mismatches": best["cds_mismatches"],
                        "context_match": ("yes" if ctx.replace("U", "X") in best["protein"].replace("U", "X") else "no") if ctx else "-",
                        "annotation_status": status, "annotated_proteins": ids, "notes": "; ".join(notes) or "-"})
            report.append(row)
            gff += ["\t".join(map(str, x)) for x in gff_lines(g, gene, contig, strand, best, (codon0, codon0 + 2))]
            prot.append(f">{gene}|{g}|{best['tx']}\n{best['protein']}\n")
            if best["utr3"]: utr.append(f">{gene}|{g}|{best['tx']}|{len(best['utr3'])} nt\n{best['utr3']}\n")
            # genomic sequence after the stop codon (or after the end of an open frame)
            if strand == "+":
                p = hi_cds + 1; seq = gs[p:p + a.after_stop]
            else:
                p = lo_cds; seq = rc(gs[max(0, p - a.after_stop):p])
            after.append(f">{gene}|{g}|after_stop|{len(seq)} nt\n{seq}\n")
        with open(f"{a.out}/{g}.gff3", "w") as o: o.write("\n".join(gff) + "\n")
        for suffix, items in (("proteins.faa", prot), ("utr3.fa", utr), ("after_stop.fa", after)):
            with open(f"{a.out}/{g}.{suffix}", "w") as o: o.write("".join(items))
        print(g, len(genes), "genes", flush=True)
    cols = ["genome", "gene", "family", "contig", "sec_codon", "strand", "baseline_status", "evidence",
            "supporting_transcripts", "protein_length", "complete", "cds_start", "cds_end", "exons", "introns",
            "non_canonical_introns", "utr5_length", "utr3_length", "utr3_max_extent", "polyA_signal", "cds_insertions",
            "cds_deletions", "cds_mismatches", "context_match", "annotation_status", "annotated_proteins", "notes"]
    os.makedirs(os.path.dirname(a.table) or ".", exist_ok=True)
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in report: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(open(a.table).read())


if __name__ == "__main__":
    main()
