#!/usr/bin/env python3
"""Step 2 of the plan (docs/plan.md, 2026-10-24): gold-standard gene models of the curated O. dioica selenoprotein
genes in OKI2018, OSKA2016 and OdB3.

For each curated selenocysteine codon (curation/baseline_sec_genes.tsv):
  1. Transcripts. The transcripts of the assembly (stage1/samplesheet.csv: the OKI2018 Trinity transcriptome, the
     OSKA2016 transcript assembly, the OdB3 cDNA collection) are aligned to it with minimap2 (-x splice:hq -un, one
     primary alignment each; cached in --work). -un turns off the GT-AG preference of minimap2, which otherwise
     creates artificial indels next to the atypical introns of O. dioica (GA-AG 9%, GC-AG 2.7% and GG-AG 0.7% of the
     introns, against 87.6% GT-AG, and no AT-AC or U12 intron: Denoeud et al. 2010, supplementary section 10.1;
     GA-AG also in Frey and Pucker 2020); junctions are then placed by the script: each intron is moved within
     +/- --window nt to the placement of lowest cost, the number of mismatches between transcript and genome in the
     window plus the cost of the splice site, log10(f(GT-AG) / f(site)) / 2 with the frequencies above (other sites
     taken as 0.1%), 2 being the log10 cost of a mismatch at about 1% divergence (GA-AG 0.49, GC-AG 0.76, GG-AG
     1.05, other sites 1.47), then the smallest move. Every exon chain is also moved, among the placements that
     give the same spliced sequence, to the most frequent class. O. dioica introns are short (peak at 40-44 nt, some
     below 35 nt; Denoeud et al. 2010, figure S22); minimap2 -x splice:hq aligns introns of 20 nt or more as
     introns and shorter ones as deletions (tested with minimap2 2.26), so deletions of 10 nt or more that have an
     equivalent placement on one of these sites are reported as possible introns. A
     transcript supports the gene when an aligned exon holds the three bases of the codon and they read TGA on the
     gene strand.
  2. Candidate models. Each candidate is an exon chain on the genome: the chain of each supporting transcript, of
     each annotated mRNA that holds the codon in an exon (exon rows, or CDS rows when the annotation has no exons),
     and of the stage 1 Selenoprofiles prediction that holds the codon (<stage1>/<genome>.selenoprofiles.gtf). The
     genome sequence is spliced along the chain, so that the coding sequence, the protein and the GFF3 rows of a
     model always agree with the assembly. In this sequence the frame of the TGA is extended upstream to the
     previous stop codon (TAA, TAG or TGA) and the coding sequence starts at the first ATG after it; downstream it
     ends at the next stop codon. The TGA is read as selenocysteine (U). A frame that reaches the end of the chain
     is open at that end (for a chain of CDS rows, the 5' end is closed when it starts with ATG); an open 5' end
     without ATG before the TGA (a fragment, for example a 3' EST) starts at the first codon of the chain.
  3. Choice. Candidates are ranked by: the protein context of the curated codon found in the protein (the part up
     to the U or the part from the U, since the curated context was read on contiguous genome); a protein at
     least 90% as long as the longest candidate; a complete frame (start and stop); at least one supporting
     transcript compatible with the chain (same introns where they overlap, allowing equivalent placements, and no
     transcript end inside an intron of the chain or the reverse, beyond 10 nt); an annotated chain (so that the
     public annotation is changed as little as possible); protein length; 3' UTR length.
  4. Assembly or transcript. The frame of the TGA is also read in the sequence of each supporting transcript, and
     the best transcript protein is chosen by the same criteria. When its transcript differs from the assembly by
     indels in the coding sequence and the two proteins differ, each protein is aligned (global, BLOSUM62, free end
     gaps) to the proteins of the same gene in the other assemblies (assembly and transcript proteins), and the
     transcript protein is used only when its best identity (identical residues over the residues facing the other
     protein, see protein_identity) exceeds that of
     the assembly protein by more than 0.02 (without another assembly: when it is more than 10% longer, or holds
     the curated context when the assembly protein does not). The transcript protein is also used when no
     candidate on the assembly gives a frame. Then the assembly does not encode the protein: the model is the exon
     chain of the transcript with the genome positions of its coding sequence, marked in the GFF3 with the
     attribute exception (as RefSeq models annotated from transcripts), and protein_source is "transcript".
Reported per gene: the protein context check; introns, the atypical ones (GA-AG, GC-AG, GG-AG) and the other splice
sites (listed with all their equivalent placements, the one used first, as "CA-TT/AG-TC/GG-CA"; such sites point to
an alignment, assembly or transcript problem); insertions, deletions and mismatches between
the assembly and the compatible transcript that covers most of the coding sequence; the 3' UTR of the model and the
longest 3' UTR of the supporting transcripts that end at the same stop codon; how the public annotation treats the
gene (annotation_status: selenocysteine declared, codon inside the coding sequence without declaration, coding
sequence ending before the codon or starting after it, no model); and the change that the annotation needs
(annotation_fix).
Outputs in --out: <genome>.gff3 (gene, mRNA, exon and CDS features; transl_except=(pos:..%2Caa:Sec) on the CDS rows),
<genome>.proteins.faa (U for selenocysteine), <genome>.utr3.fa (3' UTR of the model), <genome>.after_stop.fa
(--after-stop nucleotides of genome after the stop codon, on the gene strand), and --table.
"""
import argparse, collections, csv, gzip, math, os, re, subprocess

STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
# Intron classes of the O. dioica reference annotation (Denoeud et al. 2010, supplementary section 10.1): GT-AG 87.6%,
# GA-AG 9%, GC-AG 2.7%, GG-AG 0.7%; no AT-AC or U12 intron (no minor spliceosome). Other sites: at most 0.1%.
INTRON_CLASSES = {"GT-AG": 87.6, "GA-AG": 9.0, "GC-AG": 2.7, "GG-AG": 0.7}
OTHER_SITES = 0.1
MISMATCH_LOG10 = 2.0   # log10 cost of a mismatch at about 1% divergence between transcript and assembly
SITE_COST = {k: math.log10(INTRON_CLASSES["GT-AG"] / f) / MISMATCH_LOG10 for k, f in INTRON_CLASSES.items()}
OTHER_COST = math.log10(INTRON_CLASSES["GT-AG"] / OTHER_SITES) / MISMATCH_LOG10
END_TOLERANCE = 10


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
# splice sites and junction placement

def site(genome, a, b, strand):
    """Splice site of the intron [a, b) on the gene strand, as donor-acceptor dinucleotides."""
    s = f"{genome[a:a + 2]}-{genome[b - 2:b]}" if strand == "+" else f"{rc(genome[b - 2:b])}-{rc(genome[a:a + 2])}"
    return s.upper()


def equivalent_sites(genome, a, b, strand, limit=50):
    """Splice sites of all placements of the intron [a, b) that give the same spliced sequence (this one first)."""
    left = 0
    while left < limit and a - 1 - left >= 0 and genome[a - 1 - left].upper() == genome[b - 1 - left].upper(): left += 1
    right = 0
    while right < limit and b + right < len(genome) and genome[a + right].upper() == genome[b + right].upper(): right += 1
    out = [site(genome, a, b, strand)]
    for k in range(-left, right + 1):
        s = site(genome, a + k, b + k, strand)
        if s not in out: out.append(s)
    return out


def normalise_blocks(blocks, genome, strand):
    """Exon chain with each intron moved, among the placements that give the same spliced sequence, to the most
    frequent intron class of O. dioica (GT-AG, GA-AG, GC-AG, GG-AG), then to the smallest move."""
    b = [list(x) for x in sorted(blocks)]
    for i in range(len(b) - 1):
        a, e = b[i][1], b[i + 1][0]
        if e <= a: continue
        left = 0
        while left < b[i][1] - b[i][0] - 1 and genome[a - 1 - left] == genome[e - 1 - left]: left += 1
        right = 0
        while right < b[i + 1][1] - b[i + 1][0] - 1 and genome[a + right] == genome[e + right]: right += 1
        k = min((SITE_COST.get(site(genome, a + k, e + k, strand), OTHER_COST), abs(k), k) for k in range(-left, right + 1))[2]
        b[i][1] += k; b[i + 1][0] += k
    return [tuple(x) for x in b]


def parse_cigar(cg):
    ops = []
    for n, op in re.findall(r"(\d+)([MIDN=X])", cg):
        op = "M" if op in "=X" else op
        if ops and ops[-1][1] == op: ops[-1][0] += int(n)
        else: ops.append([int(n), op])
    return ops


def refine_junctions(aln, fwd, genome, strand, window):
    """CIGAR with each intron moved within +/- window nt to the placement of lowest cost: mismatches between the
    transcript (fwd: in the orientation of the genome forward strand) and the genome in the window, plus the cost of
    its splice site (SITE_COST, OTHER_COST), then the smallest move; and the number of moves."""
    ops = parse_cigar(aln["cigar"])
    q = aln["qs"] if aln["strand"] == "+" else aln["qlen"] - aln["qe"]
    t = aln["ts"]
    moved = 0
    for i in range(len(ops)):
        n, op = ops[i]
        if op == "N" and 0 < i < len(ops) - 1 and ops[i - 1][1] == "M" and ops[i + 1][1] == "M":
            w = min(window, ops[i - 1][0] - 1, ops[i + 1][0] - 1)
            if w > 0:
                best = None
                for k in range(-w, w + 1):
                    mm = sum(fwd[x].upper() != genome[(t + x - q) if x < q + k else (t + n + x - q)].upper()
                             for x in range(q - w, q + w))
                    key = (mm + SITE_COST.get(site(genome, t + k, t + k + n, strand), OTHER_COST), abs(k), k)
                    if best is None or key < best: best = key
                k = best[2]
                if k:
                    ops[i - 1][0] += k; ops[i + 1][0] -= k; q += k; t += k; moved += 1
        n, op = ops[i]
        if op == "M": q += n; t += n
        elif op == "I": q += n
        elif op in "DN": t += n
    return "".join(f"{n}{op}" for n, op in ops), moved


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
    (with the genome position they precede) and deleted target positions."""
    q = aln["qs"] if aln["strand"] == "+" else aln["qlen"] - aln["qe"]
    t = aln["ts"]
    g2q, q2g, ins, dels, introns, blocks, del_spans = {}, {}, [], [], [], [], []
    bstart = t
    for n, op in parse_cigar(aln["cigar"]):
        if op == "M":
            for i in range(n): g2q[t + i] = q + i; q2g[q + i] = t + i
            t += n; q += n
        elif op == "I":
            ins += [(q + i, t) for i in range(n)]; q += n
        elif op == "D":
            dels += [t + i for i in range(n)]; del_spans.append((t, t + n)); t += n
        elif op == "N":
            blocks.append((bstart, t)); introns.append((t, t + n)); t += n; bstart = t
    blocks.append((bstart, t))
    return {"g2q": g2q, "q2g": q2g, "ins": ins, "dels": dels, "introns": introns, "blocks": blocks, "del_spans": del_spans}


# ---------------------------------------------------------------------------------------------------------------
# models

def orf(mrna, s):
    """Frame through the TGA at s: (start, stop) indices (stop = first base of the stop codon), and open ends.
    The start is the first ATG after the previous in-frame stop codon or, when there is no previous stop codon
    (5' end open), the first ATG of the frame, or its first codon when it has no ATG before the TGA (a fragment)."""
    i = s - 3
    while i >= 0 and mrna[i:i + 3].upper() not in STOPS: i -= 3
    up_stop = i if i >= 0 else None
    first = (up_stop + 3) if up_stop is not None else s % 3
    start = next((k for k in range(first, s, 3) if mrna[k:k + 3].upper() == "ATG"), None)
    if start is None and up_stop is None: start = first   # 5' end open without ATG: a fragment, from its first codon
    k = s + 3
    while k + 3 <= len(mrna) and mrna[k:k + 3].upper() not in STOPS: k += 3
    stop = k if k + 3 <= len(mrna) else None
    return start, stop, up_stop is None, stop is None


def chain_introns(blocks):
    b = sorted(blocks)
    return [(x[1], y[0]) for x, y in zip(b, b[1:]) if y[0] > x[1]]


def spliced(genome, blocks, strand):
    """Sequence of the exon chain on the gene strand and the genome position of each of its bases."""
    pos = [p for a, b in sorted(blocks) for p in range(a, b)]
    if strand == "-": pos = pos[::-1]
    seq = "".join(genome[p] for p in pos).upper()
    if strand == "-": seq = seq.translate(COMP)
    return seq, pos


def model_from_chain(blocks, genome, codon0, strand, source, kind, cds_only=False):
    """Model read on the genome sequence spliced along an exon chain, or None when the chain does not hold the
    codon as TGA."""
    blocks = normalise_blocks(blocks, genome, strand)
    seq, pos = spliced(genome, blocks, strand)
    idx = {p: i for i, p in enumerate(pos)}
    s = idx.get(codon0) if strand == "+" else idx.get(codon0 + 2)
    if s is None or [idx.get(codon0 + d) for d in ((0, 1, 2) if strand == "+" else (2, 1, 0))] != [s, s + 1, s + 2]:
        return None
    if seq[s:s + 3] != "TGA": return None
    start, stop, open5, open3 = orf(seq, s)
    if cds_only and open5: open5 = seq[start:start + 3] != "ATG"   # a chain of CDS rows that starts with ATG
    m = {"source": source, "kind": kind, "blocks": sorted(blocks), "introns": chain_introns(blocks),
         "cds_only": cds_only, "open5": open5, "open3": open3, "cds": "", "protein": "", "utr3": "", "utr5_len": 0}
    if start is None: return m
    end = (stop + 3) if stop is not None else len(seq) - (len(seq) - s) % 3
    m.update({"cds": seq[start:end], "protein": translate(seq[start:end], s - start).rstrip("*"),
              "utr3": seq[end:] if stop is not None else "", "utr5_len": start, "cds_genome": pos[start:end],
              "stop_genome": pos[end - 1] if stop is not None else None})
    m["sites"] = [equivalent_sites(genome, a, b, strand) for a, b in m["introns"]]
    return m


def transcript_model(aln, tx, genome, codon0, strand):
    """Frame of the TGA read in the transcript sequence, or None when the transcript does not hold the codon as TGA
    in an aligned exon. aln carries the CIGAR after junction placement."""
    w = walk(aln)
    fwd = tx if aln["strand"] == "+" else rc(tx)
    L = len(fwd)
    qs = [w["g2q"].get(codon0 + i) for i in range(3)]
    if None in qs or qs != [qs[0], qs[0] + 1, qs[0] + 2]: return None
    if strand == "+":
        mrna, s = fwd.upper(), qs[0]
        m2g = lambda i: w["q2g"].get(i)
    else:
        mrna, s = rc(fwd).upper(), L - 1 - qs[2]
        m2g = lambda i: w["q2g"].get(L - 1 - i)
    if mrna[s:s + 3] != "TGA": return None
    start, stop, open5, open3 = orf(mrna, s)
    m = {"tx": aln["q"], "walk": w, "fwd": fwd, "blocks": w["blocks"], "open5": open5, "open3": open3,
         "protein": "", "utr3": "", "cds_indels": 0, "stop_genome": None}
    if start is None: return m
    end = (stop + 3) if stop is not None else len(mrna) - (len(mrna) - s) % 3
    m["protein"] = translate(mrna[start:end], s - start).rstrip("*")
    m["utr3"] = mrna[end:] if stop is not None else ""
    m["utr5_len"] = start
    gpos = [m2g(i) for i in range(start, end)]
    mapped = [p for p in gpos if p is not None]
    m["cds_genome"] = mapped
    if mapped:
        lo, hi = min(mapped), max(mapped)
        m["cds_inserted"] = sum(p is None for p in gpos)
        m["cds_deleted"] = sum(lo <= d <= hi for d in w["dels"])
        m["cds_indels"] = m["cds_inserted"] + m["cds_deleted"]
        if stop is not None: m["stop_genome"] = gpos[-1]
    return m


def context_check(context):
    """Test of a protein against the curated context of the codon (curation/baseline_sec_genes.tsv), which was read on
    the contiguous genome sequence and can run into an intron on either side: the protein must hold the part up to
    the U or the part from the U (up to a stop codon), each counted only when at least 4 residues long."""
    if "U" not in context: return lambda p: True
    u = context.index("U")
    sides = [x for x in (context[:u + 1], context[u:].split("*")[0]) if len(x) >= 4]
    return lambda p: not sides or any(x in p for x in sides)


def effective_length(m):
    """Protein length for ranking: from the first methionine when the 5' end of the frame is open (the frame may
    run through a 5' UTR), the whole length otherwise or when there is no methionine."""
    p = m["protein"]
    return len(p) - p.index("M") if m["open5"] and "M" in p else len(p)


def mapped_transcript_model(tm, genome, strand):
    """Model of a transcript whose protein the assembly does not encode: its exon chain and the genome positions of
    its coding sequence, as in RefSeq models annotated from transcripts (GFF3 attribute exception)."""
    exc = f"{tm['cds_inserted']} bases of the transcript coding sequence absent from the assembly, {tm['cds_deleted']} " \
          f"bases of the assembly absent from it"
    return {"source": tm["tx"], "kind": "transcript", "blocks": sorted(tm["blocks"]),
            "introns": chain_introns(tm["blocks"]), "cds_only": False, "open5": tm["open5"], "open3": tm["open3"],
            "protein": tm["protein"], "cds": "", "utr3": tm["utr3"], "utr5_len": tm["utr5_len"],
            "cds_genome": tm["cds_genome"], "stop_genome": tm["stop_genome"], "moved": tm["moved"],
            "sites": [equivalent_sites(genome, a, b, strand) for a, b in chain_introns(tm["blocks"])], "exception": exc}


SHORT_INTRON_CHECK = 10


def short_intron_candidates(tms, genome, strand, lo, hi):
    """Deletions of at least SHORT_INTRON_CHECK nt of transcripts against the assembly, between lo and hi, with an
    equivalent placement on a site of the O. dioica intron classes: possible short introns, which minimap2
    (-x splice:hq) aligns as deletions below about 20 nt (O. dioica introns peak at 40-44 nt, some below 35 nt:
    Denoeud et al. 2010, figure S22). Returns {(start, end): (site, [transcripts])}."""
    out = {}
    for tm in tms:
        for a, b in tm["walk"]["del_spans"]:
            if b - a < SHORT_INTRON_CHECK or b <= lo or a > hi: continue
            sites = [x for x in equivalent_sites(genome, a, b, strand) if x in INTRON_CLASSES]
            if sites: out.setdefault((a, b), (sites[0], []))[1].append(tm["tx"])
    return out


def same_intron(i1, i2, genome):
    """Same intron, allowing placements that give the same spliced sequence."""
    (a1, b1), (a2, b2) = sorted([i1, i2])
    if b1 - a1 != b2 - a2: return False
    d = a2 - a1
    return d == 0 or (d <= END_TOLERANCE and genome[a1:a2].upper() == genome[b1:b2].upper())


def compatible(t_blocks, c_blocks, genome):
    """A transcript chain and a candidate chain agree where they overlap: same introns, and no end of one inside an
    intron of the other (beyond END_TOLERANCE nt)."""
    t0, t1 = min(a for a, _ in t_blocks), max(b for _, b in t_blocks)
    c0, c1 = min(a for a, _ in c_blocks), max(b for _, b in c_blocks)
    o0, o1 = max(t0, c0), min(t1, c1)
    if o1 <= o0: return False
    ti = [x for x in chain_introns(t_blocks) if x[0] >= o0 and x[1] <= o1]
    ci = [x for x in chain_introns(c_blocks) if x[0] >= o0 and x[1] <= o1]
    if len(ti) != len(ci) or not all(same_intron(x, y, genome) for x, y in zip(ti, ci)): return False
    for introns, (e0, e1) in ((chain_introns(c_blocks), (t0, t1)), (chain_introns(t_blocks), (c0, c1))):
        for a, b in introns:
            if a < e0 < b and b - e0 > END_TOLERANCE: return False
            if a < e1 < b and e1 - a > END_TOLERANCE: return False
    return True


def cds_differences(model, tm, genome):
    """Bases of the model's coding sequence covered by the transcript, and the mismatches, insertions and deletions
    of the transcript against the assembly there."""
    g2q, fwd, w = tm["walk"]["g2q"], tm["fwd"], tm["walk"]
    cov = [p for p in model["cds_genome"] if p in g2q]
    if not cov: return 0, 0, 0, 0
    lo, hi = min(cov), max(cov)
    inside = lambda p: any(a <= p < b for a, b in w["blocks"])     # not in an intron of the transcript
    mism = sum(fwd[g2q[p]].upper() != genome[p].upper() for p in cov)
    dels = sum(lo <= p <= hi and p not in g2q and inside(p) for p in model["cds_genome"])
    ins = sum(lo < t <= hi for _, t in w["ins"])
    return len(cov), mism, ins, dels


def read_gtf(path):
    """Selenoprofiles GTF: contig -> prediction -> [(start0, end0 inclusive, strand)] of its CDS rows."""
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
    """Annotated transcripts on the contigs: id -> {"contig", "strand", "exons": [(s0, e0 exclusive)],
    "cds": [(s0, e0 inclusive)], "protein", "sec"}; rows are grouped by Parent."""
    out = {}
    if not path or not os.path.exists(path): return out
    with opener(path) as fh:
        for l in fh:
            if l.startswith("#"): continue
            f = l.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in contigs or f[2] not in ("CDS", "exon"): continue
            at = attributes(f[8])
            parent = at.get("Parent", at.get("ID", "?")).split(",")[0]
            r = out.setdefault(parent, {"contig": f[0], "strand": f[6], "exons": [], "cds": [], "protein": "", "sec": False})
            if f[2] == "exon":
                r["exons"].append((int(f[3]) - 1, int(f[4])))
            else:
                r["cds"].append((int(f[3]) - 1, int(f[4]) - 1))
                r["protein"] = r["protein"] or at.get("protein_id", "")
                if "Sec" in at.get("transl_except", ""): r["sec"] = True
    for k, r in out.items(): r["protein"] = r["protein"] or k
    return out


def annotation_status(ann, contig, strand, codon0, lo, hi):
    """How the annotation treats the gene: status and the overlapping annotated proteins."""
    over = []
    for mid, r in ann.items():
        if r["contig"] != contig or r["strand"] != strand or not r["cds"]: continue
        a, b = min(x for x, _ in r["cds"]), max(y for _, y in r["cds"])
        if b < lo or a > hi: continue
        over.append((mid, r, a, b))
    if not over: return "no annotated coding sequence", "-", []
    ids = ",".join(sorted({r["protein"] for _, r, _, _ in over}))
    statuses = []
    for _, r, a, b in over:
        inside = any(x <= codon0 and codon0 + 2 <= y for x, y in r["cds"])
        if inside: statuses.append("selenocysteine declared" if r["sec"] else "codon inside the coding sequence, not declared")
        elif (strand == "+" and b < codon0) or (strand == "-" and a > codon0 + 2): statuses.append("coding sequence ends before the codon")
        else: statuses.append("coding sequence starts after the codon")
    order = ["selenocysteine declared", "codon inside the coding sequence, not declared",
             "coding sequence ends before the codon", "coding sequence starts after the codon"]
    return next(s for s in order if s in statuses), ids, over


def cds_segments(model):
    """Coding exons of a model, (start0, end0 inclusive), sorted."""
    lo, hi = min(model["cds_genome"]), max(model["cds_genome"])
    return [(max(a, lo), min(b - 1, hi)) for a, b in model["blocks"] if b - 1 >= lo and a <= hi]


def annotation_fix(model, over, genome):
    """Change that the public annotation needs to give the model."""
    if not over: return "new gene model"
    segs = cds_segments(model)
    for _, r, _, _ in over:
        if sorted(r["cds"]) == segs:
            return "none (selenocysteine declared)" if r["sec"] else "declare the selenocysteine (transl_except)"
    if model["kind"] == "annotation":
        return f"change the coding sequence on the annotated exon chain of {model['source']}"
    for _, r, _, _ in over:
        if compatible(r["exons"] or [(x, y + 1) for x, y in r["cds"]], model["blocks"], genome):
            return f"change the coding sequence of {r['protein']} (compatible exon chain, to be extended)"
    return "new exon chain (overlapping " + ",".join(sorted({r["protein"] for _, r, _, _ in over})) + ")"


# ---------------------------------------------------------------------------------------------------------------

def gff_lines(g, gene, contig, strand, m, sec_genome):
    """GFF3 rows of one model (0-based inputs, 1-based output)."""
    cds = cds_segments(m)
    exons = [(a, b - 1) for a, b in m["blocks"]]
    if m["cds_only"]: exons = cds
    glo, ghi = min(a for a, _ in exons), max(b for _, b in exons)
    gid, rid = f"gene-{g}_{gene}", f"rna-{g}_{gene}"
    note = f"evidence={m['kind']}:{m['source']}"
    if m.get("exception"):
        note += f";exception=annotated by transcript or proteomic data;Note={m['exception'].replace(',', '%2C')}"
    rows = [(contig, "gold_standard", "gene", glo + 1, ghi + 1, ".", strand, ".", f"ID={gid};Name={gene}"),
            (contig, "gold_standard", "mRNA", glo + 1, ghi + 1, ".", strand, ".", f"ID={rid};Parent={gid};Name={gene};{note}")]
    for k, (a, b) in enumerate(sorted(exons), 1):
        rows.append((contig, "gold_standard", "exon", a + 1, b + 1, ".", strand, ".", f"ID=exon-{g}_{gene}-{k};Parent={rid}"))
    order = cds if strand == "+" else cds[::-1]
    done = 0
    te = f"transl_except=(pos:{sec_genome[0] + 1}..{sec_genome[1] + 1}%2Caa:Sec)" if strand == "+" else \
        f"transl_except=(pos:complement({sec_genome[0] + 1}..{sec_genome[1] + 1})%2Caa:Sec)"
    phased = []
    for a, b in order:
        phased.append((a, b, (3 - done % 3) % 3)); done += b - a + 1
    for a, b, ph in sorted(phased):
        rows.append((contig, "gold_standard", "CDS", a + 1, b + 1, ".", strand, str(ph), f"ID=cds-{g}_{gene};Parent={rid};{te}"))
    return rows


def protein_identity(q, t, aligner):
    """Identity of q to t in a global alignment with free end gaps (BLOSUM62, U read as C): identical residues over
    the residues of q that face t, that is, those in aligned blocks and internal gaps, and those of an end overhang
    of q only where t also continues (an overhang of q beyond the ends of t, for example a fragment, is not counted)."""
    q2, t2 = q.replace("U", "C"), t.replace("U", "C")
    if not q2 or not t2: return 0.0
    blocks_q, blocks_t = aligner.align(q2, t2)[0].aligned
    if len(blocks_q) == 0: return 0.0
    same = sum(q2[i] == t2[j] for (qa, qb), (ta, tb) in zip(blocks_q, blocks_t) for i, j in zip(range(qa, qb), range(ta, tb)))
    first_q, last_q = blocks_q[0][0], blocks_q[-1][1]
    facing = last_q - first_q
    if blocks_t[0][0] > 0: facing += first_q            # t continues before the first aligned block
    if blocks_t[-1][1] < len(t2): facing += len(q2) - last_q   # t continues after the last aligned block
    return same / facing


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
    ap.add_argument("--window", type=int, default=8, help="junction placement window (nt on each side)")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True); os.makedirs(a.out, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    base = [r for r in csv.DictReader(open(a.baseline), delimiter="\t") if r["genome"] in a.genomes]
    # pass 1: candidate models of each gene
    records, genomes, annotations = [], {}, {}
    for g in a.genomes:
        genes = [r for r in base if r["genome"] == g]
        if not genes: continue
        contigs = {r["contig"] for r in genes}
        genome = genomes[g] = {k: v.upper() for k, v in read_fasta(sheet[g]["fasta"], contigs)}
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
        ann = annotations[g] = read_annotation(sheet[g].get("annotation", ""), contigs)
        for r in genes:
            contig, strand = r["contig"], r["strand"]
            codon0 = int(r["sec_codon"]) - 1
            gs = genome[contig]
            has_ctx = context_check(r.get("context", ""))
            # supporting transcripts, with junctions placed
            tms = []
            for x in alns:
                if x["t"] != contig or not (x["ts"] <= codon0 and codon0 + 3 <= x["te"]): continue
                fwd = txseq[x["q"]] if x["strand"] == "+" else rc(txseq[x["q"]])
                cg, moved = refine_junctions(x, fwd, gs, strand, a.window)
                tm = transcript_model(dict(x, cigar=cg), txseq[x["q"]], gs, codon0, strand)
                if tm: tm["moved"] = moved; tms.append(tm)
            # candidate models on the assembly
            cands = []
            for tm in tms:
                m = model_from_chain(tm["blocks"], gs, codon0, strand, tm["tx"], "transcript")
                if m: m["moved"] = tm["moved"]; cands.append(m)
            for mid, ar in ann.items():
                if ar["contig"] != contig or ar["strand"] != strand: continue
                blocks = ar["exons"] or [(x, y + 1) for x, y in ar["cds"]]
                if not (min(x for x, _ in blocks) <= codon0 and codon0 + 3 <= max(y for _, y in blocks)): continue
                m = model_from_chain(blocks, gs, codon0, strand, ar["protein"], "annotation", cds_only=not ar["exons"])
                if m: cands.append(m)
            for pid, cds in gtf.get(contig, {}).items():
                if cds[0][2] != strand or not (min(x for x, _, _ in cds) <= codon0 <= max(y for _, y, _ in cds)): continue
                m = model_from_chain([(x, y + 1) for x, y, _ in cds], gs, codon0, strand, pid, "prediction", cds_only=True)
                if m: cands.append(m)
            usable = [m for m in cands if m["cds"]]
            for m in usable:
                m["ctx"] = has_ctx(m["protein"])
                m["support"] = [tm for tm in tms if compatible(tm["blocks"], m["blocks"], gs)]
            best = None
            if usable:
                pool = [m for m in usable if m["ctx"]] or usable
                longest = max(effective_length(m) for m in pool)
                rank = lambda m: (m["ctx"], effective_length(m) >= 0.9 * longest, not (m["open5"] or m["open3"]),
                                  bool(m["support"]), m["kind"] == "annotation", effective_length(m), len(m["utr3"]))
                best = max(usable, key=rank)
            tx_ok = [tm for tm in tms if tm["protein"]]
            tbest = max(tx_ok, key=lambda t: (has_ctx(t["protein"]), not (t["open5"] or t["open3"]), effective_length(t))) \
                if tx_ok else None
            records.append({"g": g, "r": r, "tms": tms, "candidates": len(usable), "best": best, "tbest": tbest,
                            "has_ctx": has_ctx})
    # pass 2: assembly or transcript protein, outputs
    from Bio import Align
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                    open_gap_score=-10, extend_gap_score=-0.5)
    aligner.target_end_gap_score = aligner.query_end_gap_score = 0.0
    proteins = collections.defaultdict(list)
    for rec in records:
        for m in (rec["best"], rec["tbest"]):
            if m and m["protein"]: proteins[rec["r"]["gene"]].append((rec["g"], m["protein"]))
    out = {g: {"gff": ["##gff-version 3"], "prot": [], "utr": [], "after": []} for g in genomes}
    report = []
    for rec in records:
        g, r, tms, best, tbest, has_ctx = rec["g"], rec["r"], rec["tms"], rec["best"], rec["tbest"], rec["has_ctx"]
        contig, strand, gene = r["contig"], r["strand"], r["gene"]
        codon0 = int(r["sec_codon"]) - 1
        gs = genomes[g][contig]
        row = {"genome": g, "gene": gene, "family": r["family"], "contig": contig, "sec_codon": r["sec_codon"],
               "strand": strand, "baseline_status": r["status"], "transcripts_with_codon": len(tms),
               "candidates": rec["candidates"]}
        if best is None and tbest is None:
            row.update({"evidence": "none", "notes": "no transcript, annotated mRNA or prediction holds the codon in a frame"})
            report.append(row); continue
        source, notes, id_g, id_t = "genome", [], "", ""
        if tbest:
            use_tx = False
            if best is None:
                use_tx = True
            elif tbest["cds_indels"] and tbest["protein"] != best["protein"]:
                others = [p for og, p in proteins[gene] if og != g]
                if others:
                    id_g = max(protein_identity(best["protein"], p, aligner) for p in others)
                    id_t = max(protein_identity(tbest["protein"], p, aligner) for p in others)
                    use_tx = id_t > id_g + 0.02
                    if not use_tx:
                        notes.append(f"{tbest['tx']} differs from the assembly by indels in the coding sequence; its protein "
                                     f"is less similar to the other assemblies ({id_t:.2f} against {id_g:.2f}): protein "
                                     f"from the assembly")
                else:
                    use_tx = effective_length(tbest) > 1.1 * effective_length(best) or \
                        (has_ctx(tbest["protein"]) and not best["ctx"])
            elif effective_length(tbest) > 1.1 * effective_length(best):
                notes.append(f"{tbest['tx']} encodes a longer protein ({len(tbest['protein'])} aa) without indel against "
                             f"the assembly (other exon chain)")
            if use_tx:
                source = "transcript"
                best = mapped_transcript_model(tbest, gs, strand)
                best["ctx"] = has_ctx(best["protein"])
                best["support"] = [tm for tm in tms if compatible(tm["blocks"], best["blocks"], gs)]
                why = f"; identity to the other assemblies {id_t:.2f} against {id_g:.2f}" if id_t != "" else ""
                notes.append(f"the assembly does not encode the protein of {tbest['tx']} ({best['exception']}{why}): "
                             f"model mapped from the transcript, protein from the transcript")
        protein = best["protein"]
        lo_cds, hi_cds = min(best["cds_genome"]), max(best["cds_genome"])
        if source == "genome":   # internal check: the GFF3 coding sequence, read on the genome, gives the protein
            cds_seq = "".join(gs[x:y + 1] for x, y in cds_segments(best))
            cds_seq = cds_seq if strand == "+" else rc(cds_seq)
            sec = best["cds_genome"].index(codon0 if strand == "+" else codon0 + 2)
            if translate(cds_seq, sec).rstrip("*") != protein:
                raise SystemExit(f"internal check failed: {g} {gene}: the GFF3 coding sequence does not give the protein")
        # differences between the assembly and the supporting transcripts
        diffs = [(cds_differences(best, tm, gs), tm["tx"]) for tm in best["support"]]
        (cov, mism, ins, dels), dtx = max(diffs, key=lambda d: (d[0][0], -sum(d[0][1:]))) if diffs else ((0, 0, 0, 0), "-")
        indel_tx = sorted({tm["tx"] for tm in best["support"] if cds_differences(best, tm, gs)[2:] != (0, 0)})
        if indel_tx and source == "genome":
            notes.append(f"{len(indel_tx)} supporting transcript(s) differ from the assembly by indels in the coding "
                         f"sequence ({','.join(indel_tx[:3])}{'...' if len(indel_tx) > 3 else ''})")
        for (x, y), (site_, txs) in sorted(short_intron_candidates(tms, gs, strand, lo_cds, hi_cds).items()):
            notes.append(f"possible intron of {y - x} nt ({site_}) aligned as a deletion at {x + 1}-{y} "
                         f"({','.join(sorted(set(txs))[:3])})")
        if best["open5"]: notes.append("5' end open (start may lie upstream)")
        if best["open3"]: notes.append("3' end open (no stop codon in the chain)")
        if not best["support"]: notes.append("no supporting transcript compatible with the exon chain")
        if not best["ctx"]: notes.append("curated context not found in the protein")
        atypical = collections.Counter(x[0] for x in best["sites"] if x[0] in INTRON_CLASSES and x[0] != "GT-AG")
        other_sites = ["/".join(x) for x in best["sites"] if x[0] not in INTRON_CLASSES]
        same_stop = [tm for tm in tms if best["stop_genome"] is not None and tm["stop_genome"] == best["stop_genome"]]
        status, ids, over = annotation_status(annotations[g], contig, strand, codon0, lo_cds, hi_cds)
        fix = annotation_fix(best, over, gs)
        if source != "genome": fix += "; the assembly does not encode the protein of the transcripts"
        row.update({"evidence": f"{best['kind']}:{best['source']}", "supporting_transcripts": len(best["support"]),
                    "protein_length": len(protein), "protein_source": source,
                    "identity_assembly_protein": f"{id_g:.3f}" if id_g != "" else "",
                    "identity_transcript_protein": f"{id_t:.3f}" if id_t != "" else "",
                    "complete": "yes" if not (best["open5"] or best["open3"]) else "no",
                    "cds_start": (lo_cds + 1) if strand == "+" else (hi_cds + 1),
                    "cds_end": (hi_cds + 1) if strand == "+" else (lo_cds + 1),
                    "exons": len(best["blocks"]), "introns": len(best["introns"]),
                    "atypical_introns": ",".join(f"{k} x{v}" if v > 1 else k for k, v in sorted(atypical.items())) or "-",
                    "other_splice_sites": ",".join(other_sites) or "-",
                    "junctions_moved": best.get("moved", "-"),
                    "utr5_length": best["utr5_len"] if not best["cds_only"] else "",
                    "utr3_length": len(best["utr3"]) if not best["cds_only"] else "",
                    "utr3_max_length": max((len(tm["utr3"]) for tm in same_stop), default=""),
                    "polyA_signal": "yes" if re.search("AATAAA|ATTAAA", best["utr3"][-60:]) else "no",
                    "difference_transcript": dtx, "cds_bases_covered": cov, "cds_mismatches": mism,
                    "cds_insertions": ins, "cds_deletions": dels,
                    "context_match": ("yes" if has_ctx(protein) else "no") if r.get("context") else "-",
                    "annotation_status": status, "annotated_proteins": ids, "annotation_fix": fix,
                    "notes": "; ".join(notes) or "-"})
        report.append(row)
        o = out[g]
        o["gff"] += ["\t".join(map(str, x)) for x in gff_lines(g, gene, contig, strand, best, (codon0, codon0 + 2))]
        o["prot"].append(f">{gene}|{g}|{best['kind']}:{best['source']}|protein_source={source}|"
                         f"complete={'yes' if not (best['open5'] or best['open3']) else 'no'}\n{protein}\n")
        if best["utr3"]:
            o["utr"].append(f">{gene}|{g}|{best['kind']}:{best['source']}|{len(best['utr3'])} nt\n{best['utr3']}\n")
        # genomic sequence after the stop codon (or after the end of an open frame)
        if strand == "+":
            p = hi_cds + 1; seq = gs[p:p + a.after_stop]
        else:
            p = lo_cds; seq = rc(gs[max(0, p - a.after_stop):p])
        o["after"].append(f">{gene}|{g}|after_{'stop' if not best['open3'] else 'open_frame'}|{len(seq)} nt\n{seq}\n")
    for g, o in out.items():
        with open(f"{a.out}/{g}.gff3", "w") as fh: fh.write("\n".join(o["gff"]) + "\n")
        for suffix, key in (("proteins.faa", "prot"), ("utr3.fa", "utr"), ("after_stop.fa", "after")):
            with open(f"{a.out}/{g}.{suffix}", "w") as fh: fh.write("".join(o[key]))
    cols = ["genome", "gene", "family", "contig", "sec_codon", "strand", "baseline_status", "evidence",
            "transcripts_with_codon", "candidates", "supporting_transcripts", "protein_length", "protein_source",
            "identity_assembly_protein", "identity_transcript_protein",
            "complete", "cds_start", "cds_end", "exons", "introns", "atypical_introns", "other_splice_sites",
            "junctions_moved",
            "utr5_length", "utr3_length", "utr3_max_length", "polyA_signal", "difference_transcript",
            "cds_bases_covered", "cds_mismatches", "cds_insertions", "cds_deletions", "context_match",
            "annotation_status", "annotated_proteins", "annotation_fix", "notes"]
    os.makedirs(os.path.dirname(a.table) or ".", exist_ok=True)
    with open(a.table, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in report: fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(open(a.table).read())


if __name__ == "__main__":
    main()
