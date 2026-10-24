#!/usr/bin/env python3
"""Step 4 of the plan (docs/plan.md, 2026-10-24): UGA/SECIS screens, independent of protein families.

scan (heavy; Infernal cmsearch):
  - genome-wide SECIS scans with the Oikopleura SECIS model of step 3 (results/models/oikopleura_secis.cm), both
    strands, hits of at least --report-bits, of the --scan genomes (six Oikopleura assemblies and four non-target
    tunicates) and of a shuffled copy of each, as background: every contig cut in windows of --shuffle-window nt,
    each shuffled with its dinucleotide counts preserved (Altschul and Erickson 1985), so that local composition
    stays as in the genome;
  - transcript sets: the annotated transcripts of the genomes that have an annotation (exon chain, or CDS rows when
    the annotation has no exon; extended on the genome by --extend3 nt after their 3' end, and by --extend5 nt before
    the 5' end of CDS-only models, so that a SECIS beyond an annotated 3' end or a selenocysteine codon before an
    annotated start is still seen), searched on their sense strand; and the transcript assemblies of the samplesheet
    (Trinity, transcripts, ESTs), searched on both strands. O. dioica sets with the Oikopleura model; the controls
    (human and C. intestinalis RefSeq transcripts) with the Oikopleura model and with Rfam RF00031.
screen (light; reads the outputs of scan):
  - threshold of genome scans: the lowest integer score at which the hits of the shuffled copies of the six
    Oikopleura assemblies, pooled, give at most --genome-fp expected false hits in the largest of them (both strands);
  - transcript screen: on each transcript strand with a SECIS hit of at least --tx-bits (20 bits, step 3), every TGA
    upstream of the SECIS is read as selenocysteine in a frame that starts at the first ATG after the previous stop
    codon of the frame (TAA, TAG or TGA; or at the first ATG of the frame when there is none: "5' open"; or at the
    first codon when there is no ATG either: "5' fragment") and ends at the next stop codon after the TGA, which must
    end before the SECIS; frames shorter than --min-codons are left out; "internal" when at least --internal codons
    follow the TGA, "C-terminal" otherwise (as SELENOK, SELENOS and TXNRD); every TGA gives at most one frame, ranked
    by length for each SECIS hit;
  - comparative screen (selenocysteine/selenocysteine pairs between species, as SGP2, Castellano et al. 2004): the
    --upstream nt before each genome-scan hit at or above the threshold, on its strand, compared by tblastx with
    those of the other species (the three O. dioica assemblies count as one species); a TGA aligned to a TGA in an
    alignment of E <= --cmp-evalue on the forward frames of both regions, with at least --cmp-flank aligned residues
    on one side of it, unless a stronger alignment of the same two regions overlaps it in another frame; its open
    stretch (codons to the previous and next stop codons) must hold at least --cmp-min-codons codons; "internal" or
    "C-terminal" as above; one candidate per TGA, from the nearest SECIS hit downstream;
  - each candidate protein (translated with U, from the frame of the transcript screen or the open stretch of the
    comparative screen) is compared by blastp (legacy BLAST) with the human RefSeq proteome (best hit, E <= --evalue, and
    the residue aligned to the U in a local alignment with the hit, BLOSUM62) and with the gold-standard proteins of
    step 2; known genes: the gold-standard models of the assembly (overlap of the TGA with their selenocysteine codon
    when genome positions are known, otherwise protein identity >= --known-identity with the U aligned to their U),
    the stage 1 predictions (overlap of the TGA with their coding sequence: selenocysteine or cysteine homologs, as
    labelled by the first round of hit validation; for the comparative screen also at the partner TGA) and, for the
    controls, the selenocysteines of the annotation
    (transl_except). O. dioica candidates are grouped across assemblies by blastp (identity >= --group-identity,
    E <= 1e-10).
Outputs: work in --work; tables stage2_secis_scan_summary.tsv (hits per genome, real and shuffled, by score),
stage2_secis_scan_hits.tsv (genome-scan hits at or above the threshold, with the nearest known gene),
stage2_uga_screen.tsv (one row per candidate) and stage2_uga_screen_summary.tsv (per genome and source: sequences,
SECIS hits, candidates, known genes recovered).
"""
import argparse, collections, concurrent.futures, csv, gzip, os, random, re, subprocess, zlib

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = {"TAA", "TAG", "TGA"}
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
        for k, c in enumerate(BASES)}
OIKOPLEURA = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan"]
DIOICA = ["OKI", "OSKA", "OdB3"]
NONTARGET = ["Cint", "Csav", "Scla", "Tdem"]
CONTROLS = ["Hsap", "Cint"]
SPECIES = {"OKI": "O. dioica", "OSKA": "O. dioica", "OdB3": "O. dioica"}
species = lambda g: SPECIES.get(g, g)
SOURCE = {"OKI": ("RefSeq", "Trinity"), "OSKA": ("ANISEED", "transcripts"), "OdB3": ("GenBank", "ESTs"),
          "Hsap": ("RefSeq", None), "Cint": ("RefSeq", None)}


def rc(s): return s.translate(COMP)[::-1]


def opener(path): return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    name, buf = None, []
    with opener(path) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = l[1:].split()[0], []
            else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def attributes(field):
    out = {}
    for x in field.strip().strip(";").split(";"):
        x = x.strip()
        if not x: continue
        if "=" in x: k, v = x.split("=", 1)
        else: k, _, v = x.partition(" ")
        out[k] = v.strip('"')
    return out


def translate(s):
    """Translation with TGA read as U except as the last codon."""
    aa = [CODE.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3)]
    return "".join("U" if (a == "*" and s[3 * k:3 * k + 3].upper() == "TGA") else a for k, a in enumerate(aa))


def run(cmd, cwd=None):
    return subprocess.run(cmd, check=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout


def dinucleotide_shuffle(seq, rng):
    """Random sequence with the dinucleotide counts and the first and last letters of seq (Altschul and Erickson
    1985: Eulerian walk with a random spanning arborescence of last exits); same as stage2/secis_model_training.py."""
    s = seq.upper()
    if len(s) < 3: return s
    edges = collections.defaultdict(list)
    for x, y in zip(s, s[1:]): edges[x].append(y)
    last = s[-1]
    while True:
        exit_ = {v: rng.choice(e) for v, e in edges.items() if v != last}
        ok = True
        for v in edges:
            seen, u = set(), v
            while u != last:
                if u in seen or u not in exit_: ok = False; break
                seen.add(u); u = exit_[u]
            if not ok: break
        if ok: break
    lists = {}
    for v, e in edges.items():
        e = list(e)
        if v in exit_:
            e.remove(exit_[v]); rng.shuffle(e); e.append(exit_[v])
        else: rng.shuffle(e)
        lists[v] = e
    out, u, idx = [s[0]], s[0], collections.Counter()
    for _ in range(len(s) - 1):
        v = lists[u][idx[u]]; idx[u] += 1; out.append(v); u = v
    return "".join(out)


# ---------------------------------------------------------------------------------------------------------------
# transcript sets

def read_gff_transcripts(path):
    """Annotated transcripts: id -> {"contig", "strand", "blocks" [(s0, e0 excl)], "cds_only", "cds" [(s0, e0 excl)],
    "sec" [genomic 0-based first base of each annotated selenocysteine codon], "gene", "protein"}; exon rows give the
    chain, or the CDS rows when the transcript has no exon row."""
    out = {}
    with opener(path) as fh:
        for l in fh:
            if l.startswith("#"): continue
            f = l.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("exon", "CDS"): continue
            at = attributes(f[8])
            parent = at.get("Parent", at.get("ID", "?")).split(",")[0]
            r = out.setdefault(parent, {"contig": f[0], "strand": f[6], "exons": [], "cds": [], "sec": set(),
                                        "gene": at.get("gene", ""), "protein": ""})
            a, b = int(f[3]) - 1, int(f[4])
            if f[2] == "exon": r["exons"].append((a, b))
            else:
                r["cds"].append((a, b))
                r["protein"] = r["protein"] or at.get("protein_id", "")
                for x, y in re.findall(r"pos:(?:complement\()?(\d+)\.\.(\d+)\)?,aa:Sec", at.get("transl_except", "").replace("%2C", ",")):
                    r["sec"].add(min(int(x), int(y)) - 1)
    for k, r in out.items():
        r["cds_only"] = not r["exons"]
        r["blocks"] = sorted(r["exons"] or r["cds"])
        r["cds"] = sorted(r["cds"])
        r["sec"] = sorted(r["sec"])
        del r["exons"]
    return out


def tx_blocks(r, e5, e3, contig_len):
    """Genomic blocks (0-based half-open) of a transcript in its orientation, with the extensions."""
    bl = sorted(r["blocks"])
    if r["strand"] == "+":
        out = [(max(0, bl[0][0] - e5), bl[0][0])] + bl + [(bl[-1][1], min(contig_len, bl[-1][1] + e3))]
    else:
        out = [(bl[-1][1], min(contig_len, bl[-1][1] + e5))] + bl[::-1] + [(max(0, bl[0][0] - e3), bl[0][0])]
    return [(a, b) for a, b in out if b > a]


def tx_to_genome(blocks, strand, t):
    """Genomic position of position t of a transcript (blocks in transcript orientation)."""
    off = 0
    for a, b in blocks:
        if t < off + b - a: return a + (t - off) if strand == "+" else b - 1 - (t - off)
        off += b - a
    return None


def genome_to_tx(blocks, strand, p):
    off = 0
    for a, b in blocks:
        if a <= p < b: return off + (p - a if strand == "+" else b - 1 - p)
        off += b - a
    return None


def build_annotation_set(genome_fa, gff, out_fa, out_tsv, ext5, ext3, mrna_only=False):
    """Transcript sequences of an annotation (in transcript orientation, with extensions) and a table of their
    coordinates: genomic blocks in transcript order, annotated CDS ends and selenocysteine codons in transcript
    coordinates (0-based)."""
    tx = read_gff_transcripts(gff)
    if mrna_only: tx = {k: v for k, v in tx.items() if v["cds"]}
    by_contig = collections.defaultdict(list)
    for k, r in tx.items(): by_contig[r["contig"]].append(k)
    n = 0
    with open(out_fa, "w") as fa, open(out_tsv, "w") as tsv:
        tsv.write("id\tcontig\tstrand\tlength\tblocks\textension5\textension3\tcds_only\tcds_start\tcds_end\tsec\tgene\tprotein\n")
        for contig, seq in read_fasta(genome_fa):
            if contig not in by_contig: continue
            seq = seq.upper()
            for k in by_contig[contig]:
                r = tx[k]
                e5 = ext5 if r["cds_only"] else 0
                bl = tx_blocks(r, e5, ext3, len(seq))
                s = "".join(seq[a:b] if r["strand"] == "+" else rc(seq[a:b]) for a, b in bl)
                ends = [genome_to_tx(bl, r["strand"], p) for a, b in r["cds"] for p in (a, b - 1)]
                ends = sorted(x for x in ends if x is not None)
                sec = [genome_to_tx(bl, r["strand"], p if r["strand"] == "+" else p + 2) for p in r["sec"]]
                sec = sorted(x for x in sec if x is not None)
                fa.write(f">{k}\n{s}\n")
                tsv.write(f"{k}\t{contig}\t{r['strand']}\t{len(s)}\t{';'.join(f'{a}-{b}' for a, b in bl)}\t{e5}\t{ext3}\t"
                          f"{'yes' if r['cds_only'] else 'no'}\t{ends[0] if ends else ''}\t{ends[-1] if ends else ''}\t"
                          f"{','.join(map(str, sec))}\t{r['gene']}\t{r['protein']}\n")
                n += 1
    return n


def copy_transcripts(src, out_fa):
    """Transcript assembly with simple names (first word), upper case."""
    n = 0
    with open(out_fa, "w") as o:
        for name, s in read_fasta(src):
            o.write(f">{name}\n{s.upper()}\n"); n += 1
    return n


# ---------------------------------------------------------------------------------------------------------------
# scan

def cmsearch(cm, fasta, tbl, cpus, bits, toponly=False):
    cmd = ["cmsearch", "--cpu", str(cpus), "-T", str(bits), "--tblout", tbl + ".tmp", "-o", os.devnull]
    if toponly: cmd.insert(1, "--toponly")
    subprocess.run(cmd + [cm, fasta], check=True, stdout=subprocess.DEVNULL)
    os.replace(tbl + ".tmp", tbl)


def write_genome(src, dest, shuffle_window=0, seed=0):
    """Uncompressed copy of a genome (Infernal reads it faster), or its window-shuffled copy; returns its size."""
    rng = random.Random(seed)
    total = 0
    with open(dest, "w") as o:
        for name, s in read_fasta(src):
            s = s.upper()
            if shuffle_window:
                s = "".join(dinucleotide_shuffle(s[i:i + shuffle_window], rng) for i in range(0, len(s), shuffle_window))
            o.write(f">{name}\n")
            for i in range(0, len(s), 100): o.write(s[i:i + 100] + "\n")
            total += len(s)
    return total


def scan_genome(g, kind, a, sheet, cpus):
    d = f"{a.work}/scan"
    tbl = f"{d}/{g}.{kind}.tbl"
    if os.path.exists(tbl): return g, kind, None
    fa = f"{d}/{g}.{kind}.fa"
    size = write_genome(sheet[g]["fasta"], fa, a.shuffle_window if kind == "shuffled" else 0, a.shuffle_seed + zlib.crc32(g.encode()))
    with open(f"{d}/{g}.{kind}.size", "w") as o: o.write(f"{size}\n")
    cmsearch(a.model, fa, tbl, cpus, a.report_bits)
    os.remove(fa)
    return g, kind, size


def scan(a):
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    os.makedirs(f"{a.work}/scan", exist_ok=True); os.makedirs(f"{a.work}/tx", exist_ok=True)
    jobs = [(g, kind) for g in a.scan for kind in ("genome", "shuffled")]
    per = max(1, a.cpus // a.jobs)
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
        for g, kind, size in ex.map(lambda j: scan_genome(j[0], j[1], a, sheet, per), jobs):
            print(f"scan {g} {kind}: {'done before' if size is None else f'{size / 1e6:.1f} Mb'}", flush=True)
    # transcript sets
    tasks = []
    for g in a.transcripts:
        ann, asm = SOURCE[g]
        ann_path = sheet[g].get("annotation", "")
        if ann and ann_path and os.path.exists(ann_path):
            fa, tsv = f"{a.work}/tx/{g}.{ann}.fa", f"{a.work}/tx/{g}.{ann}.tsv"
            if not os.path.exists(tsv):
                n = build_annotation_set(sheet[g]["fasta"], ann_path, fa + ".tmp", tsv + ".tmp", a.extend5, a.extend3,
                                         mrna_only=g in CONTROLS)
                os.replace(fa + ".tmp", fa); os.replace(tsv + ".tmp", tsv)
                print(f"transcripts {g} {ann}: {n}", flush=True)
            for m in (["oikopleura", "RF00031"] if g in CONTROLS else ["oikopleura"]):
                tasks.append((fa, f"{a.work}/tx/{g}.{ann}.{m}.tbl", m, True))
        tx_path = sheet[g].get("transcripts", "")
        if asm and tx_path and os.path.exists(tx_path):
            fa = f"{a.work}/tx/{g}.{asm}.fa"
            if not os.path.exists(fa):
                n = copy_transcripts(tx_path, fa + ".tmp"); os.replace(fa + ".tmp", fa)
                print(f"transcripts {g} {asm}: {n}", flush=True)
            tasks.append((fa, f"{a.work}/tx/{g}.{asm}.oikopleura.tbl", "oikopleura", False))
    models = {"oikopleura": a.model, "RF00031": a.rfam}
    todo = [t for t in tasks if not os.path.exists(t[1])]
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
        list(ex.map(lambda t: cmsearch(models[t[2]], t[0], t[1], per, a.report_bits, toponly=t[3]), todo))
    print(f"transcript searches: {len(todo)} run, {len(tasks) - len(todo)} done before", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# screen

def read_tbl(path, min_bits):
    """cmsearch hits: (target, strand, from, to, bits, evalue), 1-based positions as reported (from > to on -)."""
    out = []
    if not os.path.exists(path): return out
    for l in open(path):
        if l.startswith("#"): continue
        f = l.split()
        if float(f[14]) >= min_bits:
            out.append((f[0], f[9], int(f[7]), int(f[8]), float(f[14]), float(f[15])))
    return out


def strand_coords(length, strand, x, y):
    """0-based half-open interval of a hit on the strand sequence (the reverse complement for -)."""
    if strand == "+": return x - 1, y
    return length - x, length - y + 1


def sec_frames(s, h0, min_codons, internal):
    """Frames through a TGA read as selenocysteine (one per TGA), ending before position h0 of the strand sequence s,
    longest first."""
    out = []
    u = s.find("TGA", 0, max(0, h0 - 2))
    while u != -1 and u + 3 <= h0:
        i = u - 3
        while i >= 0 and s[i:i + 3] not in STOPS: i -= 3
        open5 = i < 0
        first = i + 3 if not open5 else u % 3
        start = next((k for k in range(first, u, 3) if s[k:k + 3] == "ATG"), None)
        kind = "ATG" if start is not None and not open5 else ("ATG, 5' open" if start is not None else "5' fragment")
        if start is None and open5: start = first
        if start is not None:
            k = u + 3
            while k + 3 <= len(s) and s[k:k + 3] not in STOPS: k += 3
            if k + 3 <= h0:
                nb, na = (u - start) // 3, (k - u - 3) // 3
                if nb + 1 + na >= min_codons:
                    out.append({"start": start, "tga": u, "stop": k, "before": nb, "after": na, "codons": nb + 1 + na,
                                "start_kind": kind, "class": "internal" if na >= internal else "C-terminal",
                                "distance": h0 - (k + 3), "protein": translate(s[start:k])})
        u = s.find("TGA", u + 1, max(0, h0 - 2))
    # a frame is left out when it overlaps, in another reading frame, a longer frame already kept or a longer
    # ordinary open reading frame (ATG to stop, TGA as stop) of the strand: two coding frames cannot overlap
    orfs = ordinary_orfs(s, min_codons) if out else []
    out.sort(key=lambda c: -c["codons"])
    kept = []
    for c in out:
        span, fr = (c["start"], c["stop"] + 3), c["tga"] % 3
        clash = lambda o: o[0] < span[1] and span[0] < o[1] and o[0] % 3 != fr
        if any(clash((k["start"], k["stop"] + 3)) for k in kept): continue
        if any(clash((a0, b0)) and (b0 - a0) // 3 > c["codons"] for a0, b0 in orfs): continue
        kept.append(c)
    for i, c in enumerate(kept, 1): c["rank"] = i
    return kept


def ordinary_orfs(s, min_codons):
    """Open reading frames of a strand sequence, from the first ATG after a stop codon (TGA as stop) to the next stop
    codon (included), of at least min_codons codons: [(start, end)]."""
    out = []
    for f in range(3):
        start = None
        for i in range(f, len(s) - 2, 3):
            c = s[i:i + 3]
            if start is None and c == "ATG": start = i
            elif start is not None and c in STOPS:
                if (i - start) // 3 >= min_codons: out.append((start, i + 3))
                start = None
    return out


def longest_orf(s):
    """Length in codons of the longest ATG-to-stop frame of a strand sequence (TGA as stop)."""
    best = 0
    for f in range(3):
        start = None
        for i in range(f, len(s) - 2, 3):
            c = s[i:i + 3]
            if start is None and c == "ATG": start = i
            elif start is not None and c in STOPS:
                best = max(best, (i - start) // 3); start = None
    return best


def open_stretch(s, u):
    """Codons before (to the previous stop codon or the start of s) and after (to the next stop codon, excluded) a TGA
    at u, and the position of that stop codon (None when the end of s is reached)."""
    i = u - 3
    while i >= 0 and s[i:i + 3] not in STOPS: i -= 3
    k = u + 3
    while k + 3 <= len(s) and s[k:k + 3] not in STOPS: k += 3
    return (u - i - 3) // 3, (k - u - 3) // 3, (k if k + 3 <= len(s) else None), i + 3


def comparative_screen(a, summary, T, sheet):
    """Selenocysteine/selenocysteine pairs between species: the --upstream nt before each genome-scan hit at or above
    T (on the strand of the hit) are compared by tblastx (legacy BLAST, no filter) with those of the other species;
    a TGA aligned to a TGA (stop/stop column) in an alignment of E <= --cmp-evalue on the forward frames of both
    regions, with at least --cmp-flank aligned residues on one side, is a candidate, unless a stronger alignment
    between the same regions overlaps it in another frame (the other frames of a conserved coding sequence). Returns
    the candidates (one per TGA and genome, from the nearest SECIS hit downstream, with its best partner) and the
    numbers of hits and of regions with an alignment per genome."""
    import xml.etree.ElementTree as ET
    d = f"{a.work}/comparative"
    os.makedirs(d, exist_ok=True)
    regions, stats = {}, collections.Counter()
    for g in a.genome_screen:
        rec = next((x for x in summary if x["genome"] == g and x["set"] == "genome"), None)
        if rec is None or T is None: continue
        by_contig = collections.defaultdict(list)
        for h in rec["hits"]:
            if h[4] >= T: by_contig[h[0]].append(h)
        stats[(g, "hits")] = sum(len(v) for v in by_contig.values())
        for contig, seq in read_fasta(sheet[g]["fasta"]):
            if contig not in by_contig: continue
            seq = seq.upper()
            for _, strand, x, y, bits, ev in by_contig[contig]:
                ss = seq if strand == "+" else rc(seq)
                h0, h1 = strand_coords(len(seq), strand, x, y)
                lo = max(0, h0 - a.upstream)
                regions[f"r{len(regions) + 1}"] = {"genome": g, "contig": contig, "strand": strand, "x": x, "y": y,
                                                   "bits": bits, "evalue": ev, "lo": lo, "h0": h0, "seq": ss[lo:h0],
                                                   "length": len(seq)}
    if not regions: return [], stats
    with open(f"{d}/regions.fa", "w") as o:
        for k, r in regions.items(): o.write(f">{k}\n{r['seq']}\n")
    run(["formatdb", "-i", "regions.fa", "-p", "F"], d)
    run(["blastall", "-p", "tblastx", "-d", "regions.fa", "-i", "regions.fa", "-e", str(a.cmp_evalue), "-F", "F",
         "-m", "7", "-v", "100", "-b", "100", "-a", str(a.cpus), "-o", "tblastx.xml"], d)
    best, partnered = {}, set()
    for _, it in ET.iterparse(f"{d}/tblastx.xml"):
        if it.tag != "Iteration": continue
        q = it.findtext("Iteration_query-def").split()[0]
        for hit in it.iter("Hit"):
            t = hit.findtext("Hit_def").split()[0]
            if species(regions[t]["genome"]) == species(regions[q]["genome"]): continue
            hsps = [h for h in hit.iter("Hsp") if int(h.findtext("Hsp_query-frame")) > 0 and int(h.findtext("Hsp_hit-frame")) > 0]
            spans = [(int(h.findtext("Hsp_query-frame")), int(h.findtext("Hsp_query-from")), int(h.findtext("Hsp_query-to")),
                      float(h.findtext("Hsp_bit-score"))) for h in hsps]
            for hsp, (fr, q0, q1, sc) in zip(hsps, spans):
                # an alignment in another frame of a conserved coding sequence: a stronger alignment between the same
                # regions overlaps it in another query frame
                if any(f2 != fr and a2 <= q1 and q0 <= b2 and s2 > sc for f2, a2, b2, s2 in spans): continue
                ev = float(hsp.findtext("Hsp_evalue"))
                qs, hs = hsp.findtext("Hsp_qseq"), hsp.findtext("Hsp_hseq")
                qf, hf = int(hsp.findtext("Hsp_query-from")) - 1, int(hsp.findtext("Hsp_hit-from")) - 1
                partnered.add(q)
                aligned = [c for c in range(len(qs)) if qs[c] != "-" and hs[c] != "-"]
                iq = ih = 0
                for c in range(len(qs)):
                    if qs[c] == "*" and hs[c] == "*":
                        qn, hn = qf + 3 * iq, hf + 3 * ih
                        before = sum(1 for k in aligned if k < c); after = sum(1 for k in aligned if k > c)
                        if (regions[q]["seq"][qn:qn + 3] == "TGA" and regions[t]["seq"][hn:hn + 3] == "TGA"
                                and max(before, after) >= a.cmp_flank):
                            key = (q, qn)
                            if key not in best or ev < best[key][2]: best[key] = (t, hn, ev, before, after)
                    if qs[c] != "-": iq += 1
                    if hs[c] != "-": ih += 1
        it.clear()
    for k in partnered: stats[(regions[k]["genome"], "partnered")] += 1
    # one candidate per TGA of a genome: the one of the nearest SECIS hit downstream (regions of neighbouring hits
    # overlap)
    gpos = lambda reg, n: reg["lo"] + n if reg["strand"] == "+" else reg["length"] - (reg["lo"] + n) - 3
    site = {}
    for (q, qn), v in best.items():
        r = regions[q]
        k = (r["genome"], r["contig"], r["strand"], gpos(r, qn))
        if k not in site or len(r["seq"]) - qn < len(regions[site[k][0]]["seq"]) - site[k][1]: site[k] = (q, qn)
    best = {qk: best[qk] for qk in site.values()}
    out = []
    for (q, qn), (t, hn, ev, before, after) in sorted(best.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        r, p = regions[q], regions[t]
        left, right, stop, first = open_stretch(r["seq"], qn)
        if left + 1 + right < a.cmp_min_codons: continue
        end = stop if stop is not None else first + 3 * ((len(r["seq"]) - first) // 3)
        out.append({"screen": "comparative", "genome": r["genome"], "source": "comparative", "model": "oikopleura",
                    "sequence": r["contig"], "strand": r["strand"], "secis_from": r["x"], "secis_to": r["y"],
                    "bits": r["bits"], "evalue": r["evalue"], "rank": "",
                    "class": "internal" if right >= a.internal else "C-terminal", "start_kind": "open stretch",
                    "before": left, "after": right, "codons": left + 1 + right,
                    "distance": (r["h0"] - r["lo"] - (stop + 3)) if stop is not None else "", "longest_orf": "",
                    "contig": r["contig"], "genomic_strand": r["strand"], "genomic_tga": gpos(r, qn),
                    "protein": translate(r["seq"][first:end]),
                    "partner": f"{p['genome']} {p['contig']}:{gpos(p, hn) + 1} ({p['strand']}), E={ev:.1e}, "
                               f"{before}/{after} aligned residues before/after",
                    "partner_site": (p["genome"], p["contig"], p["strand"], gpos(p, hn))})
    return out, stats


def transcript_frames(name, s, hits, meta, a):
    """Candidates of one transcript: the frames through a TGA before each SECIS hit (strongest hit first; a frame
    already found before a stronger hit is not repeated), with the genomic position of the TGA when the transcript
    comes from an annotation (meta: its row in the coordinates table)."""
    out, seen = [], set()
    for _, strand, x, y, bits, ev in sorted(hits, key=lambda h: -h[4]):
        ss = s if strand == "+" else rc(s)
        h0, h1 = strand_coords(len(s), strand, x, y)
        frames = [c for c in sec_frames(ss, h0, a.min_codons, a.internal) if (strand, c["tga"]) not in seen]
        lo = longest_orf(ss) if frames else ""
        for c in frames:
            seen.add((strand, c["tga"]))
            c.update({"screen": "transcript", "sequence": name, "strand": strand, "secis_from": x, "secis_to": y,
                      "bits": bits, "evalue": ev, "longest_orf": lo, "contig": "", "genomic_tga": ""})
            if meta:   # genomic position of the TGA (first base, lowest coordinate) and strand
                bl = [tuple(map(int, b.split("-"))) for b in meta["blocks"].split(";")]
                gp = tx_to_genome(bl, meta["strand"], c["tga"] if meta["strand"] == "+" else c["tga"] + 2)
                sec = [int(v) for v in meta["sec"].split(",") if v]
                c.update({"contig": meta["contig"], "genomic_strand": meta["strand"],
                          "genomic_tga": gp if gp is not None else "", "gene": meta["gene"],
                          "annotated_sec": ("yes" if c["tga"] in sec else "other position") if sec else "no"})
            out.append(c)
    return out


def read_stage1(a, genomes):
    """Stage 1 predictions: genome -> contig -> [(start0, end0 excl, strand, prediction, label)] (labels of the first
    round of hit validation)."""
    labels = {}
    if os.path.exists(a.hits):
        for r in csv.DictReader(open(a.hits), delimiter="\t"):
            labels[(r["genome"], r["prediction"])] = r["label"] if r["decision"] == "kept" else r["label"] + " (left out)"
    out = {}
    for g in genomes:
        p = f"{a.stage1}/{g}.selenoprofiles.gtf"
        rows = collections.defaultdict(list)
        if os.path.exists(p):
            for l in open(p):
                f = l.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS": continue
                at = attributes(f[8])
                pid = at.get("transcript_id") or at.get("gene_id")
                rows[(f[0], pid, f[6])].append((int(f[3]) - 1, int(f[4])))
        by = collections.defaultdict(list)
        for (c, pid, st), cds in rows.items():
            lab = labels.get((g, pid), pid.rsplit(".", 1)[-1])
            by[c].append((min(x for x, _ in cds), max(y for _, y in cds), st, pid, lab, cds))
        out[g] = by
    return out


def read_gold(a):
    """Gold-standard models: genome -> list of (gene, contig, strand, CDS rows (0-based half-open), sec codon first
    base 0-based genomic (lowest coordinate)), and the proteins (genome.gene -> sequence)."""
    out, prot = collections.defaultdict(list), {}
    for g in DIOICA:
        p = f"{a.gold}/{g}.gff3"
        if not os.path.exists(p): continue
        rows = collections.defaultdict(lambda: {"cds": []})
        for l in open(p):
            f = l.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS": continue
            gene = re.search(r"Parent=rna-([^;]+)", f[8]).group(1).split("_", 1)[1]
            r = rows[gene]; r["contig"], r["strand"] = f[0], f[6]
            r["cds"].append((int(f[3]) - 1, int(f[4])))
            m = re.search(r"pos:(?:complement\()?(\d+)\.\.(\d+)", f[8])
            if m: r["sec"] = min(int(m.group(1)), int(m.group(2))) - 1
        for gene, r in rows.items():
            out[g].append((gene, r["contig"], r["strand"], sorted(r["cds"]), r.get("sec")))
        fa = f"{a.gold}/{g}.proteins.faa"
        if os.path.exists(fa):
            for name, s in read_fasta(fa): prot[f"{g}.{name.split('|')[0]}"] = s
    return out, prot


class Blast:
    """blastp (legacy BLAST) of candidate proteins against a protein set; best hit and the residue aligned to U."""

    def __init__(self, seqs, work, name):
        from Bio import Align
        from Bio.Align import substitution_matrices
        self.seqs = seqs
        self.dir = f"{work}/blast_{name}"
        os.makedirs(self.dir, exist_ok=True)
        with open(f"{self.dir}/db.faa", "w") as o:
            for k, s in seqs.items(): o.write(f">{k}\n{s.replace('*', '')}\n")
        run(["formatdb", "-i", "db.faa", "-p", "T"], self.dir)
        self.al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                        open_gap_score=-11, extend_gap_score=-1)

    def search(self, queries, evalue, cpus):
        """queries: id -> protein; returns id -> (hit, evalue, identity %, aligned length, residue aligned to the
        first U, identity over the alignment computed locally)."""
        if not queries: return {}
        with open(f"{self.dir}/q.faa", "w") as o:
            for k, s in queries.items(): o.write(f">{k}\n{s.replace('U', 'X').rstrip('*')}\n")
        best = {}
        for l in run(["blastall", "-p", "blastp", "-d", "db.faa", "-i", "q.faa", "-e", str(evalue), "-m", "8",
                      "-a", str(cpus), "-v", "5", "-b", "5"], self.dir).splitlines():
            f = l.split("\t")
            if f[0] not in best: best[f[0]] = (f[1], float(f[10]), float(f[2]), int(f[3]))
        out = {}
        for k, (hit, ev, ident, alen) in best.items():
            q, t = queries[k].rstrip("*"), self.seqs[hit].replace("*", "")
            aln = self.al.align(q.replace("U", "X"), t.replace("U", "X"))[0]
            res = "-"
            if "U" in q:
                u = q.index("U")
                for (q0, q1), (t0, t1) in zip(*aln.aligned):
                    if q0 <= u < q1: res = t[t0 + u - q0]; break
            out[k] = (hit, ev, ident, alen, res)
        return out


def overlap_known(contig, strand, gpos, gold_g, stage1_g):
    """Known gene at a genomic TGA position (first base, lowest coordinate): gold-standard selenocysteine codon of the
    assembly, else a stage 1 prediction whose coding sequence holds the TGA."""
    for gene, c, st, cds, sec in gold_g:
        if c == contig and st == strand and sec is not None and sec == gpos: return f"gold standard {gene}"
    for lo, hi, st, pid, lab, cds in stage1_g.get(contig, []):
        if st == strand and lo <= gpos < hi and any(x <= gpos < y for x, y in cds): return f"stage 1 {pid} ({lab})"
    return ""


def near_known(contig, strand, h_lo, h_hi, gold_g, stage1_g, window=3000):
    """Known gene whose coding sequence ends at most window nt before a SECIS hit on the same strand."""
    best = None
    items = [(f"gold standard {gene}", st, cds) for gene, c, st, cds, sec in gold_g if c == contig]
    items += [(f"stage 1 {pid} ({lab})", st, cds) for lo, hi, st, pid, lab, cds in stage1_g.get(contig, [])]
    for name, st, cds in items:
        if st != strand: continue
        if st == "+":
            end = max(y for _, y in cds); d = h_lo - end
        else:
            end = min(x for x, _ in cds); d = end - h_hi
        if 0 <= d <= window and (best is None or d < best[1]): best = (name, d)
    return f"{best[0]} ({best[1]} nt)" if best else ""


def genome_threshold(summary, a):
    """Lowest integer score at which the pooled shuffled hits of the Oikopleura assemblies give at most a.genome_fp
    expected hits in the largest assembly (both strands)."""
    shuffled = [r for r in summary if r["set"] == "shuffled" and r["genome"] in OIKOPLEURA]
    if not shuffled: return None
    mb = sum(r["Mb"] for r in shuffled)
    largest = max(r["Mb"] for r in shuffled)
    scores = sorted(s for r in shuffled for s in r["scores"])
    for t in range(int(a.report_bits), 200):
        n = sum(1 for s in scores if s >= t)
        if n / mb * largest <= a.genome_fp: return t
    return None


def screen(a):
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    os.makedirs(a.tables, exist_ok=True)
    # genome scans
    summary = []
    cuts = [15, 20, 25, 30, 35, 40, 50]
    for g in a.scan:
        for kind in ("genome", "shuffled"):
            tbl, size = f"{a.work}/scan/{g}.{kind}.tbl", f"{a.work}/scan/{g}.{kind}.size"
            if not os.path.exists(tbl): continue
            hits = read_tbl(tbl, a.report_bits)
            mb = 2 * int(open(size).read()) / 1e6
            summary.append({"genome": g, "set": kind, "Mb": mb, "scores": [h[4] for h in hits], "hits": hits})
    T = genome_threshold(summary, a) if a.genome_bits is None else a.genome_bits
    with open(f"{a.tables}/stage2_secis_scan_summary.tsv", "w") as o:
        o.write("genome\tset\tMb (both strands)\t" + "\t".join(f">= {c} bits" for c in cuts) +
                "\thits per Mb at the genome threshold\n")
        for r in summary:
            n = sum(1 for s in r["scores"] if T is not None and s >= T)
            o.write(f"{r['genome']}\t{r['set']}\t{r['Mb']:.1f}\t" + "\t".join(str(sum(1 for s in r['scores'] if s >= c)) for c in cuts) +
                    (f"\t{n / r['Mb']:.3f}\n" if T is not None else "\t-\n"))
        o.write(f"# genome threshold: {T} bits (lowest integer score with at most {a.genome_fp} expected hit(s) of the pooled "
                f"shuffled Oikopleura assemblies in the largest one)\n")
    print(f"genome threshold: {T} bits", flush=True)
    gold, gold_prot = read_gold(a)
    stage1 = read_stage1(a, list(dict.fromkeys(a.scan + a.transcripts + a.genome_screen)))
    cands, refgenes = [], {}
    stats = collections.defaultdict(lambda: collections.Counter())
    # transcript screen
    for g in a.transcripts:
        ann, asm = SOURCE[g]
        for src, both in ((ann, False), (asm, True)):
            if not src: continue
            fa = f"{a.work}/tx/{g}.{src}.fa"
            if not os.path.exists(fa): continue
            meta = {}
            if not both:
                for r in csv.DictReader(open(f"{a.work}/tx/{g}.{src}.tsv"), delimiter="\t"): meta[r["id"]] = r
                genes = {r["gene"] or r["id"] for r in meta.values() if r["sec"]}
                for m in ("oikopleura", "RF00031"): refgenes[(g, src, m)] = len(genes)
            for m in (["oikopleura", "RF00031"] if g in CONTROLS and not both else ["oikopleura"]):
                hits = collections.defaultdict(list)
                for h in read_tbl(f"{a.work}/tx/{g}.{src}.{m}.tbl", a.tx_bits): hits[h[0]].append(h)
                key = (g, src, m)
                for name, s in read_fasta(fa):
                    stats[key]["sequences"] += 1
                    if name not in hits: continue
                    stats[key]["with SECIS"] += 1
                    cs = transcript_frames(name, s, hits[name], meta.get(name), a)
                    for c in cs: c.update({"genome": g, "source": src, "model": m})
                    cands += cs
                    stats[key]["with a frame"] += bool(cs)
    # comparative screen of the genome-scan hits
    cmp_cands, cmp_stats = comparative_screen(a, summary, T, sheet)
    cands += cmp_cands
    for g in a.genome_screen:
        stats[(g, "comparative", "oikopleura")]["SECIS hits"] = cmp_stats[(g, "hits")]
        stats[(g, "comparative", "oikopleura")]["with a frame"] = cmp_stats[(g, "partnered")]
        if g not in DIOICA:   # reference genes: the selenocysteine predictions of stage 1 kept by the hit validation
            refgenes[(g, "comparative", "oikopleura")] = sum(1 for v in stage1.get(g, {}).values() for x in v
                                                             if x[4] == "selenocysteine")
    # known genes, blastp
    for i, c in enumerate(cands): c["id"] = f"cand{i + 1}"
    human = {k: v for k, v in read_fasta(a.human)} if os.path.exists(a.human) else {}
    hb = Blast(human, a.work, "human").search({c["id"]: c["protein"] for c in cands}, a.evalue, a.cpus) if human else {}
    gb = Blast(gold_prot, a.work, "gold").search({c["id"]: c["protein"] for c in cands}, a.evalue, a.cpus) if gold_prot else {}
    for c in cands:
        g = c["genome"]
        known = ""
        if c["genomic_tga"] != "":
            known = overlap_known(c["contig"], c["genomic_strand"], c["genomic_tga"], gold.get(g, []), stage1.get(g, {}))
            if c.get("annotated_sec") == "yes": known = (known + "; " if known else "") + "annotated selenocysteine"
        h = gb.get(c["id"])
        if not known and h and h[0].split(".", 1)[0] == g and h[2] >= a.known_identity and h[4] == "U":
            known = f"gold standard {h[0].split('.', 1)[1]} (protein)"
        c["known"] = known
        if c.get("partner_site"):
            pg, pc, pst, pp = c["partner_site"]
            c["partner_known"] = overlap_known(pc, pst, pp, gold.get(pg, []), stage1.get(pg, {}))
        c["gold_hit"] = f"{h[0]} ({h[2]:.0f}%, U aligned to {h[4]})" if h else ""
        hh = hb.get(c["id"])
        c["human_hit"] = f"{hh[0]} (E={hh[1]:.1e}, {hh[2]:.0f}%, U aligned to {hh[4]})" if hh else ""
    # O. dioica candidates grouped across assemblies
    dio = [c for c in cands if c["genome"] in DIOICA and c["model"] == "oikopleura"]
    groups = group_candidates(dio, a)
    for c in cands: c["assemblies"] = groups.get(c["id"], "")
    cols = ["id", "screen", "genome", "source", "model", "sequence", "strand", "secis_from", "secis_to", "bits",
            "evalue", "rank", "class", "start_kind", "before", "after", "codons", "distance", "longest_orf",
            "contig", "genomic_tga", "known", "gold_hit", "human_hit", "assemblies", "partner", "partner_known", "protein"]
    for c in cands:
        if c["genomic_tga"] != "": c["genomic_tga"] += 1     # 1-based in the table
    with open(f"{a.tables}/stage2_uga_screen.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for c in cands: o.write("\t".join(str(c.get(k, "")) for k in cols) + "\n")
    # summary per genome, source and model
    def gene_of(c):
        k = c["known"].split("; ")[0]
        if k.startswith("gold standard"): return k.split()[2]
        if k == "annotated selenocysteine": return c.get("gene") or c["sequence"]
        return k.split()[2] if k.startswith("stage 1") else k
    with open(f"{a.tables}/stage2_uga_screen_summary.tsv", "w") as o:
        o.write("genome\tsource\tmodel\tsequences\tsequences with a SECIS hit (comparative: hits at the threshold)\t"
                "sequences with a frame (comparative: regions with an alignment to another species)\tcandidates\tinternal\t"
                "C-terminal\tknown\tknown genes recovered\treference genes\tother candidates\n")
        keys = list(stats) + [k for k in dict.fromkeys((c["genome"], c["source"], c["model"]) for c in cands) if k not in stats]
        for k in keys:
            cs = [c for c in cands if (c["genome"], c["source"], c["model"]) == k]
            kn = [c for c in cs if c["known"]]
            ref = len(gold.get(k[0], [])) if k[0] in DIOICA else refgenes.get(k, "")
            o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{stats[k].get('sequences', '')}\t"
                    f"{stats[k].get('with SECIS', stats[k].get('SECIS hits', ''))}\t{stats[k].get('with a frame', '')}\t{len(cs)}\t"
                    f"{sum(c['class'] == 'internal' for c in cs)}\t{sum(c['class'] == 'C-terminal' for c in cs)}\t"
                    f"{len(kn)}\t{len({gene_of(c) for c in kn})}\t{ref}\t{len(cs) - len(kn)}\n")
    print(f"{len(cands)} candidates", flush=True)
    # genome-scan hits at or above the threshold, with the nearest known gene upstream and their candidates
    idx = collections.defaultdict(list)
    for c in cands:
        if c["screen"] == "comparative":
            idx[(c["genome"], c["contig"], c["strand"], c["secis_from"], c["secis_to"])].append(f"{c['id']} ({c['class']})")
    with open(f"{a.tables}/stage2_secis_scan_hits.tsv", "w") as o:
        o.write("genome\tcontig\tstrand\tfrom\tto\tbits\tevalue\tknown gene upstream (distance from its stop codon)\t"
                "candidates of the comparative screen\n")
        for r in summary:
            if r["set"] != "genome" or T is None or r["genome"] not in a.genome_screen + NONTARGET: continue
            for contig, strand, x, y, bits, ev in sorted(r["hits"], key=lambda h: -h[4]):
                if bits < T: continue
                nk = near_known(contig, strand, min(x, y) - 1, max(x, y), gold.get(r["genome"], []), stage1.get(r["genome"], {}))
                o.write(f"{r['genome']}\t{contig}\t{strand}\t{x}\t{y}\t{bits}\t{ev:.2g}\t{nk}\t"
                        f"{', '.join(idx.get((r['genome'], contig, strand, x, y), []))}\n")


def group_candidates(cands, a):
    """Assemblies holding a matching candidate (blastp identity >= group_identity, E <= 1e-10), for each candidate:
    id -> "OKI,OSKA,OdB3"."""
    if not cands: return {}
    seqs = {c["id"]: c["protein"].rstrip("*") for c in cands}
    b = Blast(seqs, a.work, "candidates")
    with open(f"{b.dir}/q.faa", "w") as o:
        for k, s in seqs.items(): o.write(f">{k}\n{s.replace('U', 'X')}\n")
    pairs = collections.defaultdict(set)
    for l in run(["blastall", "-p", "blastp", "-d", "db.faa", "-i", "q.faa", "-e", "1e-10", "-m", "8", "-a",
                  str(a.cpus), "-v", "50", "-b", "50"], b.dir).splitlines():
        f = l.split("\t")
        if float(f[2]) >= a.group_identity: pairs[f[0]].add(f[1])
    genome = {c["id"]: c["genome"] for c in cands}
    return {k: ",".join(sorted({genome[x] for x in pairs[k] | {k}}, key=DIOICA.index)) for k in seqs}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["scan", "screen"])
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--model", default="results/models/oikopleura_secis.cm")
    ap.add_argument("--rfam", default="data/references/RF00031.cm")
    ap.add_argument("--scan", nargs="*", default=OIKOPLEURA + NONTARGET, help="genomes scanned whole (and shuffled)")
    ap.add_argument("--transcripts", nargs="*", default=DIOICA + CONTROLS, help="genomes with transcript sets")
    ap.add_argument("--genome-screen", nargs="*", default=OIKOPLEURA, help="genomes of the comparative screen")
    ap.add_argument("--report-bits", type=float, default=15.0)
    ap.add_argument("--shuffle-window", type=int, default=1000)
    ap.add_argument("--shuffle-seed", type=int, default=600)
    ap.add_argument("--extend3", type=int, default=1000)
    ap.add_argument("--extend5", type=int, default=300)
    ap.add_argument("--tx-bits", type=float, default=20.0)
    ap.add_argument("--genome-bits", type=float, default=None, help="threshold of genome scans (default: from the shuffled copies)")
    ap.add_argument("--genome-fp", type=float, default=1.0)
    ap.add_argument("--min-codons", type=int, default=40)
    ap.add_argument("--internal", type=int, default=10)
    ap.add_argument("--upstream", type=int, default=5000)
    ap.add_argument("--cmp-evalue", type=float, default=1e-5)
    ap.add_argument("--cmp-flank", type=int, default=10)
    ap.add_argument("--cmp-min-codons", type=int, default=20)
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--known-identity", type=float, default=90.0)
    ap.add_argument("--group-identity", type=float, default=80.0)
    ap.add_argument("--human", default="data/genomes/Hsap/protein.faa.gz")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--stage1", default="work/stage1/pipeline/selenoprofiles")
    ap.add_argument("--hits", default="results/tables/stage2_secis_hits.tsv")
    ap.add_argument("--species-map", nargs="*", default=[], help="genome=species pairs that replace the species of a genome (tests)")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--work", default="work/stage2/uga_screen")
    ap.add_argument("--tables", default="results/tables")
    a = ap.parse_args()
    SPECIES.update(dict(x.split("=", 1) for x in a.species_map))
    scan(a) if a.command == "scan" else screen(a)


if __name__ == "__main__":
    main()
