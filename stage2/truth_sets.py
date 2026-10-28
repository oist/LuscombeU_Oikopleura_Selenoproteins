#!/usr/bin/env python3
"""Reference sets of selenoprotein genes for the stage 2 benchmark (stage2/benchmark.py), one table per genome.

Three kinds of reference, by genome:
- gold: the gold-standard models of OKI2018, OSKA2016 and OdB3 (curation/gold_standard/<genome>.gff3): coding
  segments, and the selenocysteine codon of transl_except.
- lift: new assemblies of the same cryptic species (--lift TARGET=SOURCE, default KUM=OKI AOM=OSKA Bar2=OdB3). The
  genomic region of each gold-standard gene of the source (coding sequence and --flank nt on each side) is aligned to
  the target assembly with minimap2 (-x asm20 -c --cs=long); the coding segments and the selenocysteine codon are
  carried over base by base through the best alignment, completed by the other alignments of the region on the same
  target sequence and strand, in the same order and within --window nt. The lifted model is read on the target
  assembly (frame through the lifted codon read as selenocysteine) and classified: complete (codon TGA, coding
  length unchanged, no other stop codon, ends with a stop codon), complete with in-frame indels (coding length
  changed by a multiple of 3, or restored where the source assembly differs from the gold protein), disrupted
  (frameshift or premature stop codon), codon not TGA, partial (a segment end not aligned), or not found. copies:
  other alignments of the region, of at least --copy-identity, that cover the codon (for example uncollapsed
  haplotypes in the target assembly). protein_identity: identity of the lifted protein to the gold-standard protein
  of the source (curation/gold_standard/<source>.proteins.faa; global alignment, over the longer protein).
- reference: genomes without gold standard (--reference, default Olon Oalb Ovan): the stage 1 predictions labelled
  selenocysteine that step 1 kept (results/tables/stage2_secis_hits.tsv, decision kept), with the coding segments
  and selenocysteine codon of the stage 1 GTF (work/stage1/pipeline/selenoprofiles/<genome>.selenoprofiles.gtf),
  and the KKSU sites of the comparative screen (results/tables/stage2_candidates_genomes.tsv; no gene model).
Coordinates: 1-based, inclusive; the selenocysteine codon is given by its lowest coordinate, on both strands (as in
the pipeline's sec_sites tables). Outputs: --out/<genome>.truth.tsv and --table (all genomes).
"""
import argparse, collections, csv, os, re, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uga_secis_screen import read_fasta, rc, CODE, attributes   # noqa: E402

COLS = ["genome", "gene", "profile", "contig", "strand", "sec_start", "cds", "source", "status", "protein_identity",
        "copies", "note"]
GOLD = ["OKI", "OSKA", "OdB3"]


def profile_of(gene):
    return gene.split("_")[0].replace(".", "")


def gold_models(path, genome):
    """gene -> {contig, strand, cds: [(lo, hi)] in genome order, sec: lowest coordinate of the codon}"""
    out = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        gene = re.search(r"Parent=rna-([^;]+)", f[8]).group(1)[len(genome) + 1:]
        m = out.setdefault(gene, {"contig": f[0], "strand": f[6], "cds": [], "sec": None})
        m["cds"].append((int(f[3]), int(f[4])))
        s = re.search(r"pos:(?:complement\()?(\d+)\.\.(\d+)", f[8])
        if s: m["sec"] = int(s.group(1))
    for m in out.values(): m["cds"].sort()
    return out


def read_fasta_headers(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].strip(), []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def cds_text(cds): return ";".join(f"{a}-{b}" for a, b in cds)


def coding(seq, cds, strand):
    s = "".join(seq[a - 1:b] for a, b in cds)
    return s if strand == "+" else rc(s)


def translate(nt, sec_offset=None):
    aa = []
    for i in range(0, len(nt) - 2, 3):
        aa.append("U" if i == sec_offset and nt[i:i + 3].upper() == "TGA" else CODE.get(nt[i:i + 3].upper(), "X"))
    return "".join(aa)


def sec_offset(cds, strand, sec):
    """Offset of the codon in the coding sequence (transcript order)."""
    segs = cds if strand == "+" else cds[::-1]
    off = 0
    for a, b in segs:
        if a <= sec <= b:
            return off + (sec - a if strand == "+" else b - (sec + 2))
        off += b - a + 1
    return None


def identity(x, y):
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-11, extend_gap_score=-1, end_gap_score=0)
    x, y = (s.upper().replace("U", "C").replace("*", "X") for s in (x, y))
    a = al.align(x, y)[0]
    same = sum(1 for (s0, s1), (t0, _) in zip(*a.aligned) for i in range(s1 - s0) if x[s0 + i] == y[t0 + i])
    return same / max(len(x), len(y))


def paf_records(path):
    out = collections.defaultdict(list)
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        tags = dict(t.split(":", 2)[0::2] for t in f[12:])
        out[f[0]].append({"qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "t": f[5],
                          "ts": int(f[7]), "te": int(f[8]), "match": int(f[9]), "alen": int(f[10]),
                          "primary": tags.get("tp") == "P", "cs": tags.get("cs", "")})
    return out


def position_map(r):
    """query offset (0-based, forward strand of the query) -> target position (0-based)"""
    m, t = {}, r["ts"]
    q, step = (r["qs"], 1) if r["strand"] == "+" else (r["qe"] - 1, -1)
    for op, val in re.findall(r"([=:*+\-~])([A-Za-z0-9]+)", r["cs"]):
        if op == "=" or op == ":":
            n = int(val) if op == ":" else len(val)
            for i in range(n): m[q + step * i] = t + i
            q += step * n; t += n
        elif op == "*":
            m[q] = t; q += step; t += 1
        elif op == "+": q += step * len(val)
        elif op == "-": t += len(val)
    return m


def chain(records, window):
    """Best alignment, then the others on the same target sequence and strand, colinear with it and within window."""
    recs = sorted(records, key=lambda r: -r["match"])
    best = recs[0]
    chosen = [best]
    for r in recs[1:]:
        if r["t"] != best["t"] or r["strand"] != best["strand"]: continue
        if r["ts"] > best["te"] + window or r["te"] < best["ts"] - window: continue
        if any(min(r["qe"], c["qe"]) - max(r["qs"], c["qs"]) > 0.1 * (r["qe"] - r["qs"]) for c in chosen): continue
        # colinearity: target order follows query order (+ strand) or reverses it (- strand)
        ok = True
        for c in chosen:
            before_q = r["qe"] <= c["qs"]
            before_t = r["te"] <= c["ts"]
            if (before_q == before_t) != (r["strand"] == "+"): ok = False
        if ok: chosen.append(r)
    pmap = {}
    for r in chosen:
        for k, v in position_map(r).items(): pmap.setdefault(k, v)
    return best, chosen, pmap


def lift(target, source, models, src_seq, src_prot, tgt_fasta, work, flank, window, copy_identity, threads):
    rows = []
    q = f"{work}/{target}.regions.fa"
    regions = {}
    with open(q, "w") as o:
        for gene, m in models.items():
            lo = max(1, m["cds"][0][0] - flank)
            hi = min(len(src_seq[m["contig"]]), m["cds"][-1][1] + flank)
            regions[gene] = lo
            o.write(f">{gene}\n{src_seq[m['contig']][lo - 1:hi]}\n")
    paf = f"{work}/{target}.regions.paf"
    with open(paf, "w") as o:
        subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", str(threads), tgt_fasta, q], stdout=o,
                       stderr=subprocess.DEVNULL, check=True)
    recs = paf_records(paf)
    tgt = None
    for gene, m in models.items():
        row = {"genome": target, "gene": gene, "profile": profile_of(gene), "source": f"lift:{source}", "copies": "",
               "contig": "", "strand": "", "sec_start": "", "cds": "", "protein_identity": "", "note": ""}
        if gene not in recs:
            rows.append(dict(row, status="not found")); continue
        best, chosen, pmap = chain(recs[gene], window)
        lo = regions[gene]
        strand = m["strand"] if best["strand"] == "+" else ("-" if m["strand"] == "+" else "+")
        ends, partial, indel = [], False, False
        for a, b in m["cds"]:
            x, y = pmap.get(a - lo), pmap.get(b - lo)
            if x is None or y is None: partial = True; continue
            x, y = sorted((x + 1, y + 1))
            if y - x != b - a: indel = True
            ends.append((x, y))
        sec = [pmap.get(m["sec"] + i - lo) for i in range(3)]
        sec_ok = all(v is not None for v in sec)
        row.update(contig=best["t"], strand=strand, cds=cds_text(sorted(ends)))
        if sec_ok: row["sec_start"] = min(sec) + 1
        cov = [r for r in recs[gene] if r not in chosen and r["match"] / r["alen"] >= copy_identity
               and r["qs"] <= m["sec"] - lo and r["qe"] >= m["sec"] + 3 - lo]
        row["copies"] = len(cov)
        if partial or not sec_ok:
            rows.append(dict(row, status="partial", note="coding segment end or codon not aligned")); continue
        if tgt is None:
            tgt = dict(read_fasta(tgt_fasta))
        seq = tgt[best["t"]]
        cds = sorted(ends)
        nt = coding(seq, cds, strand)
        off = sec_offset(cds, strand, row["sec_start"])
        codon = nt[off:off + 3].upper() if off is not None else ""
        ref = src_prot.get(gene)
        if ref is None:                          # no gold protein: the source model read on its assembly
            src_nt = coding(src_seq[m["contig"]], m["cds"], m["strand"])
            ref = translate(src_nt, sec_offset(m["cds"], m["strand"], m["sec"]))
        prot = translate(nt, off)
        row["protein_identity"] = f"{identity(prot.rstrip('*'), ref.rstrip('*')):.3f}"
        inner = prot[:-1]
        if codon != "TGA": status = f"codon not TGA ({codon})"
        elif len(nt) % 3 or "*" in inner: status = "disrupted"
        elif not prot.endswith("*"): status = "disrupted"
        elif indel: status = "complete, in-frame indels"
        else: status = "complete"
        notes = []
        if not nt.upper().startswith("ATG"): notes.append("no ATG at the lifted start")
        if len(chosen) > 1: notes.append(f"{len(chosen)} alignments")
        row["note"] = "; ".join(notes)
        rows.append(dict(row, status=status))
    return rows


def stage1_reference(genome, gtf, validated, candidates):
    rows = []
    kept = {r["prediction"] for r in csv.DictReader(open(validated), delimiter="\t")
            if r["genome"] == genome and r["decision"] == "kept" and r["label"] == "selenocysteine"}
    cds, sec, where = collections.defaultdict(list), {}, {}
    if os.path.exists(gtf):
        for l in open(gtf):
            f = l.rstrip("\n").split("\t")
            if len(f) < 9: continue
            gid = attributes(f[8]).get("gene_id", "")
            if f[2] == "Selenocysteine":
                gid = gid.split(":", 1)[-1]
                sec[gid] = min(sec.get(gid, 10 ** 12), int(f[3]))
            elif f[2] == "CDS":
                cds[gid].append((int(f[3]), int(f[4]))); where[gid] = (f[0], f[6])
    else:
        print(f"warning: {gtf} missing; no stage 1 reference for {genome}", file=sys.stderr)
    for gid in sorted(kept):
        if gid not in sec: continue
        rows.append({"genome": genome, "gene": gid, "profile": gid.split(".")[0], "contig": where[gid][0],
                     "strand": where[gid][1], "sec_start": sec[gid], "cds": cds_text(sorted(cds[gid])),
                     "source": "stage1, validated (step 1)", "status": "reference", "protein_identity": "",
                     "copies": "", "note": ""})
    for r in csv.DictReader(open(candidates), delimiter="\t"):
        if r["candidate"] == "KKSU" and r["genome"] == genome and r["codon"] == "TGA" and r["position"]:
            rows.append({"genome": genome, "gene": f"KKSU_{genome}", "profile": "KKSU", "contig": r["contig"],
                         "strand": r["strand"], "sec_start": int(r["position"]), "cds": "",
                         "source": "comparative screen (tblastn)", "status": "reference", "protein_identity": "",
                         "copies": "", "note": "no gene model"})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--gold-genomes", nargs="*", default=GOLD)
    ap.add_argument("--lift", nargs="*", default=["KUM=OKI", "AOM=OSKA", "Bar2=OdB3"])
    ap.add_argument("--reference", nargs="*", default=["Olon", "Oalb", "Ovan"])
    ap.add_argument("--genomes", default="data/genomes", help="directory of <genome>/genome.fna.gz")
    ap.add_argument("--stage1", default="work/stage1/pipeline/selenoprofiles")
    ap.add_argument("--validated", default="results/tables/stage2_secis_hits.tsv")
    ap.add_argument("--candidates", default="results/tables/stage2_candidates_genomes.tsv")
    ap.add_argument("--flank", type=int, default=500)
    ap.add_argument("--window", type=int, default=50000)
    ap.add_argument("--copy-identity", type=float, default=0.95)
    ap.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--out", default="work/stage2/truth")
    ap.add_argument("--table", default="results/tables/stage2_truth_sets.tsv")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    allrows = []

    def write(genome, rows):
        with open(f"{a.out}/{genome}.truth.tsv", "w") as o:
            o.write("\t".join(COLS) + "\n")
            for r in rows: o.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")
        allrows.extend(rows)
        st = collections.Counter(r["status"] for r in rows)
        print(f"{genome}: {len(rows)} genes ({', '.join(f'{k} {v}' for k, v in sorted(st.items()))})", flush=True)

    models = {g: gold_models(f"{a.gold}/{g}.gff3", g) for g in a.gold_genomes}
    for g in a.gold_genomes:
        write(g, [{"genome": g, "gene": gene, "profile": profile_of(gene), "contig": m["contig"], "strand": m["strand"],
                   "sec_start": m["sec"], "cds": cds_text(m["cds"]), "source": "gold standard", "status": "gold",
                   "protein_identity": "", "copies": "", "note": ""} for gene, m in sorted(models[g].items())])
    with tempfile.TemporaryDirectory() as work:
        for pair in a.lift:
            target, source = pair.split("=")
            tf, sf = f"{a.genomes}/{target}/genome.fna.gz", f"{a.genomes}/{source}/genome.fna.gz"
            if not os.path.exists(tf) or not os.path.exists(sf):
                print(f"warning: {tf} or {sf} missing; {target} skipped", file=sys.stderr); continue
            src_seq = dict(read_fasta(sf))
            src_prot = {h.split("|")[0]: q for h, q in read_fasta_headers(f"{a.gold}/{source}.proteins.faa")}
            write(target, lift(target, source, models.get(source) or gold_models(f"{a.gold}/{source}.gff3", source),
                               src_seq, src_prot, tf, work, a.flank, a.window, a.copy_identity, a.threads))
    for g in a.reference:
        write(g, stage1_reference(g, f"{a.stage1}/{g}.selenoprofiles.gtf", a.validated, a.candidates))
    if a.table:
        with open(a.table, "w") as o:
            o.write("\t".join(COLS) + "\n")
            for r in allrows: o.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")


if __name__ == "__main__":
    main()
