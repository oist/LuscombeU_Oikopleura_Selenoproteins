#!/usr/bin/env python3
"""Step 5 of the plan (docs/plan.md, 2026-10-24): candidate selenoprotein genes outside the baseline
(curation/candidates.tsv; proteins in curation/candidates.faa, U for the candidate selenocysteine): TXNL1
(thioredoxin-like protein 1), for which the rule below was fixed in advance, and the family found by the screens of
step 4 (provisional name KKSU, after its C-terminus), to which the same rule is applied.

1. Rule in the three O. dioica assemblies: at the curated position of each assembly, the codon must be TGA, and the
   best SECIS hit of the Oikopleura model (results/models/oikopleura_secis.cm, cmsearch on the top strand) in the
   --window nt after the next stop codon of its frame must reach --secis-bits (20 bits, the threshold of 3' windows
   set in step 3). Counted when this holds in the three assemblies; otherwise reported as a candidate. The genome
   thresholds of step 4 (results/tables/stage2_secis_scan_summary.tsv) are given for comparison.
2. Other genomes (--genomes: the six Oikopleura assemblies, the three held-out appendicularians, the other tunicates
   and amphioxus): tblastn (legacy BLAST, no filter, E <= --evalue) of each candidate protein (U read as X). For each
   candidate and genome, the best locus (the HSPs within --locus-span nt of the best HSP, on its strand; E-value,
   bit score and identity of the best HSP) and, in the HSP of that locus that covers the position of the U,
   the codon aligned to it, read in the genome; when the alignment stops at most --extend codons before the U, the
   codon is read in the same frame beyond its end (flagged "extended"). Then the codons from it to the next stop codon
   of the frame, and the best SECIS hits of the Oikopleura model and of RF00031 in the --window nt after that stop
   codon, on the strand of the hit.
3. Proteomes (--proteomes: RefSeq proteins of the genomes that have them, vertebrates included): blastp (E <=
   --evalue); best hit, the residue aligned to the U in a local alignment (BLOSUM62) and the length of the hit
   protein after that residue.
Outputs: results/tables/stage2_candidates_rule.tsv (the rule), stage2_candidates_genomes.tsv (tblastn),
stage2_candidates_proteomes.tsv (blastp).
"""
import argparse, collections, concurrent.futures, csv, gzip, os, re, subprocess, sys, xml.etree.ElementTree as ET

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
STOPS = {"TAA", "TAG", "TGA"}
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
        for k, c in enumerate(BASES)}
DIOICA = ["OKI", "OSKA", "OdB3"]
GENOMES = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor", "Cint", "Crob", "Csav", "Hror",
           "Scla", "Tdem", "Bsch", "Pmam", "Pege", "Stho", "Bflo", "Blan"]
PROTEOMES = ["Hsap", "Mmus", "Drer", "Cint", "Crob", "Scla", "Tdem", "Bsch", "Pmam", "Bflo", "Blan"]


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


def run(cmd, cwd=None):
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode:
        sys.exit(f"command failed ({p.returncode}) in {cwd or '.'}: {' '.join(cmd)}\n{p.stderr[-2000:]}")
    return p.stdout


def translate(s): return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def strand_seq(seq, strand): return seq if strand == "+" else rc(seq)


def after_codon(s, u):
    """Codons from the codon at u (strand sequence) to the next stop codon (excluded), and the end of that stop codon
    (None when the sequence ends first)."""
    k = u + 3
    while k + 3 <= len(s) and s[k:k + 3] not in STOPS: k += 3
    return ((k - u - 3) // 3, k + 3) if k + 3 <= len(s) else ((len(s) - u - 3) // 3, None)


class Secis:
    """Best SECIS hits of covariance models in windows (top strand)."""

    def __init__(self, models, work):
        self.models, self.work = {k: v for k, v in models.items() if os.path.exists(v)}, work

    def best(self, windows, cpus, tag):
        """windows: id -> sequence; returns id -> {model: (bits, from, to)}."""
        out = collections.defaultdict(dict)
        if not windows: return out
        fa = f"{self.work}/{tag}.windows.fa"
        with open(fa, "w") as o:
            for k, s in windows.items():
                if s: o.write(f">{k}\n{s}\n")
        for m, cm in self.models.items():
            tbl = f"{self.work}/{tag}.{m}.tbl"
            run(["cmsearch", "--toponly", "-T", "0", "--cpu", str(cpus), "--tblout", tbl, "-o", os.devnull, cm, fa])
            for l in open(tbl):
                if l.startswith("#"): continue
                f = l.split()
                if float(f[14]) > out[f[0]].get(m, (-1,))[0]: out[f[0]][m] = (float(f[14]), int(f[7]), int(f[8]))
        return out


def read_candidates(path, faa):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    prot = dict(read_fasta(faa))
    return rows, prot


def rule(a, rows, sheet, secis):
    """The rule in the three O. dioica assemblies."""
    thr = {}
    if os.path.exists(a.scan_summary):
        for r in csv.DictReader((l for l in open(a.scan_summary) if not l.startswith("#")), delimiter="\t"):
            if r["set"] == "genome": thr[r["genome"]] = r.get("genome threshold (bits)", "")
    out, windows = [], {}
    by_genome = collections.defaultdict(list)
    for r in rows:
        if r["genome"] in DIOICA: by_genome[r["genome"]].append(r)
    for g, rs in by_genome.items():
        if g not in sheet:
            out += [dict(r, codon_found="genome absent") for r in rs]; continue
        need = {r["contig"] for r in rs}
        contigs = {k: v.upper() for k, v in read_fasta(sheet[g]["fasta"]) if k in need}
        for r in rs:
            seq, st, p = contigs.get(r["contig"], ""), r["strand"], int(r["codon"]) - 1
            if not seq: out.append(dict(r, codon_found="contig absent")); continue
            s = strand_seq(seq, st)
            u = p if st == "+" else len(seq) - p - 3
            n_after, stop_end = after_codon(s, u)
            key = f"{r['gene']}|{g}"
            if stop_end is not None: windows[key] = s[stop_end:stop_end + a.window]
            # context as in curation/candidates.tsv: 8 codons before the TGA, U, the codons to the stop codon, *, and
            # the 3 codons after it
            ctx = translate(s[max(u % 3, u - 24):u]) + "U"
            ctx += (translate(s[u + 3:stop_end - 3]) + "*" + translate(s[stop_end:stop_end + 9])) if stop_end else translate(s[u + 3:])
            out.append(dict(r, codon_found=s[u:u + 3], codons_after=n_after, context_found=ctx, key=key,
                            genome_threshold=thr.get(g, "")))
    best = secis.best(windows, a.cpus, "rule")
    for r in out:
        b = best.get(r.get("key"), {}).get("oikopleura")
        r["secis_bits"], r["secis_position"] = (b[0], f"{b[1]}-{b[2]}") if b else ("", "")
    verdict = {}
    for gene in {r["gene"] for r in out}:
        rs = {r["genome"]: r for r in out if r["gene"] == gene}
        ok = all(g in rs and rs[g].get("codon_found") == "TGA" and rs[g]["secis_bits"] != "" and
                 rs[g]["secis_bits"] >= a.secis_bits for g in DIOICA)
        verdict[gene] = "counted" if ok else "candidate"
    for r in out: r["verdict"] = verdict[r["gene"]]
    return out


def tblastn_genome(g, a, sheet, prot, secis):
    """tblastn of the candidate proteins against one genome; returns rows."""
    d = f"{a.work}/{g}"
    os.makedirs(d, exist_ok=True)
    fa = f"{d}/genome.fa"
    lengths = {}
    with open(fa, "w") as o:
        for k, s in read_fasta(sheet[g]["fasta"]):
            o.write(f">{k}\n{s.upper()}\n"); lengths[k] = len(s)
    run(["formatdb", "-i", "genome.fa", "-p", "F"], d)
    with open(f"{d}/q.faa", "w") as o:
        for k, s in prot.items(): o.write(f">{k}\n{s.replace('U', 'X')}\n")
    run(["blastall", "-p", "tblastn", "-d", "genome.fa", "-i", "q.faa", "-e", str(a.evalue), "-F", "F", "-m", "7",
         "-v", "20", "-b", "20", "-a", str(a.cpus_per_genome), "-o", "tblastn.xml"], d)
    hsps = collections.defaultdict(list)
    for _, it in ET.iterparse(f"{d}/tblastn.xml"):
        if it.tag != "Iteration": continue
        q = it.findtext("Iteration_query-def").split()[0]
        for hit in it.iter("Hit"):
            contig = hit.findtext("Hit_def").split()[0]
            for h in hit.iter("Hsp"):
                f0, f1 = int(h.findtext("Hsp_hit-from")), int(h.findtext("Hsp_hit-to"))
                hsps[q].append({"contig": contig, "frame": int(h.findtext("Hsp_hit-frame")),
                                "lo": min(f0, f1), "hi": max(f0, f1), "qfrom": int(h.findtext("Hsp_query-from")),
                                "qto": int(h.findtext("Hsp_query-to")), "qseq": h.findtext("Hsp_qseq"),
                                "hseq": h.findtext("Hsp_hseq"), "bits": float(h.findtext("Hsp_bit-score")),
                                "evalue": float(h.findtext("Hsp_evalue")), "ident": int(h.findtext("Hsp_identity")),
                                "alen": int(h.findtext("Hsp_align-len"))})
        it.clear()
    # the codon aligned to the U, in the HSPs of the best locus
    out, windows, need = [], {}, set()
    for q, hs in hsps.items():
        best = max(hs, key=lambda h: h["bits"])
        need.add(best["contig"])
    seqs = {k: s for k, s in read_fasta(fa) if k in need}
    for q, hs in hsps.items():
        best = max(hs, key=lambda h: h["bits"])
        strand = "+" if best["frame"] > 0 else "-"
        locus = [h for h in hs if h["contig"] == best["contig"] and (h["frame"] > 0) == (best["frame"] > 0) and
                 abs((h["lo"] + h["hi"]) / 2 - (best["lo"] + best["hi"]) / 2) < a.locus_span]
        upos = prot[q].find("U")   # 0-based in the query
        row = {"candidate": q, "genome": g, "contig": best["contig"], "strand": strand,
               "from": min(h["lo"] for h in locus), "to": max(h["hi"] for h in locus), "evalue": best["evalue"], "bits": best["bits"],
               "identity": f"{100 * best['ident'] / best['alen']:.0f}", "hsps": len(locus)}
        seq = seqs.get(best["contig"], "")
        s = strand_seq(seq, strand)
        codon_pos, note = None, ""
        for h in sorted(locus, key=lambda h: -h["bits"]):
            if upos < 0 or not seq: break
            start = (h["lo"] - 1) if strand == "+" else (len(seq) - h["hi"])   # first aligned base, strand coordinates
            qi, k = h["qfrom"] - 1, 0
            for c in range(len(h["qseq"])):
                if h["qseq"][c] != "-":
                    if qi == upos and h["hseq"][c] != "-": codon_pos = start + 3 * k; break
                    qi += 1
                if h["hseq"][c] != "-": k += 1
            if codon_pos is not None: break
            if 0 <= upos - (h["qto"] - 1) <= a.extend:   # the alignment stops just before the U
                n_h = sum(1 for x in h["hseq"] if x != "-")
                codon_pos, note = start + 3 * (n_h + upos - h["qto"]), "extended"
                break
        if codon_pos is not None and 0 <= codon_pos <= len(s) - 3:
            n_after, stop_end = after_codon(s, codon_pos)
            row.update({"codon": s[codon_pos:codon_pos + 3], "residue": CODE.get(s[codon_pos:codon_pos + 3], "X"),
                        "codon_note": note, "codons_after": n_after,
                        "position": (codon_pos + 1) if strand == "+" else (len(seq) - codon_pos - 2)})
            if stop_end is not None:
                row["stop"] = s[stop_end - 3:stop_end]
                windows[f"{q}|{g}"] = s[stop_end:stop_end + a.window]
        else:
            row.update({"codon": "", "residue": "", "codon_note": "U position not aligned", "codons_after": ""})
        out.append(row)
    best = secis.best(windows, a.cpus_per_genome, f"{g}/secis")
    for r in out:
        b = best.get(f"{r['candidate']}|{g}", {})
        r["secis_oikopleura"] = b["oikopleura"][0] if "oikopleura" in b else ""
        r["secis_RF00031"] = b["RF00031"][0] if "RF00031" in b else ""
    for x in os.listdir(d):   # the genome copy and its database are large
        if x.startswith("genome.fa"): os.remove(os.path.join(d, x))
    return out


def proteome_hits(g, a, prot):
    from Bio import Align
    from Bio.Align import substitution_matrices
    path = f"{a.data}/{g}/protein.faa.gz"
    if not os.path.exists(path): return []
    d = f"{a.work}/proteome_{g}"
    os.makedirs(d, exist_ok=True)
    seqs, desc = {}, {}
    with opener(path) as fh:
        name = None
        for l in fh:
            if l.startswith(">"):
                name, _, dsc = l[1:].strip().partition(" "); desc[name] = re.sub(r" \[[^]]+\]$", "", dsc); seqs[name] = []
            elif name: seqs[name].append(l.strip())
    seqs = {k: "".join(v).replace("*", "") for k, v in seqs.items()}
    with open(f"{d}/db.faa", "w") as o:
        for k, s in seqs.items(): o.write(f">{k}\n{s}\n")
    run(["formatdb", "-i", "db.faa", "-p", "T"], d)
    with open(f"{d}/q.faa", "w") as o:
        for k, s in prot.items(): o.write(f">{k}\n{s.replace('U', 'X')}\n")
    best = {}
    for l in run(["blastall", "-p", "blastp", "-d", "db.faa", "-i", "q.faa", "-e", str(a.evalue), "-m", "8", "-a",
                  str(a.cpus_per_genome)], d).splitlines():
        f = l.split("\t")
        if f[0] not in best: best[f[0]] = (f[1], float(f[10]), float(f[2]))
    al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-11, extend_gap_score=-1)
    out = []
    for q, (hit, ev, ident) in best.items():
        x, t = prot[q], seqs[hit]
        aln = al.align(x.replace("U", "C"), t.replace("U", "C"))[0]
        u, res, pos = x.find("U"), "-", None
        for (q0, q1), (t0, t1) in zip(*aln.aligned):
            if q0 <= u < q1: pos = t0 + u - q0; res = t[pos]; break
        out.append({"candidate": q, "genome": g, "protein": hit, "description": desc.get(hit, ""), "evalue": ev,
                    "identity": ident, "length": len(t), "residue": res,
                    "position": pos + 1 if pos is not None else "",
                    "residues_after": len(t) - pos - 1 if pos is not None else ""})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--candidates", default="curation/candidates.tsv")
    ap.add_argument("--proteins", default="curation/candidates.faa")
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--data", default="data/genomes")
    ap.add_argument("--model", default="results/models/oikopleura_secis.cm")
    ap.add_argument("--rfam", default="data/references/RF00031.cm")
    ap.add_argument("--scan-summary", default="results/tables/stage2_secis_scan_summary.tsv")
    ap.add_argument("--genomes", nargs="*", default=GENOMES)
    ap.add_argument("--proteomes", nargs="*", default=PROTEOMES)
    ap.add_argument("--window", type=int, default=2000)
    ap.add_argument("--secis-bits", type=float, default=20.0)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--extend", type=int, default=5)
    ap.add_argument("--locus-span", type=int, default=20000, help="HSPs within this distance of the best one form its locus")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--work", default="work/stage2/candidates")
    ap.add_argument("--tables", default="results/tables")
    a = ap.parse_args()
    a.cpus_per_genome = max(1, a.cpus // a.jobs)
    os.makedirs(a.work, exist_ok=True); os.makedirs(a.tables, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    rows, prot = read_candidates(a.candidates, a.proteins)
    secis = Secis({"oikopleura": a.model, "RF00031": a.rfam}, a.work)
    # 1. the rule
    res = rule(a, rows, sheet, secis)
    cols = ["gene", "genome", "contig", "codon", "strand", "codon_found", "codons_after", "context", "context_found",
            "secis_bits", "secis_position", "genome_threshold", "verdict", "notes"]
    with open(f"{a.tables}/stage2_candidates_rule.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in res: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    for gene in sorted({r["gene"] for r in res}):
        print(gene, next(r["verdict"] for r in res if r["gene"] == gene), flush=True)
    # 2. other genomes
    genomes = [g for g in a.genomes if g in sheet]
    out = []
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
        for rs in ex.map(lambda g: tblastn_genome(g, a, sheet, prot, secis), genomes):
            out += rs; print(f"tblastn {rs[0]['genome'] if rs else '?'}: {len(rs)} candidates with a hit", flush=True)
    cols = ["candidate", "genome", "contig", "strand", "from", "to", "evalue", "bits", "identity", "hsps", "codon",
            "residue", "codon_note", "position", "codons_after", "stop", "secis_oikopleura", "secis_RF00031"]
    order = {g: i for i, g in enumerate(a.genomes)}
    with open(f"{a.tables}/stage2_candidates_genomes.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in sorted(out, key=lambda r: (r["candidate"], order.get(r["genome"], 99))):
            o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    # 3. proteomes
    out = []
    for g in a.proteomes: out += proteome_hits(g, a, prot)
    cols = ["candidate", "genome", "protein", "description", "evalue", "identity", "length", "residue", "position",
            "residues_after"]
    with open(f"{a.tables}/stage2_candidates_proteomes.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in out: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


if __name__ == "__main__":
    main()
