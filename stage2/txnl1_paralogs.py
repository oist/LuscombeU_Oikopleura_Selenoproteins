#!/usr/bin/env python3
"""The two TXNL1-related genes of O. dioica in other genomes: the full-length TXNL1 (thioredoxin domain and PITH
region; OKI2018 XP_081633344.1, the reciprocal best hit of human TXNL1) and the selenocysteine copy (PITH region only;
gold-standard gene --sec-gene). Follows make txnl1-homologs (results/tables/stage2_txnl1_homologs_summary.tsv).

1. Queries. In each O. dioica assembly (--dioica), the selenocysteine copy (gold standard) and the full-length copy:
   the best blastp hit of human TXNL1 among the annotated proteins outside the locus of the selenocysteine copy,
   trimmed to the part aligned with human TXNL1 (and the tail after it) when the model is longer, as the OSKA2016
   model merged with a GPR107-like gene. In local alignments (BLOSUM62) with the selenocysteine copy: the residue of
   the full-length copy and of human TXNL1 aligned to the selenocysteine, and the residues after it. The copies of
   --query-genome are the queries of steps 2 and 3; the boundary between the thioredoxin part and the PITH part of
   the full-length copy is the residue aligned to human residue --split - 1.
2. Genomes (--genomes): tblastn (legacy BLAST, no filter, E <= --evalue, XML) of the two queries. Colinear HSPs within
   --locus-span nt form a locus per query (function loci of txnl1_homologs.py); overlapping loci of the two queries on
   one strand form a site. For each site: coverage of the thioredoxin part and of the PITH part of the full-length
   query, structure (full-length, PITH only, thioredoxin only, partial), the distance from its 5' end to the end of
   the sequence (a PITH-only site close to it may be an incomplete gene), the codon aligned to the selenocysteine (in
   the HSPs of the selenocysteine copy, else of the full-length copy at its aligned position; read beyond the end of
   the alignment when it stops at most --extend codons before, "extended"), the codons to the next stop codon, the best
   SECIS hits of the Oikopleura model and of RF00031 in the --window nt after it, the annotated transcripts on the same
   strand that overlap the site, the best human gene of its PITH-region sequence (subject residues at the query
   positions, from the best HSP covering each; blastp) and that sequence.
3. Proteomes (--proteomes, RefSeq; for the O. dioica assemblies, the annotated proteins translated from the genome and
   the gold-standard proteins): blastp of the two queries (E <= --evalue, --max-hits per query); every protein hit in
   the PITH part, with its length, its coverage of the two parts of the full-length query, its best human gene
   (blastp against the human proteome) and the residue aligned to the selenocysteine with the residues after it.
4. Tree of the PITH regions: the TXNL1 and PITHD1 proteins of the proteomes (best human gene TXNL1 or PITHD1, longest
   protein per gene; the O. dioica copies of --tree-dioica only), the genome sites without such a protein (genomes
   without a proteome, and sites that overlap no annotated transcript), at least --min-length residues each; MAFFT
   L-INS-i, ClipKIT (gappy, 0.9), IQ-TREE 2 (ModelFinder, 1,000 ultrafast bootstraps), rooted on the clade of the
   sequences whose best human gene is PITHD1 (by at least --root-margin bits over the next gene), else on human
   PITHD1. Tip names: genome, protein or site, best human gene, structure and the residue aligned to
   the selenocysteine.
5. Neighbourhood of the full-length copy in the O. dioica assemblies: best tblastn loci of the full-length query and of
   human --neighbour-gene (GPR107: the OSKA2016 model of the full-length copy holds a GPR107-like part): strands,
   order along the strand of the full-length copy, distance, the annotated models overlapping each (one model over
   both: merged), and the transcripts of the assembly (samplesheet) that align to both loci (blastn, at least
   --tx-identity % identity over at least --tx-length nt on each).
Outputs: <prefix>_queries.tsv (1), _sites.tsv (2), _proteins.tsv (3), _tree.tsv (4: the clade of each O. dioica copy
and its support), _neighbours.tsv (5), with --prefix results/tables/stage2_txnl1_paralogs; tree files in --tree-dir;
figure in --figure (.png and .pdf). Per-genome results are cached in --work (a rerun does not search again).
"""
import argparse, collections, concurrent.futures, csv, json, os, re, shutil, subprocess, sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uga_secis_screen import read_fasta, rc, opener, read_gff_transcripts, attributes   # noqa: E402
from candidate_genes import Secis, after_codon, CODE                          # noqa: E402
from txnl1_homologs import Db, blast, clean, annotated_proteins, gold_loci, overlaps, loci, human_genes  # noqa: E402

GENOMES = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor", "Cint", "Crob", "Csav", "Hror",
           "Scla", "Tdem", "Bsch", "Pmam", "Pege", "Stho", "Bflo", "Blan"]
PROTEOMES = ["Hsap", "Mmus", "Drer", "Cint", "Crob", "Scla", "Tdem", "Bsch", "Pmam", "Bflo", "Blan"]
DIOICA = ["OKI", "OSKA", "OdB3"]


def run(cmd, cwd=None, stdout=None):
    p = subprocess.run(cmd, cwd=cwd, stdout=stdout or subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode: sys.exit(f"command failed ({p.returncode}) in {cwd or '.'}: {' '.join(cmd)}\n{p.stderr[-2000:]}")
    return p.stdout


def aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    return Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                 open_gap_score=-11, extend_gap_score=-1)


def pair(al, x, y):
    """Aligned blocks [((x0, x1), (y0, y1))] of the best local alignment of x and y."""
    def norm(s): return s.upper().replace("U", "C").replace("*", "X")
    try:
        aln = al.align(norm(x), norm(y))[0]
    except (IndexError, ValueError):   # no alignment with a positive score
        return []
    return [tuple(map(tuple, b)) for b in zip(*aln.aligned)]


def mapped(blocks, i, extend=0, length=None):
    """Position in y aligned to position i of x (0-based), "extended" beyond the end of the alignment (within y, of
    length length), or (None, "")."""
    for (x0, x1), (y0, y1) in blocks:
        if x0 <= i < x1: return int(y0 + i - x0), ""
    if blocks:
        (_, x1), (_, y1) = blocks[-1]
        if 0 <= i - x1 < extend and (length is None or y1 + i - x1 < length): return int(y1 + i - x1), "extended"
    return None, ""


def x_coverage(blocks, lo, hi):
    """Fraction of x positions lo..hi-1 inside aligned blocks."""
    n = sum(max(0, min(x1, hi) - max(x0, lo)) for (x0, x1), _ in blocks)
    return n / (hi - lo) if hi > lo else 0


def safe(s): return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")


def tblastn_xml(dbpath, queries, work, tag, evalue, cpus):
    """HSPs of a tblastn search (XML): q, s, ss, se (subject coordinates in the orientation of the frame, as in -m 8),
    qs, qe, qseq, hseq, bits, evalue, ident."""
    fa = os.path.join(work, f"{tag}.q.faa")
    names = list(queries)
    with open(fa, "w") as o:
        for n, k in enumerate(names, 1): o.write(f">q{n}\n{queries[k]}\n")
    xml = os.path.join(work, f"{tag}.xml")
    run(["blastall", "-p", "tblastn", "-d", dbpath, "-i", fa, "-e", str(evalue), "-F", "F", "-m", "7", "-v", "50",
         "-b", "200", "-a", str(cpus), "-o", xml])
    out = []
    for _, it in ET.iterparse(xml):
        if it.tag != "Iteration": continue
        num = int(it.findtext("Iteration_iter-num") or 0)
        qdef = it.findtext("Iteration_query-def")
        q = names[int(qdef.split()[0][1:]) - 1] if qdef else names[num - 1]
        for hit in it.iter("Hit"):
            s = hit.findtext("Hit_def").split()[0]
            for h in hit.iter("Hsp"):
                f0, f1, fr = int(h.findtext("Hsp_hit-from")), int(h.findtext("Hsp_hit-to")), int(h.findtext("Hsp_hit-frame"))
                lo, hi = min(f0, f1), max(f0, f1)
                ss, se = (lo, hi) if fr > 0 else (hi, lo)
                out.append({"q": q, "s": s, "ss": ss, "se": se, "qs": int(h.findtext("Hsp_query-from")),
                            "qe": int(h.findtext("Hsp_query-to")), "qseq": h.findtext("Hsp_qseq"),
                            "hseq": h.findtext("Hsp_hseq"), "bits": float(h.findtext("Hsp_bit-score")),
                            "evalue": float(h.findtext("Hsp_evalue")),
                            "ident": 100 * int(h.findtext("Hsp_identity")) / int(h.findtext("Hsp_align-len"))})
        it.clear()
    return out


def codon_at(hsps, upos, seq, strand, extend):
    """Strand-coordinate position (0-based) of the codon aligned to query position upos (0-based), and a note."""
    for h in sorted(hsps, key=lambda h: -h["bits"]):
        lo, hi = min(h["ss"], h["se"]), max(h["ss"], h["se"])
        start = (lo - 1) if strand == "+" else (len(seq) - hi)
        qi, k = h["qs"] - 1, 0
        for c in range(len(h["qseq"])):
            if h["qseq"][c] != "-":
                if qi == upos and h["hseq"][c] != "-": return start + 3 * k, ""
                qi += 1
            if h["hseq"][c] != "-": k += 1
    for h in sorted(hsps, key=lambda h: -h["bits"]):
        if 0 <= upos - (h["qe"] - 1) <= extend:
            lo, hi = min(h["ss"], h["se"]), max(h["ss"], h["se"])
            start = (lo - 1) if strand == "+" else (len(seq) - hi)
            n_h = sum(1 for x in h["hseq"] if x != "-")
            return start + 3 * (n_h + upos - h["qe"]), "extended"
    return None, "not aligned"


def query_anchored(hsps, lo, hi):
    """Subject residues at query positions lo..hi (1-based, inclusive), from the best HSP covering each."""
    res = {}
    for h in sorted(hsps, key=lambda h: -h["bits"]):
        qi = h["qs"]
        for c, d in zip(h["qseq"], h["hseq"]):
            if c == "-": continue
            if lo <= qi <= hi and qi not in res and d != "-": res[qi] = "X" if d == "*" else d
            qi += 1
    return "".join(res[i] for i in sorted(res))


def covered(hsps, lo, hi):
    pos = set()
    for h in hsps: pos.update(range(max(lo, h["qs"]), min(hi, h["qe"]) + 1))
    return len(pos) / (hi - lo + 1) if hi >= lo else 0


def model_structure(gff, model_id):
    """Contig, strand, exons and coding segments (1-based inclusive, genome order) and the phase of the first coding
    segment in transcript order, of the annotated model whose CDS rows carry protein_id (or Parent) model_id."""
    rows = [l.rstrip("\n").split("\t") for l in opener(gff) if not l.startswith("#")]
    rows = [f for f in rows if len(f) >= 9 and f[2] in ("CDS", "exon")]
    cds, parents, contig, strand = [], set(), None, "+"
    for f in rows:
        if f[2] != "CDS": continue
        at = attributes(f[8])
        par = at.get("Parent", "").split(",")[0]
        if at.get("protein_id") == model_id or par == model_id:
            cds.append((int(f[3]), int(f[4]), int(f[7]) if f[7] in "012" else 0))
            parents.add(par); contig, strand = f[0], f[6]
    cds.sort()
    exons = sorted({(int(f[3]), int(f[4])) for f in rows if f[2] == "exon" and attributes(f[8]).get("Parent", "").split(",")[0] in parents})
    phase = (cds[0][2] if strand == "+" else cds[-1][2]) if cds else 0
    return contig, strand, exons or [(x, y) for x, y, _ in cds], [(x, y) for x, y, _ in cds], phase


def coding_positions(cds, strand):
    """Genomic positions (1-based) of the coding bases, in transcript order."""
    out = []
    for x, y in (cds if strand == "+" else cds[::-1]):
        out += list(range(x, y + 1)) if strand == "+" else list(range(y, x - 1, -1))
    return out


def splice_sites(seq, exons, strand):
    """Donor-acceptor dinucleotides of the introns between exons (genome order), read in gene orientation."""
    out = []
    for (x1, y1), (x2, y2) in zip(exons, exons[1:]):
        if x2 <= y1 + 1: continue
        d, c = seq[y1:y1 + 2], seq[x2 - 3:x2 - 1]
        out.append(f"{d}-{c}" if strand == "+" else f"{rc(c)}-{rc(d)}")
    return out


def read_proteome(path):
    seqs, desc, name = {}, {}, None
    with opener(path) as fh:
        for l in fh:
            if l.startswith(">"):
                name, _, d = l[1:].strip().partition(" ")
                desc[name], seqs[name] = re.sub(r" \[[^]]+\]$", "", d), []
            elif name: seqs[name].append(l.strip())
    return {k: clean("".join(v)) for k, v in seqs.items()}, desc


def genes_of(gff):
    out = {}
    if not gff or not os.path.exists(gff): return out
    for l in opener(gff):
        if "\tCDS\t" not in l: continue
        m, g = re.search(r"protein_id=([^;\n]+)", l), re.search(r"gene=([^;\n]+)", l)
        if m and g: out[m.group(1)] = g.group(1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--data", default="data/genomes")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--sec-gene", default="TXNL1.b", help="gold-standard name of the selenocysteine copy")
    ap.add_argument("--full-gene", default="TXNL1.a", help="name of the full-length copy")
    ap.add_argument("--figure-inputs", default="figures/inputs",
                    help="where reannotation_<sec-gene>.paralogs.{tsv,fa} go (the full-length copy, for make reannotation-pages)")
    ap.add_argument("--dioica", nargs="+", default=DIOICA)
    ap.add_argument("--query-genome", default="OKI")
    ap.add_argument("--genomes", nargs="*", default=GENOMES)
    ap.add_argument("--proteomes", nargs="*", default=PROTEOMES)
    ap.add_argument("--tree-dioica", default="OKI")
    ap.add_argument("--human-gene", default="TXNL1")
    ap.add_argument("--human-txnl1", default="NP_004777.1")
    ap.add_argument("--outgroup-gene", default="PITHD1")
    ap.add_argument("--neighbour-gene", default="GPR107")
    ap.add_argument("--split", type=int, default=118)
    ap.add_argument("--model", default="results/models/oikopleura_secis.cm")
    ap.add_argument("--rfam", default="data/references/RF00031.cm")
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--max-hits", type=int, default=20)
    ap.add_argument("--extend", type=int, default=5)
    ap.add_argument("--window", type=int, default=2000)
    ap.add_argument("--locus-span", type=int, default=20000)
    ap.add_argument("--min-length", type=int, default=60)
    ap.add_argument("--root-margin", type=float, default=10.0, help="bits between the best and the next human gene")
    ap.add_argument("--tx-identity", type=float, default=95)
    ap.add_argument("--tx-length", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--work", default="work/stage2/txnl1_paralogs")
    ap.add_argument("--prefix", default="results/tables/stage2_txnl1_paralogs")
    ap.add_argument("--tree-dir", default="results/trees/txnl1_pith")
    ap.add_argument("--figure", default="figures/txnl1_pith_tree")
    ap.add_argument("--no-tree", action="store_true")
    a = ap.parse_args()
    cpg = max(1, a.cpus // a.jobs)
    os.makedirs(a.work, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    al = aligner()

    # human proteome
    hseqs, hdesc = read_proteome(f"{a.data}/Hsap/protein.faa.gz")
    hgene = genes_of(f"{a.data}/Hsap/annotation.gff3.gz")
    hdb = Db(a.work, "human", hseqs.items(), True)
    h = hseqs[a.human_txnl1]

    def to_human(seqs, tag):
        if not seqs: return {}
        rows = blast("blastp", hdb, seqs, a.work, tag, 1e-3, a.cpus, 100)
        for r in rows: r["gene"] = hgene.get(r["s"], r["s"])
        return human_genes(rows)

    def gene_txt(genes):
        if not genes: return "", "", ""
        return genes[0][0], f"{genes[0][2]:.1f}", (f"{genes[1][0]} {genes[1][2]:.1f}" if len(genes) > 1 else "")

    neighbour = max((p for p, g in hgene.items() if g == a.neighbour_gene and p in hseqs), key=lambda p: len(hseqs[p]),
                    default=None)

    # 1. queries
    qrows, copies, fig_rows, fig_fa = [], {}, [], []
    for g in a.dioica:
        gold = {k.split("|")[0]: s.upper().rstrip("*") for k, s in read_fasta(f"{a.gold}/{g}.proteins.faa")}   # U kept
        sec = gold[a.sec_gene]
        sloc = gold_loci(f"{a.gold}/{g}.gff3").get(a.sec_gene)
        ann = annotated_proteins(sheet[g]["fasta"], sheet[g]["annotation"])
        cand = {k: r["protein"] for k, r in ann.items() if r["protein"] and
                not overlaps(sloc, r["contig"], r["strand"], r["lo"], r["hi"])}
        pdb = Db(a.work, f"{g}.annotated", cand.items(), True)
        rows = blast("blastp", pdb, {"human": h}, a.work, f"{g}.full", a.evalue, a.cpus, 5)
        shutil.rmtree(pdb.dir, ignore_errors=True)
        if not rows: sys.exit(f"{g}: no full-length copy found")
        fid = max(rows, key=lambda r: r["bits"])["s"]
        full, note = cand[fid], ""
        trim = (0, len(full))
        hb = pair(al, h, full)
        if len(full) > 1.3 * len(h):   # merged model: keep the part aligned with human TXNL1 and the tail after it
            y0, end = hb[0][1][0], min(len(full), hb[-1][1][1] + len(h) - hb[-1][0][1])
            full, note, trim = full[y0:end], f"trimmed to residues {y0 + 1}-{end} of {len(full)}", (y0, end)
            hb = pair(al, h, full)
        u = sec.find("U")
        copies[g] = {"sec": sec, "full": full, "full_id": fid, "sec_locus": sloc, "full_locus": ann[fid]}
        sb = pair(al, sec, full)
        for label, prot, pid, blocks, loc in (
                (f"{a.sec_gene} (selenocysteine copy)", sec, f"gold:{a.sec_gene}", None, sloc),
                (f"{a.full_gene} (full-length copy)", full, fid, sb, (ann[fid]["contig"], ann[fid]["strand"], ann[fid]["lo"], ann[fid]["hi"]))):
            if blocks is None: pos, pn = u, ""
            else: pos, pn = mapped(blocks, u, a.extend, len(prot))
            qrows.append([g, label, pid, len(prot), f"{loc[0]}:{loc[2] + 1}-{loc[3]}{loc[1]}" if loc else "",
                          prot[pos] if pos is not None else "", (pos + 1) if pos is not None else "", pn,
                          len(prot) - pos - 1 if pos is not None else "", prot[pos:] if pos is not None else "",
                          note if label.startswith(a.full_gene) else ""])
        # the full-length copy for Figure CANDIDATES: model, codon aligned to the selenocysteine, drawn span
        contig, strand, exons, cds, phase = model_structure(sheet[g]["annotation"], fid)
        cpos = coding_positions(cds, strand)
        pu, _ = mapped(pair(al, sec, cand[fid]), u, a.extend, len(cand[fid]))
        seq = next((s.upper() for k, s in read_fasta(sheet[g]["fasta"]) if k == contig), "")
        if pu is not None and phase + 3 * pu + 3 <= len(cpos) and seq:
            codon = sorted(cpos[phase + 3 * pu:phase + 3 * pu + 3])
            part = cpos[phase + 3 * trim[0]:min(len(cpos), phase + 3 * trim[1])]
            db = {"OKI": "RefSeq", "OSKA": "ANISEED", "OdB3": "GenBank"}.get(g, "annotation")
            src = "annotation; make txnl1-paralogs"
            fig_rows += [("paralog", g, fid, "label", a.full_gene, "make txnl1-paralogs"),
                         ("paralog", g, fid, "database", db, ""), ("paralog", g, fid, "strand", strand, src),
                         ("paralog", g, fid, "exons", ";".join(f"{x}-{y}" for x, y in exons), src),
                         ("paralog", g, fid, "cds", ";".join(f"{x}-{y}" for x, y in cds), src),
                         ("paralog", g, fid, "sites", ";".join(splice_sites(seq, exons, strand)), "genome sequence"),
                         ("paralog", g, fid, "anchor", f"{codon[0]}-{codon[-1]}",
                          f"codon aligned to the selenocysteine of {a.sec_gene} (local alignment, BLOSUM62)"),
                         ("paralog", g, fid, "span", f"{min(part)}-{max(part)}", "coding part aligned to human TXNL1")]
            if trim != (0, len(cand[fid])):
                fig_rows.append(("paralog", g, fid, "note", f"model of {len(cand[fid])} residues; drawn: residues "
                                 f"{trim[0] + 1}-{trim[1]}, aligned to human TXNL1", "make txnl1-paralogs"))
            fig_fa.append((f"paralog:{g}", full))
        if g == a.query_genome:
            trx_end = None
            for (x0, x1), (y0, y1) in hb:   # human residue split-1 -> full-length copy
                if x0 <= a.split - 2: trx_end = int(y0 + min(x1, a.split - 1) - x0)
            copies[g]["trx_end"] = trx_end or a.split - 1
    hb = pair(al, copies[a.query_genome]["sec"], h)
    pos, pn = mapped(hb, copies[a.query_genome]["sec"].find("U"), a.extend, len(h))
    qrows.append(["Hsap", "human TXNL1", a.human_txnl1, len(h), "", h[pos] if pos is not None else "",
                  (pos + 1) if pos is not None else "", pn, len(h) - pos - 1 if pos is not None else "",
                  h[pos:] if pos is not None else "", f"thioredoxin/PITH boundary of the full-length query: residue "
                  f"{copies[a.query_genome]['trx_end']} (human {a.split - 1})"])
    if fig_rows:
        os.makedirs(a.figure_inputs, exist_ok=True)
        base_in = os.path.join(a.figure_inputs, f"reannotation_{a.sec_gene}.paralogs")
        with open(f"{base_in}.tsv", "w") as o:
            o.write("kind\tassembly\tname\tfield\tvalue\tsource\n")
            for r in fig_rows: o.write("\t".join(r) + "\n")
        with open(f"{base_in}.fa", "w") as o:
            for k, s in fig_fa: o.write(f">{k}\n{s}\n")
    with open(f"{a.prefix}_queries.tsv", "w") as o:
        o.write("genome\tcopy\tprotein\tlength\tlocation\tresidue_at_sec\tposition\tnote\tresidues_after\t"
                "c_terminus\tremark\n")
        for r in qrows: o.write("\t".join(str(x) for x in r) + "\n")
    Q = copies[a.query_genome]
    sec_q, full_q, trx_end = Q["sec"], Q["full"], Q["trx_end"]
    upos = sec_q.find("U")
    fpos, _ = mapped(pair(al, sec_q, full_q), upos, a.extend, len(full_q))
    queries = {"sec": sec_q.replace("U", "X"), "full": full_q}
    print(f"queries: {a.query_genome} {a.sec_gene} ({len(sec_q)} aa), {a.full_gene} {Q['full_id']} "
          f"({len(full_q)} aa, thioredoxin part 1-{trx_end})", flush=True)

    # 2 and 5. genomes
    secis = Secis({"oikopleura": a.model, "RF00031": a.rfam}, a.work)

    def genome(g):
        d = os.path.join(a.work, g)
        os.makedirs(d, exist_ok=True)
        done = os.path.join(d, "result.json")
        if os.path.exists(done): return json.load(open(done))
        lengths = {}
        with open(os.path.join(d, "genome.fa"), "w") as o:
            for k, s in read_fasta(sheet[g]["fasta"]):
                o.write(f">{k}\n{s.upper()}\n"); lengths[k] = len(s)
        run(["formatdb", "-i", "genome.fa", "-p", "F", "-o", "F"], d)
        qs = dict(queries)
        if g in a.dioica and neighbour: qs["neighbour"] = hseqs[neighbour]
        hs = tblastn_xml(os.path.join(d, "genome.fa"), qs, d, "tblastn", a.evalue, cpg)
        L = loci(hs, a.locus_span)
        sites = []
        for q in ("sec", "full"):
            for x in L.get(q, []):
                for S in sites:
                    if S["contig"] == x["contig"] and S["strand"] == x["strand"] and x["lo"] < S["hi"] and S["lo"] < x["hi"]:
                        S["lo"], S["hi"] = min(S["lo"], x["lo"]), max(S["hi"], x["hi"])
                        S.setdefault(q, x)
                        break
                else: sites.append({"contig": x["contig"], "strand": x["strand"], "lo": x["lo"], "hi": x["hi"], q: x})
        need = {S["contig"] for S in sites}
        nb = {q: (L.get(q) or [None])[0] for q in ("full", "neighbour")}
        need |= {x["contig"] for x in nb.values() if x}
        seqs = {k: s.upper() for k, s in read_fasta(os.path.join(d, "genome.fa")) if k in need}
        txs = read_gff_transcripts(sheet[g]["annotation"]) if sheet[g].get("annotation") and os.path.exists(sheet[g]["annotation"]) else {}

        def models_at(contig, strand, lo, hi):
            return sorted(f"{k}{' (' + t['gene'] + ')' if t['gene'] else ''}" for k, t in txs.items()
                          if t["contig"] == contig and t["strand"] == strand and t["blocks"][0][0] < hi and lo < t["blocks"][-1][1])

        out, windows = [], {}
        for n, S in enumerate(sites, 1):
            seq, strand = seqs[S["contig"]], S["strand"]
            s = seq if strand == "+" else rc(seq)
            fh = S["full"]["hsps"] if "full" in S else []
            sh = S["sec"]["hsps"] if "sec" in S else []
            trx_cov = covered(fh, 1, trx_end)
            pith_cov = max(covered(fh, trx_end + 1, len(full_q)), covered(sh, 1, len(sec_q)))
            structure = ("full-length" if trx_cov >= 0.3 and pith_cov >= 0.3 else "PITH only" if pith_cov >= 0.3
                         else "thioredoxin only" if trx_cov >= 0.3 else "partial")
            up = S["lo"] if strand == "+" else lengths[S["contig"]] - S["hi"]
            cp, note = codon_at(sh, upos, seq, strand, a.extend) if sh else (None, "not aligned")
            src = "sec"   # names given at output (the cached results do not depend on them)
            if cp is None and fh and fpos is not None:
                cp, note = codon_at(fh, fpos, seq, strand, a.extend)
                src = "full"
            row = {"genome": g, "site": n, "contig": S["contig"], "strand": strand, "from": S["lo"] + 1, "to": S["hi"],
                   "structure": structure, "thioredoxin_coverage": f"{trx_cov:.2f}", "pith_coverage": f"{pith_cov:.2f}",
                   "bits_sec": f"{S['sec']['bits']:.1f}" if "sec" in S else "",
                   "identity_sec": f"{max(x['ident'] for x in sh):.0f}" if sh else "",
                   "bits_full": f"{S['full']['bits']:.1f}" if "full" in S else "",
                   "identity_full": f"{max(x['ident'] for x in fh):.0f}" if fh else "",
                   "upstream_to_sequence_end": up, "annotation": ", ".join(models_at(S["contig"], strand, S["lo"], S["hi"])),
                   "codon_source": src if cp is not None else "", "codon_note": note}
            if cp is not None and 0 <= cp <= len(s) - 3:
                n_after, stop_end = after_codon(s, cp)
                row.update({"codon": s[cp:cp + 3], "residue": CODE.get(s[cp:cp + 3], "X"), "codons_after": n_after,
                            "stop": s[stop_end - 3:stop_end] if stop_end else ""})
                if stop_end: windows[f"site{n}"] = s[stop_end:stop_end + a.window]
            pith = query_anchored(sh, 1, len(sec_q)) if sh else query_anchored(fh, trx_end + 1, len(full_q))
            row["pith_sequence"] = pith
            out.append(row)
        best = secis.best(windows, cpg, f"{g}/secis") if windows else {}
        for n, row in enumerate(out, 1):
            b = best.get(f"site{n}", {})
            row["secis_oikopleura"] = b["oikopleura"][0] if "oikopleura" in b else ""
            row["secis_RF00031"] = b["RF00031"][0] if "RF00031" in b else ""
        # 5. neighbourhood (O. dioica)
        nbrow = None
        if g in a.dioica and nb["full"] and nb["neighbour"]:
            F, N = nb["full"], nb["neighbour"]
            same = F["contig"] == N["contig"] and F["strand"] == N["strand"]
            if same:
                gap = (F["lo"] - N["hi"]) if F["strand"] == "+" else (N["lo"] - F["hi"])
                order = f"{a.neighbour_gene}-like upstream" if gap >= 0 else f"{a.neighbour_gene}-like downstream"
                if gap < 0: gap = (N["lo"] - F["hi"]) if F["strand"] == "+" else (F["lo"] - N["hi"])
            else:
                order, gap = ("other strand" if F["contig"] == N["contig"] else "other sequence"), ""
            mf, mn = models_at(F["contig"], F["strand"], F["lo"], F["hi"]), models_at(N["contig"], N["strand"], N["lo"], N["hi"])
            span_tx = []
            tx = sheet[g].get("transcripts", "")
            if tx and os.path.exists(tx):
                td = os.path.join(d, "tx")
                os.makedirs(td, exist_ok=True)
                tdb = Db(td, "db", ((k, s.upper()) for k, s in read_fasta(tx)), False)
                reg = {}
                for lab, X in (("full", F), ("neighbour", N)):
                    r = seqs[X["contig"]][X["lo"]:X["hi"]]
                    reg[lab] = r if X["strand"] == "+" else rc(r)
                hit = collections.defaultdict(set)
                for r in blast("blastn", tdb, reg, td, "tx", 1e-20, cpg, 500):
                    if r["ident"] >= a.tx_identity and r["len"] >= a.tx_length: hit[r["s"]].add(r["q"])
                span_tx = sorted(k for k, v in hit.items() if len(v) == 2)
                shutil.rmtree(td, ignore_errors=True)
            nbrow = {"genome": g, "full_locus": f"{F['contig']}:{F['lo'] + 1}-{F['hi']}{F['strand']}",
                     "neighbour_locus": f"{N['contig']}:{N['lo'] + 1}-{N['hi']}{N['strand']}",
                     "neighbour_bits": f"{N['bits']:.1f}", "order": order, "distance": gap,
                     "models_full": ", ".join(mf), "models_neighbour": ", ".join(mn),
                     "merged_model": ", ".join(sorted(set(mf) & set(mn))),
                     "transcripts_spanning": len(span_tx) if tx and os.path.exists(tx) else "no transcripts",
                     "transcript_examples": ", ".join(span_tx[:5])}
        for x in os.listdir(d):
            if x.startswith("genome.fa"): os.remove(os.path.join(d, x))
        res = {"sites": out, "neighbour": nbrow}
        json.dump(res, open(done + ".tmp", "w")); os.replace(done + ".tmp", done)
        return res

    genomes = [g for g in a.genomes if g in sheet and os.path.exists(sheet[g]["fasta"])]
    missing = [g for g in a.genomes if g not in genomes]
    if missing: print(f"genomes not found, left out: {' '.join(missing)}", file=sys.stderr)
    results = {}
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
        for g, r in zip(genomes, ex.map(genome, genomes)):
            results[g] = r
            print(f"{g}: {len(r['sites'])} sites", flush=True)
    sites = [r for g in genomes for r in results[g]["sites"]]
    names = {"sec": a.sec_gene, "full": a.full_gene, "selenocysteine copy": a.sec_gene, "full-length copy": a.full_gene}
    for r in sites:
        if r.get("codon_source"): r["codon_source"] = names.get(r["codon_source"], r["codon_source"])
    hg = to_human({f"{r['genome']}|{r['site']}": r["pith_sequence"] for r in sites if len(r["pith_sequence"]) >= 20}, "sites_to_human")
    for r in sites:
        gl = hg.get(f"{r['genome']}|{r['site']}", [])
        r["best_human_gene"], r["best_human_bits"], r["next_human_gene"] = gene_txt(gl)
        r["human_margin"] = gl[0][2] - gl[1][2] if len(gl) > 1 else 1e9
    cols = ["genome", "site", "contig", "strand", "from", "to", "structure", "thioredoxin_coverage", "pith_coverage",
            "bits_sec", "identity_sec", "bits_full", "identity_full", "best_human_gene", "best_human_bits",
            "next_human_gene", "codon", "residue", "codon_source", "codon_note", "codons_after", "stop",
            "secis_oikopleura", "secis_RF00031", "upstream_to_sequence_end", "annotation", "pith_sequence"]
    with open(f"{a.prefix}_sites.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in sites: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    cols = ["genome", "full_locus", "neighbour_locus", "neighbour_bits", "order", "distance", "models_full",
            "models_neighbour", "merged_model", "transcripts_spanning", "transcript_examples"]
    with open(f"{a.prefix}_neighbours.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for g in genomes:
            r = results[g]["neighbour"]
            if r: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # 3. proteomes
    prows, tree_seqs = [], {}
    for g in a.proteomes + [x for x in a.dioica if x not in a.proteomes]:
        if g in a.dioica:
            if g not in sheet: continue
            ann = annotated_proteins(sheet[g]["fasta"], sheet[g]["annotation"])
            sloc = copies[g]["sec_locus"] if g in copies else None   # the annotated model of the selenocysteine copy
            seqs = {k: r["protein"] for k, r in ann.items() if r["protein"] and
                    not overlaps(sloc, r["contig"], r["strand"], r["lo"], r["hi"])}
            desc = {k: r["gene"] for k, r in ann.items()}
            for k, s in read_fasta(f"{a.gold}/{g}.proteins.faa"):
                n = k.split("|")[0]
                if n == a.sec_gene: seqs[f"gold:{n}"], desc[f"gold:{n}"] = s.upper().rstrip("*"), "selenocysteine copy (gold standard)"
            gene = {k: desc[k] for k in seqs}
        else:
            path = f"{a.data}/{g}/protein.faa.gz"
            if not os.path.exists(path): continue
            seqs, desc = read_proteome(path)
            gene = genes_of(f"{a.data}/{g}/annotation.gff3.gz")
        pdb = Db(a.work, f"{g}.proteome", ((k, clean(s)) for k, s in seqs.items()), True)
        rows = blast("blastp", pdb, queries, a.work, f"{g}.proteome", a.evalue, a.cpus, a.max_hits)
        shutil.rmtree(pdb.dir, ignore_errors=True)
        hits = {}
        for r in rows:
            if r["q"] == "full" and r["qe"] <= trx_end: continue   # thioredoxin part only
            hits[r["s"]] = max(hits.get(r["s"], 0), r["bits"])
        hg = to_human({k: clean(seqs[k]) for k in hits}, f"{g}.proteome_to_human")
        for k in sorted(hits, key=lambda k: -hits[k]):
            t = seqs[k]
            fb, sb = pair(al, full_q, t), pair(al, sec_q, t)
            pos, pn = mapped(sb, upos, a.extend, len(t))
            bg, bb, ng = gene_txt(hg.get(k, []))
            trx_cov, pith_cov = x_coverage(fb, 0, trx_end), x_coverage(fb, trx_end, len(full_q))
            structure = ("full-length" if trx_cov >= 0.3 and pith_cov >= 0.3 else "PITH only" if pith_cov >= 0.3
                         else "thioredoxin only" if trx_cov >= 0.3 else "partial")
            prows.append([g, k, desc.get(k, ""), gene.get(k, ""), len(t), f"{hits[k]:.1f}", structure, f"{trx_cov:.2f}",
                          f"{pith_cov:.2f}", bg, bb, ng, t[pos] if pos is not None else "", (pos + 1) if pos is not None else "",
                          pn, len(t) - pos - 1 if pos is not None else ""])
            keep = bg in (a.human_gene, a.outgroup_gene) or k.startswith("gold:")
            if g in a.dioica and g != a.tree_dioica: keep = False
            if keep and sb:
                region = t[sb[0][1][0]:]
                if len(region) < a.min_length: continue
                label = gene.get(k, "") if g not in a.dioica else (a.sec_gene if k.startswith("gold:") else
                                                                    a.full_gene if k == copies[g]["full_id"] else bg + "-like")
                key = (g, label or k) if g not in a.dioica else (g, k)
                res = t[pos] if pos is not None else "-"
                gl = hg.get(k, [])
                margin = gl[0][2] - gl[1][2] if len(gl) > 1 else 1e9
                rank = (k.startswith("NP_"), len(region))   # one protein per gene: curated RefSeq first, then the longest
                if key not in tree_seqs or rank > tree_seqs[key][4]:
                    tree_seqs[key] = (safe(f"{g}|{label or k}|{k}|{bg or 'none'}|{structure}|{res}"), clean(region), bg, margin, rank)
        print(f"proteome {g}: {len(hits)} proteins", flush=True)
    cols = ["genome", "protein", "description", "gene", "length", "bits", "structure", "thioredoxin_coverage",
            "pith_coverage", "best_human_gene", "best_human_bits", "next_human_gene", "residue_at_sec", "position",
            "residue_note", "residues_after"]
    with open(f"{a.prefix}_proteins.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in prows: o.write("\t".join(str(x) for x in r) + "\n")

    # 4. tree
    if a.no_tree: return
    with_proteome = set(a.proteomes) | set(a.dioica)
    for r in sites:
        if r["genome"] in with_proteome and (r["annotation"] or r["genome"] in a.dioica): continue
        if len(r["pith_sequence"]) < a.min_length or r["best_human_gene"] not in (a.human_gene, a.outgroup_gene): continue
        name = safe(f"{r['genome']}|site{r['site']}|{r['contig']}|{r['best_human_gene']}|{r['structure']}|{(r.get('codon') if r.get('residue') == '*' else r.get('residue')) or '-'}")
        tree_seqs[(r["genome"], f"site{r['site']}")] = (name, r["pith_sequence"], r["best_human_gene"], r["human_margin"], None)
    os.makedirs(a.tree_dir, exist_ok=True)
    base = os.path.join(a.tree_dir, "txnl1_pith")
    faa = "".join(f">{v[0]}\n{v[1]}\n" for v in tree_seqs.values())
    if os.path.exists(f"{base}.treefile") and os.path.exists(f"{base}.faa") and open(f"{base}.faa").read() == faa:
        print("tree: same sequences as the existing tree, which is reused", flush=True)
    else:
        for ext in ("faa", "aln", "trimmed", "treefile", "iqtree", "log", "ckp.gz", "rooted.nwk"):
            if os.path.exists(f"{base}.{ext}"): os.remove(f"{base}.{ext}")
        with open(f"{base}.faa", "w") as o: o.write(faa)
        with open(f"{base}.aln", "w") as o:
            run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", str(a.cpus), "--quiet", f"{base}.faa"], stdout=o)
        run(["clipkit", f"{base}.aln", "-o", f"{base}.trimmed", "-m", "gappy", "-g", "0.9"])
        run(["iqtree2", "-s", f"{base}.trimmed", "-m", "MFP", "-B", "1000", "-T", str(a.cpus), "--prefix", base,
             "--quiet", "--seed", "1"])
        for ext in ("bionj", "ckp.gz", "mldist", "model.gz", "splits.nex", "uniqueseq.phy", "contree"):
            if os.path.exists(f"{base}.{ext}"): os.remove(f"{base}.{ext}")
    from Bio import Phylo
    tree = Phylo.read(f"{base}.treefile", "newick")
    for c in tree.get_nonterminals():
        if c.name and re.fullmatch(r"\d+(\.\d+)?", c.name): c.confidence, c.name = float(c.name), None
    # supports belong to bipartitions; Bio.Phylo keeps them on nodes, which rerooting reverses along the new root path,
    # so they are read back by bipartition after rerooting
    all_tips = frozenset(x.name for x in tree.get_terminals())
    support = {}
    for c in tree.get_nonterminals():
        if c.confidence is None: continue
        s = frozenset(x.name for x in c.get_terminals())
        support[s] = support[all_tips - s] = c.confidence
    gene_of_tip = {v[0]: v[2] for v in tree_seqs.values()}
    margin_of_tip = {v[0]: v[3] for v in tree_seqs.values()}
    # root: the clade of the sequences whose best human gene is --outgroup-gene (by at least --root-margin bits), found
    # after a provisional rooting on a sequence of the ingroup; else the human outgroup sequence alone
    tips = tree.get_terminals()
    out = [c for c in tips if gene_of_tip.get(c.name) == a.outgroup_gene and margin_of_tip.get(c.name, 0) >= a.root_margin]
    inn = [c for c in tips if gene_of_tip.get(c.name) == a.human_gene and margin_of_tip.get(c.name, 0) >= a.root_margin]
    rooted = ""
    if out:
        if inn: tree.root_with_outgroup(inn[0])
        ca = tree.common_ancestor(out) if len(out) > 1 else out[0]
        if ca is not tree.root and all(gene_of_tip.get(x.name) == a.outgroup_gene for x in ca.get_terminals()):
            rooted = f"{a.outgroup_gene} clade ({len(ca.get_terminals())} sequences)"
            tree.root_with_outgroup(ca)
        else:
            o1 = next((c for c in out if c.name.startswith("Hsap")), out[0])
            tree.root_with_outgroup(o1)
            rooted = f"{o1.name} ({a.outgroup_gene} sequences not monophyletic)"
    for c in tree.get_nonterminals():
        c.confidence = support.get(frozenset(x.name for x in c.get_terminals()))
    Phylo.write(tree, f"{base}.rooted.nwk", "newick")
    # the clade of each O. dioica copy: its sister group and support
    trows = []
    for c in tree.get_terminals():
        if not c.name.startswith(safe(a.tree_dioica) + "_"): continue
        path = [tree.root] + tree.get_path(c)
        parent = path[-2]
        grand = path[-3] if len(path) > 2 else None
        sisters = [t.name for t in parent.get_terminals() if t is not c]
        nxt = [t.name for t in grand.get_terminals() if t not in parent.get_terminals()] if grand is not None else []
        trows.append([c.name, parent.confidence if parent.confidence is not None else "", len(sisters), ", ".join(sisters),
                      grand.confidence if grand is not None and grand.confidence is not None else "", len(nxt), ", ".join(nxt)])
    with open(f"{a.prefix}_tree.tsv", "w") as o:
        o.write(f"tip\tsupport_of_smallest_clade\tsisters\tsister_tips\tsupport_of_next_clade\tnext_clade_added\t"
                f"next_clade_added_tips (tree rooted on {rooted or 'nothing'}; supports read by bipartition)\n")
        for r in trows: o.write("\t".join(str(x) for x in r) + "\n")
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(tree.get_terminals())
    fig = plt.figure(figsize=(11, max(4, 0.22 * n)))
    ax = fig.add_subplot(1, 1, 1)
    colours = {c.name: ("#D55E00" if safe(a.sec_gene) in c.name else "#0072B2" if c.name.startswith(safe(a.tree_dioica)) else
                        "#777777" if gene_of_tip.get(c.name) == a.outgroup_gene and margin_of_tip.get(c.name, 0) >= a.root_margin
                        else "black") for c in tree.get_terminals()}
    Phylo.draw(tree, axes=ax, do_show=False, label_colors=colours, label_func=lambda c: c.name if c.is_terminal() else "",
               branch_labels=lambda cl: f"{cl.confidence:.0f}" if cl.confidence is not None and not cl.is_terminal() else "")
    ax.set_title(f"PITH regions of TXNL1 and {a.outgroup_gene}; rooted on {rooted or 'none'}; IQ-TREE, ultrafast bootstrap "
                 f"support; orange: {a.sec_gene}, blue: other {a.tree_dioica} proteins, grey: {a.outgroup_gene}",
                 fontsize=8, loc="left")
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.set_yticks([]); ax.set_ylabel("")
    fig.tight_layout()
    os.makedirs(os.path.dirname(a.figure) or ".", exist_ok=True)
    for ext in ("png", "pdf"): fig.savefig(f"{a.figure}.{ext}", dpi=200)
    print(f"tree: {n} sequences -> {base}.rooted.nwk, {a.figure}.png", flush=True)


if __name__ == "__main__":
    main()
