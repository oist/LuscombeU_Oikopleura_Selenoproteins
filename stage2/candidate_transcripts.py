#!/usr/bin/env python3
"""Codons of candidate selenoproteins located through transcripts (step 9, absence tests, 2026-10-30).

For SELENOI, the tblastn search of candidate_genes.py did not place the candidate codon in OSKA2016 (its last coding
exon is short). Here the transcripts of the three O. dioica assemblies that encode the candidate (rows of type
"transcriptome" of the absence test of the gene, results/tables/stage2_absence_<gene>.candidates.tsv) are aligned to
their assembly: the minimap2 alignments of curation/gold_standard.py (-x splice:hq -un, primary alignments; read from
--paf-dir, or computed for these transcripts only when that file is missing), with the junctions placed as there
(refine_junctions). For each transcript:
  - the open reading frame read through TGA (only TAA and TAG end it) that aligns best to the candidate protein
    (curation/candidates.faa; local alignment, BLOSUM62, U scored as C), and the codon aligned to the first U of the
    candidate protein, the "candidate codon" (for SELENOI, the U aligned to the selenocysteine of vertebrate
    SELENOI); "outside" when the alignment does not reach it (an isoform with another C-terminus);
  - the C-terminal residues of the protein when TGA is read as a stop codon;
  - the genome position of the candidate codon (1-based start on the forward strand, the convention of
    curation/candidates.tsv), the codon read in the assembly on the gene strand, and the exon of the alignment that
    holds it (counted from the 5' end of the transcript) with its span;
  - the stop codon (TGA, TAA or TAG) that follows the candidate codon in the transcript, the length of transcript
    after it, and the best SECIS hits of the Oikopleura model and of RF00031 there (cmsearch, top strand).
Outputs: --table (one row per transcript) and --positions (one row per gene, assembly and genome position of a
candidate codon read as TGA in at least one transcript, with the transcripts), which candidate_genes.py
(--transcripts) uses for the rule in place of the tblastn positions.
"""
import argparse, collections, csv, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "curation"))
from gold_standard import parse_cigar, refine_junctions, walk          # noqa: E402
from candidate_genes import Secis, read_fasta, rc, CODE, DIOICA        # noqa: E402

STOPS_RT = {"TAA", "TAG"}


def tr(s, u=False):
    return "".join("U" if (u and s[i:i + 3] == "TGA") else CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    return Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                 open_gap_score=-10, extend_gap_score=-0.5)


def clean(x): return "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else ("C" if c == "U" else "X") for c in x)


def best_frame(nt, prot, al):
    """Best read-through ORF (TGA as U) of the sequence against the protein: (score, nucleotide offset of the ORF's
    first codon, ORF protein, alignment)."""
    best = None
    for fr in range(3):
        rt = tr(nt[fr:], u=True)
        start = 0
        for seg in rt.split("*"):
            if len(seg) >= 30:
                try: aln = al.align(clean(seg), clean(prot))[0]
                except (IndexError, ValueError): aln = None
                if aln is not None and (best is None or aln.score > best[0]):
                    best = (aln.score, fr + 3 * start, seg, aln)
            start += len(seg) + 1
    return best


def residue_at(aln, i):
    """Index in the first sequence of the residue aligned to position i of the second, or None."""
    for (c0, c1), (r0, r1) in zip(*aln.aligned):
        if r0 <= i < r1: return c0 + i - r0
    return None


def read_paf(paf, names):
    out = {}
    with open(paf) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if f[0] not in names: continue
            tags = {x.split(":", 2)[0]: x.split(":", 2)[2] for x in f[12:] if x.count(":") >= 2}
            if tags.get("tp", "P") != "P" or "cg" not in tags: continue
            out[f[0]] = {"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "t": f[5],
                         "ts": int(f[7]), "te": int(f[8]), "cigar": tags["cg"]}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", nargs="+", default=["SELENOI"])
    ap.add_argument("--absence", default="results/tables/stage2_absence_{gene}.candidates.tsv")
    ap.add_argument("--proteins", default="curation/candidates.faa")
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--paf-dir", default="work/stage2/gold_standard")
    ap.add_argument("--model", default="results/models/oikopleura_secis.cm")
    ap.add_argument("--rfam", default="data/references/RF00031.cm")
    ap.add_argument("--window", type=int, default=8, help="junction placement window (as curation/gold_standard.py)")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--work", default="work/stage2/candidate_transcripts")
    ap.add_argument("--table", default="results/tables/stage2_candidate_transcripts.tsv")
    ap.add_argument("--positions", default="results/tables/stage2_candidate_transcripts_positions.tsv")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    prot = dict(read_fasta(a.proteins))
    al = aligner()
    secis = Secis({"oikopleura": a.model, "RF00031": a.rfam}, a.work)
    rows, windows = [], {}
    for gene in a.genes:
        path = a.absence.format(gene=gene)
        if gene not in prot or not os.path.exists(path):
            print(f"{gene}: no candidate protein or no absence table ({path})", flush=True); continue
        upos = prot[gene].find("U")
        tx = collections.defaultdict(dict)   # genome -> transcript -> strand of the hits on the transcript
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["type"] != "transcriptome" or r["target"] not in DIOICA: continue
            _, name, _, strand = r["locus"].split(":")
            tx[r["target"]][name] = strand
        for g in DIOICA:
            if not tx[g] or g not in sheet: continue
            seqs = {k: v.upper() for k, v in read_fasta(sheet[g]["transcripts"]) if k in tx[g]}
            paf = f"{a.paf_dir}/{g}.transcripts.paf"
            if not os.path.exists(paf):   # these transcripts only
                paf = f"{a.work}/{g}.{gene}.paf"
                with open(f"{a.work}/{g}.{gene}.fa", "w") as o:
                    for k, v in seqs.items(): o.write(f">{k}\n{v}\n")
                with open(paf, "w") as o:
                    subprocess.run(["minimap2", "-x", "splice:hq", "-un", "--secondary=no", "-c", "-t", str(a.cpus),
                                    sheet[g]["fasta"], f"{a.work}/{g}.{gene}.fa"], stdout=o, check=True,
                                   stderr=subprocess.DEVNULL)
            alns = read_paf(paf, set(seqs))
            contigs = {k: v.upper() for k, v in read_fasta(sheet[g]["fasta"]) if k in {x["t"] for x in alns.values()}}
            for name in sorted(seqs):
                raw, st = seqs[name], tx[g][name]
                nt = raw if st == "+" else rc(raw)
                row = {"gene": gene, "genome": g, "transcript": name, "transcript_length": len(raw),
                       "transcript_strand": st}
                b = best_frame(nt, prot[gene], al)
                if b is None:
                    rows.append(dict(row, codon_transcript="no frame")); continue
                score, orf0, seg, aln = b
                i = residue_at(aln, upos)
                # protein with TGA as a stop codon: the standard frame from the first aligned residue
                first = orf0 + 3 * aln.aligned[0][0][0]
                std = tr(nt[first:]).split("*")[0]
                row.update({"alignment_score": f"{score:.0f}", "c_terminus_tga_as_stop": std[-15:]})
                if i is None:
                    rows.append(dict(row, codon_transcript="outside")); continue
                off = orf0 + 3 * i
                codon = nt[off:off + 3]
                k = off + 3
                while k + 3 <= len(nt) and nt[k:k + 3] not in {"TGA", "TAA", "TAG"}: k += 3
                stop_end = k + 3 if k + 3 <= len(nt) else None
                row.update({"codon_transcript": codon, "context": tr(nt[max(off % 3, off - 24):off]) + "[" +
                            tr(nt[off:off + 3], u=True) + "]" + tr(nt[off + 3:off + 45], u=True),
                            "codons_to_stop": (k - off - 3) // 3,
                            "stop": nt[k:k + 3] if stop_end else "",
                            "after_stop_nt": len(nt) - stop_end if stop_end else ""})
                if stop_end and len(nt) > stop_end: windows[f"{gene}|{g}|{name}"] = nt[stop_end:]
                x = alns.get(name)
                if x is None:
                    rows.append(dict(row, contig="not aligned")); continue
                gene_strand = x["strand"] if st == "+" else ("-" if x["strand"] == "+" else "+")
                fwd = raw if x["strand"] == "+" else rc(raw)
                cg, moved = refine_junctions(x, fwd, contigs[x["t"]], gene_strand, a.window)
                w = walk(dict(x, cigar=cg))
                L = len(raw)
                gpos = []
                for j in range(off, off + 3):
                    rawi = j if st == "+" else L - 1 - j
                    fi = rawi if x["strand"] == "+" else L - 1 - rawi
                    gpos.append(w["q2g"].get(fi))
                blocks = w["blocks"] if gene_strand == "+" else w["blocks"][::-1]
                row.update({"contig": x["t"], "gene_strand": gene_strand,
                            "aligned_transcript": f"{x['qs'] + 1}-{x['qe']}", "exons": len(blocks),
                            "junctions_moved": moved})
                if None in gpos:
                    rows.append(dict(row, codon_genome="not aligned")); continue
                lo = min(gpos)
                gc = contigs[x["t"]][lo:lo + 3]
                gc = gc if gene_strand == "+" else rc(gc)
                ex = next((n for n, (s, e) in enumerate(blocks, 1) if s <= lo and lo + 3 <= e), "split")
                row.update({"codon": lo + 1, "codon_genome": gc if max(gpos) - lo == 2 else gc + " (split)",
                            "codon_exon": ex,
                            "exon_span": f"{blocks[ex - 1][0] + 1}-{blocks[ex - 1][1]}" if ex != "split" else ""})
                rows.append(row)
    best = secis.best(windows, a.cpus, "transcripts")
    for r in rows:
        b = best.get(f"{r['gene']}|{r['genome']}|{r['transcript']}", {})
        r["secis_oikopleura"] = b["oikopleura"][0] if "oikopleura" in b else ""
        r["secis_RF00031"] = b["RF00031"][0] if "RF00031" in b else ""
    cols = ["gene", "genome", "transcript", "transcript_length", "transcript_strand", "alignment_score",
            "codon_transcript", "context", "codons_to_stop", "stop", "after_stop_nt", "secis_oikopleura",
            "secis_RF00031", "c_terminus_tga_as_stop", "contig", "gene_strand", "aligned_transcript", "exons",
            "junctions_moved", "codon", "codon_genome", "codon_exon", "exon_span"]
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    pos = collections.defaultdict(list)
    for r in rows:
        if r.get("codon_transcript") == "TGA" and isinstance(r.get("codon"), int):
            pos[(r["gene"], r["genome"], r["contig"], r["codon"], r["gene_strand"])].append(r["transcript"])
    with open(a.positions, "w") as o:
        o.write("gene\tgenome\tcontig\tcodon\tstrand\ttranscripts\n")
        for k, v in sorted(pos.items(), key=lambda kv: (kv[0][0], DIOICA.index(kv[0][1]), -len(kv[1]))):
            o.write("\t".join(map(str, k)) + "\t" + ",".join(v) + "\n")
    for r in rows:
        print(r["genome"], r["transcript"], r.get("codon_transcript", ""), r.get("contig", ""), r.get("codon", ""),
              r.get("codon_genome", ""), r.get("secis_oikopleura", ""), flush=True)


if __name__ == "__main__":
    main()
