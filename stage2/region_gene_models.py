#!/usr/bin/env python3
"""Gene models of a selenoprotein in genome regions, by protein-to-genome alignment (2026-10-31; first use: the
SELENOK genes of the other Oikopleura species, which the final SELENOK profile does not find, to become members of
that profile).

Input: regions written by stage2/extract_regions.py (--regions; header <name>|<genome>|<contig>:<start>-<end>|
<strand>|offset=<n>, sequence on the gene strand) and query proteins with U for the selenocysteine (--queries; for
SELENOK, the gold-standard proteins of the three O. dioica assemblies). For each region:
  1. exonerate --model protein2genome (U read as X) of each query against the region, sense strand only; the
     alignment of best score is kept (its coding exons, from the GFF of exonerate);
  2. the frame is extended from the end of the alignment to the next TAA or TAG, TGA being read as selenocysteine
     (U): exonerate, like miniprot, stops before an in-frame TGA near the C-terminus; the extension is kept only if it
     is at most --max-extension codons;
  3. the protein is the translation of the exons and the extension (TGA as U); its 5' end is "complete" when the
     first codon is ATG; the codon aligned to the U of the best query (local alignment, BLOSUM62, U scored as C) and
     whether it is the TGA of the screen (within 2 nt of the reference offset of the region, so that the leftmost
     coordinate of a codon on the minus strand also matches) are reported, with the identity to the query.
Outputs: --faa (proteins, header <family>|<genome>|<source>, the format of build_profiles.py --curated) and
--table.
"""
import argparse, os, shutil, subprocess, tempfile

CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def tr_u(s): return "".join("U" if s[i:i + 3] == "TGA" else CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].strip(), []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def parse(out):
    """Alignments of exonerate output (GFF block, then the ryo line of the same alignment)."""
    alns, exons = [], []
    for l in out.splitlines():
        if l.startswith("##ryo "):
            _, q, score, strand = l.split()
            alns.append((int(score), q, strand, sorted(exons)))
            exons = []
        elif not l.startswith("#") and l.count("\t") >= 8:
            f = l.split("\t")
            if f[2] == "cds": exons.append((int(f[3]) - 1, int(f[4])))
    return alns


def align_u(prot, ref):
    """Residue of prot aligned to the first U of ref, its index, and the identity over the aligned columns."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    clean = lambda x: "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else ("C" if c == "U" else "X") for c in x)
    al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-10, extend_gap_score=-0.5)
    aln = al.align(clean(prot), clean(ref))[0]
    i = ref.find("U")
    res, idx = "outside", None
    same = tot = 0
    for (c0, c1), (r0, r1) in zip(*aln.aligned):
        for k in range(c1 - c0):
            tot += 1; same += prot[c0 + k] == ref[r0 + k]
        if r0 <= i < r1: idx = c0 + i - r0; res = prot[idx]
    return res, idx, (100 * same / tot if tot else 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--queries", nargs="+", required=True, help="FASTA files; entries of --family only (header "
                    "starting with <family>|) when --family is given")
    ap.add_argument("--family", default="")
    ap.add_argument("--max-extension", type=int, default=30)
    ap.add_argument("--faa", required=True)
    ap.add_argument("--table", required=True)
    a = ap.parse_args()
    queries = {}
    for f in a.queries:
        for k, s in read_fasta(f):
            if a.family and k.split()[0].split("|")[0] != a.family: continue
            queries[k.split()[0].replace("|", "_")] = s
    if not queries: raise SystemExit(f"no query of family {a.family!r} in {', '.join(a.queries)}")
    work = tempfile.mkdtemp(prefix="region_models_")
    qf = f"{work}/q.faa"
    with open(qf, "w") as o:
        for k, s in queries.items(): o.write(f">{k}\n{s.replace('U', 'X')}\n")
    rows, faa = [], []
    for name, seq in read_fasta(a.regions):
        seq = seq.upper()
        rname, genome, loc, strand = name.split("|")[:4]
        off = int(name.split("offset=")[1].split("|")[0]) if "offset=" in name else None
        tf = f"{work}/t.fa"
        open(tf, "w").write(f">{rname}\n{seq}\n")
        out = subprocess.run(["exonerate", "--model", "protein2genome", "--showalignment", "no", "--showvulgar", "no",
                              "--showtargetgff", "yes", "--ryo", "##ryo %qi %s %tS\n", "--bestn", "1", qf, tf],
                             capture_output=True, text=True, check=True).stdout
        alns = [x for x in parse(out) if x[2] == "+" and x[3]]
        r = {"region": rname, "genome": genome, "location": loc, "strand": strand}
        if not alns:
            rows.append(dict(r, status="no alignment")); continue
        score, q, _, exons = max(alns)
        cds = "".join(seq[s:e] for s, e in exons)
        end = exons[-1][1]
        k = end
        while k + 3 <= len(seq) and seq[k:k + 3] not in ("TAA", "TAG"): k += 3
        ext = (k - end) // 3
        stop = seq[k:k + 3] if k + 3 <= len(seq) else "none"
        if ext <= a.max_extension and stop != "none": cds += seq[end:k]
        else: ext = 0
        prot = tr_u(cds)
        res, idx, ident = align_u(prot, queries[q])
        # genome position of the residue aligned to the U
        tga = ""
        if idx is not None:
            nt = 3 * idx
            for s, e in exons + ([(end, k)] if ext else []):
                if nt < e - s: pos = s + nt; break
                nt -= e - s
            else: pos = None
            if pos is not None:
                tga = f"{seq[pos:pos + 3]} at {pos}" + (" (screen codon)" if off is not None and abs(pos - off) <= 2
                                                       else "")
        r.update({"status": "model", "query": q, "exonerate_score": score, "exons": len(exons) + (1 if ext else 0),
                  "coding_exons": ";".join(f"{s + 1}-{e}" for s, e in exons), "extension_codons": ext,
                  "stop": stop, "length": len(prot), "complete_5prime": "yes" if cds[:3] == "ATG" else "no",
                  "residue_at_query_U": res, "codon_at_query_U": tga, "identity_to_query": f"{ident:.0f}",
                  "protein": prot})
        rows.append(r)
        if res == "U":
            faa.append((f"{a.family or q.split('_')[0]}|{genome}|{rname} {loc} ({strand}), exonerate protein2genome "
                        f"of {q} and extension through TGA", prot))
    os.makedirs(os.path.dirname(a.faa) or ".", exist_ok=True)
    with open(a.faa, "w") as o:
        for h, s in faa: o.write(f">{h}\n{s}\n")
    cols = ["region", "genome", "location", "strand", "status", "query", "exonerate_score", "exons", "coding_exons",
            "extension_codons", "stop", "length", "complete_5prime", "residue_at_query_U", "codon_at_query_U",
            "identity_to_query", "protein"]
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    for r in rows: print("\t".join(str(r.get(c, "")) for c in cols))
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
