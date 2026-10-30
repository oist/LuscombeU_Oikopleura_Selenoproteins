#!/usr/bin/env python3
"""Gene models of a selenoprotein in genome regions, by protein-to-genome alignment (2026-10-31; first use: the
SELENOK genes of the other Oikopleura species, which the final SELENOK profile does not find, to become members of
that profile).

Input: regions written by stage2/extract_regions.py (--regions; header <name>|<genome>|<contig>:<start>-<end>|
<strand>|offset=<n>, sequence on the gene strand, the offset pointing at the selenocysteine codon found by a screen,
within 2 nt) and query proteins with U for the selenocysteine (--queries; for SELENOK, the gold-standard proteins of
the three O. dioica assemblies). For each region:
  1. genewise (wise2; selenocysteine read as cysteine) of each query against the sense strand of the region; the
     model of best score is kept (its coding exons). exonerate protein2genome found fewer exons in these regions
     (no alignment in two of four), and like miniprot, these tools stop before an in-frame TGA near the C-terminus.
  2. C-terminus. When the model ends before the screen codon (--screen-window nt or more upstream of it), the exon
     holding that codon is joined to the model: a donor (GT, GC or GA, the O. dioica classes; Denoeud et al. 2010)
     within --donor-window nt of the end of the last exon, an acceptor (AG) at most --acceptor-span nt upstream of the
     codon, at least 30 nt after the donor, such that the codon is in frame and no TAA or TAG comes before it; among
     these, the protein of best global alignment score to the query (BLOSUM62, end gaps free; +2 for GT). Otherwise
     the frame is extended from the end of the model to the next TAA or TAG, TGA read as selenocysteine, if at most
     --max-extension codons away.
  3. N-terminus: when the model does not start with ATG, an ATG at most --start-window codons upstream in frame,
     without stop codon in between, starts it, or else an ATG among its first --start-downstream codons (genewise
     can begin one codon early); otherwise the 5' end is "incomplete" (the first exon was not found).
  4. The protein (TGA as U), the residue aligned to the U of the query (local alignment), whether it is the screen
     codon, and the identity to the query.
Outputs: --faa (proteins whose residue at the query U is U; header <family>|<genome>|<source>, the format of
build_profiles.py --curated) and --table (exons in region and contig coordinates).
"""
import argparse, os, shutil, subprocess, sys, tempfile

CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
STOPS = ("TAA", "TAG")


def tr_u(s): return "".join("U" if s[i:i + 3] == "TGA" else CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].strip(), []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def clean(x): return "".join(c if c in "ACDEFGHIKLMNPQRSTVWY" else ("C" if c == "U" else "X") for c in x)


def aligners():
    from Bio import Align
    from Bio.Align import substitution_matrices
    m = substitution_matrices.load("BLOSUM62")
    loc = Align.PairwiseAligner(mode="local", substitution_matrix=m, open_gap_score=-10, extend_gap_score=-0.5)
    glo = Align.PairwiseAligner(mode="global", substitution_matrix=m, open_gap_score=-10, extend_gap_score=-0.5)
    glo.end_gap_score = 0
    return loc, glo


def wise_env():
    env = dict(os.environ)
    if "WISECONFIGDIR" not in env:
        for d in (os.path.join(sys.prefix, "share", "wise2", "wisecfg"), "/usr/share/wise",
                  os.path.join(os.environ.get("CONDA_PREFIX", sys.prefix), "share", "wise2", "wisecfg")):
            if os.path.exists(os.path.join(d, "BLOSUM62.bla")) or os.path.exists(os.path.join(d, "gene.stat")):
                env["WISECONFIGDIR"] = d; break
    return env


def genewise(query, target, env):
    """(score, [(start0, end0) of the coding exons]) of genewise on the sense strand."""
    out = subprocess.run(["genewise", "-gff", "-tfor", "-silent", "-quiet", query, target], capture_output=True,
                         text=True, env=env).stdout
    score, exons = None, []
    for l in out.splitlines():
        f = l.split("\t")
        if len(f) < 9: continue
        if f[2] == "match" and score is None: score = float(f[5])
        elif f[2] == "cds" and score is not None and len(exons) < 1000:
            if f[8].endswith("prediction-1"): exons.append((int(f[3]) - 1, int(f[4])))
    return score, sorted(exons)


def join_screen_exon(seq, exons, tga, ref, glo, a):
    """Best (score, donor, acceptor, splice, cds_end) joining the last exon to the exon holding the TGA at tga."""
    head = "".join(seq[s:e] for s, e in exons[:-1])
    s_last, ge = exons[-1]
    best = None
    for d in range(max(s_last + 3, ge - a.donor_window), min(ge + a.donor_window, tga - 30)):
        dn = seq[d:d + 2]
        if dn not in ("GT", "GC", "GA"): continue
        for acc in range(max(d + 30, tga - a.acceptor_span), tga - 2):
            if seq[acc - 2:acc] != "AG": continue
            body = head + seq[s_last:d]
            tpos = len(body) + tga - acc
            if tpos % 3: continue
            cds = body + seq[acc:]
            k = tpos + 3
            while k + 3 <= len(cds) and cds[k:k + 3] not in STOPS: k += 3
            if k + 3 > len(cds): continue
            p = tr_u(cds[:k])
            if "*" in p or p[tpos // 3] != "U": continue
            sc = glo.score(clean(p), clean(ref)) + (2 if dn == "GT" else 0)
            if best is None or sc > best[0]: best = (sc, d, acc, f"{dn}-{seq[acc - 2:acc]}", acc + (k - len(body)))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--queries", nargs="+", required=True, help="FASTA files; with --family, only the entries whose "
                    "first word, up to '|', is the family")
    ap.add_argument("--family", default="")
    ap.add_argument("--max-extension", type=int, default=30)
    ap.add_argument("--screen-window", type=int, default=6)
    ap.add_argument("--donor-window", type=int, default=40)
    ap.add_argument("--acceptor-span", type=int, default=400)
    ap.add_argument("--start-window", type=int, default=10)
    ap.add_argument("--start-downstream", type=int, default=2)
    ap.add_argument("--faa", required=True)
    ap.add_argument("--table", required=True)
    a = ap.parse_args()
    queries = {}
    for f in a.queries:
        for k, s in read_fasta(f):
            if a.family and k.split()[0].split("|")[0] != a.family: continue
            queries["_".join(k.split()[0].split("|")[:2])] = s
    if not queries: raise SystemExit(f"no query of family {a.family!r} in {', '.join(a.queries)}")
    loc, glo = aligners()
    env = wise_env()
    work = tempfile.mkdtemp(prefix="region_models_")
    rows, faa = [], []
    for name, seq in read_fasta(a.regions):
        seq = seq.upper()
        rname, genome, place, strand = name.split("|")[:4]
        contig, span = place.rsplit(":", 1)
        s0, e0 = map(int, span.split("-"))
        off = int(name.split("offset=")[1].split("|")[0]) if "offset=" in name else None
        tf = f"{work}/t.fa"
        open(tf, "w").write(f">{rname}\n{seq}\n")
        best = None
        for q, p in queries.items():
            qf = f"{work}/q.faa"
            open(qf, "w").write(f">{q}\n{clean(p)}\n")
            score, exons = genewise(qf, tf, env)
            if exons and (best is None or score > best[0]): best = (score, q, exons)
        r = {"region": rname, "genome": genome, "location": place, "strand": strand}
        if best is None:
            rows.append(dict(r, status="no model")); continue
        score, q, exons = best
        ref = queries[q]
        # screen codon: the TGA within 2 nt of the offset
        tga = next((x for x in ((off - 2, off - 1, off, off + 1, off + 2) if off is not None else ())
                    if 0 <= x and seq[x:x + 3] == "TGA"), None)
        how, splice = "", ""
        if tga is not None and exons[-1][1] <= tga - a.screen_window:
            j = join_screen_exon(seq, exons, tga, ref, glo, a)
            if j:
                _, d, acc, splice, end = j
                exons = exons[:-1] + [(exons[-1][0], d), (acc, end)]
                how = "joined to the screen exon"
        if not how:
            end = exons[-1][1]
            body = "".join(seq[s:e] for s, e in exons)
            k = end + (-len(body)) % 3     # complete the last codon
            while k + 3 <= len(seq) and seq[k:k + 3] not in STOPS: k += 3
            if k + 3 <= len(seq) and (k - end) // 3 <= a.max_extension:
                exons[-1] = (exons[-1][0], k + 3); how = f"extended {(k - end) // 3} codons through TGA"
            else: how = "not extended"
        # N-terminus
        s1 = exons[0][0]
        start = "complete" if seq[s1:s1 + 3] == "ATG" else "incomplete"
        if start == "incomplete":
            k = s1
            for _ in range(a.start_window):
                k -= 3
                if k < 0 or seq[k:k + 3] in STOPS + ("TGA",): break
                if seq[k:k + 3] == "ATG": s1, start = k, "complete (start upstream)"
            if start == "incomplete":   # an ATG among the first codons of the model (not an internal methionine)
                for k in range(exons[0][0], min(exons[0][1], exons[0][0] + 3 * a.start_downstream), 3):
                    if seq[k:k + 3] == "ATG": s1, start = k, "complete (start downstream)"; break
            exons[0] = (s1, exons[0][1])
        cds = "".join(seq[s:e] for s, e in exons)
        prot = tr_u(cds)
        if prot.endswith("*"): prot = prot[:-1]
        aln = loc.align(clean(prot), clean(ref))[0]
        ui, res, idx, same, tot = ref.find("U"), "outside", None, 0, 0
        for (c0, c1), (q0, q1) in zip(*aln.aligned):
            for k in range(c1 - c0): tot += 1; same += prot[c0 + k] == ref[q0 + k]
            if q0 <= ui < q1: idx = c0 + ui - q0; res = prot[idx]
        codon = ""
        if idx is not None:
            nt = 3 * idx
            for s, e in exons:
                if nt < e - s:
                    pos = s + nt
                    codon = f"{seq[pos:pos + 3]} at {pos}" + (" (screen codon)" if tga is not None and pos == tga else "")
                    break
                nt -= e - s
        to_contig = (lambda x: s0 + x) if strand == "+" else (lambda x: e0 - x)
        r.update({"status": "model", "query": q, "genewise_score": score, "exons": len(exons),
                  "coding_exons": ";".join(f"{s + 1}-{e}" for s, e in exons),
                  "coding_exons_contig": ";".join(f"{min(to_contig(s), to_contig(e - 1))}-{max(to_contig(s), to_contig(e - 1))}"
                                                  for s, e in exons),
                  "c_terminus": how, "splice_to_screen_exon": splice, "n_terminus": start, "length": len(prot),
                  "residue_at_query_U": res, "codon_at_query_U": codon,
                  "identity_to_query": f"{100 * same / tot:.0f}" if tot else "", "protein": prot})
        rows.append(r)
        if res == "U":
            faa.append((f"{a.family or q.split('_')[0]}|{genome}|{rname} {contig} ({strand}), genewise of {q}, "
                        f"{how}, N-terminus {start}", prot))
    os.makedirs(os.path.dirname(a.faa) or ".", exist_ok=True)
    with open(a.faa, "w") as o:
        for h, s in faa: o.write(f">{h}\n{s}\n")
    cols = ["region", "genome", "location", "strand", "status", "query", "genewise_score", "exons", "coding_exons",
            "coding_exons_contig", "c_terminus", "splice_to_screen_exon", "n_terminus", "length",
            "residue_at_query_U", "codon_at_query_U", "identity_to_query", "protein"]
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    for r in rows: print("\t".join(str(r.get(c, "")) for c in cols))
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
