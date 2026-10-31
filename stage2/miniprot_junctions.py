#!/usr/bin/env python3
"""Why miniprot does not call the selenocysteine codon of SELENOT and SEPHS2 (manuscript: "their selenocysteine codons
lie 2 and 1 nucleotides from an intron, where miniprot may move the junction to avoid an in-frame stop codon (not
tested)").

Test: the gold-standard proteins of the three O. dioica assemblies (curation/gold_standard/<genome>.proteins.faa) are
aligned by miniprot to each assembly as in the pipeline (miniprot -I --gff on the indexed genome), with the
selenocysteine as U (as in the pipeline's profile queries) and, as a control, replaced by C. For each gold-standard
gene and query, the best alignment (alignment score) overlapping the gene on its strand is compared with the
gold-standard model:
  - sec_support: the pipeline's rule (bin/sec_check.py of the pipeline, copied here): the query U aligned to a TGA
    ("aligned"), or U the first residue after the end of the alignment with TGA as next codon ("next_codon"), else
    "none"; for the C queries, a TGA aligned to that C;
  - codon_in_cds: whether the three bases of the selenocysteine codon lie in one CDS feature of the alignment, and its
    frame there (0 = in frame);
  - the gold-standard intron nearest to the codon (distance in nt from the codon), the alignment's intron that
    overlaps it, the shift of its donor and acceptor (nt, transcript orientation, alignment minus gold standard) and
    the splice sites of both;
  - the genome codon the alignment puts against the query's selenocysteine (sec_aligned_to);
  - frameshifts and stop codons reported by miniprot (fs, st tags; for miniprot 0.11, which does not write them,
    counted from the cg and cs tags).
Default genes: all gold-standard genes, so that SELENOT and SEPHS2 can be compared with the genes miniprot finds.
Output: --table; miniprot outputs in --work.

Other uses (both read existing miniprot GFF files instead of aligning):
  - --from-gff DIR: the alignments in DIR/<genome>.gff, for example made with another miniprot version from
    <work>/queries.faa (--queries-only writes it and stops; scripts/miniprot_v011.sh does this with miniprot 0.11,
    the version of the pipeline);
  - --from-gff DIR --queries FASTA: the pipeline's own alignments (DIR = <outdir>/miniprot of a pipeline run, FASTA =
    its <outdir>/summary/selenoprotein_queries.faa): one row per query whose best alignment overlaps a gold-standard
    gene on its strand, residue U.
"""
import argparse, collections, csv, os, re, subprocess, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "curation"))
from gold_standard import attributes, read_fasta, rc   # noqa: E402

DIOICA = ["OKI", "OSKA", "OdB3"]
CS_TOKEN = re.compile(r"(:\d+)|(\*[acgtn]+[A-Z*])|(\+[A-Z*]+)|(-[acgtn]+)|(~[acgtn]{2}\d+[acgtn]{2})")


def gold_models(path):
    """name -> {contig, strand, cds: sorted [(s0, e0 inclusive)], span}."""
    out = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        name = attributes(f[8])["Parent"].split("-", 1)[1].split("_", 1)[1]
        out.setdefault(name, {"contig": f[0], "strand": f[6], "cds": []})["cds"].append((int(f[3]) - 1, int(f[4]) - 1))
    for m in out.values():
        m["cds"].sort(); m["span"] = (m["cds"][0][0], m["cds"][-1][1])
    return out


def read_miniprot_gff(path):
    """Alignments: list of {paf fields, cs, AS, fs, st, cds: sorted [(s0, e0 inclusive)]}."""
    alns, cur = [], None
    for l in open(path):
        if l.startswith("##PAF\t"):
            f = l.rstrip("\n").split("\t")[1:]
            tags = {x.split(":", 2)[0]: x.split(":", 2)[2] for x in f[12:] if x.count(":") >= 2}
            cs = tags.get("cs", "")
            # miniprot 0.11 writes no fs and st tags: count frameshifts from the cg tag (F and G operations) and
            # in-frame stop codons from the cs tag (stop codons aligned to a residue)
            fs = int(tags["fs"]) if "fs" in tags else len(re.findall(r"\d+[FG]", tags.get("cg", "")))
            st = int(tags["st"]) if "st" in tags else len(re.findall(r"\*(?:taa|tag|tga)[A-Z]", cs))
            cur = {"query": f[0], "qlen": int(f[1]), "qstart": int(f[2]), "qend": int(f[3]), "strand": f[4],
                   "contig": f[5], "tstart": int(f[7]), "tend": int(f[8]), "AS": int(tags.get("AS", 0)),
                   "fs": fs, "st": st, "cs": cs, "cds": []}
            alns.append(cur)
        elif not l.startswith("#") and cur is not None:
            f = l.rstrip("\n").split("\t")
            if len(f) >= 9 and f[2] == "CDS": cur["cds"].append((int(f[3]) - 1, int(f[4]) - 1))
    for a in alns: a["cds"].sort()
    return alns


def aligned_codons(rec, query_seq, residue, uidx=None):
    """Genomic intervals (0-based half-open) of TGA codons aligned to `residue` at the query's selenocysteine position
    (pipeline bin/sec_check.py: aligned_sec_codons), the query position after the alignment and, if uidx is given,
    what the alignment puts against query position uidx (a codon, 'gap' for an insertion, '' if outside)."""
    sites, qpos, at = [], rec["qstart"], ""
    fwd = rec["strand"] == "+"
    gpos = rec["tstart"] if fwd else rec["tend"]

    def consume(n):
        nonlocal gpos
        iv = (gpos, gpos + n) if fwd else (gpos - n, gpos)
        gpos = gpos + n if fwd else gpos - n
        return iv
    for m in CS_TOKEN.finditer(rec["cs"]):
        t = m.group(0)
        if t[0] == ":":
            n = int(t[1:])
            if uidx is not None and qpos <= uidx < qpos + n: at = "identical"
            consume(3 * n); qpos += n
        elif t[0] == "*":
            codon = t[1:-1]; iv = consume(len(codon))
            if codon == "tga" and qpos < len(query_seq) and query_seq[qpos] == residue: sites.append(iv)
            if qpos == uidx: at = codon
            qpos += 1
        elif t[0] == "+":
            if uidx is not None and qpos <= uidx < qpos + len(t) - 1: at = "gap"
            qpos += len(t) - 1
        elif t[0] == "-":
            consume(len(t) - 1)
        else:
            consume(int(t[3:-2]))
    return (sites, qpos, at) if uidx is not None else (sites, qpos)


def introns(cds):
    """Introns between consecutive CDS features: [(first intron base, last intron base)] (0-based, inclusive)."""
    return [(cds[i][1] + 1, cds[i + 1][0] - 1) for i in range(len(cds) - 1)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--gold-table", default="results/tables/gold_standard.tsv")
    ap.add_argument("--genomes", nargs="+", default=DIOICA)
    ap.add_argument("--genes", nargs="*", default=None, help="default: all gold-standard genes")
    ap.add_argument("--miniprot", default="miniprot")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--work", default="work/stage2/miniprot_junctions")
    ap.add_argument("--table", default="results/tables/stage2_miniprot_junctions.tsv")
    ap.add_argument("--queries-only", action="store_true", help="write <work>/queries.faa and stop")
    ap.add_argument("--from-gff", help="read DIR/<genome>.gff instead of running miniprot")
    ap.add_argument("--queries", help="with --from-gff: the query proteins of those alignments (pipeline queries)")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    sec = {(r["genome"], r["gene"]): int(r["sec_codon"]) - 1 for r in csv.DictReader(open(a.gold_table), delimiter="\t")
           if r.get("sec_codon", "").isdigit()}
    # queries: gold proteins of all assemblies, U as is and U -> C
    prots = {}
    for g in a.genomes:
        for name, s in read_fasta(f"{a.gold}/{g}.proteins.faa"):
            gene = name.split("|")[0]
            if a.genes and gene not in a.genes: continue
            if "U" not in s: continue
            prots[(gene, g)] = s.rstrip("*")
    # queries: name -> {seq, gene (None: any gold-standard gene), label, residue, uidx}
    queries = {}
    if a.queries:
        for name, s in read_fasta(a.queries):
            s = s.upper().rstrip("*")
            if "U" in s:
                queries[name.split()[0]] = {"seq": s, "gene": None, "label": name.split()[0], "residue": "U",
                                            "uidx": s.index("U")}
    else:
        text = []
        for (gene, src), s in sorted(prots.items()):
            for var, seq in (("U", s), ("C", s.replace("U", "C"))):
                q = f"{gene}|{src}|{var}"; text.append(f">{q}\n{seq}\n")
                queries[q] = {"seq": seq, "gene": gene, "label": f"{src} protein, U as {var}", "residue": var,
                              "uidx": s.index("U")}
        qfa = os.path.join(a.work, "queries.faa")
        text = "".join(text)
        if not os.path.exists(qfa) or open(qfa).read() != text:   # unchanged file: alignments made from it stay valid
            with open(qfa, "w") as o: o.write(text)
        if a.queries_only:
            print(f"wrote {qfa}: {len(queries)} queries", flush=True); return
    rows = []
    for g in a.genomes:
        golds = gold_models(f"{a.gold}/{g}.gff3")
        if a.from_gff:
            gff = os.path.join(a.from_gff, f"{g}.gff")
        else:
            idx, gff = os.path.join(a.work, f"{g}.mpi"), os.path.join(a.work, f"{g}.gff")
            if not os.path.exists(idx):
                subprocess.run([a.miniprot, "-t", str(a.cpus), "-d", idx, sheet[g]["fasta"]], check=True)
            if not os.path.exists(gff) or os.path.getmtime(gff) < os.path.getmtime(qfa):
                with open(gff + ".part", "w") as o:
                    subprocess.run([a.miniprot, "-t", str(a.cpus), "-I", "--gff", idx, qfa], check=True, stdout=o)
                os.replace(gff + ".part", gff)
        alns = read_miniprot_gff(gff)
        contigs = {k: v.upper() for k, v in read_fasta(sheet[g]["fasta"], {m["contig"] for m in golds.values()})}
        by_q = collections.defaultdict(list)
        for al in alns: by_q[al["query"]].append(al)
        for gene, q in sorted((gn, q) for gn in golds for q in queries if queries[q]["gene"] in (None, gn)):
            if a.genes and gene not in a.genes: continue
            Q = queries[q]
            m = golds.get(gene)
            if m is None or (g, gene) not in sec: continue
            plus = m["strand"] == "+"
            c0 = sec[(g, gene)]                       # first base (lowest coordinate) of the TGA
            qseq, var, uidx = Q["seq"], Q["residue"], Q["uidx"]
            hits = [al for al in by_q.get(q, []) if al["contig"] == m["contig"] and al["strand"] == m["strand"]
                    and al["tstart"] <= m["span"][1] and m["span"][0] <= al["tend"] - 1]
            if Q["gene"] is None and not hits: continue
            row = {"genome": g, "gene": gene, "query": Q["label"], "alignment": "none"}
            gi = introns(m["cds"])
            near = min(gi, key=lambda iv: min(abs(iv[0] - c0), abs(iv[1] - c0))) if gi else None
            if near:
                dist = (near[0] - (c0 + 3)) if near[0] >= c0 + 3 else (c0 - near[1] - 1) if near[1] < c0 else 0
                row["gold_intron_distance"] = dist
                row["gold_intron_sites"] = (contigs[m["contig"]][near[0]:near[0] + 2] + "-" +
                                            contigs[m["contig"]][near[1] - 1:near[1] + 1]) if plus else \
                    (rc(contigs[m["contig"]][near[1] - 1:near[1] + 1]) + "-" + rc(contigs[m["contig"]][near[0]:near[0] + 2]))
            if hits:
                al = max(hits, key=lambda x: x["AS"])
                sites, qpos, at = aligned_codons(al, qseq, var, uidx)
                support = "aligned" if any(s == c0 for s, _ in sites) else "none"
                if support == "none" and var == "U" and qpos == al["qend"] and qpos < len(qseq) and qseq[qpos] == "U":
                    ctg = contigs[al["contig"]]
                    nxt = ctg[al["tend"]:al["tend"] + 3] if plus else rc(ctg[al["tstart"] - 3:al["tstart"]])
                    nstart = al["tend"] if plus else al["tstart"] - 3
                    if nxt == "TGA" and nstart == c0: support = "next_codon"
                inside = [cd for cd in al["cds"] if cd[0] <= c0 and c0 + 2 <= cd[1]]
                frame = ""
                if inside:
                    rows_ = al["cds"] if plus else al["cds"][::-1]
                    off = 0
                    for x, y in rows_:
                        if (x, y) == inside[0]:
                            off += (c0 - x) if plus else (y - (c0 + 2)); break
                        off += y - x + 1
                    frame = off % 3
                mi = introns(al["cds"])
                over = [iv for iv in mi if near and iv[0] <= near[1] and near[0] <= iv[1]]
                row.update({"alignment": f"{al['contig']}:{al['tstart'] + 1}-{al['tend']}", "score": al["AS"],
                            "query_covered": f"{al['qstart'] + 1}-{al['qend']} of {al['qlen']}",
                            "sec_position": uidx + 1, "sec_aligned_to": at, "sec_support": support,
                            "codon_in_cds": ("yes, frame " + str(frame)) if inside else "no",
                            "frameshifts": al["fs"], "stop_codons": al["st"]})
                if near and over:
                    iv = over[0]
                    d_shift = (iv[0] - near[0]) if plus else (near[1] - iv[1])
                    a_shift = (iv[1] - near[1]) if plus else (near[0] - iv[0])
                    ctg = contigs[m["contig"]]
                    sites_ = (ctg[iv[0]:iv[0] + 2] + "-" + ctg[iv[1] - 1:iv[1] + 1]) if plus else \
                        (rc(ctg[iv[1] - 1:iv[1] + 1]) + "-" + rc(ctg[iv[0]:iv[0] + 2]))
                    row.update({"intron_donor_shift": d_shift, "intron_acceptor_shift": a_shift,
                                "alignment_intron_sites": sites_})
                elif near:
                    row["alignment_intron_sites"] = "no intron there"
            rows.append(row)
    cols = ["genome", "gene", "query", "alignment", "score", "query_covered", "sec_position", "sec_aligned_to",
            "sec_support",
            "codon_in_cds", "gold_intron_distance", "gold_intron_sites", "intron_donor_shift", "intron_acceptor_shift",
            "alignment_intron_sites", "frameshifts", "stop_codons"]
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"wrote {a.table}: {len(rows)} rows", flush=True)
    # summary: per assembly and gene, alignments with U at the selenocysteine position, with the codon called, and
    # with the gold-standard intron nearest to the codon (no shift), without it, or with frameshifts
    for gene in sorted({r["gene"] for r in rows}):
        for g in a.genomes:
            rs = [r for r in rows if r["gene"] == gene and r["genome"] == g and r["alignment"] != "none"
                  and (a.queries or r["query"].endswith("U as U"))]
            if not rs: continue
            called = sum(r.get("sec_support") in ("aligned", "next_codon") for r in rs)
            exact = sum(r.get("intron_donor_shift") == 0 and r.get("intron_acceptor_shift") == 0 for r in rs)
            nointr = sum(r.get("alignment_intron_sites") == "no intron there" for r in rs)
            fs = sum(int(r.get("frameshifts") or 0) > 0 for r in rs)
            at = collections.Counter(r.get("sec_aligned_to") or "outside" for r in rs)
            print(f"{gene} {g}: {len(rs)} U alignments, codon called {called}, nearest intron exact {exact}, "
                  f"absent {nointr}, with frameshifts {fs}; U against " +
                  ", ".join(f"{k} {v}" for k, v in at.most_common()), flush=True)


if __name__ == "__main__":
    main()
