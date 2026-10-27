#!/usr/bin/env python3
"""Genome-scan SECIS hits without a selenoprotein gene upstream: what are they?

Input: results/tables/stage2_secis_scan_classes.tsv (hits at the genome threshold, classed by
figures/secis_model_figures.py); by default the classes "nothing found upstream" and "stage 1 prediction without
selenocysteine". For each hit (genome coordinates, strand of the hit):
1. annotation, when the genome has one (stage1/samplesheet.csv): the annotated transcripts on the same strand that
   overlap the hit, with the part they cover (coding exon, untranslated exon, intron), and the nearest annotated
   transcript that ends upstream on the same strand within --upstream nt (distance from the end of its coding
   sequence and from its 3' end);
2. blastx (legacy BLAST, top strand, no filter, E <= --evalue) of the --upstream nt before the hit against the human
   RefSeq proteome and the gold-standard proteins: best human hit and best gold-standard hit, each with its distance
   to the hit and the stop codons inside the
   aligned part of the upstream region, each as codon>aligned residue (TGA>U, a selenocysteine of a gold-standard
   protein, or a TGA facing a residue of a human protein would point to a missed gene; stops facing gaps usually
   mark introns crossed by the ungapped extension);
3. copies in its own genome: blastn of the hit and --flank nt on each side, HSPs of at least 80% identity over at least
   80% of the query (a count well above 1 points to a repeated element);
4. conservation: the same query against the other Oikopleura assemblies (best HSP: identity and coverage of the
   query), and whether a genome-scan hit (results/tables/stage2_secis_scan_hits.tsv, all scores listed there) lies on
   the same strand within --near nt of the aligned position.
Output: results/tables/stage2_secis_orphan_hits.tsv. Genome databases are built in --work (one per genome, removed at
the end unless --keep).
"""
import argparse, collections, csv, gzip, os, re, shutil, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from uga_secis_screen import read_gff_transcripts, read_fasta, rc, opener   # noqa: E402

OIKOPLEURA = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan"]
CODE = {a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[16 * i + 4 * j + k]
        for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}


def run(cmd, cwd=None):
    p = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode: sys.exit(f"command failed ({p.returncode}) in {cwd or '.'}: {' '.join(cmd)}\n{p.stderr[-2000:]}")
    return p.stdout


def tsv(path): return list(csv.DictReader(open(path), delimiter="\t"))


def genome_db(g, path, work):
    """Legacy BLAST nucleotide database of a genome in work/<g>/ (built once); returns its path."""
    d = os.path.join(work, g)
    os.makedirs(d, exist_ok=True)
    db = os.path.join(d, "genome.fa")
    if not (os.path.exists(db + ".nsq") or os.path.exists(db + ".nal")):
        with open(db, "w") as o:
            for k, s in read_fasta(path): o.write(f">{k}\n{s.upper()}\n")
        run(["formatdb", "-i", "genome.fa", "-p", "F", "-o", "F"], d)
    return db


def blast_tab(cmd):
    """Rows of a -m 8 BLAST table."""
    out = []
    for l in run(cmd).splitlines():
        f = l.split("\t")
        if len(f) < 12: continue
        out.append({"q": f[0], "s": f[1], "ident": float(f[2]), "len": int(f[3]), "qs": int(f[6]), "qe": int(f[7]),
                    "ss": int(f[8]), "se": int(f[9]), "evalue": float(f[10]), "bits": float(f[11])})
    return out


def annotation(txs, contig, strand, lo, hi, upstream):
    """Transcripts on the strand overlapping [lo, hi) (0-based), and the nearest one ending upstream."""
    over, best = [], None
    for k, t in txs.items():
        if t["contig"] != contig or t["strand"] != strand: continue
        a, b = t["blocks"][0][0], t["blocks"][-1][1]
        if a < hi and lo < b:
            if any(x < hi and lo < y for x, y in t["cds"]): part = "coding exon"
            elif any(x < hi and lo < y for x, y in t["blocks"]): part = "untranslated exon"
            else: part = "intron"
            over.append(f"{k}{' (' + t['gene'] + ')' if t['gene'] else ''}: {part}")
            continue
        # ends upstream: on +, its 3' end b <= lo; on -, its 3' end a >= hi
        d3 = lo - b if strand == "+" else a - hi
        if 0 <= d3 <= upstream and (best is None or d3 < best[1]):
            cds_end = (max(y for _, y in t["cds"]) if strand == "+" else min(x for x, _ in t["cds"])) if t["cds"] else None
            dc = (lo - cds_end if strand == "+" else cds_end - hi) if cds_end is not None else None
            best = (k, d3, dc, t["gene"])
    up = "" if best is None else (f"{best[0]}{' (' + best[3] + ')' if best[3] else ''}: 3' end {best[1]} nt upstream" +
                                  (f", coding sequence end {best[2]} nt" if best[2] is not None else ""))
    return "; ".join(over), up


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--classes-table", default="results/tables/stage2_secis_scan_classes.tsv")
    ap.add_argument("--classes", nargs="+", default=["nothing found upstream", "stage 1 prediction without selenocysteine"])
    ap.add_argument("--scan-hits", default="results/tables/stage2_secis_scan_hits.tsv")
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--human-proteome", default="data/genomes/Hsap/protein.faa.gz")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--upstream", type=int, default=5000)
    ap.add_argument("--flank", type=int, default=50, help="nt on each side of the hit for the copy and conservation searches")
    ap.add_argument("--near", type=int, default=300)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--work", default="work/stage2/secis_orphans")
    ap.add_argument("--keep", action="store_true", help="keep the genome databases")
    ap.add_argument("--out", default="results/tables/stage2_secis_orphan_hits.tsv")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    hits = [r for r in tsv(a.classes_table) if r["class"] in a.classes]
    if not hits: sys.exit("no hit of the requested classes")
    scan = tsv(a.scan_hits)
    genomes = sorted({h["genome"] for h in hits}, key=lambda g: list(sheet).index(g))
    # sequences: upstream region and hit with flanks, in the orientation of the hit
    need = collections.defaultdict(set)
    for h in hits: need[h["genome"]].add(h["contig"])
    seqs = {}
    for g in genomes:
        for k, s in read_fasta(sheet[g]["fasta"]):
            if k in need[g]: seqs[(g, k)] = s.upper()
    q_up, q_hit = {}, {}
    for n, h in enumerate(hits):
        s = seqs[(h["genome"], h["contig"])]
        lo, hi = sorted((int(h["from"]), int(h["to"])))
        lo -= 1   # 0-based half-open
        h["lo"], h["hi"], h["id"] = lo, hi, f"hit{n + 1}"
        if h["strand"] == "+":
            q_up[h["id"]] = s[max(0, lo - a.upstream):lo]
            q_hit[h["id"]] = s[max(0, lo - a.flank):hi + a.flank]
        else:
            q_up[h["id"]] = rc(s[hi:hi + a.upstream])
            q_hit[h["id"]] = rc(s[max(0, lo - a.flank):hi + a.flank])
    fa_up, fa_hit = os.path.join(a.work, "upstream.fa"), os.path.join(a.work, "hits.fa")
    with open(fa_up, "w") as o:
        for k, s in q_up.items():
            if s: o.write(f">{k}\n{s}\n")
    with open(fa_hit, "w") as o:
        for k, s in q_hit.items(): o.write(f">{k}\n{s}\n")
    # 1. annotation
    ann = {}
    for g in genomes:
        path = sheet[g].get("annotation", "")
        if not path or not os.path.exists(path): continue
        txs = read_gff_transcripts(path)
        for h in hits:
            if h["genome"] == g: ann[h["id"]] = annotation(txs, h["contig"], h["strand"], h["lo"], h["hi"], a.upstream)
    # 2. blastx against human and gold-standard proteins
    pdir = os.path.join(a.work, "proteins")
    os.makedirs(pdir, exist_ok=True)
    pdb = os.path.join(pdir, "proteins.faa")
    desc = {}
    if not os.path.exists(pdb + ".psq"):
        with open(pdb, "w") as o:
            with opener(a.human_proteome) as fh:
                name = None
                for l in fh:
                    if l.startswith(">"):
                        name, _, d = l[1:].strip().partition(" ")
                        o.write(f">{name}\n")
                    else: o.write(l.replace("*", ""))
            for f in sorted(os.listdir(a.gold)):
                if f.endswith(".proteins.faa"):
                    for k, s in read_fasta(os.path.join(a.gold, f)): o.write(f">gold|{f.split('.')[0]}|{k}\n{s}\n")
        run(["formatdb", "-i", "proteins.faa", "-p", "T", "-o", "F"], pdir)
    with opener(a.human_proteome) as fh:
        for l in fh:
            if l.startswith(">"):
                name, _, d = l[1:].strip().partition(" ")
                desc[name] = re.sub(r" \[[^]]+\]$", "", d)
    bx = collections.defaultdict(dict)   # query -> {"human": best human hit, "gold": best gold-standard hit}
    xml = run(["blastall", "-p", "blastx", "-d", pdb, "-i", fa_up, "-e", str(a.evalue), "-F", "F", "-S", "1", "-m", "7",
               "-v", "100", "-b", "100", "-a", str(a.cpus)])
    import xml.etree.ElementTree as ET
    names = list(q_up)
    for it in ET.fromstring(xml).iter("Iteration"):
        num = int(it.findtext("Iteration_iter-num") or 0)
        qdef = it.findtext("Iteration_query-def")
        q = qdef.split()[0] if qdef else names[num - 1]
        for hit in it.iter("Hit"):   # in decreasing order of score
            sid = hit.findtext("Hit_def").split()[0]
            kind = "gold" if sid.startswith("gold|") else "human"
            if kind in bx[q]: continue
            hsp = hit.find("Hit_hsps/Hsp")
            qf, qt = int(hsp.findtext("Hsp_query-from")), int(hsp.findtext("Hsp_query-to"))
            qseq, hseq = hsp.findtext("Hsp_qseq"), hsp.findtext("Hsp_hseq")
            # stop codons inside the HSP: the codon read in the upstream region, and the residue aligned to it
            stops, pos = [], min(qf, qt) - 1
            for c, d in zip(qseq, hseq):
                if c == "-": continue
                if c == "*": stops.append(f"{q_up[q][pos:pos + 3]}>{d}")
                pos += 3
            bx[q][kind] = {"subject": sid if kind == "human" else sid.split("|")[2] + " (" + sid.split("|")[1] + ")",
                           "description": desc.get(sid, ""), "evalue": float(hsp.findtext("Hsp_evalue")),
                           "identity": f"{100 * int(hsp.findtext('Hsp_identity')) / int(hsp.findtext('Hsp_align-len')):.0f}",
                           "to_hit": len(q_up[q]) - max(qf, qt), "stops": len(stops), "aligned_to_stops": ",".join(stops)}
    # 3-4. copies in the own genome, conservation in the other Oikopleura assemblies
    targets = [g for g in sheet if g in set(genomes) | set(OIKOPLEURA)]
    missing = [g for g in targets if not os.path.exists(sheet[g]["fasta"])]
    if missing: print(f"genomes not found, left out of the copy and conservation searches: {' '.join(missing)}", file=sys.stderr)
    dbs = {g: genome_db(g, sheet[g]["fasta"], a.work) for g in targets if g not in missing}
    qlen = {k: len(s) for k, s in q_hit.items()}
    copies, cons = collections.Counter(), collections.defaultdict(dict)
    scan_by = collections.defaultdict(list)
    for r in scan: scan_by[(r["genome"], r["contig"], r["strand"])].append((min(int(r["from"]), int(r["to"])), float(r["bits"])))
    for g, db in dbs.items():
        rows = blast_tab(["blastall", "-p", "blastn", "-d", db, "-i", fa_hit, "-e", "1e-5", "-F", "F", "-m", "8", "-a", str(a.cpus)])
        best = {}
        for r in rows:
            h = next(x for x in hits if x["id"] == r["q"])
            cov = (abs(r["qe"] - r["qs"]) + 1) / qlen[r["q"]]
            if g == h["genome"]:
                if r["ident"] >= 80 and cov >= 0.8: copies[r["q"]] += 1
            elif r["q"] not in best or r["bits"] > best[r["q"]]["bits"]:
                best[r["q"]] = dict(r, cov=cov)
        for q, r in best.items():
            h = next(x for x in hits if x["id"] == q)
            strand = "+" if r["se"] > r["ss"] else "-"   # the query is in the orientation of the hit
            pos = min(r["ss"], r["se"])
            near = [b for p, b in scan_by[(g, r["s"], strand)] if abs(p - pos) <= a.near + qlen[q]]
            cons[q][g] = f"{r['ident']:.0f}%/{100 * r['cov']:.0f}%" + (f" SECIS {max(near):.1f}" if near else "")
    # output
    cols = ["genome", "contig", "strand", "from", "to", "bits", "class", "annotation_at_hit", "annotated_upstream"]
    for kind in ("human", "gold"):
        cols += [f"blastx_{kind}", f"blastx_{kind}_evalue", f"blastx_{kind}_identity", f"blastx_{kind}_end_to_hit",
                 f"blastx_{kind}_stops", f"blastx_{kind}_aligned_to_stops"]
    cols += ["copies_in_genome", "conservation (identity/coverage; scan hit nearby)"]
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for h in hits:
            an = ann.get(h["id"], ("", "")) if sheet[h["genome"]].get("annotation") else ("no annotation", "no annotation")
            row = [h["genome"], h["contig"], h["strand"], h["from"], h["to"], h["bits"], h["class"], an[0], an[1]]
            for kind in ("human", "gold"):
                b = bx[h["id"]].get(kind, {})
                row += [f"{b['subject']} {b['description']}".strip() if b else "", f"{b['evalue']:.1e}" if b else "",
                        b.get("identity", ""), b.get("to_hit", ""), b.get("stops", ""), b.get("aligned_to_stops", "")]
            row += [copies.get(h["id"], 0), "; ".join(f"{g} {v}" for g, v in cons[h["id"]].items())]
            o.write("\t".join(str(x) for x in row) + "\n")
    print(f"{len(hits)} hits -> {a.out}")
    if not a.keep:
        for g in dbs: shutil.rmtree(os.path.join(a.work, g), ignore_errors=True)


if __name__ == "__main__":
    main()
