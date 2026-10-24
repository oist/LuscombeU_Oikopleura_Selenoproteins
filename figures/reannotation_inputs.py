#!/usr/bin/env python3
"""Inputs of the reannotation figures (figures/reannotation_pilot.py): the public gene models that overlap the
gold-standard model of a gene in each O. dioica assembly, and reference proteins with their Pfam domains.

Public models: the annotated transcripts of the candidate table of step 2 (curation/gold_standard/<genome>
.candidate_models.tsv) whose exons overlap the coding sequence of the gold-standard model, on its strand. Their coding sequence is taken, in this
order, from the candidate table (annotated CDS rows), from the EMBL record of the
contig when given (--embl, GenBank models of OdB3: CDS and mRNA locations, translation), or by finding the public
protein (EBI dbfetch, RefSeq proteins) in the translation of the exon chain read on the genome (--genome; TGA read
as U). A public protein with less than 30% identity to the gold-standard protein is marked as read in another
frame. Models of CDS rows only (OSKA2016) use their exons as coding sequence.
Splice sites of the drawn introns of the public and gold-standard models: read on the genome when it is given,
otherwise taken from the rows of the candidate table that have the same intron ("?" when none has it).
Reference proteins (--reference label=UniProt accession, the first one giving the coordinates of panel B): sequence,
selenocysteine positions and Pfam domains (short names) from the InterPro API.
Needs network access to www.ebi.ac.uk. Outputs figures/inputs/reannotation_<gene>.tsv (with the source of each
value) and .fa (public and reference proteins).
"""
import argparse, csv, json, os, re, sys, urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "curation"))
from gold_standard import read_fasta, site, translate  # noqa: E402

DB = {"OKI": "RefSeq", "OSKA": "ANISEED", "OdB3": "GenBank"}


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as r: return r.read().decode()


def blocks(s): return [tuple(map(int, x.split("-"))) for x in s.split(";")] if s else []
def fmt(b): return ";".join(f"{a}-{c}" for a, c in sorted(b))


def gold(path, gene):
    """Strand, exons and CDS rows (1-based inclusive) of a gold-standard model."""
    exons, cds, strand = [], [], "+"
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("exon", "CDS"): continue
        if re.search(r"Parent=rna-([^;]+)", f[8]).group(1).split("_", 1)[1] == gene:
            (exons if f[2] == "exon" else cds).append((int(f[3]), int(f[4]))); strand = f[6]
    return strand, exons, cds


def embl_models(text):
    """protein_id -> {"CDS": location, "mRNA": location of the mRNA before it, "translation"} from an EMBL file."""
    feats = []
    for l in text.splitlines():
        if not l.startswith("FT"): continue
        key, body = l[5:21].strip(), l[21:].rstrip()
        if key: feats.append({"key": key, "loc": body, "q": []})
        elif not feats: continue
        elif body.startswith("/"): feats[-1]["q"].append(body)
        elif feats[-1]["q"] and feats[-1]["q"][-1].count('"') % 2 == 1: feats[-1]["q"][-1] += body
        else: feats[-1]["loc"] += body
    out, last_mrna = {}, None
    for f in feats:
        q = {}
        for x in f["q"]:
            k, _, v = x[1:].partition("=")
            q[k] = v.strip('"')
        if f["key"] == "mRNA": last_mrna = f["loc"]
        if f["key"] == "CDS" and "protein_id" in q:
            out[q["protein_id"]] = {"CDS": f["loc"], "mRNA": last_mrna, "translation": q.get("translation", "")}
    return out


def location(s):
    return [(int(a), int(b)) for a, b in re.findall(r"(\d+)\.\.(\d+)", s.replace("<", "").replace(">", ""))]


def locate_protein(protein, exons, genome, strand):
    """CDS (1-based segments) of a protein found in the translation of an exon chain, TGA read as U; None if absent."""
    pos = [p for a, b in sorted(exons) for p in range(a, b + 1)]
    if strand == "-": pos = pos[::-1]
    nt = "".join(genome[p - 1] for p in pos).upper()
    if strand == "-": nt = nt.translate(str.maketrans("ACGT", "TGCA"))
    for f in (0, 1, 2):
        aa = "".join("U" if nt[i:i + 3] == "TGA" else translate(nt[i:i + 3]) for i in range(f, len(nt) - 2, 3))
        k = aa.replace("U", "*").find(protein.replace("U", "*"))
        if k < 0: continue
        start = f + 3 * k
        end = start + 3 * len(protein) + (3 if aa[k + len(protein):k + len(protein) + 1] == "*" else 0)
        keep = sorted(pos[start:end])
        segs, s0 = [], keep[0]
        for x, y in zip(keep, keep[1:] + [None]):
            if y != x + 1: segs.append((s0, x)); s0 = y
        return segs
    return None


def intron_sites(exons, genome, strand, table):
    """Splice site of each intron of an exon chain (1-based inclusive exons), with its source."""
    introns = [(b, a2 - 1) for (_, b), (a2, _) in zip(sorted(exons), sorted(exons)[1:]) if a2 > b + 1]
    if genome: return ";".join(site(genome, a, b, strand) for a, b in introns), "genome sequence"
    return ";".join(table.get(i, "?") for i in introns), "candidate table of step 2 (same intron)"


def table_sites(rows):
    """Intron (0-based, end exclusive) -> splice site, from the rows of a candidate table."""
    out = {}
    for r in rows:
        if not r["exons"] or not r["splice_sites"] or r["splice_sites"].startswith(("CDS", "no CDS")): continue
        ex = blocks(r["exons"])
        introns = [(b, a2 - 1) for (_, b), (a2, _) in zip(ex, ex[1:]) if a2 > b + 1]
        sites = r["splice_sites"].split(";")
        if len(sites) == len(introns): out.update(zip(introns, sites))
    return out


def identity(a, b):
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                               open_gap_score=-10, extend_gap_score=-0.5)
    x, y = a.replace("U", "C").replace("*", ""), b.replace("U", "C").replace("*", "")
    bq, br = al.align(x, y)[0].aligned
    same = sum(x[i] == y[j] for (q0, q1), (r0, r1) in zip(bq, br) for i, j in zip(range(q0, q1), range(r0, r1)))
    return same / min(len(x), len(y))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--genome", nargs="*", default=[], help="genome=FASTA of the contigs (to locate RefSeq proteins)")
    ap.add_argument("--embl", nargs="*", default=[], help="EMBL accessions of OdB3 contigs (GenBank models)")
    ap.add_argument("--reference", nargs="+", required=True, help="label=UniProt accession; the first gives the axis")
    ap.add_argument("--reference-fasta", nargs="*", default=[],
                    help="label=FASTA:sequence name:accession, for a reference protein without UniProt entry (gaps removed)")
    ap.add_argument("--note", nargs="*", default=[], help="genome:model=text shown next to a public model")
    ap.add_argument("--reference-note", nargs="*", default=[], help="label=text shown under a reference protein")
    ap.add_argument("--gold-note", nargs="*", default=[], help="genome=text shown next to a gold-standard model")
    ap.add_argument("--out", default="figures/inputs")
    a = ap.parse_args()
    gene = a.gene
    genomes = {k: v for k, v in (x.split("=", 1) for x in a.genome)}
    embl = {}
    for acc in a.embl:
        embl.update(embl_models(fetch(f"https://www.ebi.ac.uk/ena/browser/api/embl/{acc}?download=false")))
    rows, fasta = [("kind", "assembly", "name", "field", "value", "source")], []
    gprot = {}
    for g in ("OKI", "OSKA", "OdB3"):
        name = None
        for l in open(f"{a.gold}/{g}.proteins.faa"):
            if l.startswith(">"): name = l[1:].split("|")[0]
            elif name == gene: gprot[g] = l.strip()
    for g in ("OKI", "OSKA", "OdB3"):
        strand, gex, gcds = gold(f"{a.gold}/{g}.gff3", gene)
        lo, hi = min(x for x, _ in gcds), max(y for _, y in gcds)   # public models must overlap the coding sequence
        allrows = [r for r in csv.DictReader(open(f"{a.gold}/{g}.candidate_models.tsv"), delimiter="\t") if r["gene"] == gene]
        cand = [r for r in allrows if r["kind"] == "annotation" and r["item"] != "transcript frame"]
        tsites = table_sites(allrows)
        seen = set()
        contig = next(r2["contig"] for r2 in csv.DictReader(open("results/tables/gold_standard.tsv"), delimiter="\t")
                      if r2["genome"] == g and r2["gene"] == gene)
        genome = dict(read_fasta(genomes[g], {contig})) if g in genomes else {}
        cseq = next((v for k, v in genome.items() if k.split()[0] == contig), "")
        for r in cand:
            pid = r["source"]
            ex = blocks(r["exons"])
            if pid in seen or not any(x <= hi and y >= lo for x, y in ex): continue
            seen.add(pid)
            src = f"candidate table of step 2 ({g}.candidate_models.tsv)"
            cds, prot, csrc = None, "", ""
            if r.get("annotated_cds") and r["annotated_cds"] != "no CDS":
                cds, csrc = blocks(r["annotated_cds"]), src
            elif r["item"].startswith("annotated") and r["splice_sites"].startswith("CDS"):   # tables before annotated_cds
                cds, csrc = blocks(r["splice_sites"][4:]), src
            if pid in embl:
                cds, csrc = location(embl[pid]["CDS"]), "ENA EMBL record (CDS location)"
                if embl[pid]["mRNA"]: ex = location(embl[pid]["mRNA"])
                prot = embl[pid]["translation"]
            elif pid.startswith("XP_"):
                prot = "".join(l.strip() for l in fetch(f"https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=refseqp&id={pid}&format=fasta&style=raw").splitlines()[1:])
                if cds is None and cseq:
                    cds = locate_protein(prot, ex, cseq, strand)
                    csrc = "RefSeq protein (EBI dbfetch) located in the translation of the RefSeq exon chain"
            if cds is None and g == "OSKA": cds, csrc = ex, "CDS rows of the ANISEED model (no UTR in this annotation)"
            rows.append(("public", g, pid, "database", DB[g], ""))
            rows.append(("public", g, pid, "exons", fmt(ex), src if pid not in embl else "ENA EMBL record (mRNA location)"))
            if cds: rows.append(("public", g, pid, "cds", fmt(cds), csrc))
            if len(ex) > 1: rows.append(("public", g, pid, "sites", *intron_sites(ex, cseq, strand, tsites)))
            if prot:
                fasta.append(f">{g}_{pid}\n{prot}\n")
                if identity(prot, gprot[g]) < 0.3: rows.append(("public", g, pid, "frame", "other", "identity to the gold-standard protein below 30%"))
        for n in a.note:
            k, _, text = n.partition("=")
            ng, _, nm = k.partition(":")
            if ng == g and nm in seen: rows.append(("public", g, nm, "note", text, "curator"))
        rows.append(("gold", g, gene, "sites", *intron_sites(gex, cseq, strand, tsites)))
        for n in a.gold_note:
            ng, _, text = n.partition("=")
            if ng == g: rows.append(("gold", g, gene, "note", text, "curator"))
    for i, ref in enumerate(a.reference):
        label, acc = ref.split("=", 1)
        md = json.loads(fetch(f"https://www.ebi.ac.uk/interpro/api/protein/uniprot/{acc}/"))["metadata"]
        dom = json.loads(fetch(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{acc}/"))
        s = md["sequence"]
        fasta.append(f">{label}\n{s}\n")
        src = f"UniProt {acc} (InterPro API)"
        rows.append(("reference", "-", label, "accession", f"UniProt {acc}", src))
        rows.append(("reference", "-", label, "organism", md["source_organism"]["scientificName"], src))
        if i == 0: rows.append(("reference", "-", label, "axis", "yes", ""))
        rows.append(("reference", "-", label, "selenocysteine", ",".join(str(k + 1) for k, c in enumerate(s) if c == "U"), src))
        for r in dom["results"]:
            f = r["proteins"][0]["entry_protein_locations"][0]["fragments"][0]
            pf = r["metadata"]["accession"]
            short = json.loads(fetch(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pf}/"))["metadata"]["name"].get("short", pf)
            rows.append(("reference", "-", label, "domain", f"{pf}:{f['start']}-{f['end']}:{short}", src + ", Pfam short name"))
    for ref in a.reference_fasta:
        label, spec = ref.split("=", 1)
        path, name, acc = spec.split(":", 2)
        seq = next(v for k, v in read_fasta(path) if name in k).replace("-", "").replace(".", "").upper()
        fasta.append(f">{label}\n{seq}\n")
        rows.append(("reference", "-", label, "accession", acc, path))
    for n in a.reference_note:
        label, _, text = n.partition("=")
        rows.append(("reference", "-", label, "note", text, "curator"))
    os.makedirs(a.out, exist_ok=True)
    with open(f"{a.out}/reannotation_{gene}.tsv", "w") as o:
        for r in rows: o.write("\t".join(r) + "\n")
    with open(f"{a.out}/reannotation_{gene}.fa", "w") as o: o.write("".join(fasta))
    print(open(f"{a.out}/reannotation_{gene}.tsv").read())


if __name__ == "__main__":
    main()
