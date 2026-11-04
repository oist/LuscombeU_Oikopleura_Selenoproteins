#!/usr/bin/env python3
"""PSTK trees with more eukaryotic PSTK and with archaeal PSTK: a test of the position of the appendicularian clade.

The PSTK tree (stage2/machinery_phylogeny.py, results/trees/pstk.faa) samples few eukaryotes outside the animals with
genomes in this study, and none of the insects and protists whose PSTK are the closest to the appendicularian
proteins in a BLASTP search of ClusteredNR (curation/pstk_blast). This script adds them and builds the trees aside
from the main one.
Selection, from the downloads of the searches limited to Eukaryota and to Opisthokonta and of a PSI-BLAST search
(--eukaryotes), and of the search limited to Archaea (--archaea), each sequence once (accession), O. dioica and
bacteria left out (archaeal sequences of the PSI-BLAST search are treated as archaeal):
  1. family: for eukaryotes, the best BLASTP hit (legacy BLAST, bit score) among the sequences of the PSTK tree
     outside the appendicularians: a PSTK of at least --min-bits (kept) or a KTI12 of the outgroup (left out); for archaea, whose PSTK are
     about as close to the KTI12 as to the eukaryotic PSTK (the best hit of the characterised Methanocaldococcus
     jannaschii PSTK is a KTI12), the name given by NCBI: L-seryl-tRNA(Sec) kinase, O-phosphoseryl-tRNA(Sec) kinase or
     PSTK (the Asgard "AAA family ATPase" and adenylyl-sulfate kinases are left out); the best hits of both kinds are
     reported for all;
  2. complete proteins: --min-length to --max-length residues (leaves out fragments and gene models joined to other
     genes), no "partial" in the name, at most 2% X;
  3. new taxa: genera absent from the PSTK tree, and no vertebrates (birds, reptiles, mammals, already sampled);
     the taxonomic group of the genus from config/taxon_groups.tsv (genera without a group are left out);
  4. one sequence per genus, the longest.
Every candidate and the reason it was kept or left out: <out>/selection.tsv.
Trees, as in stage2/machinery_phylogeny.py: the sequences of results/trees/pstk.faa and the added ones (named
Genus_species_accession, group "added" in pstk.groups.tsv), MAFFT L-INS-i; added sequences that cover less than
--min-core of the core columns are removed and the rest realigned; ClipKIT gappy mode (columns with at least 90%
gaps removed); IQ-TREE 2 with ModelFinder (LG, WAG, JTT) and 1000 ultrafast bootstraps, then LG+C60+F+G with the
PMSF approximation on that tree. Two sets: the eukaryotes added (<out>/), and the eukaryotes and archaea added
(<out>_archaea/); KTI12 stays the outgroup in both, so that the archaeal PSTK show where they join the tree.
Outputs: <dir>/pstk.{faa,groups.tsv,removed.tsv,aln,trimmed,treefile,iqtree} and <dir>/pstk.pmsf.*.
"""
import argparse, collections, os, re, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from machinery_lba import taxon, iqtree, finished           # noqa: E402
from machinery_phylogeny import align, core_filter, read_fasta   # noqa: E402
from taxonomy import load_taxa                              # noqa: E402

VERTEBRATES = {"Mammals", "Birds", "Reptiles", "Amphibians", "Ray-finned fishes", "Lobe-finned fishes",
               "Jawless fishes", "Cartilaginous fishes"}
ARCHAEAL_PSTK = r"L-seryl-tRNA\(Sec\) kinase|O-phosphoseryl-tRNA\(Sec\) kinase|\bPSTK\b"


def downloads(paths):
    """accession -> (organism, description, sequence, files), from NCBI FASTA downloads (>ACC description [Organism])."""
    out = collections.OrderedDict()
    for p in paths:
        for name, seq in read_fasta(p):
            m = re.match(r"(\S+)\s+(.*?)\s*\[([^\]]+)\]\s*$", name)
            if not m: continue
            acc, desc, org = m.groups()
            if acc in out: out[acc][3].append(os.path.basename(p))
            else: out[acc] = [org, desc, seq.upper().rstrip("*"), [os.path.basename(p)]]
    return out


def genus_of(org):
    return next(t for t in org.split() if t not in ("Candidatus", "unclassified"))


def leaf_name(acc, org):
    words = [w for w in org.split() if w not in ("Candidatus", "unclassified")][:2]
    return "_".join(re.sub(r"[^A-Za-z0-9]", "", w) for w in words) + "_" + acc


def best_hits(queries, refs, work):
    """query -> {reference: best bit score}, legacy BLASTP (no filter)."""
    with tempfile.TemporaryDirectory(dir=work) as d:
        with open(f"{d}/db.faa", "w") as o:
            for n, s in refs.items(): o.write(f">{n}\n{s.replace('-', '')}\n")
        subprocess.run(["formatdb", "-i", "db.faa", "-p", "T", "-o", "F"], cwd=d, check=True, stdout=subprocess.DEVNULL)
        with open(f"{d}/q.faa", "w") as o:
            for n, s in queries.items(): o.write(f">{n}\n{s}\n")
        out = subprocess.run(["blastall", "-p", "blastp", "-d", f"{d}/db.faa", "-i", f"{d}/q.faa", "-e", "10", "-F", "F",
                              "-m", "8", "-b", "500", "-v", "500"], capture_output=True, text=True, check=True).stdout
    hits = collections.defaultdict(dict)
    for l in out.splitlines():
        f = l.split("\t")
        hits[f[0]][f[1]] = max(hits[f[0]].get(f[1], 0.0), float(f[11]))
    return hits


def select(a, base, groups, taxa):
    tax = {n: taxon(n, groups, taxa) for n in base}
    refs = {n: s for n, s in base.items() if tax[n] != "Appendicularians"}
    in_tree = {re.split(r"[_.]", n)[0] for n in base}
    cands = downloads(a.eukaryotes)
    arch = downloads([a.archaea]) if a.archaea else {}
    for k, v in arch.items(): cands.setdefault(k, v + ["archaea"])
    hits = best_hits({k: v[2] for k, v in cands.items()}, refs, a.work)
    rows, chosen = [], {}
    for acc, v in cands.items():
        org, desc, seq, files = v[:4]
        g = genus_of(org)
        grp = taxa.get(g, "unclassified")
        archaeal = len(v) > 4 or grp == "Archaea"       # also the archaea of the PSI-BLAST search
        h = hits.get(acc, {})
        best_p = max(((s, n) for n, s in h.items() if tax[n] != "outgroup"), default=(0.0, ""))
        best_k = max(((s, n) for n, s in h.items() if tax[n] == "outgroup"), default=(0.0, ""))
        family = "PSTK" if best_p[0] > best_k[0] else ("KTI12" if best_k[0] > 0 else "no hit")
        reason = ""
        if g == "Oikopleura": reason = "O. dioica (in the tree)"
        elif grp == "Bacteria": reason = "bacterial sequence"
        elif archaeal and not re.search(ARCHAEAL_PSTK, desc): reason = "archaeal sequence not named PSTK"
        elif not archaeal and family != "PSTK": reason = f"best hit {family}"
        elif not archaeal and best_p[0] < a.min_bits: reason = f"best PSTK hit below {a.min_bits:g} bits"
        elif not a.min_length <= len(seq) <= a.max_length: reason = f"length {len(seq)}"
        elif "partial" in desc.lower(): reason = "partial"
        elif seq.count("X") > 0.02 * len(seq): reason = "more than 2% X"
        elif g in in_tree: reason = "genus in the tree"
        elif grp in VERTEBRATES: reason = f"{grp}: already sampled"
        elif grp == "unclassified": reason = "genus without group in config/taxon_groups.tsv"
        rows.append(dict(accession=acc, organism=org, description=desc, length=len(seq), searches=",".join(files),
                         group=grp, best_pstk=f"{best_p[0]:.1f} {best_p[1]}".strip(),
                         best_kti12=f"{best_k[0]:.1f} {best_k[1]}".strip(), family=family, decision="",
                         reason=reason, archaeal=archaeal))
        if not reason:
            prev = chosen.get(g)
            if prev is None or len(seq) > len(cands[prev][2]): chosen[g] = acc
    keep = set(chosen.values())
    for r in rows:
        if not r["reason"] and r["accession"] not in keep: r["reason"] = "another sequence of the genus is longer"
        r["decision"] = "added" if r["accession"] in keep else "left out"
    arch_acc = {r["accession"] for r in rows if r["archaeal"]}
    added = {leaf_name(acc, cands[acc][0]): (cands[acc][2], acc in arch_acc) for acc in keep}
    return rows, added


def build(d, base, groups, added, a):
    os.makedirs(d, exist_ok=True)
    b = f"{d}/pstk"
    seqs = dict(base); grp = dict(groups)
    for n, (s, _) in added.items(): seqs[n], grp[n] = s, "added"
    removed = []

    def write():
        with open(f"{b}.faa", "w") as o:
            for n, s in seqs.items(): o.write(f">{n}\n{s.replace('-', '')}\n")
        with open(f"{b}.groups.tsv", "w") as o:
            for n in seqs: o.write(f"{n}\t{grp[n]}\n")
        with open(f"{b}.removed.tsv", "w") as o:
            o.write("sequence\tgroup\treason\n" + "".join("\t".join(r) + "\n" for r in removed))
    if not os.path.exists(f"{b}.trimmed"):
        write(); align(b, a.cpus)
        partial = {n: c for n, c in core_filter(f"{b}.aln", grp, a.min_core).items() if grp[n] == "added"}
        if partial:
            for n, c in partial.items():
                removed.append((n, "added", f"covers {c:.0%} of the core columns")); seqs.pop(n); grp.pop(n)
            write(); align(b, a.cpus)
        subprocess.run(["clipkit", f"{b}.aln", "-o", f"{b}.trimmed", "-m", "gappy", "-g", "0.9"], check=True,
                       stdout=subprocess.DEVNULL)
    if not finished(b):
        iqtree(["-s", f"{b}.trimmed", "-m", "MFP", "-mset", "LG,WAG,JTT", "-B", "1000", "--prefix", b], a.cpus)
    if not finished(f"{b}.pmsf"):
        iqtree(["-s", f"{b}.trimmed", "-m", "LG+C60+F+G", "-ft", f"{b}.treefile", "-B", "1000", "--prefix", f"{b}.pmsf"],
               a.cpus)
    print(d, sum(1 for v in grp.values() if v == "added"), "added sequences in the tree", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eukaryotes", nargs="+", default=["curation/pstk_blast/eukaryota.fa",
                                                        "curation/pstk_blast/opisthokonta.fa",
                                                        "curation/pstk_blast/psiblast.fa"])
    ap.add_argument("--archaea", default="curation/pstk_blast/archaea.fa")
    ap.add_argument("--trees", default="results/trees")
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--min-bits", type=float, default=50, help="minimum bit score of the best PSTK hit (eukaryotes)")
    ap.add_argument("--min-length", type=int, default=200)
    ap.add_argument("--max-length", type=int, default=450)
    ap.add_argument("--min-core", type=float, default=0.5)
    ap.add_argument("--select-only", action="store_true", help="write the selection table and stop")
    ap.add_argument("--work", default="work")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", default="results/trees/pstk_enriched")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); os.makedirs(a.work, exist_ok=True)
    taxa = load_taxa(a.taxa)
    base = dict(read_fasta(f"{a.trees}/pstk.faa"))
    groups = dict(l.rstrip("\n").split("\t")[:2] for l in open(f"{a.trees}/pstk.groups.tsv"))
    rows, added = select(a, base, groups, taxa)
    cols = ["accession", "organism", "description", "length", "searches", "group", "best_pstk", "best_kti12", "family",
            "decision", "reason"]
    with open(f"{a.out}/selection.tsv", "w") as o:
        o.write("\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows))
    euk = {n: v for n, v in added.items() if not v[1]}
    print(f"{len(rows)} candidates; added: {len(euk)} eukaryotic, {len(added) - len(euk)} archaeal", flush=True)
    for n in sorted(added): print(" ", n, taxon(n, {}, taxa), flush=True)
    if a.select_only: return
    build(a.out, base, groups, euk, a)
    build(a.out + "_archaea", base, groups, added, a)


if __name__ == "__main__":
    main()
