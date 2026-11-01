#!/usr/bin/env python3
"""Protein-to-genome alignments (miniprot) of the machinery reference proteins, as candidate models of
curation/gold_standard.py --anchor sense --projections.

Queries (--queries, FASTA files whose headers start with the gene name: curation/machinery_references.faa, the human
proteins of each gene, and the <genome>.proteins.faa of a first gold-standard pass on the three O. dioica assemblies);
selenocysteine (U) and stop codons (*) are written X. Each genome of --genomes (stage1/samplesheet.csv) is searched
with miniprot (--gff, maximum intron size set from the genome size with -I). The alignments of a gene give exon
chains on the genome that hold the complete protein of a related assembly or species, where Selenoprofiles predicts
only the region of its profile and the transcripts of a species may be partial or absent.
Output: <out>/<genome>.projections.gff (miniprot GFF3: mRNA, CDS and stop_codon rows; Target = query name).
"""
import argparse, csv, gzip, os, subprocess


def read_fasta(path):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = l[1:].strip(), []
            else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--queries", nargs="+", required=True)
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan"])
    ap.add_argument("--samplesheet", default="stage1/samplesheet.csv")
    ap.add_argument("--out", default="work/stage2/gold_standard_machinery/projections")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sheet = {r["sample"]: r for r in csv.DictReader(open(a.samplesheet))}
    q = f"{a.out}/queries.faa"
    with open(q, "w") as o:
        seen = set()
        for f in a.queries:
            for n, s in read_fasta(f):
                # gold-standard headers: <gene>|<genome>|...; references: <gene>|Human_<symbol>|<id>
                parts = n.split("|")
                name = "|".join(parts[:3]).replace(" ", "_")
                if name in seen: continue
                seen.add(name)
                o.write(f">{name}\n{s.replace('U', 'X').replace('*', '')}\n")
    print(f"{len(seen)} queries", flush=True)
    for g in a.genomes:
        out = f"{a.out}/{g}.projections.gff"
        with open(out + ".tmp", "w") as fh:
            subprocess.run(["miniprot", "-t", str(a.cpus), "-I", "--gff", sheet[g]["fasta"], q], stdout=fh, check=True)
        os.replace(out + ".tmp", out)
        n = sum(1 for l in open(out) if "\tmRNA\t" in l)
        print(f"{g}: {n} alignments", flush=True)


if __name__ == "__main__":
    main()
