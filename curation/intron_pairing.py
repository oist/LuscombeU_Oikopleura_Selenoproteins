#!/usr/bin/env python3
"""Paired intron ends and splice sites in the O. dioica gold-standard introns, with and without the insertions.

Reads results/tables/intron_content.tsv (--content: ends_z, sites, repeated_bases, loci) and
results/tables/intron_presence.tsv (--presence: status in the assembly of another individual of the species) and
prints, for the introns of --genomes with an ends z-score: the number whose ends pair with z at most --z, how many of
them are not GT-AG, the same for the other introns, and the medians of z for the GT-AG and other introns (all
assemblies and each one); first for all introns, then without those that are absent from the other individual's
assembly (insertions that the transcripts span), which are listed. Output: stdout, and --out (TSV of the counts).
"""
import argparse, csv, statistics


def summary(rows, z):
    strong = [r for r in rows if float(r["ends_z"]) <= z]
    other = [r for r in rows if float(r["ends_z"]) > z]
    nong = lambda L: sum(r["sites"] != "GT-AG" for r in L)
    med = lambda L: statistics.median(float(r["ends_z"]) for r in L) if L else float("nan")
    return dict(introns=len(rows), strong=len(strong), strong_not_gtag=nong(strong), other=len(other),
                other_not_gtag=nong(other), median_not_gtag=med([r for r in rows if r["sites"] != "GT-AG"]),
                median_gtag=med([r for r in rows if r["sites"] == "GT-AG"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--content", default="results/tables/intron_content.tsv")
    ap.add_argument("--presence", default="results/tables/intron_presence.tsv")
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3"])
    ap.add_argument("--z", type=float, default=-3.0)
    ap.add_argument("--out", default="results/tables/intron_pairing.tsv")
    a = ap.parse_args()
    status = {(r["genome"], r["gene"], r["intron"]): r["status"]
              for r in csv.DictReader(open(a.presence), delimiter="\t")}
    content = [r for r in csv.DictReader(open(a.content), delimiter="\t") if r["genome"] in a.genomes]
    scored = [r for r in content if r["ends_z"] not in ("", None)]
    absent = [r for r in content if status.get((r["genome"], r["gene"], r["intron"])) == "absent"]
    print("introns absent from the other individual's assembly:")
    for r in absent:
        print(f"  {r['genome']} {r['gene']} intron {r['intron']}: {r['length']} nt, {r['sites']}, "
              f"{r['repeated_bases']} repeated bases, {r['loci'] or 0} loci, ends z {r['ends_z'] or '-'}, fold z {r['mfe_z'] or '-'}")
    keep = [r for r in scored if status.get((r["genome"], r["gene"], r["intron"])) != "absent"]
    out = []
    for label, rows in [("all", scored), ("without insertions", keep)]:
        for g in ["all"] + a.genomes:
            sub = rows if g == "all" else [r for r in rows if r["genome"] == g]
            s = summary(sub, a.z)
            out.append(dict(set=label, genome=g, **s))
            print(f"{label:18s} {g:5s} {s['introns']} introns: z <= {a.z}: {s['strong']} ({s['strong_not_gtag']} not GT-AG); "
                  f"others {s['other']} ({s['other_not_gtag']} not GT-AG); median z not GT-AG {s['median_not_gtag']:.2f}, "
                  f"GT-AG {s['median_gtag']:.2f}")
    with open(a.out, "w") as o:
        keys = list(out[0])
        o.write("\t".join(keys) + "\n")
        for r in out: o.write("\t".join(f"{r[k]:.2f}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")


if __name__ == "__main__":
    main()
