#!/usr/bin/env python3
"""Presence of each gold-standard intron in the assembly of another individual of the same species.

An intron of a gold-standard model can be an insertion that the transcripts span: a transposon inserted inside an exon
in the assembled haplotype, absent from the alleles that the transcripts come from (or spliced out). Such an insertion
has the marks of the introns of interest here: sites outside GT-AG, bases repeated at its two ends (the duplicated
target site) and paired ends (terminal inverted repeats). The assembly of another individual tells the two apart: if
the intron is a polymorphic insertion, the other assembly can have the two exons contiguous.

For each intron of --introns (results/tables/intron_content.tsv: genome, gene, intron, contig, start, end, strand), three
queries are built from the genome of its assembly, on the gene strand: the junction without the intron (--flank nt of
each exon, joined), the donor junction (exon and intron start) and the acceptor junction (intron end and exon); each is
searched in the genome of the other individual (--pairs; BLASTN, E at most 1e-5). A query is found across its junction
when an HSP of at least --identity % covers three quarters of the flank on both sides of the junction point (so that
the bases repeated at the two ends of an intron do not make a joined query match an intron junction). Status: "present" (a
donor or acceptor junction found, the joined exons not), "absent" (the joined exons found contiguous, no intron
junction), "both" (both kinds: two haplotypes or copies in the other assembly), "not found" (neither: diverged or
missing locus). Exons shorter than the flank contribute the bases of the next intron or exon: such queries can miss.
Output: --out (TSV: genome, gene, intron, length, sites, other, status, joined_hit, donor_hit, acceptor_hit; a hit is
contig:start-end, identity).
Genome FASTA paths: --samplesheet files (sample, fasta columns), or --genome NAME=FASTA.
"""
import argparse, csv, gzip, os, subprocess, sys, tempfile


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))


def read_fasta(path, keep=None):
    out, name, seq = {}, None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for l in f:
            if l.startswith(">"):
                if name is not None and (keep is None or name in keep): out[name] = "".join(seq).upper()
                name, seq = l[1:].split()[0], []
            else: seq.append(l.strip())
    if name is not None and (keep is None or name in keep): out[name] = "".join(seq).upper()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--introns", default="results/tables/intron_content.tsv")
    ap.add_argument("--pairs", nargs="+", default=["OKI=KUM", "OSKA=AOM", "OdB3=Bar2"])
    ap.add_argument("--samplesheet", nargs="+", default=["stage1/samplesheet.csv", "stage2/run/samplesheet_final.csv"])
    ap.add_argument("--genome", nargs="*", default=[], help="NAME=FASTA, overriding the samplesheets")
    ap.add_argument("--flank", type=int, default=40)
    ap.add_argument("--identity", type=float, default=85.0)
    ap.add_argument("--work", default="work/intron_presence")
    ap.add_argument("--out", default="results/tables/intron_presence.tsv")
    a = ap.parse_args()
    paths = {}
    for s in a.samplesheet:
        if os.path.exists(s):
            for r in csv.DictReader(open(s)): paths.setdefault(r["sample"], r["fasta"])
    paths.update(dict(x.split("=", 1) for x in a.genome))
    pairs = dict(x.split("=") for x in a.pairs)
    introns = [r for r in csv.DictReader(open(a.introns), delimiter="\t") if r["genome"] in pairs]
    os.makedirs(a.work, exist_ok=True)
    F, q4 = a.flank, (3 * a.flank) // 4
    rows = []
    for g, other in pairs.items():
        mine = [r for r in introns if r["genome"] == g]
        if not mine: continue
        seqs = read_fasta(paths[g], keep={r["contig"] for r in mine})
        db = os.path.join(a.work, other)
        if not os.path.exists(db + ".nsq") and not os.path.exists(db + ".00.nsq"):
            plain = db + ".fa"
            with open(plain, "w") as o:
                for n, s in read_fasta(paths[other]).items(): o.write(f">{n}\n{s}\n")
            subprocess.run(["makeblastdb", "-dbtype", "nucl", "-in", plain, "-out", db], check=True, capture_output=True)
            os.remove(plain)
        queries = {}
        for r in mine:
            c, s0, e0 = r["contig"], int(r["start"]), int(r["end"])   # 1-based intron coordinates
            sq = seqs[c]
            left, intron, right = sq[s0 - 1 - F:s0 - 1], sq[s0 - 1:e0], sq[e0:e0 + F]
            if r["strand"] == "-": left, intron, right = rc(right), rc(intron), rc(left)
            key = f"{r['gene']}|{r['intron']}"
            queries[key + "|joined"] = left + right
            queries[key + "|donor"] = left + intron[:F]
            queries[key + "|acceptor"] = intron[-F:] + right
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/q.fa", "w").write("".join(f">{k}\n{v}\n" for k, v in queries.items()))
            out = subprocess.run(["blastn", "-task", "blastn", "-query", f"{d}/q.fa", "-db", db, "-evalue", "1e-5",
                                  "-max_target_seqs", "10", "-max_hsps", "5", "-outfmt",
                                  "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"],
                                 capture_output=True, text=True, check=True).stdout
        found = {}
        for l in out.splitlines():
            q, sid, pid, ln, qs, qe, ss, se = l.split("\t")[:8]
            if float(pid) < a.identity or int(qs) > F - q4 + 1 or int(qe) < F + q4: continue
            if q not in found or float(pid) > found[q][1]:
                found[q] = (f"{sid}:{min(int(ss), int(se))}-{max(int(ss), int(se))}", float(pid))
        for r in mine:
            key = f"{r['gene']}|{r['intron']}"
            j, dn, ac = (found.get(f"{key}|{k}") for k in ("joined", "donor", "acceptor"))
            status = ("both" if j and (dn or ac) else "absent" if j else "present" if (dn or ac) else "not found")
            fmt = lambda h: f"{h[0]} {h[1]:.1f}" if h else "-"
            rows.append([g, r["gene"], r["intron"], r["length"], r["sites"], other, status, fmt(j), fmt(dn), fmt(ac)])
    with open(a.out, "w") as o:
        o.write("genome\tgene\tintron\tlength\tsites\tother\tstatus\tjoined_hit\tdonor_hit\tacceptor_hit\n")
        for r in rows: o.write("\t".join(map(str, r)) + "\n")
    counts = {}
    for r in rows: counts[(r[0], r[6])] = counts.get((r[0], r[6]), 0) + 1
    print(a.out, "; ".join(f"{g} {s}: {n}" for (g, s), n in sorted(counts.items())))


if __name__ == "__main__":
    main()
