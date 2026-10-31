#!/usr/bin/env python3
"""SELENOI: the junction that removes the TGA segment in some O. dioica transcripts (docs/plan.md, 2026-10-31: "The
only placement of that junction consistent with the genome has AG...AG ends, not one of the O. dioica intron classes:
to check").

Inputs: the codon regions and transcripts of make selenoi-regions (results/sequences/selenoi_regions.fa; region in
gene orientation, TGA at offset 5,000). Each transcript of OKI2018, OSKA2016 and OdB3 is aligned to the codon region of
its assembly with minimap2 (-x splice:hq -un: no preference for GT-AG, as in curation/gold_standard.py, -c --cs=long).
For each junction from 450 nt before the TGA to the end of the region:
  - placements: all positions of the junction that give the same transcript sequence (sliding over identical bases,
    up to 12 nt each way), with their splice sites;
  - canonical: placements within 15 nt with GT-AG, GC-AG, GA-AG or GG-AG splice sites if 1 or 2 mismatches between
    transcript and genome are accepted at the junction;
  - the 12 bases before the 5' end and before the 3' end of the removed segment (to show repeats).
Output: --table (one row per transcript junction), coordinates relative to the first base of the TGA.
"""
import argparse, re, subprocess, tempfile, os

CANON = {"GT-AG", "GC-AG", "GA-AG", "GG-AG"}
TOKEN = re.compile(r"(=[A-Za-z]+)|(\*[a-z][a-z])|(\+[a-z]+)|(-[a-z]+)|(~[a-z]{2}\d+[a-z]{2})")


def read_fasta(path):
    out, name = {}, None
    for l in open(path):
        if l.startswith(">"): name = l[1:].strip(); out[name] = []
        else: out[name].append(l.strip())
    return {k: "".join(v) for k, v in out.items()}


def junctions(paf_line):
    """(transcript, [(intron start, intron end exclusive)]) on the target, from a PAF line with a long cs tag."""
    f = paf_line.split("\t")
    pos = int(f[7]); cs = next(x for x in f if x.startswith("cs:Z:"))[5:]
    out = []
    for m in TOKEN.finditer(cs):
        t = m.group(0)
        if t[0] in "=-": pos += len(t) - 1
        elif t[0] == "*": pos += 1
        elif t[0] == "~":
            n = int(re.search(r"\d+", t).group()); out.append((pos, pos + n)); pos += n
    return f[0], f[4], int(f[2]), int(f[3]), out


def placements(r, e, s, span=12):
    """Equivalent placements of intron [e, s) by sliding over identical bases."""
    out = []
    for k in range(-span, span + 1):
        if k > 0: ok = all(r[e + i] == r[s + i] for i in range(k))
        else: ok = all(r[e + i] == r[s + i] for i in range(k, 0))
        if ok: out.append((k, r[e + k:e + k + 2] + "-" + r[s + k - 2:s + k]))
    return out


def canonical_near(r, e, s, span=15, max_mm=2):
    out = []
    for k in range(-span, span + 1):
        mm = sum(r[e + i] != r[s + i] for i in (range(k) if k > 0 else range(k, 0)))
        d = r[e + k:e + k + 2] + "-" + r[s + k - 2:s + k]
        if 0 < mm <= max_mm and d in CANON: out.append(f"{k:+d}:{d} ({mm} mismatch{'es' if mm > 1 else ''})")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", default="results/sequences/selenoi_regions.fa")
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3"])
    ap.add_argument("--offset", type=int, default=5000, help="position of the TGA in the codon regions")
    ap.add_argument("--table", default="results/tables/stage2_selenoi_junction.tsv")
    a = ap.parse_args()
    seqs = read_fasta(a.regions)
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for g in a.genomes:
            reg = next(k for k in seqs if k.startswith(f"{g}|codon|"))
            r = seqs[reg].upper()
            assert r[a.offset:a.offset + 3] == "TGA", (g, r[a.offset:a.offset + 3])
            open(f"{tmp}/{g}.reg.fa", "w").write(f">{g}\n{r}\n")
            with open(f"{tmp}/{g}.tx.fa", "w") as o:
                for k, s in seqs.items():
                    f = k.split("|")
                    if f[0] == g and f[1] == "transcript": o.write(f">{f[2]}\n{s}\n")
            out = subprocess.run(["minimap2", "-x", "splice:hq", "-un", "-c", "--cs=long", f"{tmp}/{g}.reg.fa",
                                  f"{tmp}/{g}.tx.fa"], capture_output=True, text=True, check=True).stdout
            for l in out.splitlines():
                tx, strand, qs, qe, introns = junctions(l)
                f = l.split("\t"); ts, te = int(f[7]), int(f[8])
                if te < a.offset - 450 or strand != "+": continue
                keeps = ts <= a.offset and a.offset + 3 <= te and not any(e <= a.offset < s for e, s in introns)
                for e, s in introns:
                    if s < a.offset - 450: continue
                    pl = placements(r, e, s)
                    canon = [f"{k:+d}:{d}" for k, d in pl if d in CANON]
                    rows.append({"genome": g, "transcript": tx, "codon": "retained" if keeps else "removed",
                                 "intron": f"{e - a.offset}..{s - a.offset}", "length": s - e,
                                 "contains_codon": "yes" if e <= a.offset < s else "no",
                                 "sites": r[e:e + 2] + "-" + r[s - 2:s],
                                 "equivalent_placements": ";".join(f"{k:+d}:{d}" for k, d in pl),
                                 "canonical_placement": ";".join(canon) or "none",
                                 "canonical_with_mismatches": ";".join(canonical_near(r, e, s)) or "none",
                                 "before_5prime_end": r[e - 12:e].lower() + "|" + r[e:e + 2],
                                 "before_3prime_end": r[s - 12:s] + "|" + r[s:s + 2].lower()})
    cols = ["genome", "transcript", "codon", "intron", "length", "contains_codon", "sites", "equivalent_placements",
            "canonical_placement", "canonical_with_mismatches", "before_5prime_end", "before_3prime_end"]
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for x in rows: o.write("\t".join(str(x[c]) for c in cols) + "\n")
    print(f"wrote {a.table}: {len(rows)} junctions")
    for x in rows:
        if x["contains_codon"] == "yes":
            print(x["genome"], x["transcript"], x["intron"], x["length"], x["sites"], x["canonical_placement"],
                  x["canonical_with_mismatches"], x["before_5prime_end"], x["before_3prime_end"])


if __name__ == "__main__":
    main()
