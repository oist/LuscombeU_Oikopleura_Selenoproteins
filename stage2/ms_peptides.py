#!/usr/bin/env python3
"""Peptides for a mass-spectrometry validation of the gold-standard selenoproteins (planning; no MS data).

For each gold-standard protein with a selenocysteine (curation/gold_standard/<genome>.proteins.faa), in silico
digestion with trypsin (after K or R, not before P), Arg-C (after R, not before P), Glu-C (after E), chymotrypsin
(after F, Y or W, not before P; high-specificity rule) and Asp-N (before D), up to one missed cleavage:
  - sec_peptide_<protease>: the shortest peptide of 7 to 30 residues that contains the selenocysteine (peptides of 26
    to 30 residues only when no shorter one exists), or "-";
  - tryptic_after_sec: number of distinct tryptic peptides (7-30 residues, no missed cleavage) entirely after the
    selenocysteine: they exist only if the UGA codon is read through, so their detection shows read-through without
    the selenocysteine peptide itself;
  - mass_kda: average mass of the protein (selenocysteine 150.0379 Da), to place it on a gel lane;
  - sec_peptide_mono_da: monoisotopic mass of the tryptic selenocysteine peptide with carbamidomethyl on C and U.
Output: --table.
"""
import argparse

AVG = {"G": 57.0513, "A": 71.0779, "S": 87.0773, "P": 97.1152, "V": 99.1311, "T": 101.1039, "C": 103.1429,
       "L": 113.1576, "I": 113.1576, "N": 114.1026, "D": 115.0874, "Q": 128.1292, "K": 128.1723, "E": 129.1140,
       "M": 131.1961, "H": 137.1393, "F": 147.1739, "R": 156.1857, "Y": 163.1733, "W": 186.2099, "U": 150.0379}
MONO = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768, "C": 103.00919,
        "L": 113.08406, "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496,
        "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333,
        "W": 186.07931, "U": 150.95364}
WATER_AVG, WATER_MONO, CAM = 18.0153, 18.01056, 57.02146

RULES = {
    "trypsin": lambda s, i: s[i] in "KR" and (i + 1 >= len(s) or s[i + 1] != "P"),
    "argc": lambda s, i: s[i] == "R" and (i + 1 >= len(s) or s[i + 1] != "P"),
    "gluc": lambda s, i: s[i] == "E",
    "chymotrypsin": lambda s, i: s[i] in "FYW" and (i + 1 >= len(s) or s[i + 1] != "P"),
    "aspn": lambda s, i: i + 1 < len(s) and s[i + 1] == "D",
}


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split("|")[0].strip(); out[name] = ""
        elif name:
            out[name] += line.strip()
    return out


def peptides(seq, rule, missed=1):
    cuts = sorted({0, len(seq)} | {i + 1 for i in range(len(seq) - 1) if rule(seq, i)})
    return [(cuts[a], cuts[b], b - a - 1) for a in range(len(cuts) - 1)
            for b in range(a + 1, min(a + 2 + missed, len(cuts)))]


def sec_peptide(seq, u, rule):
    cands = [seq[a:b] for a, b, _ in peptides(seq, rule) if a <= u < b and 7 <= b - a <= 30]
    cands.sort(key=lambda p: (len(p) > 25, len(p)))
    return cands[0] if cands else "-"


def mono_cam(pep):
    return sum(MONO[x] for x in pep) + WATER_MONO + CAM * sum(x in "CU" for x in pep)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--genomes", nargs="+", default=["OKI", "OSKA", "OdB3"])
    ap.add_argument("--table", default="results/tables/stage2_ms_peptides.tsv")
    a = ap.parse_args()
    cols = ["genome", "gene", "length", "mass_kda", "sec_position", "residues_after_sec", "tryptic_after_sec"] + \
        [f"sec_peptide_{p}" for p in RULES] + ["sec_peptide_mono_da"]
    rows = []
    for g in a.genomes:
        for gene, seq in read_fasta(f"{a.gold}/{g}.proteins.faa").items():
            seq = seq.rstrip("*")
            if "U" not in seq:
                continue
            u = seq.index("U")
            r = {"genome": g, "gene": gene, "length": len(seq), "sec_position": u + 1,
                 "residues_after_sec": len(seq) - u - 1,
                 "mass_kda": f"{(sum(AVG[x] for x in seq if x in AVG) + WATER_AVG) / 1000:.1f}",
                 "tryptic_after_sec": len({seq[s:e] for s, e, m in peptides(seq, RULES["trypsin"])
                                           if s > u and m == 0 and 7 <= e - s <= 30})}
            for p, rule in RULES.items():
                r[f"sec_peptide_{p}"] = sec_peptide(seq, u, rule)
            t = r["sec_peptide_trypsin"]
            r["sec_peptide_mono_da"] = f"{mono_cam(t):.4f}" if t != "-" else ""
            rows.append(r)
    with open(a.table, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")
    for g in a.genomes:
        rs = [r for r in rows if r["genome"] == g]
        tryp = sum(r["sec_peptide_trypsin"] != "-" for r in rs)
        anyp = sum(any(r[f"sec_peptide_{p}"] != "-" for p in RULES) for r in rs)
        down = sum(r["tryptic_after_sec"] >= 2 for r in rs)
        cterm = [r["gene"] for r in rs if r["residues_after_sec"] <= 2]
        print(f"{g}: {len(rs)} genes; tryptic Sec peptide {tryp}; any protease {anyp}; >=2 tryptic peptides after "
              f"the Sec {down}; Sec within 2 residues of the C-terminus: {', '.join(cterm)}")
    print(f"wrote {a.table}")


if __name__ == "__main__":
    main()
