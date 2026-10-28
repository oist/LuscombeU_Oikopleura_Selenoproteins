#!/usr/bin/env python3
"""Write a copy of an Infernal covariance model with gathering (GA), trusted (TC) and noise (NC) cutoffs, so that
programs that search with --cut_ga (the pipeline's default for --secis_cm) use the chosen threshold.

The lines are written in the CM header after NULL, where Infernal and Rfam put them; existing GA, TC and NC lines of
the CM header are replaced. The parameters of the model are not changed: scores and E-values stay the same.
usage: cm_thresholds.py --ga 20 [--tc 20] [--nc 19.99] in.cm out.cm
"""
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ga", type=float, required=True)
    ap.add_argument("--tc", type=float)
    ap.add_argument("--nc", type=float)
    ap.add_argument("cm")
    ap.add_argument("out")
    a = ap.parse_args()
    tc = a.ga if a.tc is None else a.tc
    nc = a.ga - 0.01 if a.nc is None else a.nc
    lines = open(a.cm).read().split("\n")
    out, in_cm_header, done = [], False, False
    for l in lines:
        if l.startswith("INFERNAL1/"): in_cm_header, done = True, False
        if in_cm_header and l == "CM": in_cm_header = False
        if in_cm_header and l.split(" ")[0] in ("GA", "TC", "NC"): continue
        out.append(l)
        if in_cm_header and l.startswith("NULL ") and not done:
            out += [f"GA       {a.ga:.2f}", f"TC       {tc:.2f}", f"NC       {nc:.2f}"]
            done = True
    if not done: raise SystemExit(f"{a.cm}: no NULL line in a CM header")
    open(a.out, "w").write("\n".join(out))


if __name__ == "__main__":
    main()
