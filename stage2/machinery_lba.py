#!/usr/bin/env python3
"""Long-branch analyses of the machinery trees.

Appendicularian proteins evolve fast (O. dioica about twice as fast as Ciona; Berna et al. 2012), so their position in
the trees of stage2/machinery_phylogeny.py is exposed to long-branch attraction. Four standard checks, each on the
trimmed alignment of results/trees/<family>.trimmed:
  pmsf      LG+C60+F+G, a 60-class site-heterogeneous model, with the posterior mean site frequency approximation
            (PMSF) on the maximum-likelihood tree as guide
  slowfast  the fastest sites (by IQ-TREE site rates on the maximum-likelihood tree) removed in steps of
            --step (default 10 %, up to --max-removed), each alignment re-analysed with ModelFinder
  taxa      sequences matching --exclude (non-target long branches, for example the insect-gut fungi in PSTK) removed
  au        approximately unbiased test of the maximum-likelihood tree against the best tree in which the
            appendicularian sequences group with the other tunicates
For each tree, the summary records the support of the clades that the conclusions rest on: the appendicularian clade,
the family against the outgroup (PSTK) and the appendicularians with the other tunicates.
Outputs: <out>/<family>.* and <out>/<family>.summary.tsv
"""
import argparse, os, re, shutil, subprocess


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].strip(), []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def iqtree(args, cpus):
    subprocess.run(["iqtree2"] + args + ["-T", str(cpus), "--quiet", "--seed", "1"], check=True)


def best_model(iqtree_report):
    m = re.search(r"Best-fit model according to BIC: (\S+)", open(iqtree_report).read())
    return m.group(1) if m else "LG+G4"


def supports(treefile, groups):
    """Support of the clades the conclusions rest on, in the unrooted tree (as bipartitions)."""
    from Bio import Phylo
    t = Phylo.read(treefile, "newick")
    for c in t.get_nonterminals():
        if c.name and re.fullmatch(r"\d+(\.\d+)?", c.name): c.confidence, c.name = float(c.name), None
    names = {x.name for x in t.get_terminals()}
    app = {n for n in names if groups.get(n) in ("appendicularian", "curated") and not n.startswith("Fbor")}
    tun = {n for n in names if groups.get(n) == "tunicate"}
    out = {n for n in names if groups.get(n) == "outgroup"}

    def split_support(side):
        if not side or side == names: return "NA"
        for c in t.get_nonterminals():
            leaves = {x.name for x in c.get_terminals()}
            if leaves == side or leaves == names - side:
                return f"{c.confidence:.0f}" if c.confidence is not None else "yes"
        return "no"
    return {"appendicularians": split_support(app), "appendicularians+tunicates": split_support(app | tun),
            "family vs outgroup": split_support(names - out) if out else "NA"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", nargs="+", default=["pstk", "SEPHS2", "SBP2"])
    ap.add_argument("--trees", default="results/trees")
    ap.add_argument("--analyses", nargs="+", default=["pmsf", "slowfast", "taxa", "au"])
    ap.add_argument("--step", type=float, default=0.1)
    ap.add_argument("--max-removed", type=float, default=0.4)
    ap.add_argument("--exclude", default="Smittium|Zancudomyces", help="regular expression of sequences to remove")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", default="results/trees/lba")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for fam in a.families:
        aln, ml = f"{a.trees}/{fam}.trimmed", f"{a.trees}/{fam}.treefile"
        groups = dict(l.rstrip("\n").split("\t") for l in open(f"{a.trees}/{fam}.groups.tsv"))
        # the inputs must come from the same run of machinery_phylogeny.py (not one being rewritten)
        from Bio import Phylo
        aln_names = {n for n, _ in read_fasta(aln)}
        tree_names = {x.name for x in Phylo.read(ml, "newick").get_terminals()}
        if aln_names != tree_names or not aln_names <= set(groups):
            raise SystemExit(f"{fam}: alignment, tree and groups do not match (was machinery-phylogeny still running?)")
        model = best_model(f"{a.trees}/{fam}.iqtree")
        base = f"{a.out}/{fam}"
        rows = [("analysis", "model", "sites", "sequences") + tuple(supports(ml, groups).keys())]
        seqs = list(read_fasta(aln)); nsites = len(seqs[0][1])
        rows.append(("reference", model, nsites, len(seqs)) + tuple(supports(ml, groups).values()))
        if "pmsf" in a.analyses:
            if not os.path.exists(f"{base}.pmsf.treefile"):
                iqtree(["-s", aln, "-m", "LG+C60+F+G", "-ft", ml, "-B", "1000", "--prefix", f"{base}.pmsf"], a.cpus)
            rows.append(("PMSF LG+C60+F+G", "LG+C60+F+G (PMSF)", nsites, len(seqs)) +
                        tuple(supports(f"{base}.pmsf.treefile", groups).values()))
        if "slowfast" in a.analyses:
            if not os.path.exists(f"{base}.rates.rate"):
                iqtree(["-s", aln, "-te", ml, "-m", model, "-wsr", "--prefix", f"{base}.rates"], a.cpus)
            rates = [float(l.split()[1]) for l in open(f"{base}.rates.rate") if l[:1].isdigit()]
            order = sorted(range(len(rates)), key=lambda i: rates[i])   # slowest first
            k = a.step
            while k <= a.max_removed + 1e-9:
                keep = sorted(order[:int(round(len(order) * (1 - k)))])
                tag = f"slowfast{int(round(k * 100))}"
                f = f"{base}.{tag}.fa"
                with open(f, "w") as o:
                    for n, s in seqs: o.write(f">{n}\n{''.join(s[i] for i in keep)}\n")
                if not os.path.exists(f"{base}.{tag}.treefile"):
                    iqtree(["-s", f, "-m", "MFP", "-mset", "LG,WAG,JTT", "-B", "1000", "--prefix", f"{base}.{tag}"], a.cpus)
                rows.append((f"fastest {int(round(k * 100))} % of sites removed", best_model(f"{base}.{tag}.iqtree"),
                             len(keep), len(seqs)) + tuple(supports(f"{base}.{tag}.treefile", groups).values()))
                k += a.step
        if "taxa" in a.analyses and a.exclude:
            kept = [(n, s) for n, s in seqs if not re.search(a.exclude, n)]
            if len(kept) < len(seqs):
                f = f"{base}.taxa.fa"
                with open(f, "w") as o:
                    for n, s in kept: o.write(f">{n}\n{s}\n")
                if not os.path.exists(f"{base}.taxa.treefile"):
                    iqtree(["-s", f, "-m", "MFP", "-mset", "LG,WAG,JTT", "-B", "1000", "--prefix", f"{base}.taxa"], a.cpus)
                rows.append((f"without /{a.exclude}/", best_model(f"{base}.taxa.iqtree"), nsites, len(kept)) +
                            tuple(supports(f"{base}.taxa.treefile", groups).values()))
        if "au" in a.analyses:
            names = [n for n, _ in seqs]
            ingroup = [n for n in names if groups.get(n) in ("appendicularian", "curated", "tunicate")
                       and not n.startswith("Fbor")]
            rest = [n for n in names if n not in ingroup]
            with open(f"{base}.constraint.nwk", "w") as o: o.write(f"(({','.join(ingroup)}),{','.join(rest)});\n")
            if not os.path.exists(f"{base}.constrained.treefile"):
                iqtree(["-s", aln, "-m", model, "-g", f"{base}.constraint.nwk", "--prefix", f"{base}.constrained"], a.cpus)
            with open(f"{base}.candidates.treels", "w") as o:
                o.write(open(ml).read().strip() + "\n" + open(f"{base}.constrained.treefile").read().strip() + "\n")
            if not os.path.exists(f"{base}.au.iqtree"):
                iqtree(["-s", aln, "-m", model, "-z", f"{base}.candidates.treels", "-te", ml, "-n", "0", "-zb", "10000",
                        "-au", "--prefix", f"{base}.au"], a.cpus)
            text = open(f"{base}.au.iqtree").read()
            m = re.findall(r"^\s+2\s+.*$", text, re.M)
            p_au = m[0].split()[-2] if m else "see report"
            rows.append(("AU test: appendicularians with the other tunicates", model, nsites, len(seqs),
                         f"p-AU of the constrained tree = {p_au}", "", ""))
        with open(f"{base}.summary.tsv", "w") as o:
            for r in rows: o.write("\t".join(map(str, r)) + "\n")
        print(open(f"{base}.summary.tsv").read(), flush=True)


if __name__ == "__main__":
    main()
