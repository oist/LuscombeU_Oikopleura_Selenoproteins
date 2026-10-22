#!/usr/bin/env python3
"""Domain analysis of the machinery trees: is the tree of each domain the same as that of the other domain?

For each family, the trimmed alignment of results/trees/<family>.trimmed is split into two domains defined on a human
protein (curation/machinery_domains.tsv: UniProt coordinates and their source; sequences in
curation/machinery_domain_references.faa). The UniProt coordinates are carried to the human sequence of the alignment
by a global pairwise alignment, and each alignment column gets the domain of the human residue it holds (or of the
last human residue before it).

Only the sequences with residues in at least --min-coverage of the columns of both domains are kept, so that the two
domains are compared on the same sequences. Then, for each domain:
  - a maximum-likelihood tree (IQ-TREE, ModelFinder restricted to LG, WAG and JTT, 1,000 ultrafast bootstraps);
  - the supports of the clades that the conclusions rest on (as in stage2/machinery_lba.py): the appendicularian
    clade, the appendicularian clade with the other tunicates of the same form, and for PSTK the family against KTI12;
  - the best trees under the constraints of stage2/machinery_lba.py (appendicularians with the other tunicates; for
    SEPHS, inside the clade of all tunicate SEPHS and at the position of an ortholog);
  - AU tests (10,000 RELL replicates) of the trees of both domains, of the full-length tree (pruned to the same
    sequences) and of the constrained trees, against the data of this domain.
A domain whose data reject the tree of the other domain (p-AU < 0.05) indicates that the two domains do not share one
history, or that one of them is misled; a hypothesis rejected by one domain only points to that domain. Because the
tree of a short domain is poorly resolved in places, such rejections can come from weakly supported parts of the
trees; --conflicts therefore lists, for each threshold of --strong, the clades supported by at least that ultrafast
bootstrap in one domain tree that are incompatible with a clade so supported in the other (each clade described by
the taxonomic groups of its smaller side, with its support).
Display trees of the domains: make machinery-domain-figures (figures/phylogeny_<family>_collapsed_domain_<domain>.pdf).
Outputs: <out>/<family>.<domain>.*, --table (one row per family, domain and tested tree) and --conflicts.
"""
import argparse, os, re, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import machinery_lba as lba

GENOMES = ("OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor")


def iqtree(args, cpus):
    subprocess.run(["iqtree2"] + args + ["-T", str(cpus), "--quiet", "--seed", "1", "-redo"], check=True)


def column_map(aln, trimmed):
    """Index in the untrimmed alignment of each column of the trimmed alignment (in-order matching of columns)."""
    names = list(trimmed)
    cols = ["".join(aln[n][i] for n in names) for i in range(len(aln[names[0]]))]
    out, j = [], 0
    for k in range(len(trimmed[names[0]])):
        col = "".join(trimmed[n][k] for n in names)
        while cols[j] != col: j += 1
        out.append(j); j += 1
    return out


def residue_index(seq):
    """For each alignment column, the index (1-based) of the last residue of seq at or before it (0 before the first)."""
    out, k = [], 0
    for c in seq:
        if c != "-": k += 1
        out.append(k)
    return out


def map_coordinates(uniprot_seq, ref_seq):
    """UniProt position (1-based) -> position in ref_seq (1-based), by a global alignment; positions in gaps take the
    next aligned position."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a = al.align(uniprot_seq.replace("U", "C"), ref_seq.replace("X", "C").replace("U", "C"))[0]
    m = {}
    for (us, ue), (rs, re_) in zip(*a.aligned):
        for k in range(ue - us): m[us + k + 1] = rs + k + 1
    out, nxt = {}, None
    for p in range(len(uniprot_seq), 0, -1):
        nxt = m.get(p, nxt); out[p] = nxt if nxt is not None else len(ref_seq)
    return out


def write_fasta(path, seqs, cols):
    with open(path, "w") as o:
        for n, s in seqs.items(): o.write(f">{n}\n{''.join(s[i] for i in cols)}\n")


def prune(newick_file, keep, out):
    t = lba.read_tree(newick_file)
    for x in list(t.get_terminals()):
        if x.name not in keep: t.prune(x)
    for c in t.get_nonterminals(): c.confidence = None
    from Bio import Phylo
    Phylo.write(t, out, "newick")


def p_values(report, n):
    table = open(report).read().split("USER TREES", 1)[-1]
    out = []
    for i in range(1, n + 1):
        m = re.findall(rf"^\s+{i}\s+-[\d.]+\s.*$", table, re.M)
        out.append(m[0].split()[-2] if m else "NA")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", nargs="+", default=["pstk", "SEPHS2", "SBP2"])
    ap.add_argument("--trees", default="results/trees")
    ap.add_argument("--domains", default="curation/machinery_domains.tsv")
    ap.add_argument("--references", default="curation/machinery_domain_references.faa")
    ap.add_argument("--min-coverage", type=float, default=0.3)
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", default="results/trees/domains")
    ap.add_argument("--table", default="results/tables/stage2_machinery_domains.tsv")
    ap.add_argument("--conflicts", default="results/tables/stage2_machinery_domain_conflicts.tsv")
    ap.add_argument("--strong", type=float, nargs="+", default=[95, 90],
                    help="ultrafast bootstrap of a strongly supported clade (one conflict row per threshold)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    refs = {n.split()[0]: s for n, s in lba.read_fasta(a.references)}
    dom = {}
    for l in open(a.domains):
        f = l.rstrip("\n").split("\t")
        if f[0] == "family": continue
        dom.setdefault(f[0], []).append((f[1], f[2], f[3], int(f[4]), int(f[5])))
    taxa = lba.load_taxa(a.taxa)
    conflicts = [("family", "domains", "threshold", "strong clades in the first", "strong clades in the second",
                  "conflicting pairs", "conflicts (smaller side of each clade, by taxonomic group, and its support)")]
    table = [("family", "domain", "columns", "sequences", "model", "appendicularian clade",
              "appendicularians+tunicates", "family vs outgroup", "tested tree", "p-AU")]
    for fam in a.families:
        aln = dict(lba.read_fasta(f"{a.trees}/{fam}.aln"))
        trimmed = dict(lba.read_fasta(f"{a.trees}/{fam}.trimmed"))
        groups = dict(l.rstrip("\n").split("\t") for l in open(f"{a.trees}/{fam}.groups.tsv"))
        ml_full = f"{a.trees}/{fam}.treefile"
        pattern, uniprot = dom[fam][0][0], dom[fam][0][1]
        ref = next(n for n in trimmed if re.search(pattern, n))
        pos = map_coordinates(refs[uniprot], aln[ref].replace("-", ""))
        cmap = column_map(aln, trimmed)
        ridx = residue_index(aln[ref])
        cols = {}
        for k, j in enumerate(cmap):
            r = max(ridx[j], 1)
            for _, _, name, start, end in dom[fam]:
                if pos[start] <= r <= (pos[end] if end < len(refs[uniprot]) else 10 ** 9):
                    cols.setdefault(name, []).append(k); break
        names = [d[2] for d in dom[fam]]
        cover = {n: {d: sum(trimmed[n][k] != "-" for k in cols[d]) / len(cols[d]) for d in names} for n in trimmed}
        keep = [n for n in trimmed if all(cover[n][d] >= a.min_coverage for d in names)]
        dropped = [n for n in trimmed if n not in keep]
        print(f"{fam}: domains " + ", ".join(f"{d} {len(cols[d])} columns" for d in names) +
              f"; {len(keep)} sequences in both, {len(dropped)} left out: {' '.join(dropped)}", flush=True)
        seqs = {n: trimmed[n] for n in keep}
        tax = {n: lba.taxon(n, groups, taxa) for n in keep}
        excl = lba.ORTHOLOG_EXCLUDE.get(fam)
        app_full = lba.appendicularian_clade(ml_full, {n: lba.taxon(n, groups, taxa) for n in trimmed})
        sets = {"app": app_full & set(keep),
                "tun": {n for n in keep if tax[n] == "Other tunicates" and not (excl and re.search(excl, n))},
                "out": {n for n in keep if tax[n] == "outgroup"}}
        # constraints, as in machinery_lba.py, restricted to the kept sequences
        j = ",".join
        cons = {"appendicularians with the other tunicates":
                f"(({j(sorted(sets['app'] | sets['tun']))}),{j(sorted(set(keep) - sets['app'] - sets['tun']))});"}
        if excl:
            tunpar = {n for n in keep if tax[n] == "Other tunicates" and re.search(excl, n)}
            apppar = {n for n in keep if tax[n] == "Appendicularians" and n not in sets["app"]}
            alltun = sets["app"] | sets["tun"] | tunpar | apppar
            cons["inside the clade of all tunicate sequences"] = f"(({j(sorted(alltun))}),{j(sorted(set(keep) - alltun))});"
            cons["at the position of an ortholog"] = (f"((({j(sorted(sets['app'] | sets['tun']))}),({j(sorted(tunpar | apppar))})),"
                                                      f"{j(sorted(set(keep) - alltun))});")
        base = f"{a.out}/{fam}"
        prune(ml_full, set(keep), f"{base}.full.pruned.treefile")
        for d in names:
            p = f"{base}.{re.sub(r'[^A-Za-z0-9]+', '_', d)}"
            write_fasta(f"{p}.fa", seqs, cols[d])
            if not (os.path.exists(f"{p}.iqtree") and os.path.exists(f"{p}.treefile")):
                iqtree(["-s", f"{p}.fa", "-m", "MFP", "-mset", "LG,WAG,JTT", "-B", "1000", "--prefix", p], a.cpus)
        models = {d: lba.best_model(f"{base}.{re.sub(r'[^A-Za-z0-9]+', '_', d)}.iqtree") for d in names}
        for d in names:
            p = f"{base}.{re.sub(r'[^A-Za-z0-9]+', '_', d)}"
            candidates = [(f"tree of the {e} domain", f"{base}.{re.sub(r'[^A-Za-z0-9]+', '_', e)}.treefile") for e in names]
            candidates.append(("full-length tree", f"{base}.full.pruned.treefile"))
            for k, (label, nwk) in enumerate(cons.items()):
                cp = f"{p}.constraint{k}"
                with open(f"{cp}.nwk", "w") as o: o.write(nwk + "\n")
                if not (os.path.exists(f"{cp}.iqtree") and os.path.exists(f"{cp}.treefile")):
                    iqtree(["-s", f"{p}.fa", "-m", models[d], "-g", f"{cp}.nwk", "--prefix", cp], a.cpus)
                candidates.append((label, f"{cp}.treefile"))
            with open(f"{p}.candidates.treels", "w") as o:
                for _, f in candidates: o.write(open(f).read().strip() + "\n")
            if not os.path.exists(f"{p}.au.iqtree"):
                iqtree(["-s", f"{p}.fa", "-m", models[d], "-z", f"{p}.candidates.treels", "-te", f"{p}.treefile",
                        "-n", "0", "-zb", "10000", "-au", "--prefix", f"{p}.au"], a.cpus)
            sup = lba.supports(f"{p}.treefile", sets)
            pv = p_values(f"{p}.au.iqtree", len(candidates))
            for (label, _), pval in zip(candidates, pv):
                table.append((fam, d, len(cols[d]), len(keep), models[d]) + tuple(sup.values()) + (label, pval))
        f1, f2 = (f"{base}.{re.sub(r'[^A-Za-z0-9]+', '_', d)}.treefile" for d in names)
        for threshold in a.strong:
            s1, s2 = strong_splits(f1, threshold), strong_splits(f2, threshold)
            bad = [(x, sx, y, sy) for x, sx in s1 for y, sy in s2 if not compatible(x, y, set(keep))]
            desc = "; ".join(f"{describe(x, groups, taxa)} {sx:.0f} vs {describe(y, groups, taxa)} {sy:.0f}"
                             for x, sx, y, sy in bad)
            conflicts.append((fam, " / ".join(names), f"{threshold:g}", len(s1), len(s2), len(bad), desc or "none"))
    with open(a.table, "w") as o:
        for r in table: o.write("\t".join(map(str, r)) + "\n")
    with open(a.conflicts, "w") as o:
        for r in conflicts: o.write("\t".join(map(str, r)) + "\n")
    print(open(a.table).read()); print(open(a.conflicts).read())


def strong_splits(treefile, threshold):
    t = lba.read_tree(treefile)
    names = frozenset(x.name for x in t.get_terminals())
    out = []
    for c in t.get_nonterminals():
        side = frozenset(x.name for x in c.get_terminals())
        if c.confidence is not None and c.confidence >= threshold and 1 < len(side) < len(names) - 1:
            out.append((min(side, names - side, key=len), c.confidence))
    return out


def compatible(a, b, names):
    a2, b2 = set(names) - a, set(names) - b
    return not (a & b) or not (a & b2) or not (a2 & b) or not (a2 & b2)


def describe(side, groups, taxa):
    import collections
    c = collections.Counter(lba.taxon(n, groups, taxa) for n in side)
    return "[" + ", ".join(f"{k} {v}" for k, v in c.most_common()) + "]"


if __name__ == "__main__":
    main()
