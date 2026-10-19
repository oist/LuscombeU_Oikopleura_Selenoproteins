#!/usr/bin/env python3
"""Maximum-likelihood trees of the machinery families whose orthology best hits cannot settle.

Families and outgroups:
  pstk    PSTK, rooted on KTI12 (human RefSeq, O. dioica, the KTI12 of the RefSeq proteomes, and KTI12 entries of the
          Selenoprofiles UniRef50 database)
  SEPHS2  selenophosphate synthetases, with the SPS1/SPS2-labelled sequences of Mariotti et al. 2015 (--references)
          and human SEPHS1 and SEPHS2; rooted on the prokaryotic SelD when they form a clade, else at the midpoint
  SBP2    SECIS-binding proteins: the profile contains SECISBP2, SECISBP2-like and non-animal sequences; human
          SECISBP2 and SECISBP2L are added; midpoint rooting
Sequences: the profile sequences (optionally at most --max-per-group per taxonomic group, the longest kept), the stage 1
predictions of the other tunicates, cephalochordates and vertebrates, the stage 2 predictions of the appendicularians
(appendicularian profiles), and the curated transcript proteins. For PSTK and the SECIS-binding proteins, which have no
selenocysteine, the RefSeq proteins of a genome (best human hit PSTK or KTI12, SECISBP2 or SECISBP2L; longest
isoform per gene) replace its predictions, which cover only the region of the profile. Predictions
shorter than --min-fraction of the median profile length are left out (fragments make trees unstable).
Steps: MAFFT (L-INS-i), ClipKIT (smart-gap trimming), IQ-TREE 2 (ModelFinder restricted to --mset and --mrate, 1000
ultrafast bootstraps). Outputs: <out>/<family>.{faa,aln,trimmed,treefile,iqtree} and <figures>/phylogeny_<family>.*
"""
import argparse, glob, gzip, os, re, statistics, subprocess

APPENDICULARIANS = {"OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"}
TUNICATES = {"Cint", "Crob", "Csav", "Hror", "Scla", "Tdem", "Bsch", "Pmam", "Pege", "Stho"}
VERTEBRATES = {"Hsap", "Mmus", "Drer"}
CEPHALOCHORDATES = {"Bflo", "Blan"}
# families without selenocysteine: their RefSeq proteins, when a genome has them, replace the Selenoprofiles
# predictions, which cover only the region of the profile (for example the C-terminal domain of SECIS-binding proteins)
REFSEQ_FAMILIES = {"pstk": ["PSTK", "KTI12"], "SBP2": ["SECISBP2", "SECISBP2L"]}
HUMAN_GENES = {"pstk": ["PSTK", "KTI12"], "SEPHS2": ["SEPHS1", "SEPHS2"], "SBP2": ["SECISBP2", "SECISBP2L"]}
OUTGROUP_PATTERN = {"pstk": r"KTI12|[Kk]iller toxin.insensitive"}
COLOURS = {"appendicularian": "#c0392b", "tunicate": "#e67e22", "vertebrate": "#2471a3", "outgroup": "#7f8c8d",
           "cephalochordate": "#117a65", "reference": "#2c3e50",
           "curated": "#8e44ad", "profile": "#2c3e50"}


def read_fasta(path):
    name, buf = None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for l in fh:
            if l.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = l[1:].rstrip("\n"), []
            else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def clean(s): return re.sub(r"[^A-Za-z]", "", s).upper().replace("U", "X")


def label(s): return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")[:60]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", nargs="+", default=["pstk", "SEPHS2", "SBP2"])
    ap.add_argument("--profiles", default="data/references/selenoprofiles_data/selenoprotein_profiles")
    ap.add_argument("--uniref", default="data/references/selenoprofiles_data/uniref50_seleno_sim.fasta")
    ap.add_argument("--stage1", default="work/stage1/pipeline/selenoprofiles")
    ap.add_argument("--stage2", default="work/stage2/machinery_test")
    ap.add_argument("--curated", default="curation/machinery_proteins.faa")
    ap.add_argument("--references", default="curation/mariotti2015_sps.faa",
                    help="labelled reference proteins, header <family>|<label>|<source>")
    ap.add_argument("--genomes", default="data/genomes", help="RefSeq proteins and annotations of the genomes")
    ap.add_argument("--max-per-group", type=int, default=0,
                    help="keep at most this many profile sequences per taxonomic group (config/taxon_groups.tsv); 0: all")
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--collect-only", action="store_true", help="write the sequence sets and stop (no alignment, no tree)")
    ap.add_argument("--outgroups", default="curation/machinery_outgroups.faa")
    ap.add_argument("--human", default="data/genomes/Hsap")
    ap.add_argument("--min-fraction", type=float, default=0.5)
    ap.add_argument("--max-uniref-outgroups", type=int, default=12)
    ap.add_argument("--clipkit-mode", default="smart-gap", help="ClipKIT trimming mode (for example gappy)")
    ap.add_argument("--clipkit-gaps", default="", help="ClipKIT gap threshold, for the gappy modes (for example 0.5)")
    ap.add_argument("--model", default="", help="fixed model instead of ModelFinder (for example LG+C20+F+G)")
    ap.add_argument("--mset", default="LG,WAG,JTT", help="substitution matrices tested by ModelFinder")
    ap.add_argument("--mrate", default="", help="rate models tested by ModelFinder (for example G4,R4,R5; default: all)")
    ap.add_argument("--cpus", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)))
    ap.add_argument("--out", default="results/trees")
    ap.add_argument("--figures", default="figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    gene_of = {}
    with gzip.open(f"{a.human}/annotation.gff3.gz", "rt") as fh:
        for l in fh:
            if "\tCDS\t" not in l: continue
            m, g = re.search(r"protein_id=([^;\n]+)", l), re.search(r"gene=([^;\n]+)", l)
            if m and g: gene_of[m.group(1)] = g.group(1)
    human = {n.split()[0]: s for n, s in read_fasta(f"{a.human}/protein.faa.gz")}
    for fam in a.families:
        seqs, groups = {}, {}

        def add(name, seq, group):
            name = label(name)
            if name in seqs: return
            if len(seq) and sum(seq.upper().count(c) for c in "ACGTN-") / len(seq) > 0.9:
                print(f"WARNING: {name} looks like a nucleotide sequence; left out", flush=True); return
            seqs[name] = clean(seq); groups[name] = group

        prof = list(read_fasta(f"{a.profiles}/{fam}.fa"))
        med = statistics.median(len(clean(s.replace("-", ""))) for _, s in prof)
        taxa = load_taxa(a.taxa)
        per_group = {}
        for n, s in sorted(prof, key=lambda x: -len(clean(x[1].replace("-", "")))):   # longest first
            sp = re.search(r"\[([^\]]+)\]", n)
            name = f"{(sp.group(1) if sp else n.split()[0])}|{n.split()[0]}"
            grp = taxa.get(re.split(r"[ _|]", re.sub(r"^(sps1|sps2|sps)_", "", name))[0], "unclassified")
            if a.max_per_group and grp != "unclassified":
                if per_group.get(grp, 0) >= a.max_per_group: continue
                per_group[grp] = per_group.get(grp, 0) + 1
            add(name, s.replace("-", ""), "profile")
        if os.path.exists(a.references):
            per_type = {}
            for n, s in sorted(read_fasta(a.references), key=lambda x: -len(x[1])):
                f = n.split("|")
                if f[0] != fam: continue
                key = (taxa.get(re.split(r"[ _|.]", f[1])[0], "unclassified"), f[1].rsplit(".", 1)[-1])
                if a.max_per_group and key[0] != "unclassified":
                    if per_type.get(key, 0) >= a.max_per_group: continue
                    per_type[key] = per_type.get(key, 0) + 1
                add(f[1], s, "reference")
        kind = lambda g: "tunicate" if g in TUNICATES else "cephalochordate" if g in CEPHALOCHORDATES else "vertebrate"
        replaced = set()
        if fam in REFSEQ_FAMILIES:
            for g in sorted((TUNICATES | VERTEBRATES | CEPHALOCHORDATES) - {"Hsap"}):
                prots = refseq_family_proteins(f"{a.genomes}/{g}", REFSEQ_FAMILIES[fam], gene_of, human, a.cpus)
                for pid, (gene, seq) in prots.items():
                    add(f"{g}|RefSeq_{gene}_{pid}", seq, "outgroup" if gene == "KTI12" else kind(g))
                if prots: replaced.add(g)
            if replaced: print(fam, "RefSeq proteins instead of predictions for", " ".join(sorted(replaced)), flush=True)
        for faa in sorted(glob.glob(f"{a.stage1}/*.selenoprofiles.faa")):
            g = os.path.basename(faa).split(".")[0]
            if g not in TUNICATES | VERTEBRATES | CEPHALOCHORDATES or g in replaced: continue
            for n, s in read_fasta(faa):
                pid = n.split()[0]
                if pid.split(".")[0] == fam and len(clean(s)) >= a.min_fraction * med:
                    add(f"{g}|{pid}", s, kind(g))
        for faa in sorted(glob.glob(f"{a.stage2}/*.faa")):
            g = os.path.basename(faa).split(".")[0]
            for n, s in read_fasta(faa):
                pid = n.split()[0]
                if pid.split(".")[0].removesuffix("_app") == fam and len(clean(s)) >= a.min_fraction * med:
                    add(f"{g}|{pid.replace('_app', '')}", s, "appendicularian")
        if os.path.exists(a.curated):
            for n, s in read_fasta(a.curated):
                f = n.split("|")
                if f[0] == fam: add(f"{f[1]}|curated_{f[0]}", s, "curated")
        for gene in HUMAN_GENES.get(fam, []):
            prots = [p for p, gg in gene_of.items() if gg == gene and p in human]
            if prots:
                p = max(prots, key=lambda x: (x.startswith("NP_"), len(human[x])))
                add(f"Human_{gene}|{p}", human[p], "outgroup" if gene == "KTI12" else "vertebrate")
        if fam in OUTGROUP_PATTERN:
            if os.path.exists(a.outgroups):
                for n, s in read_fasta(a.outgroups):
                    f = n.split("|")
                    if f[0] == fam: add(f"KTI12|{f[1]}", s, "outgroup")
            k = 0
            for n, s in read_fasta(a.uniref):
                if k >= a.max_uniref_outgroups: break
                if re.search(OUTGROUP_PATTERN[fam], n):
                    sp = re.search(r"Tax=([^=]+?) TaxID", n)
                    add(f"KTI12|{(sp.group(1) if sp else '')}|{n.split()[0]}", s, "outgroup"); k += 1
        base = f"{a.out}/{fam}"
        with open(f"{base}.faa", "w") as o:
            for n, s in seqs.items(): o.write(f">{n}\n{s}\n")
        with open(f"{base}.groups.tsv", "w") as o:
            for n, g in groups.items(): o.write(f"{n}\t{g}\n")
        print(fam, len(seqs), "sequences", flush=True)
        if a.collect_only: continue
        if not os.path.exists(f"{base}.trimmed"):
          with open(f"{base}.aln", "w") as o:
            subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", str(a.cpus), "--quiet",
                            f"{base}.faa"], stdout=o, check=True)
          subprocess.run(["clipkit", f"{base}.aln", "-o", f"{base}.trimmed", "-m", a.clipkit_mode] +
                         (["-g", a.clipkit_gaps] if a.clipkit_gaps else []), check=True, stdout=subprocess.DEVNULL)
        if not os.path.exists(f"{base}.iqtree"):  # resumes from the IQ-TREE checkpoint after an interruption
            models = ["-m", a.model] if a.model else ["-m", "MFP", "-mset", a.mset] + (["-mrate", a.mrate] if a.mrate else [])
            subprocess.run(["iqtree2", "-s", f"{base}.trimmed"] + models +
                           ["-B", "1000", "-T", str(a.cpus), "--prefix", base, "--quiet", "--seed", "1"], check=True)
        draw(f"{base}.treefile", groups, fam, a.figures)
        for ext in ("bionj", "ckp.gz", "mldist", "model.gz", "splits.nex", "uniqueseq.phy", "contree"):
            if os.path.exists(f"{base}.{ext}"): os.remove(f"{base}.{ext}")


def draw(treefile, groups, fam, figdir):
    from Bio import Phylo
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tree = Phylo.read(treefile, "newick")
    for c in tree.get_nonterminals():  # IQ-TREE writes bootstraps as clade names
        if c.name and re.fullmatch(r"\d+(\.\d+)?", c.name): c.confidence, c.name = float(c.name), None
    out = [c for c in tree.get_terminals() if groups.get(c.name) == "outgroup"]
    taxa = load_taxa("config/taxon_groups.tsv")
    prok = [c for c in tree.get_terminals() if groups.get(c.name) == "profile" and
            taxa.get(re.sub(r"^(sps1|sps2|sps)_", "", c.name).split("_")[0]) in ("Bacteria", "Archaea")]

    def reroot(t):
        try:
            if out: t.root_with_outgroup(out if len(out) > 1 and t.is_monophyletic(out) else out[0])
            elif prok and t.is_monophyletic(prok): t.root_with_outgroup(prok)
            else: t.root_at_midpoint()
        except Exception:
            t.root_at_midpoint()
    reroot_keeping_supports(tree, reroot)
    tree.ladderize()
    n = len(tree.get_terminals())
    fig, ax = plt.subplots(figsize=(10, max(6, 0.16 * n)))
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: c.name[:45] if c.is_terminal() else None,
               label_colors=lambda lab: COLOURS.get(next((groups[k] for k in groups if k[:45] == lab), "profile"), "#2c3e50"),
               branch_labels=lambda c: f"{c.confidence:.0f}" if c.confidence is not None and c.confidence >= 70 and
               not c.is_terminal() else "")
    for t in ax.texts: t.set_fontsize(5)
    length_column(ax, tree, sequence_lengths(treefile), fontsize=5)
    ax.set_title(f"{fam}: maximum-likelihood tree (IQ-TREE, ultrafast bootstraps >= 70 shown); red: appendicularians, "
                 f"orange: other tunicates, blue: vertebrates, purple: curated transcripts, grey: outgroup", fontsize=7)
    ax.set_ylabel(""); ax.set_yticks([])
    fig.tight_layout(rect=(0, 0, 0.9, 1))
    for ext in ("png", "pdf"): fig.savefig(f"{figdir}/phylogeny_{fam}.{ext}", dpi=200)
    plt.close(fig)


def load_taxa(path):
    """genus (or genome code) -> taxonomic group, from config/taxon_groups.tsv"""
    taxa = {}
    if os.path.exists(path):
        for l in open(path):
            if l.startswith("#") or l.startswith("genus\t") or "\t" not in l: continue
            g, grp = l.rstrip("\n").split("\t"); taxa[g] = grp
    return taxa


def refseq_family_proteins(genome_dir, genes, gene_of, human, cpus):
    """Proteins of a RefSeq proteome whose best human hit is one of genes (the family gene and its paralog or
    related family), longest isoform per gene; {protein_id: (best human gene, sequence)}. Empty if the genome has no
    RefSeq proteins."""
    import shutil, tempfile
    faa, gff = f"{genome_dir}/protein.faa.gz", f"{genome_dir}/annotation.gff3.gz"
    if not (os.path.exists(faa) and os.path.exists(gff)): return {}
    prot_gene = {}
    with gzip.open(gff, "rt") as fh:
        for l in fh:
            if "\tCDS\t" not in l: continue
            m, g = re.search(r"protein_id=([^;\n]+)", l), re.search(r"gene=([^;\n]+)", l)
            if m: prot_gene[m.group(1)] = g.group(1) if g else m.group(1)
    proteome = {n.split()[0]: s for n, s in read_fasta(faa)}
    d = tempfile.mkdtemp()
    try:
        with open(f"{d}/p.faa", "w") as o:
            for n, s in proteome.items(): o.write(f">{n}\n{s}\n")
        with open(f"{d}/q.faa", "w") as o:
            for p, gg in gene_of.items():
                if gg in genes and p in human: o.write(f">{p}\n{human[p]}\n")
        with open(f"{d}/h.faa", "w") as o:
            for n, s in human.items(): o.write(f">{n}\n{s}\n")
        for db in ("p.faa", "h.faa"):
            subprocess.run(["formatdb", "-i", db, "-p", "T"], cwd=d, check=True)
        cand = set()
        out = subprocess.run(["blastall", "-p", "blastp", "-d", "p.faa", "-i", "q.faa", "-e", "1e-5", "-m", "8",
                              "-a", str(cpus)], cwd=d, capture_output=True, text=True, check=True).stdout
        for l in out.splitlines(): cand.add(l.split("\t")[1])
        if not cand: return {}
        with open(f"{d}/c.faa", "w") as o:
            for c in cand: o.write(f">{c}\n{proteome[c]}\n")
        best = {}
        out = subprocess.run(["blastall", "-p", "blastp", "-d", "h.faa", "-i", "c.faa", "-e", "1e-5", "-m", "8",
                              "-a", str(cpus)], cwd=d, capture_output=True, text=True, check=True).stdout
        for l in out.splitlines():
            f = l.split("\t"); bits = float(f[11])
            if bits > best.get(f[0], (0, ""))[0]: best[f[0]] = (bits, gene_of.get(f[1], f[1]))
    finally:
        shutil.rmtree(d)
    chosen = {}
    for c, (bits, hg) in best.items():
        if hg not in genes: continue
        g = prot_gene.get(c, c)
        if g not in chosen or len(proteome[c]) > len(proteome[chosen[g]]): chosen[g] = c
    return {c: (best[c][1], proteome[c]) for c in chosen.values()}


def reroot_keeping_supports(tree, reroot):
    """Reroot with reroot(tree), then put each support back on its bipartition: Bio.Phylo keeps supports on node
    objects, so that rerooting moves them by one branch along the path between the old and the new root."""
    names = frozenset(x.name for x in tree.get_terminals())

    def key(c):
        leaves = frozenset(x.name for x in c.get_terminals())
        return min(leaves, names - leaves, key=lambda s: (len(s), sorted(s)))
    before = {key(c): c.confidence for c in tree.get_nonterminals() if c is not tree.root}
    reroot(tree)
    for c in tree.get_nonterminals():
        c.confidence = None if c is tree.root else before.get(key(c))


def sequence_lengths(treefile):
    """Unaligned lengths of the sequences of a tree, from <family>.faa next to the tree or in the parent directory."""
    import os
    base = treefile[:-len(".treefile")]
    for faa in (f"{base}.faa", os.path.join(os.path.dirname(os.path.dirname(base)), os.path.basename(base) + ".faa")):
        if os.path.exists(faa):
            lengths, name = {}, None
            for l in open(faa):
                if l.startswith(">"): name = l[1:].strip(); lengths[name] = 0
                else: lengths[name] += len(l.strip().replace("-", ""))
            return lengths
    return {}


def length_column(ax, tree, lengths, members=None, fontsize=6):
    """Protein length to the right of each tip; for a collapsed tip, the median and the range of its sequences."""
    import statistics
    from matplotlib.transforms import blended_transform_factory
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(1.01, 0.2, "length (aa)", transform=tr, fontsize=fontsize, fontweight="bold", va="center")
    for i, tip in enumerate(tree.get_terminals()):
        names = (members or {}).get(tip.name, [tip.name])
        ls = [lengths[n] for n in names if n in lengths]
        if not ls: continue
        text = f"{ls[0]}" if len(ls) == 1 else f"{statistics.median(ls):.0f} ({min(ls)}-{max(ls)})"
        ax.text(1.01, i + 1, text, transform=tr, fontsize=fontsize, va="center")


if __name__ == "__main__":
    main()
