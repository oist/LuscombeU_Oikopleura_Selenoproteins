#!/usr/bin/env python3
"""Maximum-likelihood trees of the machinery families whose orthology best hits cannot settle.

Families and outgroups:
  pstk    PSTK, rooted on KTI12 (human RefSeq, O. dioica, and KTI12 entries of the Selenoprofiles UniRef50 database)
  SEPHS2  selenophosphate synthetases: the Selenoprofiles profile names its sequences sps1, sps2 or sps after Mariotti
          et al. 2015; human SEPHS1 and SEPHS2 are added; midpoint rooting
  SBP2    SECIS-binding proteins: the profile contains SECISBP2, SECISBP2-like and non-animal sequences; human
          SECISBP2 and SECISBP2L are added; midpoint rooting
Sequences: the profile sequences, the stage 1 predictions of the other tunicates and the vertebrates, the stage 2
predictions of the appendicularians (appendicularian profiles), and the curated transcript proteins. Predictions
shorter than --min-fraction of the median profile length are left out (fragments make trees unstable).
Steps: MAFFT (L-INS-i), ClipKIT (smart-gap trimming), IQ-TREE 2 (ModelFinder restricted to --mset and --mrate, 1000
ultrafast bootstraps). Outputs: <out>/<family>.{faa,aln,trimmed,treefile,iqtree} and <figures>/phylogeny_<family>.*
"""
import argparse, glob, gzip, os, re, statistics, subprocess

APPENDICULARIANS = {"OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"}
TUNICATES = {"Cint", "Crob", "Csav", "Hror", "Scla", "Tdem"}
VERTEBRATES = {"Hsap", "Mmus", "Drer"}
HUMAN_GENES = {"pstk": ["PSTK", "KTI12"], "SEPHS2": ["SEPHS1", "SEPHS2"], "SBP2": ["SECISBP2", "SECISBP2L"]}
OUTGROUP_PATTERN = {"pstk": r"KTI12|[Kk]iller toxin.insensitive"}
COLOURS = {"appendicularian": "#c0392b", "tunicate": "#e67e22", "vertebrate": "#2471a3", "outgroup": "#7f8c8d",
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
            seqs[name] = clean(seq); groups[name] = group

        prof = list(read_fasta(f"{a.profiles}/{fam}.fa"))
        med = statistics.median(len(clean(s.replace("-", ""))) for _, s in prof)
        for n, s in prof:
            sp = re.search(r"\[([^\]]+)\]", n)
            add(f"{(sp.group(1) if sp else n.split()[0])}|{n.split()[0]}", s.replace("-", ""), "profile")
        for faa in sorted(glob.glob(f"{a.stage1}/*.selenoprofiles.faa")):
            g = os.path.basename(faa).split(".")[0]
            if g not in TUNICATES | VERTEBRATES: continue
            for n, s in read_fasta(faa):
                pid = n.split()[0]
                if pid.split(".")[0] == fam and len(clean(s)) >= a.min_fraction * med:
                    add(f"{g}|{pid}", s, "tunicate" if g in TUNICATES else "vertebrate")
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
    out = [c for c in tree.get_terminals() if groups.get(c.name) == "outgroup"]
    if out:
        try: tree.root_with_outgroup(out if len(out) > 1 and tree.is_monophyletic(out) else out[0])
        except Exception: tree.root_at_midpoint()
    else:
        tree.root_at_midpoint()
    tree.ladderize()
    for c in tree.get_nonterminals():  # IQ-TREE writes bootstraps as clade names
        if c.name and re.fullmatch(r"\d+(\.\d+)?", c.name): c.confidence, c.name = float(c.name), None
    n = len(tree.get_terminals())
    fig, ax = plt.subplots(figsize=(9, max(6, 0.16 * n)))
    Phylo.draw(tree, axes=ax, do_show=False,
               label_colors=lambda name: COLOURS.get(groups.get(name, "profile"), "#2c3e50"),
               branch_labels=lambda c: f"{c.confidence:.0f}" if c.confidence is not None and c.confidence >= 70 and
               not c.is_terminal() else "")
    for t in ax.texts: t.set_fontsize(5)
    ax.set_title(f"{fam}: maximum-likelihood tree (IQ-TREE, ultrafast bootstraps >= 70 shown); red: appendicularians, "
                 f"orange: other tunicates, blue: vertebrates, purple: curated transcripts, grey: outgroup", fontsize=7)
    ax.set_ylabel(""); ax.set_yticks([])
    fig.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(f"{figdir}/phylogeny_{fam}.{ext}", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
