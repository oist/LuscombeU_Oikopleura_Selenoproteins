#!/usr/bin/env python3
"""Appendicularian Selenoprofiles profiles: standard profiles extended with appendicularian proteins, and new profiles
for families without a standard profile.

Members of a family, in this order (every candidate is listed in members.tsv, with the reason it was kept or not):
  1. gold-standard proteins (curation/gold_standard/<genome>.proteins.faa) of --gold-genomes; the profile of a gene
     is its name up to the first "_", without dots (DIO_PAR -> DIO, TXNL1.b -> TXNL1b: Selenoprofiles reads the
     family of a prediction up to the first dot);
  2. stage 1 Selenoprofiles predictions of --stage1-genomes, for the families of --stage1-families (default: the
     machinery), labelled "homologue" or "selenocysteine" and at least --min-fraction of the median length of the
     standard profile sequences; in a genome of --gold-genomes, the predictions labelled "selenocysteine" of a family
     that has gold-standard proteins are left out (the gold standard is their curated form);
  3. stage 1 predictions labelled "selenocysteine" of --validated-genomes that step 1 kept (--validated, column
     decision = kept), for the other families, with the same length rule;
  4. curated proteins: --curated (header family|genome|source) and, from --candidates, the entries named
     <gene>_<genome> (proteins of other species; the O. dioica entries repeat the gold standard).
Proteins of --held-out genomes are never used.

Standard families: the members are added to the standard alignment with "mafft --add --keeplength", which keeps its
columns, so that the column-based options of the standard configuration stay valid; a member whose selenocysteine
(U) is not in a column where the standard alignment has a U, or is lost in an insertion, is removed (reason in
members.tsv). New families (--new-families, no standard profile): the members are aligned with MAFFT L-INS-i, and the
configuration sets p2g_refiltering to x.label != 'pseudo' instead of the default AWSI filter, which needs more
sequences. Each alignment is then built with "selenoprofiles build" (all sequences as BLAST queries), loading the
standard configuration when there is one; keys of the standard configuration that the build does not write back
(go_terms, used by the GPX filter) are appended to the new configuration. Profiles of standard families keep the
standard names, so that they replace the standard ones in a copy of the Selenoprofiles data (pipeline option
--selenoprofiles_overrides); new names must be requested in the profiles column of the samplesheet.

Outputs in --out: <family>.fa, .fa.config and .fa.profile_data (only profiles, so that the directory can be passed to
--selenoprofiles_overrides), members.tsv; alignments before the build in --out/alignments/; with --archive, a .tar.gz
of the profiles. --no-build stops after the alignments (no Selenoprofiles needed).
"""
import argparse, csv, os, re, shutil, statistics, subprocess, sys, tarfile, tempfile

MACHINERY = ["SEPHS2", "pstk", "SecS", "eEFsec", "SBP2", "secp43"]
NEW = ["KKSU", "TXNL1b"]
GOLD_GENOMES = ["OKI", "OSKA", "OdB3"]
STAGE1_GENOMES = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]
HELD_OUT = ["Bsty", "Mery", "Fbor"]          # decision of 2026-10-24: not used to train stage 2 tools
NEW_CONFIG = "p2g_refiltering = x.label!='pseudo'\n"


def read_fasta(path):
    name, buf = None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].rstrip("\n"), []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


def profile_of(gene):
    return gene.split("_")[0].replace(".", "")


def safe(name):
    return re.sub(r"[^A-Za-z0-9_.]", "_", name)


def gold_families(gold, genomes):
    fams = []
    for g in genomes:
        path = f"{gold}/{g}.proteins.faa"
        if os.path.exists(path):
            for h, _ in read_fasta(path):
                f = profile_of(h.split("|")[0])
                if f not in fams: fams.append(f)
    return fams


def u_columns(aln):
    return {i for _, s in aln for i, c in enumerate(s) if c in "Uu"}


def config_keys(path):
    """Keys of a Selenoprofiles profile configuration (lines "key = value"; the space before "=" is optional)."""
    out = {}
    for l in open(path):
        if "=" in l and l.split("=", 1)[0].strip():
            out.setdefault(l.split("=", 1)[0].strip(), l.rstrip("\n"))
    return out


def restore_config_keys(loaded, built):
    """Append to the built configuration the lines of the loaded one whose key it lacks. "selenoprofiles build"
    writes only the keys of its parameter list (filters, options, tags), so that, for example, go_terms is lost and
    a filter calling x.go_score() then fails (go_score returns None without GO terms: GPX)."""
    have, want = config_keys(built), config_keys(loaded)
    missing = [k for k in want if k not in have and k != "queries"]
    if missing:
        with open(built, "a") as o:
            for k in missing: o.write(want[k] + "\n")
    return missing


def run(cmd, **kw):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kw)
    if p.returncode:
        sys.exit(f"command failed ({p.returncode}): {' '.join(cmd)}\n{p.stderr[-2000:]}")
    return p.stdout


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--profiles", default="data/references/selenoprofiles_data/selenoprotein_profiles",
                    help="standard Selenoprofiles profiles")
    ap.add_argument("--families", nargs="+", help="profiles to build (default: the families of the gold standard, "
                    "then the machinery)")
    ap.add_argument("--new-families", nargs="*", default=NEW, help="families without a standard profile")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--gold-genomes", nargs="*", default=GOLD_GENOMES)
    ap.add_argument("--stage1", default="work/stage1/pipeline/selenoprofiles")
    ap.add_argument("--stage1-genomes", nargs="*", default=STAGE1_GENOMES)
    ap.add_argument("--stage1-families", nargs="*", default=MACHINERY)
    ap.add_argument("--validated", default="results/tables/stage2_secis_hits.tsv")
    ap.add_argument("--validated-genomes", nargs="*", default=[])
    ap.add_argument("--curated", nargs="*", default=["curation/machinery_proteins.faa"])
    ap.add_argument("--candidates", default="")
    ap.add_argument("--held-out", nargs="*", default=HELD_OUT)
    ap.add_argument("--min-fraction", type=float, default=0.5)
    ap.add_argument("--out", default="work/stage2/profiles")
    ap.add_argument("--archive", default="")
    ap.add_argument("--no-build", action="store_true")
    a = ap.parse_args()

    families = a.families or gold_families(a.gold, [g for g in a.gold_genomes if g not in a.held_out]) + \
        [f for f in MACHINERY if f not in gold_families(a.gold, a.gold_genomes)]
    out = os.path.abspath(a.out)
    os.makedirs(f"{out}/alignments", exist_ok=True)
    for f in os.listdir(out):                    # profiles of an earlier build (the directory feeds the pipeline)
        if re.search(r"\.fa(\.config|\.profile_data)?$", f): os.remove(f"{out}/{f}")

    gold = {}
    for g in a.gold_genomes:
        for h, s in read_fasta(f"{a.gold}/{g}.proteins.faa"):
            gene = h.split("|")[0]
            gold.setdefault(profile_of(gene), []).append((g, gene, s.replace("*", "")))
    validated = set()
    if a.validated_genomes:
        for r in csv.DictReader(open(a.validated), delimiter="\t"):
            if r["genome"] in a.validated_genomes and r["decision"] == "kept" and r["label"] == "selenocysteine":
                validated.add((r["genome"], r["prediction"]))
    curated = []
    for path in a.curated:
        for h, s in read_fasta(path):
            f, g, src = h.split("|")[:3]
            curated.append((f, g, src, s))
    if a.candidates:
        for h, s in read_fasta(a.candidates):
            m = re.fullmatch(r"(.+)_([A-Za-z0-9]+)", h.split()[0])
            if m and m.group(2) in STAGE1_GENOMES:
                curated.append((profile_of(m.group(1)), m.group(2), h.split(None, 1)[1] if " " in h else "", s))
    stage1 = {}
    for g in sorted(set(a.stage1_genomes) | set(a.validated_genomes), key=STAGE1_GENOMES.index):
        path = f"{a.stage1}/{g}.selenoprofiles.faa"
        if os.path.exists(path):
            stage1[g] = [(h.split()[0], s.replace("*", "")) for h, s in read_fasta(path)]
        elif not a.no_build or g in a.validated_genomes:
            print(f"warning: no stage 1 predictions for {g} ({path})", file=sys.stderr)

    log = [("family", "source", "genome", "name", "length", "kept", "reason", "residues_in_profile")]
    env = None
    if not a.no_build:
        home = tempfile.mkdtemp()                # selenoprofiles reads its configuration from $HOME
        env = dict(os.environ, HOME=home)
        subprocess.run(["selenoprofiles", "-setup"], env=env, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    for fam in families:
        new = fam in a.new_families
        std = f"{a.profiles}/{fam}.fa"
        if not new and not os.path.exists(std): sys.exit(f"{fam}: no standard profile {std}; add it to --new-families")
        prof = [] if new else list(read_fasta(std))
        median = statistics.median(len(s.replace("-", "")) for _, s in prof) if prof else 0
        members = []                             # (name in the alignment, sequence, index of its row in log)

        def note(src, g, name, s, why, aligned_name):
            log.append((fam, src, g, name, len(s), "yes" if not why else "no", why or "kept", ""))
            if not why: members.append((safe(aligned_name), s, len(log) - 1))

        for g, gene, s in gold.get(fam, []):
            note("gold", g, gene, s, "held-out genome" if g in a.held_out else "", f"{g}_{gene}")
        for g, preds in stage1.items():
            for pid, s in preds:
                if pid.split(".")[0] != fam: continue
                label = pid.rsplit(".", 1)[1]
                from_stage1 = fam in a.stage1_families and g in a.stage1_genomes
                from_valid = (g, pid) in validated and fam not in a.stage1_families
                if not (from_stage1 or from_valid): continue
                if g in a.held_out: why = "held-out genome"
                elif from_stage1 and label not in ("homologue", "selenocysteine"): why = f"label {label}"
                elif from_stage1 and label == "selenocysteine" and any(x[0] == g for x in gold.get(fam, [])):
                    why = "selenocysteine prediction of a gold-standard genome (gold protein used)"
                elif median and len(s) < a.min_fraction * median:
                    why = f"shorter than {a.min_fraction} x median ({median:.0f})"
                else: why = ""
                note("stage1" if from_stage1 else "validated", g, pid, s, why, f"{g}_{pid}")
        for f, g, src, s in curated:
            if f != fam: continue
            in_gold = any(x[0] == g for x in gold.get(fam, []))
            why = ("held-out genome" if g in a.held_out else
                   "curated protein of a gold-standard genome (gold protein used)" if in_gold else "")
            note("curated", g, src, s, why, f"{g}_curated_{len(members)}")
        if not members:
            if new: sys.exit(f"{fam}: no sequence for a new family")
            print(f"{fam}: no appendicularian sequence; standard profile kept (not written)", flush=True)
            continue

        with tempfile.TemporaryDirectory() as d:
            add = f"{d}/add.fa"
            with open(add, "w") as o:
                for n, s, _ in members: o.write(f">{n}\n{s}\n")
            if new:
                cmd = ["mafft", "--quiet", "--anysymbol", "--localpair", "--maxiterate", "1000", add]
            else:
                cmd = ["mafft", "--quiet", "--anysymbol", "--keeplength", "--add", add, std]
            aln = list(read_fasta_text(run(cmd)))
        added = {n: (s, i) for n, s, i in members}
        std_u = u_columns(prof) if prof else set()
        keep, removed = [], 0
        for n, s in aln:
            if n in added and "U" in added[n][0].upper():
                cols = {i for i, c in enumerate(s) if c in "Uu"}
                why = ("selenocysteine lost in an insertion (mafft --keeplength)" if not cols else
                       "selenocysteine outside the selenocysteine columns of the standard profile"
                       if std_u and not cols <= std_u else "")
                if why:
                    log[added[n][1]] = log[added[n][1]][:5] + ("no", why, ""); removed += 1; continue
            if n in added:                       # residues left after mafft --keeplength
                log[added[n][1]] = log[added[n][1]][:7] + (len(re.sub(r"[-.]", "", s)),)
            keep.append((n, s))
        if new and len(u_columns(keep)) > 1:
            print(f"warning: {fam}: selenocysteine in {len(u_columns(keep))} columns", file=sys.stderr)
        merged = f"{out}/alignments/{fam}.fa"
        with open(merged, "w") as o:
            for n, s in keep: o.write(f">{n}\n{s}\n")
        n_added = sum(1 for n, _ in keep if n in added)
        msg = (f"{fam}: {len(prof)} standard + {n_added} appendicularian sequences "
               f"({removed} removed by the selenocysteine check)")
        if not a.no_build:
            dest = f"{out}/{fam}.fa"
            if new:
                cfg = f"{out}/alignments/{fam}.new.config"
                with open(cfg, "w") as o: o.write(f"name = {fam}\n{NEW_CONFIG}")
            else:
                cfg = f"{std}.config"
            run(["selenoprofiles", "build", "-i", merged, "-o", dest, "-load", cfg, "-queries", "all", "-name", fam,
                 "-y"], env=env)
            restored = restore_config_keys(cfg, f"{dest}.config")
            msg += f" -> {sum(1 for _ in read_fasta(dest))} in the profile"
            if restored: msg += f"; configuration keys restored: {', '.join(restored)}"
        print(msg, flush=True)
    with open(f"{out}/members.tsv", "w") as o:
        for row in log: o.write("\t".join(map(str, row)) + "\n")
    if env: shutil.rmtree(env["HOME"], ignore_errors=True)
    if a.archive and not a.no_build:
        with tarfile.open(a.archive, "w:gz") as t:
            for f in sorted(os.listdir(out)):
                if re.search(r"\.fa(\.config|\.profile_data)?$", f): t.add(f"{out}/{f}", arcname=f)
        print(f"archive: {a.archive}")


def read_fasta_text(text):
    name, buf = None, []
    for l in text.splitlines():
        if l.startswith(">"):
            if name is not None: yield name, "".join(buf)
            name, buf = l[1:].split()[0], []
        else: buf.append(l.strip())
    if name is not None: yield name, "".join(buf)


if __name__ == "__main__":
    main()
