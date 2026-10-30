#!/usr/bin/env python3
"""Census of the selenoprotein genes of the appendicularians (step 9), from the stage 2 final run, with the
contamination screen and the list of candidates for the strict screen of horizontal gene transfer.

Inputs (tables of the repository; Deigo steps in the Makefile: stage2-collect, census-contigs):
  --sites        selenocysteine codons of the final run (stage1/collect.py --prefix stage2_final)
  --predictions  Selenoprofiles predictions of the final run, all labels (same)
  --contigs      contig checks (stage2/census_contigs.py): composition robust z-scores, Tiara classes
  --rbh          best human hits of the final-run predictions (make selenoprotein-rbh)
  --trees        round 2 trees (make selenoprotein-phylogeny, selenoprotein-pmsf): ML and PMSF trees, groups,
                 removed.tsv (predictions replaced by gold-standard proteins, identical sequences)
  --gold         gold-standard gene models of OKI2018, OSKA2016 and OdB3 (curation/gold_standard/<genome>.gff3)
  --stage1       selenocysteine codons of stage 1 (other tunicates, cephalochordates, vertebrates), for comparison

Loci. The selenocysteine codons of a genome (Selenoprofiles predictions labelled selenocysteine, and miniprot calls;
the sites that the pipeline reports) are merged into loci when they lie on the same contig and strand within
MERGE_NT nucleotides, or belong to the same Selenoprofiles prediction (miniprot calls a few nucleotides from a codon
found in the same gene: stage 2 benchmark). Each locus takes the family, label and SECIS score of its Selenoprofiles
prediction when it has one, else those of its miniprot call.

Decision for each locus (rules fixed before the census; contamination is the default explanation of foreign-looking
sequences, but a locus is left out only on evidence against it):
  counted (gold standard)  the locus lies in a gold-standard gene of OKI2018, OSKA2016 or OdB3, or in a gold-standard
                           model lifted onto KUM-M3-7f, AOM-5-5f or Bar2_p4 (always counted);
  left out                 Tiara classifies its contig as bacteria, archaea, prokarya or organelle; or the best human
                           hit of its prediction lies in another family (RELATED genes count as the family: MIEN1 for
                           SELENOW, SEPHS1 for SEPHS2); or its contig is a composition outlier (robust z-score of the
                           tetranucleotide distance above COMPOSITION_Z, on the largest window or the 2-kb window) or
                           ORF-dense (at least ORF_DENSE of the contig in open reading frames of 100 codons or more,
                           and ORF_DENSE_MARGIN above the median of the 2-kb windows of its genome, at GC up to
                           ORF_DENSE_MAX_GC, contig of ORF_DENSE_MIN_LENGTH nt or more), and its nearest
                           relatives in the ML tree of its family are all non-animal;
  counted                  otherwise.
Flags reported with each locus (not decisive alone): composition outlier; ORF-dense contig; GC outlier (robust z-score beyond GC_Z);
nearest relatives in the ML and PMSF trees (prokaryotic, bacterial or fungal, other non-animal, animal, not in the
tree); no SECIS hit (the model of the final run, the Oikopleura model at 20 bits); fewer than END_NT nucleotides
after the codon on the contig (no room for a SECIS element in the search window); miniprot only; contig shorter than
3,000 nt (no Tiara class).

Strict HGT screen (candidates only; the tests come after): every appendicularian leaf of the round 2 ML trees whose
nearest relatives outside the appendicularians are all bacterial, archaeal or fungal (a transfer after the split
from the other tunicates, or contamination), or whose appendicularian clade lies in a clade of tunicate sequences
whose nearest relatives outside the tunicates are all bacterial, archaeal or fungal (a transfer before that split,
like that of cellulose synthase; the slower ascidian sequences then share the placement, which long-branch
attraction of the appendicularian sequences alone would not explain), whatever its label. An old transfer is not
expected to differ from its genome in composition (amelioration), nor to sit on a contig that Tiara calls
prokaryotic; the evidence for it is its presence at the same place of the tree in several species and assemblies.
A placement means something only for sequences of the family of the tree: a profile also finds homologs of other
families (thioredoxins with the DIO profile, for example), whose nearest relatives in the tree are then the few
sequences of their own family that the profile holds, often bacterial, because the animal members of that family are
not in the tree. Status: "contamination (Tiara ...)"; "contamination likely" (composition outlier, or ORF-dense contig); "other family"
(best human hit in another family: to be placed in a tree of that family); "family not confirmed" (no hit in the
human proteome although the family has human members); "placement not stable (PMSF)"
when the PMSF tree does not give the same kind of nearest relatives; "candidate, no intron" (bacterial genes have
no intron: such a candidate needs presence in other assemblies or species, or spliced transcripts); "candidate"
otherwise, to be tested (AU tests against vertical descent, other assemblies or individuals, transcripts).

Gene of the clade: for a locus in a gold-standard gene, that gene; otherwise the gold-standard genes of the smallest
clade of appendicularian sequences of the ML tree that holds the prediction and a gold-standard protein (several
when that clade holds several genes); for the families without tree of one O. dioica gene (SEPHS2, TXNL1.b, KKSU),
that gene.

Outputs: --loci (one row per locus), --table and --markdown (counted loci per family and genome, appendicularians
from this census, other genomes from stage 1 with the same locus definition and no screen), --genes-table and
--genes-markdown (counted loci per gene of the clade and appendicularian genome; in brackets, those with a SECIS
hit; "+n", loci found only by the screens of step 4), --hgt.

Loci of the step 4 screens (--screen, --scan) that the final run did not predict are added as "counted (screen)":
comparative-screen TGAs classified as homologs of a known selenoprotein, candidates whose partner in another genome
leads to a gene, and SECIS-scan hits classified as copies of a known selenoprotein; their gene comes from their best
gold-standard hit, else from their partner. They have no gene model and stay out of the family table.
"""
import argparse, collections, csv, glob, os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import machinery_lba as lba                       # tree reading, taxonomic groups
from selenoprotein_orthology import sides         # bipartition sides with supports
from taxonomy import ANIMALS

APPENDICULARIANS = ["OKI", "OSKA", "OdB3", "Olon", "Oalb", "Ovan", "Bsty", "Mery", "Fbor"]
ODIOICA_TESTS = ["KUM", "AOM", "Bar2"]
OTHERS = ["Cint", "Crob", "Csav", "Hror", "Scla", "Bsch", "Pmam", "Tdem", "Pege", "Stho", "Bflo", "Blan", "Hsap",
          "Mmus", "Drer"]
NAMES = {"OKI": "*O. dioica* OKI2018", "OSKA": "*O. dioica* OSKA2016", "OdB3": "*O. dioica* OdB3",
         "KUM": "*O. dioica* KUM-M3-7f", "AOM": "*O. dioica* AOM-5-5f", "Bar2": "*O. dioica* Bar2_p4",
         "Olon": "*O. longicauda*", "Oalb": "*O. albicans*", "Ovan": "*O. vanhoeffeni*", "Bsty": "*B. stygius*",
         "Mery": "*M. erythrocephalus*", "Fbor": "*F. borealis*", "Cint": "*C. intestinalis*", "Crob": "*C. robusta*",
         "Csav": "*C. savignyi*", "Hror": "*H. roretzi*", "Scla": "*S. clava*", "Bsch": "*B. schlosseri*",
         "Pmam": "*P. mammillata*", "Tdem": "*T. democratica*", "Pege": "*Pegea* sp.", "Stho": "*S. thompsoni*",
         "Bflo": "*B. floridae*", "Blan": "*B. lanceolatum*", "Hsap": "human", "Mmus": "mouse", "Drer": "zebrafish"}
MERGE_NT = 30
END_NT = 500
COMPOSITION_Z = 5.0
ORF_DENSE = 0.6        # fraction of a contig in long open reading frames (stage2/census_contigs.py) ...
ORF_DENSE_MAX_GC = 0.5 # ... counted only at GC up to 50%, where random sequence gives at most about 0.3 ...
ORF_DENSE_MIN_LENGTH = 5000
ORF_DENSE_MARGIN = 0.3 # ... and at least this much above the median of the 2-kb windows of its genome (added
                       # 2026-10-30 after the second contig run: the median of F. borealis is 0.80)
GENOME_ORF = {}        # genome -> median open reading frame density of its 2-kb windows (read in main)
GC_Z = 5.0
RELATED = {"SELENOW": {"MIEN1"}, "SEPHS2": {"SEPHS1"}}
SINGLE_GENE = {"KKSU": "KKSU", "TXNL1b": "TXNL1.b", "SEPHS2": "SEPHS2"}   # families of one O. dioica gene, no tree
GENE_ORDER = ["DIO_PAR", "DIO_chr1", "GPX_PAR", "GPX_chr1a", "GPX_chr1b", "GPX_chr2", "SELENOF", "SELENOK", "SELENOO",
              "SELENOS", "SELENOT", "SELENOW_chr1a", "SELENOW_pair1", "SELENOW_pair2", "SELENOW_PAR", "SELENOW_chr2",
              "SEPHS2", "TXNRD", "TXNL1.b", "KKSU"]
NONEUKARYOTIC = {"bacteria", "archaea", "prokarya", "organelle", "mitochondrion", "plastid"}
PROKARYOTES = {"Bacteria", "Archaea"}


def rows(path):
    if not path or not os.path.exists(path): return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def is_sec(r):
    return r["selenoprofiles_label"] == "selenocysteine" or r["sources"] in ("miniprot", "both")


def loci(sites):
    """Merge the selenocysteine codons of one genome into loci (lists of site rows)."""
    by = collections.defaultdict(list)
    for r in sites:
        if is_sec(r): by[(r["contig"], r["strand"])].append(r)
    out = []
    for key in sorted(by):
        cur = []
        for r in sorted(by[key], key=lambda r: int(r["sec_codon_start"])):
            pid = r.get("selenoprofiles_id", "") or "."
            if cur and (int(r["sec_codon_start"]) - int(cur[-1]["sec_codon_end"]) <= MERGE_NT or
                        (pid != "." and pid in {x.get("selenoprofiles_id") for x in cur})):
                cur.append(r)
            else:
                if cur: out.append(cur)
                cur = [r]
        if cur: out.append(cur)
    return out


def gold_genes(gold, genomes, truth=None):
    """genome -> [(contig, strand, start, end, gene)] from the gene features of the gold-standard GFF3 files, and,
    for the new O. dioica assemblies (ODIOICA_TESTS), from the gold-standard models lifted onto them
    (results/tables/stage2_truth_sets.tsv; span of the lifted coding sequence, or the codon when there is none)."""
    out = collections.defaultdict(list)
    for r in rows(truth):
        if r["genome"] not in ODIOICA_TESTS or not r["contig"]: continue
        pos = [int(x) for x in re.findall(r"\d+", r["cds"])] or ([int(r["sec_start"])] if r["sec_start"] else [])
        if pos: out[r["genome"]].append((r["contig"], r["strand"], min(pos), max(pos) + 2, r["gene"]))
    for g in genomes:
        path = f"{gold}/{g}.gff3"
        if not os.path.exists(path): continue
        for l in open(path):
            f = l.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene": continue
            m = re.search(r"Name=([^;]+)", f[8])
            out[g].append((f[0], f[6], int(f[3]), int(f[4]), m.group(1) if m else f[8]))
    return out


def kind_of(others):
    kinds = set(others) - {"unclassified"}
    if kinds and kinds <= PROKARYOTES: return "prokaryotic"
    if kinds and kinds <= PROKARYOTES | {"Fungi"}: return "bacterial or fungal"
    if kinds and not kinds & ANIMALS: return "non-animal"
    return "animal" if kinds else "unclassified"


def placements(trees, taxa):
    """(genome, prediction) -> {tree: (kind, groups, support, family)} for the appendicularian leaves of the round 2
    trees, at two levels: the nearest relatives outside the appendicularians (keys ML, PMSF) and outside the
    tunicates (keys ML_tunicates, PMSF_tunicates: the smallest side holding the leaf and sequences other than
    appendicularian or other tunicate ones, with the number of other tunicate sequences in it); kind: prokaryotic,
    bacterial or fungal, non-animal, animal. Key "genomes": the genomes of the largest clade of appendicularian
    sequences holding the leaf (ML tree). Gold-standard leaves stand for the predictions they replaced, and kept
    sequences for the identical ones removed."""
    out = collections.defaultdict(dict)
    for tf in sorted(glob.glob(f"{trees}/*.treefile")):
        fam = os.path.basename(tf)[: -len(".treefile")]
        gpath = f"{trees}/{fam}.groups.tsv"
        if not os.path.exists(gpath): continue
        groups = dict(l.rstrip("\n").split("\t") for l in open(gpath))
        stand_for = collections.defaultdict(list)
        for r in rows(f"{trees}/{fam}.removed.tsv"):
            m = re.fullmatch(r"gold-standard gene (\S+) used", r["reason"])
            if m: stand_for[f"{r['sequence'].partition('_')[0]}_gold_{m.group(1)}"].append(r["sequence"])
            m = re.fullmatch(r"identical to (\S+)", r["reason"])
            if m: stand_for[m.group(1)].append(r["sequence"])
        for label, path in (("ML", tf), ("PMSF", f"{trees}/lba/{fam}.pmsf.treefile")):
            if not os.path.exists(path): continue
            names, sd = sides(lba.read_tree(path))
            tax = {n: lba.taxon(n, groups, taxa) for n in names}
            tunicates = {"Appendicularians", "Other tunicates", "unclassified"}
            for n in names:
                if tax[n] != "Appendicularians": continue
                res = {}
                for key, inside in ((label, {"Appendicularians"}), (f"{label}_tunicates", tunicates)):
                    cand = [s for s in sd if n in s and any(tax[x] not in inside for x in s)]
                    if not cand: continue
                    near = min(cand, key=len)
                    others = collections.Counter(tax[x] for x in near if tax[x] not in inside)
                    support = sd.get(near)
                    res[key] = (kind_of(others), ", ".join(f"{k} {v}" for k, v in others.most_common()),
                                "NA" if support is None else f"{support:.0f}", fam,
                                sum(tax[x] == "Other tunicates" for x in near),
                                sorted(x for x in near if tax[x] not in inside))
                if label == "ML":
                    own = [s for s in sd if n in s and all(tax[x] == "Appendicularians" for x in s)]
                    clade = max(own, key=len) if own else {n}
                    res["genomes"] = sorted({x.partition("_")[0] for x in clade})
                    # gold-standard genes of the smallest appendicularian clade holding the leaf and a gold leaf
                    gold = [s for s in own + [frozenset([n])] if any("_gold_" in x for x in s)]
                    res["gold_genes"] = sorted({x.split("_gold_", 1)[1] for x in min(gold, key=len) if "_gold_" in x}) \
                        if gold else []
                for leaf in [n] + stand_for.get(n, []):
                    g, _, pid = leaf.partition("_")
                    out[(g, pid)].update(res)
    return out


def orf_dense(ctg):
    """Contig mostly made of long open reading frames, as bacterial DNA is, and denser than its genome (see the
    constants)."""
    orf, gc, length = fnum(ctg.get("contig_orf_fraction")), fnum(ctg.get("gc")), fnum(ctg.get("length"))
    base = GENOME_ORF.get(ctg.get("genome"))
    return (orf is not None and gc is not None and length is not None and base is not None and
            orf >= max(ORF_DENSE, base + ORF_DENSE_MARGIN) and gc <= ORF_DENSE_MAX_GC and length >= ORF_DENSE_MIN_LENGTH)


def composition_z(ctg):
    """Largest composition robust z-score of the contig (window of the largest reference length, and 2-kb window)."""
    zs = [fnum(ctg.get(k)) for k in ("composition_robust_z", "composition_robust_z_2kb")]
    zs = [z for z in zs if z is not None]
    return max(zs) if zs else None


SCREEN_STATUSES = ("homolog of a known selenoprotein", "candidate")


def parse_partner(text):
    """'OKI NC_137294.1:7274834 (+), E=...' -> (genome, contig, position, strand), or None."""
    m = re.match(r"(\S+) (\S+):(\d+) \(([+-])\)", text or "")
    return (m.group(1), m.group(2), int(m.group(3)), m.group(4)) if m else None


def screen_loci(screen, scan, genomes):
    """Selenoprotein loci found by the UGA/SECIS screens of step 4 (stage2/uga_secis_screen.py): comparative-screen
    TGAs classified as homologs of a known selenoprotein, or candidates (kept when their partner in another genome
    leads to a gene), and SECIS-scan hits classified as copies of a known selenoprotein. Returns dictionaries with
    genome, contig, strand, position (TGA, or SECIS hit), gene (from the best gold-standard hit, else to resolve
    through the partner), partner, evidence, SECIS bits."""
    out = []
    for r in rows(screen):
        if r["screen"] != "comparative" or r["genome"] not in genomes or r["status"] not in SCREEN_STATUSES: continue
        if not r.get("genomic_tga", "").isdigit(): continue
        m = re.match(r"\S+?\.(\S+) \(\d+%, U aligned to U\)", r.get("gold_hit", ""))
        out.append({"genome": r["genome"], "contig": r["contig"], "strand": r["strand"], "pos": int(r["genomic_tga"]),
                    "gene": m.group(1) if m else "", "partner": parse_partner(r.get("partner")), "kind": "TGA",
                    "evidence": f"UGA/SECIS screen ({r['id']}, {r['status']})", "secis": r["bits"],
                    "candidate": r["status"] == "candidate"})
    for r in rows(scan):
        if r.get("final_class") != "copy of a known selenoprotein" or r["genome"] not in genomes: continue
        m = re.match(r"blastx: (\S+) \(", r.get("evidence", "")) or re.search(r"gold standard (\S+) in", r.get("evidence", ""))
        out.append({"genome": r["genome"], "contig": r["contig"], "strand": r["strand"], "pos": int(r["from"]),
                    "gene": m.group(1) if m else "", "partner": None, "kind": "SECIS",
                    "evidence": f"SECIS scan ({r['evidence']})", "secis": r["bits"], "candidate": False})
    return out


def near(a_, genome, contig, strand, pos, kind):
    """Census or screen loci of the same contig and strand that hold the same gene: a TGA within 300 nt, or a SECIS hit
    within 5,000 nt downstream of the TGA (on the strand of the gene)."""
    out = []
    for r in a_:
        if (r["genome"], r["contig"], r["strand"]) != (genome, contig, strand): continue
        tga = int(r["sec_codon_start"]) if "sec_codon_start" in r else r["pos"]
        rkind = "TGA" if "sec_codon_start" in r else r["kind"]
        if kind == rkind == "TGA" or kind == rkind == "SECIS": ok = abs(tga - pos) <= 300 if kind == "TGA" else abs(tga - pos) <= 3000
        else:
            t, h = (tga, pos) if rkind == "TGA" else (pos, tga)
            ok = 0 <= (h - t if strand == "+" else t - h) <= 5000
        if ok: out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", default="results/tables/stage2_final_sec_sites.tsv")
    ap.add_argument("--predictions", default="results/tables/stage2_final_predictions.tsv")
    ap.add_argument("--contigs", default="results/tables/stage2_census_contigs.tsv")
    ap.add_argument("--genome-table", default="results/tables/stage2_census_genomes.tsv")
    ap.add_argument("--rbh", default="results/tables/stage2_selenoprotein_rbh.tsv")
    ap.add_argument("--trees", default="results/trees/selenoproteins")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--stage1", default="results/tables/stage1_sec_sites.tsv")
    ap.add_argument("--screen", default="results/tables/stage2_uga_screen.tsv")
    ap.add_argument("--scan", default="results/tables/stage2_secis_scan_classes.tsv")
    ap.add_argument("--truth", default="results/tables/stage2_truth_sets.tsv",
                    help="reference genes of the stage 2 benchmark: gold-standard models lifted onto KUM, AOM, Bar2")
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--loci", default="results/tables/stage2_census_loci.tsv")
    ap.add_argument("--table", default="results/tables/stage2_census.tsv")
    ap.add_argument("--markdown", default="results/tables/stage2_census.md")
    ap.add_argument("--hgt", default="results/tables/stage2_hgt_screen.tsv")
    ap.add_argument("--genes-table", default="results/tables/stage2_census_genes.tsv")
    ap.add_argument("--genes-markdown", default="results/tables/stage2_census_genes.md")
    a = ap.parse_args()
    sites = collections.defaultdict(list)
    for r in rows(a.sites): sites[r["genome"]].append(r)
    if not sites: raise SystemExit(f"no selenocysteine codon in {a.sites}")
    preds = {(r["genome"], r["prediction"]): r for r in rows(a.predictions)}
    contigs = {(r["genome"], r["contig"]): r for r in rows(a.contigs)}
    for r in rows(a.genome_table):
        if fnum(r.get("median_orf_fraction_2kb")) is not None: GENOME_ORF[r["genome"]] = fnum(r["median_orf_fraction_2kb"])
    if not contigs: print(f"WARNING: {a.contigs} not found; no contig checks", flush=True)
    rbh = {(r["genome"], r["prediction"]): r for r in rows(a.rbh)}
    gold = gold_genes(a.gold, ["OKI", "OSKA", "OdB3"], a.truth)
    place = placements(a.trees, lba.load_taxa(a.taxa))
    out = []
    for g in APPENDICULARIANS + ODIOICA_TESTS:
        for loc in loci(sites.get(g, [])):
            sp = [r for r in loc if r["selenoprofiles_label"] == "selenocysteine"]
            main_ = sp[0] if sp else loc[0]
            pid = main_.get("selenoprofiles_id", "") if sp else ""
            pid = "" if pid == "." else pid
            contig, strand = main_["contig"], main_["strand"]
            start = min(int(r["sec_codon_start"]) for r in loc)
            end = max(int(r["sec_codon_end"]) for r in loc)
            fam = main_["family"]
            secis = max((fnum(r["secis_bits"]) for r in loc if fnum(r["secis_bits"]) is not None), default=None)
            tx = max((int(r["transcripts"]) for r in loc if (r.get("transcripts") or "").isdigit()), default=None)
            ctg = contigs.get((g, contig), {})
            length = int(ctg["length"]) if ctg.get("length") else None
            after = (length - end if strand == "+" else start - 1) if length else None
            gname = next((n for c, s, b, e, n in gold.get(g, []) if c == contig and s == strand and b <= end and
                          e >= start), "")
            hit = rbh.get((g, pid), {})
            pl = place.get((g, pid), {})
            ml, pm = pl.get("ML"), pl.get("PMSF")
            comp_z, gc_z = composition_z(ctg), fnum(ctg.get("gc_robust_z"))
            tiara = (ctg.get("tiara_first_stage") or "").lower()
            tiara2 = (ctg.get("tiara_second_stage") or "").lower()
            flags = []
            if comp_z is not None and comp_z > COMPOSITION_Z: flags.append("composition outlier")
            if orf_dense(ctg): flags.append(f"ORF-dense contig ({ctg.get('contig_orf_fraction')})")
            if gc_z is not None and abs(gc_z) > GC_Z: flags.append("GC outlier")
            if ml and ml[0] != "animal": flags.append(f"nearest relatives {ml[0]} (ML)")
            if pm and pm[0] != "animal": flags.append(f"nearest relatives {pm[0]} (PMSF)")
            if pid and not ml: flags.append("not in the tree")
            if secis is None: flags.append("no SECIS hit")
            if after is not None and after < END_NT: flags.append(f"{after} nt after the codon")
            if not sp: flags.append("miniprot only")
            if length is not None and length < 3000: flags.append("contig < 3 kb")
            reason = ""
            best = hit.get("best_human_gene", "")
            if gname: decision = "counted (gold standard)"
            elif tiara in NONEUKARYOTIC or tiara2 in NONEUKARYOTIC - {"organelle"}:
                decision, reason = "left out", f"Tiara: {tiara2 if tiara2 in NONEUKARYOTIC else tiara}"
            elif hit.get("status", "").startswith("other gene") and best not in RELATED.get(fam, set()):
                decision, reason = "left out", f"best human hit in another family ({best})"
            elif "composition outlier" in flags and ml and ml[0] in ("prokaryotic", "bacterial or fungal", "non-animal"):
                decision, reason = "left out", f"composition outlier with {ml[0]} nearest relatives"
            elif orf_dense(ctg) and ml and ml[0] in ("prokaryotic", "bacterial or fungal", "non-animal"):
                decision, reason = "left out", f"ORF-dense contig with {ml[0]} nearest relatives"
            else: decision = "counted"
            p = preds.get((g, pid), {})
            out.append({"genome": g, "family": fam, "contig": contig, "strand": strand, "sec_codon_start": start,
                        "sec_codon_end": end, "codons": len(loc), "selenoprofiles_id": pid or "-",
                        "evidence": "Selenoprofiles and miniprot" if sp and any(r["sources"] in ("miniprot", "both")
                                                                                 for r in loc)
                        else "Selenoprofiles" if sp else "miniprot",
                        "coding_exons": p.get("coding_exons", ""), "gold_standard_gene": gname or "-",
                        "clade_gene": gname or "/".join(pl.get("gold_genes", [])) or SINGLE_GENE.get(fam, "-"),
                        "secis_bits": "" if secis is None else f"{secis:g}", "transcripts": "" if tx is None else tx,
                        "best_human_hit": f"{best} ({hit['status']})" if hit else "-",
                        "nearest_relatives_ML": f"{ml[0]}: {ml[1]} ({ml[2]})" if ml else "-",
                        "nearest_relatives_PMSF": f"{pm[0]}: {pm[1]} ({pm[2]})" if pm else "-",
                        "contig_length": length or "", "composition_z": ctg.get("composition_robust_z", ""),
                        "gc_z": ctg.get("gc_robust_z", ""), "tiara": tiara or "", "decision": decision,
                        "reason": reason, "flags": "; ".join(flags)})
    # loci of the step 4 screens that the final run did not predict
    extra = []
    for x in screen_loci(a.screen, a.scan, APPENDICULARIANS):
        if near(out, x["genome"], x["contig"], x["strand"], x["pos"], x["kind"]): continue
        if near(extra, x["genome"], x["contig"], x["strand"], x["pos"], x["kind"]): continue
        extra.append(x)
    gene_at = lambda g, c, st, pos: next((n_ for c_, s_, b, e, n_ in gold.get(g, []) if c_ == c and s_ == st and
                                          b - 300 <= pos <= e + 300), "")
    for _ in range(3):   # resolve genes through partners (a partner may itself be a screen locus)
        for x in extra:
            if x["gene"] or not x["partner"]: continue
            pg, pc, pp, ps = x["partner"]
            hit = gene_at(pg, pc, ps, pp) or next((r["clade_gene"] for r in near(out, pg, pc, ps, pp, "TGA")
                                                   if r["clade_gene"] != "-"), "") or \
                next((y["gene"] for y in near(extra, pg, pc, ps, pp, "TGA") if y["gene"]), "")
            x["gene"] = hit
    for x in extra:
        if x["candidate"] and not x["gene"]: continue          # unsupported candidates (repeats) left aside
        fam = re.split(r"[_.]", x["gene"])[0] if x["gene"] else "-"
        out.append({k: "" for k in out[0]} | {
            "genome": x["genome"], "family": fam, "contig": x["contig"], "strand": x["strand"],
            "sec_codon_start": x["pos"] if x["kind"] == "TGA" else "", "sec_codon_end": "", "codons": 1 if x["kind"] == "TGA" else 0,
            "selenoprofiles_id": "-", "evidence": x["evidence"], "gold_standard_gene": "-",
            "clade_gene": x["gene"] or "-", "secis_bits": x["secis"], "best_human_hit": "-",
            "nearest_relatives_ML": "-", "nearest_relatives_PMSF": "-", "decision": "counted (screen)",
            "flags": "no gene model" + ("; SECIS hit, codon not placed" if x["kind"] == "SECIS" else "")})
    with open(a.loci, "w") as o:
        w = csv.DictWriter(o, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(out)
    # counts: appendicularians from the census, other genomes from stage 1 (same locus definition, no screen)
    counts = collections.defaultdict(collections.Counter)
    for r in out:
        if r["decision"] in ("counted", "counted (gold standard)"): counts[r["genome"]][r["family"]] += 1
    st1 = collections.defaultdict(list)
    for r in rows(a.stage1): st1[r["genome"]].append(r)
    for g in OTHERS:
        for loc in loci(st1.get(g, [])):
            sp = [r for r in loc if r["selenoprofiles_label"] == "selenocysteine"]
            counts[g][(sp[0] if sp else loc[0])["family"]] += 1
    genomes = [g for g in APPENDICULARIANS + ODIOICA_TESTS + OTHERS if counts.get(g)]
    fams = sorted({f for g in genomes for f in counts[g]})
    with open(a.table, "w") as o:
        o.write("family\t" + "\t".join(genomes) + "\n")
        for f in fams: o.write(f + "\t" + "\t".join(str(counts[g][f]) for g in genomes) + "\n")
        o.write("total\t" + "\t".join(str(sum(counts[g].values())) for g in genomes) + "\n")
    left = collections.Counter(r["genome"] for r in out if r["decision"] == "left out")
    with open(a.markdown, "w") as o:
        o.write("| Family | " + " | ".join(NAMES.get(g, g) for g in genomes) + " |\n|---|" + "---|" * len(genomes) + "\n")
        for f in fams:
            o.write(f"| {f} | " + " | ".join(str(counts[g][f]) if counts[g][f] else "" for g in genomes) + " |\n")
        o.write("| Total | " + " | ".join(str(sum(counts[g].values())) for g in genomes) + " |\n")
        o.write("| Left out | " + " | ".join(str(left[g]) if g in APPENDICULARIANS + ODIOICA_TESTS else "not screened"
                                             for g in genomes) + " |\n")
    # counted loci by the O. dioica gene of their clade (appendicularians only)
    by_gene = collections.defaultdict(collections.Counter)
    with_secis = collections.defaultdict(collections.Counter)
    by_screen = collections.defaultdict(collections.Counter)
    for r in out:
        if not r["decision"].startswith("counted"): continue
        gene = r["clade_gene"] if r["clade_gene"] != "-" else f"{r['family']} (no O. dioica gene in its clade)"
        if r["decision"] == "counted (screen)": by_screen[gene][r["genome"]] += 1; continue
        by_gene[gene][r["genome"]] += 1
        with_secis[gene][r["genome"]] += r["secis_bits"] != ""
    for x in by_screen: by_gene[x] += collections.Counter()
    app = [g for g in APPENDICULARIANS + ODIOICA_TESTS if any(by_gene[x][g] or by_screen[x][g] for x in by_gene)]
    order = sorted(by_gene, key=lambda x: (GENE_ORDER.index(x) if x in GENE_ORDER else len(GENE_ORDER), x))
    with open(a.genes_table, "w") as o:
        o.write("gene\t" + "\t".join(app) + "\n")
        o.write("# counts: loci predicted by the final run; '+n': loci found only by the step 4 screens\n")
        for x in order: o.write(x + "\t" + "\t".join(f"{by_gene[x][g]}" + (f"+{by_screen[x][g]}" if by_screen[x][g] else "")
                                                     for g in app) + "\n")
    app_md = [g for g in app if g in APPENDICULARIANS]          # the new O. dioica assemblies: in the TSV only
    with open(a.genes_markdown, "w") as o:
        o.write("| Gene of the clade | " + " | ".join(NAMES.get(g, g) for g in app_md) + " |\n|---|" + "---|" * len(app_md) + "\n")
        for x in order:
            if not any(by_gene[x][g] or by_screen[x][g] for g in app_md): continue
            label = x.replace(" (no O. dioica gene in its clade)", ", other clades")
            o.write(f"| {label} | " + " | ".join(
                ((f"{by_gene[x][g]} ({with_secis[x][g]})" if by_gene[x][g] else "") +
                 (f" +{by_screen[x][g]}" if by_screen[x][g] else "")).strip() for g in app_md) + " |\n")
    # strict HGT screen
    hgt = []
    by_contig = collections.defaultdict(list)
    for (g, pid), p in preds.items(): by_contig[(g, p["contig"])].append(pid)
    foreign = ("prokaryotic", "bacterial or fungal")
    for (g, pid), pl in sorted(place.items()):
        if pid.startswith("gold_"): continue
        mlt, pmt = pl.get("ML_tunicates"), pl.get("PMSF_tunicates")
        if pl.get("ML") and pl["ML"][0] in foreign:          # foreign relatives next to the appendicularian sequences
            ml, pm, level = pl.get("ML"), pl.get("PMSF"), "appendicularians"
        elif mlt and mlt[0] in foreign and mlt[4]:            # ... next to a clade of tunicate sequences
            ml, pm, level = mlt, pmt, "tunicates"
        else: continue
        p = preds.get((g, pid), {})
        ctg = contigs.get((g, p.get("contig", "")), {})
        tiara = (ctg.get("tiara_first_stage") or "").lower()
        comp_z = composition_z(ctg)
        animal_neighbours = sorted(x for x in by_contig.get((g, p.get("contig", "")), []) if x != pid and
                                   place.get((g, x), {}).get("ML", ("",))[0] == "animal")
        hit = rbh.get((g, pid), {})
        best = hit.get("best_human_gene", "")
        if tiara in NONEUKARYOTIC: status = f"contamination (Tiara: {tiara})"
        elif comp_z is not None and comp_z > COMPOSITION_Z: status = "contamination likely (composition outlier)"
        elif orf_dense(ctg): status = f"contamination likely (ORF-dense contig, {ctg.get('contig_orf_fraction')})"
        elif hit.get("status", "").startswith("other gene") and best not in RELATED.get(ml[3], set()):
            status = f"other family ({best}): its animal homologs are not in this tree"
        elif hit.get("status") == "no hit" and hit.get("expected_gene", "-") != "-":
            status = "family not confirmed (no human hit)"
        elif not pm or pm[0] not in foreign: status = "placement not stable (PMSF)"
        elif p.get("coding_exons") == "1": status = "candidate, no intron (needs presence in other assemblies or species)"
        else: status = "candidate"
        hgt.append({"genome": g, "prediction": pid, "family": ml[3], "label": pid.rsplit(".", 1)[-1],
                    "foreign relatives next to": level,
                    "other tunicate sequences in that clade": ml[4] if level == "tunicates" else 0,
                    "genomes of its appendicularian clade": ",".join(pl.get("genomes", [])),
                    "best_human_hit": f"{best} ({hit['status']})" if hit else "-",
                    "nearest_relatives_ML": f"{ml[0]}: {ml[1]} ({ml[2]})",
                    "nearest_relatives_PMSF": f"{pm[0]}: {pm[1]} ({pm[2]})" if pm else "-",
                    "foreign sequences (ML)": " ".join(ml[5]),
                    "contig": p.get("contig", ""), "contig_length": ctg.get("length", ""),
                    "coding_exons": p.get("coding_exons", ""), "composition_z": ctg.get("composition_robust_z", ""),
                    "composition_z_2kb": ctg.get("composition_robust_z_2kb", ""),
                    "contig_orf_fraction": ctg.get("contig_orf_fraction", ""), "contig_gc": ctg.get("gc", ""),
                    "gc_z": ctg.get("gc_robust_z", ""), "tiara": tiara,
                    "animal-placed predictions on the contig": " ".join(animal_neighbours) or "-", "status": status})
    if hgt:
        with open(a.hgt, "w") as o:
            w = csv.DictWriter(o, fieldnames=list(hgt[0]), delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(hgt)
    print(open(a.markdown).read())
    print(collections.Counter(r["decision"] + (f": {r['reason']}" if r["reason"] else "") for r in out))
    print(collections.Counter(h["status"] for h in hgt))


if __name__ == "__main__":
    main()
