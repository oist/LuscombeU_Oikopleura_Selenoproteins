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
  counted (gold standard)  the locus lies in a gold-standard gene of OKI2018, OSKA2016 or OdB3 (always counted);
  left out                 Tiara classifies its contig as bacteria, archaea, prokarya or organelle; or the best human
                           hit of its prediction lies in another family (RELATED genes count as the family: MIEN1 for
                           SELENOW, SEPHS1 for SEPHS2); or its contig is a composition outlier (robust z-score of the
                           tetranucleotide distance above COMPOSITION_Z) and its nearest relatives in the ML tree of
                           its family are all non-animal;
  counted                  otherwise.
Flags reported with each locus (not decisive alone): composition outlier; GC outlier (robust z-score beyond GC_Z);
nearest relatives in the ML and PMSF trees (prokaryotic, bacterial or fungal, other non-animal, animal, not in the
tree); no SECIS hit (the model of the final run, the Oikopleura model at 20 bits); fewer than END_NT nucleotides
after the codon on the contig (no room for a SECIS element in the search window); miniprot only; contig shorter than
3,000 nt (no Tiara class).

Strict HGT screen (candidates only; the tests come after): every appendicularian leaf of the round 2 ML trees whose
nearest relatives are all bacterial, archaeal or fungal, whatever its label. Status: "contamination (Tiara ...)";
"contamination likely (composition outlier)"; "not testable (one coding exon)" (bacterial genes have no intron, and
the strict standard requires introns or spliced transcripts); "placement not stable (PMSF)" when the PMSF tree does
not give the same kind of nearest relatives; "candidate" otherwise, to be tested (AU tests against vertical descent,
other assemblies or individuals, transcripts).

Outputs: --loci (one row per locus), --table and --markdown (counted loci per family and genome, appendicularians
from this census, other genomes from stage 1 with the same locus definition and no screen), --hgt.
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
GC_Z = 5.0
RELATED = {"SELENOW": {"MIEN1"}, "SEPHS2": {"SEPHS1"}}
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


def gold_genes(gold, genomes):
    """genome -> [(contig, strand, start, end, gene)] from the gene features of the gold-standard GFF3 files."""
    out = collections.defaultdict(list)
    for g in genomes:
        path = f"{gold}/{g}.gff3"
        if not os.path.exists(path): continue
        for l in open(path):
            f = l.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene": continue
            m = re.search(r"Name=([^;]+)", f[8])
            out[g].append((f[0], f[6], int(f[3]), int(f[4]), m.group(1) if m else f[8]))
    return out


def placements(trees, taxa):
    """(genome, prediction) -> {tree: (kind, groups, support)} for the appendicularian leaves of the round 2 trees;
    kind: prokaryotic, bacterial or fungal, non-animal, animal. Gold-standard leaves stand for the predictions they
    replaced, and kept sequences for the identical ones removed."""
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
            for n in names:
                if tax[n] != "Appendicularians": continue
                cand = [s for s in sd if n in s and any(tax[x] != "Appendicularians" for x in s)]
                if not cand: continue
                near = min(cand, key=len)
                others = collections.Counter(tax[x] for x in near if tax[x] != "Appendicularians")
                kinds = set(others) - {"unclassified"}
                if kinds and kinds <= PROKARYOTES: kind = "prokaryotic"
                elif kinds and kinds <= PROKARYOTES | {"Fungi"}: kind = "bacterial or fungal"
                elif kinds and not kinds & ANIMALS: kind = "non-animal"
                elif kinds: kind = "animal"
                else: kind = "unclassified"
                text = ", ".join(f"{k} {v}" for k, v in others.most_common())
                support = sd.get(near)
                for leaf in [n] + stand_for.get(n, []):
                    g, _, pid = leaf.partition("_")
                    out[(g, pid)][label] = (kind, text, "NA" if support is None else f"{support:.0f}", fam)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", default="results/tables/stage2_final_sec_sites.tsv")
    ap.add_argument("--predictions", default="results/tables/stage2_final_predictions.tsv")
    ap.add_argument("--contigs", default="results/tables/stage2_census_contigs.tsv")
    ap.add_argument("--rbh", default="results/tables/stage2_selenoprotein_rbh.tsv")
    ap.add_argument("--trees", default="results/trees/selenoproteins")
    ap.add_argument("--gold", default="curation/gold_standard")
    ap.add_argument("--stage1", default="results/tables/stage1_sec_sites.tsv")
    ap.add_argument("--taxa", default="config/taxon_groups.tsv")
    ap.add_argument("--loci", default="results/tables/stage2_census_loci.tsv")
    ap.add_argument("--table", default="results/tables/stage2_census.tsv")
    ap.add_argument("--markdown", default="results/tables/stage2_census.md")
    ap.add_argument("--hgt", default="results/tables/stage2_hgt_screen.tsv")
    a = ap.parse_args()
    sites = collections.defaultdict(list)
    for r in rows(a.sites): sites[r["genome"]].append(r)
    if not sites: raise SystemExit(f"no selenocysteine codon in {a.sites}")
    preds = {(r["genome"], r["prediction"]): r for r in rows(a.predictions)}
    contigs = {(r["genome"], r["contig"]): r for r in rows(a.contigs)}
    if not contigs: print(f"WARNING: {a.contigs} not found; no contig checks", flush=True)
    rbh = {(r["genome"], r["prediction"]): r for r in rows(a.rbh)}
    gold = gold_genes(a.gold, ["OKI", "OSKA", "OdB3"])
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
            comp_z, gc_z = fnum(ctg.get("composition_robust_z")), fnum(ctg.get("gc_robust_z"))
            tiara = (ctg.get("tiara_first_stage") or "").lower()
            tiara2 = (ctg.get("tiara_second_stage") or "").lower()
            flags = []
            if comp_z is not None and comp_z > COMPOSITION_Z: flags.append("composition outlier")
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
            else: decision = "counted"
            p = preds.get((g, pid), {})
            out.append({"genome": g, "family": fam, "contig": contig, "strand": strand, "sec_codon_start": start,
                        "sec_codon_end": end, "codons": len(loc), "selenoprofiles_id": pid or "-",
                        "evidence": "Selenoprofiles and miniprot" if sp and any(r["sources"] in ("miniprot", "both")
                                                                                 for r in loc)
                        else "Selenoprofiles" if sp else "miniprot",
                        "coding_exons": p.get("coding_exons", ""), "gold_standard_gene": gname or "-",
                        "secis_bits": "" if secis is None else f"{secis:g}", "transcripts": "" if tx is None else tx,
                        "best_human_hit": f"{best} ({hit['status']})" if hit else "-",
                        "nearest_relatives_ML": f"{ml[0]}: {ml[1]} ({ml[2]})" if ml else "-",
                        "nearest_relatives_PMSF": f"{pm[0]}: {pm[1]} ({pm[2]})" if pm else "-",
                        "contig_length": length or "", "composition_z": ctg.get("composition_robust_z", ""),
                        "gc_z": ctg.get("gc_robust_z", ""), "tiara": tiara or "", "decision": decision,
                        "reason": reason, "flags": "; ".join(flags)})
    with open(a.loci, "w") as o:
        w = csv.DictWriter(o, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(out)
    # counts: appendicularians from the census, other genomes from stage 1 (same locus definition, no screen)
    counts = collections.defaultdict(collections.Counter)
    for r in out:
        if r["decision"].startswith("counted"): counts[r["genome"]][r["family"]] += 1
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
    # strict HGT screen
    hgt = []
    by_contig = collections.defaultdict(list)
    for (g, pid), p in preds.items(): by_contig[(g, p["contig"])].append(pid)
    for (g, pid), pl in sorted(place.items()):
        ml, pm = pl.get("ML"), pl.get("PMSF")
        if not ml or ml[0] not in ("prokaryotic", "bacterial or fungal") or pid.startswith("gold_"): continue
        p = preds.get((g, pid), {})
        ctg = contigs.get((g, p.get("contig", "")), {})
        tiara = (ctg.get("tiara_first_stage") or "").lower()
        comp_z = fnum(ctg.get("composition_robust_z"))
        animal_neighbours = sorted(x for x in by_contig.get((g, p.get("contig", "")), []) if x != pid and
                                   place.get((g, x), {}).get("ML", ("",))[0] == "animal")
        if tiara in NONEUKARYOTIC: status = f"contamination (Tiara: {tiara})"
        elif comp_z is not None and comp_z > COMPOSITION_Z: status = "contamination likely (composition outlier)"
        elif p.get("coding_exons") == "1": status = "not testable (one coding exon)"
        elif not pm or pm[0] not in ("prokaryotic", "bacterial or fungal"): status = "placement not stable (PMSF)"
        else: status = "candidate"
        hgt.append({"genome": g, "prediction": pid, "family": ml[3], "label": pid.rsplit(".", 1)[-1],
                    "nearest_relatives_ML": f"{ml[0]}: {ml[1]} ({ml[2]})",
                    "nearest_relatives_PMSF": f"{pm[0]}: {pm[1]} ({pm[2]})" if pm else "-",
                    "contig": p.get("contig", ""), "contig_length": ctg.get("length", ""),
                    "coding_exons": p.get("coding_exons", ""), "composition_z": ctg.get("composition_robust_z", ""),
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
