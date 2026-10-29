"""Taxonomic group of a tree leaf from its name (shared by the tree, orthology and figure scripts).

Leaf names come from several sources: genome codes of this study (OKI_GPX.1.selenocysteine), human references
(Human_GPX4_NP_...), and the sequences of the Selenoprofiles profiles, whose names take many forms:
  Genus_species_...                       Gallus_gallus_gi_143771645_ref_NP_001006
  <family>.Genus_species / <family>_Genus  Sel15.Gallus_gallus_Sel15..., SelL_salmo_salar_..., TXNRD1_XM_..._Sparus
  Se_user_<name>_Genus_species             Se_user_SelK_Strongylocentrotus_purpuratus
  UniProt entry names (MNEMONIC_SPECIES)   gi_182702126_sp_Q8UVX8.3_IOD2_NEOFS_...
  Ensembl protein identifiers              ENSGMOP00000019532_Gmor_ENSGMOP...
  descriptive bacterial names              marine_gamma_proteobacterium_HTCC2080_...
The group is looked up, in this order: genome code or genus at the start of the name; Ensembl species prefix
(config/taxon_codes.tsv, kind ensembl); UniProt species mnemonic (kind uniprot); the first token anywhere in the
name that is a known genus (config/taxon_groups.tsv, case-insensitive); descriptive bacterial names. Otherwise
"unclassified".
"""
import functools, os, re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# groups of config/taxon_groups.tsv that are animals (Metazoa); the others, except "unclassified", are not
ANIMALS = {"Mammals", "Birds", "Reptiles", "Amphibians", "Ray-finned fishes", "Lobe-finned fishes",
           "Cartilaginous fishes", "Jawless fishes", "Cephalochordates", "Other tunicates", "Appendicularians",
           "Echinoderms", "Hemichordates", "Other insects", "Flies", "Mosquitoes", "Crustaceans", "Chelicerates",
           "Myriapods", "Nematodes", "Molluscs", "Annelids", "Brachiopods", "Flatworms", "Priapulids", "Cnidarians",
           "Placozoans"}


def load_taxa(path):
    """genus (or genome code) -> taxonomic group, from config/taxon_groups.tsv"""
    taxa = {}
    if os.path.exists(path):
        for l in open(path):
            if l.startswith("#") or l.startswith("genus\t") or "\t" not in l: continue
            g, grp = l.rstrip("\n").split("\t")[:2]; taxa[g] = grp
    return taxa


@functools.lru_cache(maxsize=None)
def load_codes(path=os.path.join(ROOT, "config", "taxon_codes.tsv")):
    """(kind, code) -> group, from config/taxon_codes.tsv (kind: ensembl or uniprot)."""
    codes = {}
    if os.path.exists(path):
        for l in open(path):
            if l.startswith("#") or l.startswith("kind\t") or "\t" not in l: continue
            kind, code, _species, grp = l.rstrip("\n").split("\t")[:4]
            codes[(kind, code)] = grp
    return codes


def classify(name, taxa):
    """Taxonomic group of a leaf name (see module documentation)."""
    first = re.sub(r"^(sps1|sps2|sps)_", "", name.split("|")[0])
    if first.startswith("Human_"): return "Mammals"
    tokens = [t for t in re.split(r"[_.|\s]+", first) if t]
    if tokens and tokens[0] in taxa: return taxa[tokens[0]]
    codes = load_codes()
    m = re.match(r"ENS([A-Z]*)P\d", first)
    if m and ("ensembl", m.group(1)) in codes: return codes[("ensembl", m.group(1))]
    for m in re.finditer(r"(?=(?:^|_)[A-Z0-9]{2,6}_([A-Z0-9]{3,5})(?:_|$))", first):   # overlapping matches
        if ("uniprot", m.group(1)) in codes: return codes[("uniprot", m.group(1))]
    lower = {k.lower(): v for k, v in taxa.items()}
    for t in tokens:
        if len(t) > 2 and t[0].isalpha() and t.lower() in lower: return lower[t.lower()]
    if any(re.fullmatch(r"(proteo)?bacterium|bacteria", t.lower()) for t in tokens): return "Bacteria"
    return "unclassified"
