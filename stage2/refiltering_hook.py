"""Selenoprofiles action hook for stage2/diagnose_refiltering.py (not a script).

Called by a pre_filtering action of Selenoprofiles (-ACTION.pre_filtering...) with the candidate prediction x, before
the filters of the profile are applied. Writes one line per prediction to $DIAG_OUT: prediction, position, label,
measures used by the filters, and the value of the p2g_filtering and p2g_refiltering expressions of the profile
configuration ($DIAG_CONFIG) and of each top-level "and" clause of p2g_refiltering, so that the clause that removes a
prediction can be read.
"""
import ast, os

MEASURES = ["x.profile.nseq()", "x.label", "x.coverage()", "x.sequence_identity_with_profile()", "x.awsi_z_score()",
            "x.awsi_z_score(with_coverage=False)", "x.sec_is_aligned()"]


DEFAULTS = {"p2g_filtering": "len(x.protein()) >60 or x.coverage()> 0.4",   # Selenoprofiles 4.5.8 configuration
            "p2g_refiltering": "x.awsi_filter()"}


def expression(config, key):
    for l in open(config):
        k, _, v = l.partition("=")
        if k.strip() == key and v.strip(): return v.strip()
    return DEFAULTS[key]


def clauses(expr):
    tree = ast.parse(expr, mode="eval").body
    if isinstance(tree, ast.BoolOp) and isinstance(tree.op, ast.And):
        return [ast.unparse(v) for v in tree.values]
    return []


def value(code, x):
    try:
        v = eval(code, {"x": x, "int": int, "len": len})
    except Exception as e:                   # a measure that cannot be computed is reported, not fatal
        return f"error: {type(e).__name__}: {e}"
    return f"{v:.3f}" if isinstance(v, float) else str(v)


def report(x):
    config, out = os.environ["DIAG_CONFIG"], os.environ["DIAG_OUT"]
    filt, refilt = expression(config, "p2g_filtering"), expression(config, "p2g_refiltering")
    exprs = MEASURES + [filt, refilt] + clauses(refilt)
    if not os.path.exists(out):
        with open(out, "w") as o:
            o.write("\t".join(["prediction", "contig", "strand", "start", "end"] + exprs) + "\n")
    try:
        lo, hi = x.boundaries()
    except Exception:
        lo, hi = "", ""
    row = [value("x.output_id()", x), str(getattr(x, "chromosome", "")), str(getattr(x, "strand", "")), str(lo),
           str(hi)] + [value(e, x) for e in exprs]
    with open(out, "a") as o:
        o.write("\t".join(row) + "\n")
