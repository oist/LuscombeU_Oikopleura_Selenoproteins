#!/usr/bin/env python3
"""Add BibTeX entries to the public bibliography, keeping only bibliographic metadata.

Abstracts and other text of the cited works (abstract, annotations, notes, keywords, file links) are not copied:
the public file holds only the fields listed in KEEP. Entries are merged by citation key (a new entry replaces an
older one with the same key), and the file is written sorted by key.

usage: scripts/bib_public.py [--bib manuscript/references.bib] < entries.bib
"""
import argparse, re, sys

KEEP = ["author", "editor", "title", "journal", "journaltitle", "booktitle", "year", "date", "month", "volume",
        "number", "issue", "pages", "eid", "publisher", "organization", "institution", "school", "series",
        "edition", "address", "doi", "pmid", "pmcid", "eprint", "eprinttype", "archiveprefix", "primaryclass",
        "url", "urldate", "issn", "isbn", "language", "howpublished", "type", "version"]


def entries(text):
    """Yield (type, key, {field: raw value}) for each @type{key, ...} entry, with balanced braces."""
    i = 0
    while True:
        m = re.compile(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", re.S).search(text, i)
        if not m:
            return
        depth, j = 1, m.end()
        while j < len(text) and depth:
            depth += {"{": 1, "}": -1}.get(text[j], 0)
            j += 1
        body = text[m.end(): j - 1]
        fields, k = {}, 0
        while True:
            f = re.compile(r"\s*(\w[\w-]*)\s*=\s*", re.S).match(body, k)
            if not f:
                break
            k = f.end()
            if body[k] == "{":
                d, s = 1, k + 1
                k += 1
                while d:
                    d += {"{": 1, "}": -1}.get(body[k], 0)
                    k += 1
                value = "{" + body[s:k - 1] + "}"
            elif body[k] == '"':
                e = body.index('"', k + 1)
                value, k = body[k:e + 1], e + 1
            else:
                e = re.compile(r"[,\s]").search(body, k)
                e = e.start() if e else len(body)
                value, k = body[k:e], e
            fields[f.group(1).lower()] = value
            c = re.compile(r"\s*,?").match(body, k)
            k = c.end()
        yield m.group(1).lower(), m.group(2), fields
        i = j


def render(etype, key, fields):
    kept = [(f, fields[f]) for f in KEEP if f in fields]
    width = max([len(f) for f, _ in kept] or [0])
    return f"@{etype}{{{key},\n" + "".join(f"  {f.ljust(width)} = {v},\n" for f, v in kept) + "}\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bib", default="manuscript/references.bib")
    a = ap.parse_args()
    try:
        existing = {k: (t, f) for t, k, f in entries(open(a.bib).read())}
    except FileNotFoundError:
        existing = {}
    new = list(entries(sys.stdin.read()))
    dropped = set()
    for t, k, f in new:
        dropped |= {x for x in f if x not in KEEP}
        existing[k] = (t, f)
    with open(a.bib, "w") as out:
        out.write("\n".join(render(t, k, f) for k, (t, f) in sorted(existing.items(), key=lambda x: x[0].lower())))
    print(f"{len(new)} entries added or updated; {len(existing)} in {a.bib}; fields not kept: {', '.join(sorted(dropped)) or 'none'}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
