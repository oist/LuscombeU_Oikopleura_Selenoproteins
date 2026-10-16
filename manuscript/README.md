# Manuscript

## Bibliography

`references.bib` holds bibliographic metadata only: authors, title, journal, volume, pages, year, and the
identifiers (DOI, PMID, PMCID, URL). Abstracts and other text of the cited works are not stored in this public
repository, because they are copyrighted by their authors or publishers.

To add references, pass full BibTeX entries (for example a Zotero export with abstracts) through the filter, which
keeps the metadata fields only and merges by citation key:

```bash
scripts/bib_public.py --bib manuscript/references.bib < new_entries.bib
```

The identifiers let anyone, or an AI assistant with access to PubMed or bioRxiv, retrieve an abstract when needed.
Entries should therefore carry a DOI, and a PMID when the article is indexed in PubMed.
