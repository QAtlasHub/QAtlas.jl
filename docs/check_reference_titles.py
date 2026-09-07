#!/usr/bin/env python3
"""Cross-check every bib entry's TITLE against the one its DOI actually resolves to.

`doiget verify` asks only whether a reference resolves, and that is not enough: the most
likely way to get a citation wrong is a one-character slip in the DOI, which usually lands on
a real neighbouring paper. Measured on this bibliography: 2 of 18 sampled DOIs have a resolving neighbour one digit
away, so a slip lands on somebody else's real paper about one time in nine and the
resolution check passes it. Only comparing the title catches that.

Crossref REST directly, so the runner needs no tooling beyond Python.
"""
import json
import re
import sys
import urllib.parse
import urllib.request

BIB = sys.argv[1] if len(sys.argv) > 1 else "docs/references.bib"
ALLOW = sys.argv[2] if len(sys.argv) > 2 else "docs/references.title-allow"


# Crossref returns APS/Springer titles with the maths as inline MathML or HTML, where the bib
# has LaTeX. Stripping the markup first is what makes the words comparable; without it every
# title containing a symbol reads as a mismatch (measured: 25 of 231 on this bibliography).
_TAG = re.compile(r"<[^>]*>")
_ENTITY = re.compile(r"&[a-zA-Z]+;|&#\d+;")
# `\pm`, `\infty`, `\mathrm` … survive as words on the bib side while the same symbol is a
# glyph on the Crossref side and gets stripped, so the two disagree on a difference that is
# only notation.
_LATEX = re.compile(r"\\[a-zA-Z]+")


def normalise(s):
    # Compare on words only: publishers vary punctuation, case and LaTeX braces freely, and a
    # mismatch in those is not a wrong citation.
    s = _LATEX.sub(" ", _ENTITY.sub(" ", _TAG.sub(" ", s)))
    return " ".join(re.sub(r"[^a-z0-9 ]", " ", s.lower()).split())


def allowed():
    try:
        with open(ALLOW) as f:
            return {l.split("#")[0].strip() for l in f if l.split("#")[0].strip()}
    except FileNotFoundError:
        return set()


def entries(path):
    with open(path) as f:
        text = f.read()
    # Strip `%` comments before parsing: a literal `@…` inside one would read as an entry.
    text = "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("%"))
    for block in re.findall(r"@\w+\s*\{(.*?)\n\}", text, re.S):
        key = block.split(",", 1)[0].strip()
        fields = dict(re.findall(r"(\w+)\s*=\s*\{(.*?)\}\s*,?\s*\n", block, re.S))
        yield key, fields


def crossref_title(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    req = urllib.request.Request(url, headers={"User-Agent": "reference-title-check"})
    with urllib.request.urlopen(req, timeout=30) as r:
        msg = json.load(r)["message"]
    return msg.get("title", [""])[0]


skip = allowed()
bad = 0
for key, f in entries(BIB):
    doi = f.get("doi", "").strip()
    if not doi:
        print(f"SKIP  {key}: no doi field")
        continue
    if doi in skip or key in skip:
        print(f"ALLOW {key}: {doi}")
        continue
    want = f.get("title", "").strip()
    try:
        got = crossref_title(doi)
    except Exception as e:  # noqa: BLE001 — a network failure must not read as a pass
        print(f"ERROR {key}: {doi}: {e}")
        bad += 1
        continue
    if normalise(want) == normalise(got):
        print(f"OK    {key}: {doi}")
    else:
        print(f"FAIL  {key}: {doi}\n        bib says : {want}\n        DOI is   : {got}")
        bad += 1

print(f"\n{bad} mismatch(es)")
sys.exit(1 if bad else 0)
