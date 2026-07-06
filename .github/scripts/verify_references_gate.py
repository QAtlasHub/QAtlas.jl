#!/usr/bin/env python3
"""Gate `doiget verify` output against an allowlist and emit a Markdown report.

Reads doiget's JSON-Lines `verify` output on stdin (one record per reference,
with a `status` field), plus an optional allowlist file. Classifies each entry:

    valid                        -> OK, resolves on Crossref / arXiv
    illegal | absent             -> BROKEN: malformed id, or does not resolve
                                    (fabricated / mistyped / not indexed)
    unreachable | unverifiable   -> TRANSIENT: network / 429 / no-id / coverage
                                    gap; NOT gated (would make CI flaky)

A BROKEN entry fails the gate (exit 1) UNLESS its DOI / arXiv id *or* its bibkey
appears in the allowlist — i.e. a human has explicitly vouched for it. The
allowlist (`docs/references.allow` by default) is one id/key per line, `#` starts
a comment; put the reason in the comment so every exception is self-documenting.

Completeness (the "no verification left out" guarantee): when `--bib` is given,
every reference entry in the bibliography must have produced a verify record. If
doiget aborted, timed out, or skipped an entry, that entry is UNVERIFIED and
fails the gate — a truncated run can never masquerade as a clean one. Unverified
entries are NOT allowlist-downgradable (the check simply did not run).

So a reference that does not exist — or that was never checked at all — can never
pass silently. Always writes a Markdown report (for the PR comment / job summary)
regardless of outcome.

Usage (also runnable locally):
    doiget verify docs/references.bib --mode json \\
      | python3 .github/scripts/verify_references_gate.py \\
          --bib docs/references.bib --allow docs/references.allow
"""
import argparse
import json
import os
import re
import sys

BROKEN = {"illegal", "absent"}
TRANSIENT = {"unreachable", "unverifiable"}

# `@type{key,` — a reference entry. @comment / @string / @preamble / @set are
# BibTeX machinery, not references, so they are excluded from the completeness count.
_ENTRY_RE = re.compile(r"@(\w+)\s*\{\s*([^,\s}]+)", re.IGNORECASE)
_NON_REF_TYPES = {"comment", "string", "preamble", "set"}


def load_allow(path):
    allow = set()
    if path and os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                token = line.split("#", 1)[0].strip()
                if token:
                    allow.add(token.lower())
    return allow


def bib_keys(path):
    keys = []
    if path and os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for m in _ENTRY_RE.finditer(fh.read()):
                if m.group(1).lower() not in _NON_REF_TYPES:
                    keys.append(m.group(2))
    return keys


def detail_of(entry):
    err = entry.get("error")
    if isinstance(err, dict):
        return str(err.get("message", ""))[:90]
    return ""


def render(ok, transient, excepted, broken, unverified):
    out = ["<!-- verify-references-gate -->", "## Reference check — `doiget verify`", ""]
    out.append(f"- ✅ resolved: **{len(ok)}**")
    if transient:
        out.append(f"- ⚠️ transient (network / no-id / not-yet-indexed — not gated): **{len(transient)}**")
    if excepted:
        out.append(f"- 🟡 allowlisted exceptions: **{len(excepted)}**")
    out.append(f"- {'❌' if broken else '☑️'} broken (unresolved / malformed): **{len(broken)}**")
    if unverified:
        out.append(f"- ❌ unverified (no record — check did not run): **{len(unverified)}**")
    out.append("")

    def table(title, rows):
        if not rows:
            return []
        block = [f"### {title}", "", "| bibkey | ref | status | detail |", "|---|---|---|---|"]
        for e in rows:
            block.append(
                f"| `{e.get('entry_key') or ''}` | `{e.get('ref') or ''}` | "
                f"{e.get('status') or ''} | {detail_of(e)} |"
            )
        block.append("")
        return block

    out += table("❌ Broken references — fix the id, or add an explicit exception", broken)
    if unverified:
        out += ["### ❌ Unverified references — the check did not run on these", ""]
        out += ["| bibkey |", "|---|"]
        out += [f"| `{k}` |" for k in unverified]
        out += [""]
    out += table("⚠️ Transient (not gated)", transient)
    out += table("🟡 Allowlisted exceptions", excepted)

    if broken:
        out += [
            "A broken reference does not resolve on Crossref / arXiv — it is "
            "**fabricated**, mistyped, or genuinely not indexed. Fix the id in the "
            "bibliography, or, if it is real but unresolvable, add its DOI / arXiv id "
            "(or bibkey) to `docs/references.allow` **with a comment saying why**.",
            "",
        ]
    if unverified:
        out += [
            "An unverified reference produced no result — doiget aborted, timed out, "
            "or could not parse it, so it was **never checked**. Re-run the job; if it "
            "persists, the bibliography entry is malformed. This is not allowlistable — "
            "the point is that every reference is actually verified.",
            "",
        ]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bib", default="", help="bibliography path, for the completeness check")
    ap.add_argument("--allow", default="docs/references.allow", help="allowlist file")
    ap.add_argument("--report", default="", help="write the Markdown report to this path")
    ap.add_argument("--github-output", default="", help="write broken=/unverified=/fail=/ok= here")
    args = ap.parse_args()

    allow = load_allow(args.allow)

    ok, transient, excepted, broken = [], [], [], []
    seen = set()
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            entry = json.loads(line)
        except json.JSONDecodeError:
            continue  # skip any non-JSON human-summary line
        if not isinstance(entry, dict) or "status" not in entry:
            continue
        # `ref` is JSON null for id-less entries, so `.get(k, "")` returns None,
        # not "" — coerce with `or ""` before touching .lower().
        key = (entry.get("entry_key") or "").lower()
        ref = (entry.get("ref") or "").lower()
        seen.add(key)
        status = entry.get("status", "")
        allowed = ref in allow or key in allow
        if status in BROKEN:
            (excepted if allowed else broken).append(entry)
        elif status in TRANSIENT:
            transient.append(entry)
        else:
            ok.append(entry)

    # Completeness: every reference entry in the bib must have a verify record.
    unverified = [k for k in bib_keys(args.bib) if k.lower() not in seen]

    report = render(ok, transient, excepted, broken, unverified)
    if args.report:
        with open(args.report, "w", encoding="utf-8") as fh:
            fh.write(report + "\n")
    fail = len(broken) + len(unverified)
    if args.github_output:
        with open(args.github_output, "a", encoding="utf-8") as fh:
            fh.write(
                f"broken={len(broken)}\nunverified={len(unverified)}\n"
                f"fail={fail}\nok={len(ok)}\ntransient={len(transient)}\n"
            )
    print(report)

    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
