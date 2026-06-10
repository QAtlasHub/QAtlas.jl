# Citations and the DOI cross-check

Every value QAtlas stores is a claim about the published literature. The claim
is only as good as its citation — so the citation must be **precise** and
**checked against the paper itself**.

## 1. Cite a precise DOI, not "Author (Year)"

For each `fetch` method and each `@register` row:

- The docstring cites the **DOI** and the **specific equation / table / theorem**:

  ```julia
  # Author, *Journal* Vol, page (Year), Eq. (3.6).  doi:10.1103/PhysRevB.89.014415
  ```

- `references = ["BibKey", ...]` in `@register` points to a key in
  `docs/references.bib`. Add the entry there if it is missing — full,
  `doi2bib`-quality (title, authors, journal, volume, year, **doi**).
- "Author (Year)" with no equation number is not acceptable: the reader must be
  able to open the paper and land on the exact result.

A value with no traceable DOI does not go into `src/`.

## 2. Download the paper and cross-check — do not trust your own derivation alone

A self-derivation (a hand calculation, or re-deriving the same integral that
`fetch` evaluates) can be **internally consistent and still wrong** by a
convention factor — spin normalisation `S` vs `s(s+1)`, a sign, per-site vs
per-bond, a rescaled coupling. A self-consistent check (autodiff identity, two
of your own routes) *cannot* catch this, because both routes share the same
convention mistake.

The only thing that settles it is the paper's own numbers. So:

1. **Fetch the source** with `doiget`:

   ```bash
   doiget fetch arXiv:1309.0940              # arXiv id
   doiget fetch 10.1103/PhysRevB.89.014415   # DOI (OA-first; falls back to metadata-only)
   ```

   If the DOI is paywalled with no OA copy (`metadata-only`), find the arXiv
   version or an equivalent OA paper; if none exists, say so explicitly in the
   PR and verify by an independent route instead (see [verification.md](verification.md)).

2. **Read the published value in the paper's conventions.** PDFs are often
   protected — rewrite with `pypdf` then extract text:

   ```bash
   python -c "import pypdf; w=pypdf.PdfWriter(); [w.add_page(p) for p in pypdf.PdfReader('paper.pdf').pages]; w.write('out.pdf')"
   pdftotext -layout out.pdf out.txt
   ```

3. **Anchor the convention.** Match one clean, unambiguous published coefficient
   to your implementation. Example (AKLT HTSE, #506): in the bilinear limit the
   per-bond specific-heat coefficient is `d₂ = r²/3` with `r = s(s+1)` — exactly
   Lohmann–Schmidt–Richter 2014's value. Reproducing that one number from the
   code pins the spin normalisation, the per-site convention, and the β power
   all at once.

4. **Beware papers with rescaled conventions.** A paper that rescales spins to
   unit norm `√(s(s+1))`, or whose PDF maths is garbled by extraction, is *not*
   a clean anchor — use it only after converting, or pick a cleaner source.

This is the literature counterpart of "tests must hit real code paths": the
reference value must come from outside your own head.

## 3. References live in `docs/references.bib`

The same `references.bib` that `doiget verify` checks in CI is rendered into the
atlas bibliography. Keep keys stable; reuse an existing key rather than adding a
near-duplicate.
