# JCIM resubmission: what to edit

There is one editable LaTeX source directory:

```text
resubmission_JCIM/Files/
```

Edit only the files in that directory:

- `GPX6_manuscript.tex` and its section files;
- `GPX6_suppinfo.tex`;
- `Cover_letter_and_response_JCIM.tex`;
- `references.bib`.

Figures are in `Figures_resub/`. Do not edit anything in `build/`: it contains temporary
files and the automatically generated `latexdiff` sources used for the marked PDFs.
`submission_JCIM/` is the immutable original submission used as the comparison baseline.

From this directory, run:

```sh
make
```

This produces the five files to submit directly in `resubmission_JCIM/`:

- `GPX6_manuscript_clean.pdf`
- `GPX6_manuscript_marked.pdf`
- `GPX6_suppinfo_clean.pdf`
- `GPX6_suppinfo_marked.pdf`
- `Cover_letter_and_response_JCIM.pdf`

`make clean` removes only generated working files. `make distclean` also removes the five
final PDFs.
