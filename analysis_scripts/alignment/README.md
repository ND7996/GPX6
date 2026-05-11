# GPX6 human/mouse alignment input

UniProt FASTA files downloaded for GPX6_HUMAN (P59796) and GPX6_MOUSE (Q91WR8).
Clustal Omega was run from the conda environment `evb`.
`render_gpx6_alignment.py` renders the alignment PNG used as Figure S11 in the
current Supporting Information and maps
the EVB mature/PDB residue numbering onto the full UniProt alignment using a
+24 offset.

## Reproduction commands

Run from the repository root:

```sh
cat alignment/P59796_GPX6_HUMAN.fasta \
    alignment/Q91WR8_GPX6_MOUSE.fasta \
    > alignment/gpx6_human_mouse_input.fasta

conda run -n evb clustalo \
    -i alignment/gpx6_human_mouse_input.fasta \
    -o alignment/gpx6_human_mouse_clustalo.aln \
    --outfmt=clu \
    --force

conda run -n evb clustalo \
    -i alignment/gpx6_human_mouse_input.fasta \
    -o alignment/gpx6_human_mouse_clustalo.fasta \
    --outfmt=fasta \
    --force

python3 alignment/render_gpx6_alignment.py
```
