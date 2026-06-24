# SEM data-provenance audit

## Recomputed values from the displayed SD and N

The current Table S8 uses `SEM = SD / sqrt(N)` and rounds to two decimals.

| System | Quantity | SD | N | SEM (unrounded) | SEM (table) |
|---|---:|---:|---:|---:|---:|
| Human Sec-GPX6 | activation barrier | 3.21 | 22 | 0.6844 | 0.68 |
| Human Sec-GPX6 | reaction free energy | 3.48 | 22 | 0.7419 | 0.74 |
| Human Cys-GPX6 | activation barrier | 3.63 | 14 | 0.9702 | 0.97 |
| Human Cys-GPX6 | reaction free energy | 5.00 | 14 | 1.3363 | 1.34 |
| Mouse Cys-GPX6 | activation barrier | 4.18 | 21 | 0.9122 | 0.91 |
| Mouse Cys-GPX6 | reaction free energy | 13.87 | 21 | 3.0267 | 3.03 |
| Mouse Sec-GPX6 | activation barrier | 2.43 | 22 | 0.5181 | 0.52 |
| Mouse Sec-GPX6 | reaction free energy | 2.75 | 22 | 0.5863 | 0.59 |

## Provenance status

- The path supplied as `../../RESEARCH/GPX6` does not exist from the article workspace.
- The simulation repository found at `../../STUDENTS/GPX6` is the local clone of
  `github.com/ND7996/GPX6`.
- Its README states that raw `.dcd`, `.en`, and `.log` files are stored externally and
  are not included in Git.
- The repository contains old Qtools logs for 100-replica analyses, whose statistics do
  not reproduce the four-system values in Table S8. These old and current analyses must
  not be mixed.
- The current individual replica vectors and exclusion manifest were not found locally.
  Therefore the arithmetic of Table S8 is verified, but its means, SDs, medians, N values,
  and exclusions cannot yet be independently reconstructed from the available files.
- DOI `10.34810/DATA3230` resolves to an unrelated amyloid-beta dataset and must
  not be cited for this study. The authors have identified
  `https://doi.org/10.34810/DATA3330` as the correct GPX6 Dataverse record; this
  is the DOI used in the revised manuscript.

## Required files for full certification

For each of Human Sec, Human Cys, Mouse Cys, and Mouse Sec, provide one row per attempted
replica containing: replica ID, status, objective exclusion reason, activation barrier,
reaction free energy, and mapping free energy. The final statistics and Table S8 should
then be generated from this manifest without manual selection.
