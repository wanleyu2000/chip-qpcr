# Candidate-design reference

## Inputs and defaults

| Input | Default | Meaning |
|---|---:|---|
| Species | `human` | Ensembl species alias, e.g. `mouse` |
| Assembly | Ensembl current assembly | Record in every output |
| Promoter | `-2000/+500 bp` | Relative to canonical transcript TSS |
| Threshold | `0.85` | JASPAR PWM relative score threshold |
| Selected hits | 2 | Highest-scoring non-overlapping hits |

Use `--upstream` and `--downstream` to change the promoter definition. For a transcript-specific target, use `--transcript-id` rather than relying on the gene's canonical transcript.

## Data sources

- JASPAR REST `matrix/?search=` discovers candidate matrices; the script retains only exact name matches.
- JASPAR REST `matrix/<ID>/` provides the count PFM and versioned matrix metadata.
- Ensembl REST `lookup/symbol/<species>/<gene>?expand=1` resolves the gene and canonical transcript.
- Ensembl REST `lookup/id/<transcript>` provides the transcript TSS and strand.
- Ensembl REST `sequence/region/<species>/<region>` returns the oriented promoter sequence.

All source URLs and retrieval time are written to `candidate_report.md`.

## Score and coordinates

For each PWM position, the script computes a log-odds score using uniform background frequencies and a small pseudocount. It normalizes the raw score between the PWM's theoretical minimum and maximum:

`relative_score = (raw_score - minimum_PWM_score) / (maximum_PWM_score - minimum_PWM_score)`

`relative_score >= threshold` selects a sequence match. It is not a calibrated probability, affinity, occupancy measurement, or p-value.

The FASTA sequence is always oriented 5′→3′ in the transcript direction. Thus its first base is `-upstream` relative to TSS, including for genes on the minus strand. BED files use 0-based half-open genomic coordinates. The TSV reports both genomic and TSS-relative coordinates.

## Decision rules

- Select the newest version of the exact JASPAR protein-name match by default.
- If several biologically distinct exact matrices remain, report them and select only after the user chooses or use `--matrix-id`.
- Select two highest-scoring **non-overlapping** hits for primer-target candidates; retain all threshold-passing hits for transparency.
- When no hit passes the threshold, lower the threshold only at user direction and clearly report the change.
- If protein names are aliases, complexes, antibody epitopes, or post-translational states, resolve a JASPAR matrix ID before scanning.

## Example

```bash
python scripts/design_chip_qpcr_candidates.py HES1 SYT7 \
  --species human --threshold 0.85 --upstream 2000 --downstream 500 \
  --output-dir hes1_syt7_candidates
```
