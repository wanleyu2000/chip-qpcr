# ChIP-qPCR Site Designer

`chip-qpcr` is a Codex skill and command-line toolkit for designing motif-based ChIP-qPCR validation targets and analysing ChIP-qPCR percent-input data.

Given a transcription factor and a target gene, it can retrieve an exact JASPAR motif, download a strand-aware promoter sequence from Ensembl, scan both strands with a configurable PWM relative-score threshold, and nominate the top two non-overlapping candidate sites. It also includes an auditable Ct/Cq percent-input workflow.

## What it does

- Resolve an exact, versioned JASPAR transcription-factor matrix.
- Resolve the Ensembl canonical transcript and retrieve its promoter in transcriptional 5'-to-3' orientation.
- Scan both strands of a configurable promoter window with a PWM/PSSM.
- Rank all threshold-passing motif matches and select the top two non-overlapping candidates.
- Export motif sequences, genomic and BED coordinates, TSS-relative positions, FASTA, TSV/BED tables, and an editable SVG promoter schematic.
- Calculate ChIP-qPCR percent input from common Excel/CSV Ct exports, including input dilution correction, QC auditing, Excel/Prism exports, and plots.

## Installation as a Codex skill

Install the repository into your local Codex skills directory:

```bash
python3 "${CODEX_HOME:-$HOME/.codex}/skills/.system/skill-installer/scripts/install-skill-from-github.py" \
  --repo wanleyu2000/chip-qpcr \
  --path . \
  --name chip-qpcr
```

Start a new Codex turn after installation, then use a prompt such as:

```text
Use $chip-qpcr to identify the top two HES1 motif candidates in the human SYT7 promoter at a relative-score threshold of 0.85. Export the sequences, coordinates, and promoter schematic.
```

## Motif-based candidate design

Run the bundled script directly from a clone of this repository:

```bash
python scripts/design_chip_qpcr_candidates.py HES1 SYT7 \
  --species human \
  --threshold 0.85 \
  --output-dir hes1_syt7_candidates
```

Defaults:

- Species: `human`
- Promoter: canonical-transcript TSS `-2000/+500 bp`
- Threshold: `0.85`

Use `--matrix-id` to select a specific JASPAR matrix, `--transcript-id` for transcript-specific designs, and `--upstream` / `--downstream` to redefine the promoter window.

### Candidate-design outputs

| File | Contents |
|---|---|
| `candidate_report.md` | Method, source URLs, selected motif, top candidates, and caveats |
| `motif_selection.json` | Exact JASPAR-name candidates and selected version |
| `promoter.fasta` | Strand-aware promoter sequence in transcriptional orientation |
| `all_hits.tsv` / `all_hits.bed` | All threshold-passing motif matches |
| `top_two_sites.tsv` / `top_two_sites.bed` | Top two non-overlapping candidates |
| `promoter_schematic.svg` | Editable promoter diagram with TSS and highlighted sites |

## ChIP-qPCR percent-input analysis

```bash
python scripts/analyze_chip_qpcr.py run.xlsx \
  --input-label input \
  --input-percent 2 \
  --additional-dilution 10 \
  --exclude-flag OUTLIERRG \
  --output-dir chip_qpcr_results
```

This workflow writes an Excel audit workbook, a compact CSV summary, a Prism-ready table, and a plot of biological replicate values.

The Ct/Cq workflow requires Python 3.10+, `openpyxl`, and `matplotlib`. The motif-design script uses only the Python standard library and internet access to JASPAR and Ensembl REST APIs.

## Scientific note

A PWM relative score is a sequence-match score. It is **not** a calibrated probability of transcription-factor binding, an affinity measurement, or evidence of in vivo occupancy. Treat the selected sites as ChIP-qPCR candidates and validate them with appropriate primer design, biological replication, input/IgG controls, and assay QC.

## Citation and provenance

Each candidate-design run records the selected JASPAR matrix/version, Ensembl assembly, canonical transcript, promoter definition, source URLs, threshold, and retrieval time in `candidate_report.md`.
