---
name: chip-qpcr
description: Design and analyze ChIP-qPCR experiments. Use when Codex needs to turn a transcription factor/protein and gene symbol into JASPAR PWM motif candidates, download a strand-aware promoter sequence, score both strands at a user-selected JASPAR relative-score threshold, select the top two candidate binding sites, produce sequence/coordinate tables and promoter schematics, design ChIP-qPCR validation targets, or calculate auditable ChIP-qPCR percent input from Ct/Cq exports.
---

# ChIP-qPCR candidate design and analysis

Use the candidate-design workflow before wet-lab ChIP-qPCR whenever the user provides a protein/TF and a gene. Keep computational prediction separate from evidence of in vivo occupancy.

## Protein + gene workflow

1. Confirm species, gene symbol, protein/TF name, genome assembly, promoter window, and JASPAR relative-score threshold. Default to human, GRCh38, canonical-transcript promoter `-2000/+500 bp`, and threshold `0.85` only when the user has not specified them.
2. Run the bundled script. It queries JASPAR, resolves the gene/canonical transcript and sequence through Ensembl REST, scans both strands, and ranks hits.

```bash
python scripts/design_chip_qpcr_candidates.py HES1 SYT7 \
  --species human --threshold 0.85 --output-dir chip_qpcr_candidates
```

3. Inspect `motif_selection.json`. Use only exact case-insensitive JASPAR-name matches. If there are multiple exact matrices, default to the highest matrix version and disclose the selected ID; use `--matrix-id` when the user specifies one. If no exact matrix exists, stop and ask for an approved matrix or TF alias—do not substitute a related protein silently.
4. Treat `relative_score` as a normalized PWM match score, **not a probability of protein binding**. Higher thresholds reduce false-positive sequence matches but can miss degenerate sites.
5. Use the top two non-overlapping hits in `top_two_sites.tsv` as motif-based candidates. Report their motif-oriented sequences, genomic/BED coordinates, TSS-relative positions, strand, matrix ID/version, and threshold.
6. Open `promoter_schematic.svg` to check that the TSS, gene direction, all hits, and highlighted top two hits are coherent. Deliver it with `candidate_report.md`, `all_hits.tsv`, and `top_two_sites.tsv`.
7. Suggest qPCR amplicons that centre on each site where feasible (typically 70–150 bp), plus a motif-negative control locus. Primer design, uniqueness, repeat masking, and genome build/assembly verification remain separate validation steps.

## Outputs

The candidate script creates only a new output directory:

- `candidate_report.md`: methods, source URLs, selected motif, and caveats.
- `motif_selection.json`: all exact JASPAR matrix candidates and the selected version.
- `promoter.fasta`: transcription-oriented promoter sequence, with precise interval metadata.
- `all_hits.tsv` and `all_hits.bed`: every threshold-passing scan result.
- `top_two_sites.tsv` and `top_two_sites.bed`: ranked non-overlapping motif candidates.
- `promoter_schematic.svg`: publication-ready, editable promoter diagram.

Read [references/candidate-design.md](references/candidate-design.md) for the coordinate model, threshold interpretation, and troubleshooting API identifiers.

## ChIP-qPCR Ct analysis

After obtaining Ct/Cq results, use `scripts/analyze_chip_qpcr.py` to compute percent input with the correct input-fraction and dilution correction. Read [references/input-schema.md](references/input-schema.md) for expected columns and QC rules.

```bash
python scripts/analyze_chip_qpcr.py run.xlsx \
  --input-label input --input-percent 2 --additional-dilution 10 \
  --exclude-flag OUTLIERRG --output-dir chip_qpcr_results
```

## Interpretation guardrails

- A motif hit is a sequence hypothesis, not a binding probability or ChIP-positive result.
- Do not call occupancy without biological evidence; include IgG/input controls and biological replication.
- Preserve JASPAR matrix version, Ensembl assembly, promoter definition, timestamp, threshold, and every exclusion in the final report.
- Do not mix human and mouse coordinates or use a gene-level TSS when a transcript-specific design is required.
