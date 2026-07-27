---
name: chip-qpcr
description: Analyze ChIP-qPCR Ct/Cq exports with the percent-input method, including input-dilution correction, technical-replicate aggregation, instrument-flag auditing, outlier-aware summaries, plots, and Prism-ready tables. Use when Codex receives ChIP-qPCR Excel/CSV results, needs to calculate or verify percent input, compare IP/IgG/control groups, troubleshoot plate layouts or QC flags, or produce reproducible analysis deliverables.
---

# ChIP-qPCR analysis

Analyze results reproducibly and preserve an audit trail. Treat instrument flags as review signals, not automatic biological conclusions.

## Workflow

1. Inspect the source without modifying it. Identify the results sheet, well, sample, target, Ct/Cq, omit, and QC-flag columns.
2. Confirm the experimental mapping:
   - Input sample label.
   - IP, IgG, and other comparison labels.
   - Input fraction saved before any additional dilution.
   - Additional dilution applied to the input aliquot before qPCR.
   - Biological replicate identifier and technical-replicate layout.
3. If any of these materially affect the result and cannot be inferred safely, ask the user. Never silently assume 1%, 2%, or 10% input.
4. Run `scripts/analyze_chip_qpcr.py` for standard plate exports. Use `--help` to inspect all options.
5. Review the generated well audit before reporting the filtered result. Report exclusions with well IDs, flags, and the rule used.
6. Verify the formula and at least one replicate manually:

   `adjusted_input_Ct = mean_input_Ct - log2(1 / effective_input_fraction)`

   `%Input = 100 × 2^(adjusted_input_Ct - mean_IP_Ct)`

   `effective_input_fraction = saved_input_fraction / additional_dilution`

7. Summarize biological replicate values with mean, sample SD, SEM, and `n`. Do not treat technical wells as independent biological replicates.
8. State limitations: percent input is enrichment, not proof of direct binding; primer efficiency, specificity, chromatin quality, antibody quality, and background controls affect interpretation.

## Standard command

```bash
python scripts/analyze_chip_qpcr.py run.xlsx \
  --input-label input \
  --input-percent 2 \
  --additional-dilution 10 \
  --exclude-flag OUTLIERRG \
  --output-dir chip_qpcr_results
```

The script reads `.xlsx` or `.csv`, auto-detects common QuantStudio-style headers, groups wells by plate row by default, and writes:

- `chip_qpcr_analysis.xlsx`: summary, replicate calculations, well-level audit, and Prism-ready data.
- `chip_qpcr_summary.csv`: compact group statistics.
- `chip_qpcr_plot.png`: individual biological replicate points with mean ± SD.

For a tidy CSV, provide columns equivalent to `Well Position`, `Sample Name`, `Target Name`, and `CT`. Use `--replicate-regex` when the replicate ID is encoded in the well or sample name instead of the plate row.

## QC and interpretation rules

- Exclude wells explicitly marked by the instrument/user omit column.
- Exclude only flags named with `--exclude-flag`; keep all flagged wells in the audit.
- Prefer predeclared rules. If exploring exclusions after seeing the data, show both all-well and filtered analyses.
- Investigate technical replicate spread, non-amplification, multi-peak melt curves, implausible Ct, and inconsistent input Ct before aggregation.
- Keep IgG and positive/negative loci as separate sample groups; never merge them into the target IP.
- Use sample SD (`n-1`) for biological replicates. Leave SD/SEM blank when `n < 2`.
- Do not invent significance tests. Paired biological designs require paired tests; multiple loci/groups require an explicit multiplicity plan.

## Nonstandard inputs

Read [references/input-schema.md](references/input-schema.md) when column detection fails or the plate layout is not row-based. Normalize a copy to the documented tidy schema, then run the script on that copy. Preserve the original file unchanged.

## Dependencies

Require Python 3.10+, `openpyxl`, and `matplotlib`. If a dependency is unavailable, explain the missing package and provide the exact installation command; do not alter the user's environment without permission.
