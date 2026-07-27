# Input schema and mapping

## Minimum tidy schema

| Canonical field | Common aliases | Meaning |
|---|---|---|
| `well` | Well Position, Well | Plate well such as `C5` |
| `sample` | Sample Name, Sample | `input`, target IP, IgG, or control |
| `target` | Target Name, Assay | Primer/amplicon identity |
| `ct` | CT, Ct, Cq | Quantification cycle |

Optional fields:

| Canonical field | Common aliases | Handling |
|---|---|---|
| `omit` | Omit, Exclude | Truthy values are always excluded |
| `OUTLIERRG` | OUTLIERRG | Retain in audit; exclude only when requested |
| `HIGHSD` | HIGHSD | Retain in audit; exclude only when requested |
| `THOLDFAIL` | THOLDFAIL | Retain in audit; exclude only when requested |
| `replicate` | Replicate, Bio Rep | Use directly when present |

For CSV input, the header must be the first non-empty row. For Excel input, the script scans sheets and rows for the required fields.

## Default plate-layout rule

The default replicate ID is the alphabetic plate-row prefix of `well`. Within each replicate and target:

- All wells whose sample exactly matches `--input-label` are technical input wells.
- Every other sample label becomes a separate IP/control group.
- Ct values are averaged within the input wells and within each sample group's wells.

Use `--replicate-regex` with one capture group to override this. The expression is matched against `well`, then `sample`; for example, `--replicate-regex 'Bio(\\d+)'`.

## Input dilution

Pass `--input-percent P` for the percentage of total chromatin retained as input. If that aliquot was then diluted `D`-fold before qPCR, the effective fraction is:

`(P / 100) / D`

Examples:

- 10% input, no further dilution: effective fraction `0.10`, correction `log2(10)`.
- 2% input diluted 10-fold: effective fraction `0.002`, correction `log2(500)`.

Do not encode pipetted qPCR reaction volume as `additional-dilution` when input and IP received the same reaction-volume treatment.

## Pairing and aggregation

The primary calculation uses the mean Ct of technical input wells and mean Ct of technical IP wells for each biological replicate. A well-level `%Input` column is diagnostic only and pairs sorted input/IP wells by position where possible; do not use those paired technical values as biological `n`.
