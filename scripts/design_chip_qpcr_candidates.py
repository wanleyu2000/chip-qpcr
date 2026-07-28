#!/usr/bin/env python3
"""Find JASPAR motif candidates in an Ensembl canonical-transcript promoter."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import re
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen
from xml.sax.saxutils import escape

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
JASPAR = "https://jaspar.elixir.no/api/v1"
ENSEMBL = "https://rest.ensembl.org"


def api_json(url: str) -> dict:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": "chip-qpcr-skill/1.0"})
    try:
        with urlopen(request, timeout=45) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")[:300]
        raise RuntimeError(f"HTTP {exc.code} from {url}: {detail}") from exc
    except URLError as exc:
        raise RuntimeError(f"Cannot reach {url}: {exc.reason}") from exc


def api_text(url: str) -> str:
    request = Request(url, headers={"Accept": "text/plain", "User-Agent": "chip-qpcr-skill/1.0"})
    try:
        with urlopen(request, timeout=45) as response:
            return response.read().decode("utf-8").strip().upper()
    except (HTTPError, URLError) as exc:
        raise RuntimeError(f"Cannot retrieve sequence from {url}: {exc}") from exc


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def exact_jaspar_candidates(protein: str) -> list[dict]:
    url = f"{JASPAR}/matrix/?search={quote(protein)}&page_size=100"
    payload = api_json(url)
    desired = protein.casefold()
    candidates = [item for item in payload.get("results", []) if str(item.get("name", "")).casefold() == desired]
    if not candidates:
        names = sorted({str(item.get("name", "")) for item in payload.get("results", [])})
        hint = f" Nearby search results: {', '.join(names[:8])}." if names else ""
        raise ValueError(
            f"No exact JASPAR matrix name for '{protein}'. Supply --matrix-id after resolving the TF/alias.{hint}"
        )
    return candidates


def select_matrix(protein: str, matrix_id: str | None) -> tuple[dict, list[dict], str]:
    if matrix_id:
        selected_url = f"{JASPAR}/matrix/{quote(matrix_id)}/"
        selected = api_json(selected_url)
        return selected, [selected], selected_url
    candidates = exact_jaspar_candidates(protein)
    candidates.sort(key=lambda item: (int(item.get("version", 0)), item["matrix_id"]), reverse=True)
    selected_url = candidates[0]["url"]
    return api_json(selected_url), candidates, selected_url


def resolve_transcript(species: str, gene: str, transcript_id: str | None) -> tuple[dict, dict, str, str]:
    if transcript_id:
        transcript_url = f"{ENSEMBL}/lookup/id/{quote(transcript_id)}?expand=0"
        transcript = api_json(transcript_url)
        gene_url = f"{ENSEMBL}/lookup/id/{quote(str(transcript.get('Parent', '')))}?expand=0"
        gene_info = api_json(gene_url) if transcript.get("Parent") else {}
        return gene_info, transcript, gene_url, transcript_url
    gene_url = f"{ENSEMBL}/lookup/symbol/{quote(species)}/{quote(gene)}?expand=1"
    gene_info = api_json(gene_url)
    canonical = str(gene_info.get("canonical_transcript") or "").split(".")[0]
    if not canonical:
        transcripts = gene_info.get("Transcript", [])
        if not transcripts:
            raise ValueError(f"No transcript returned for {gene} in {species}")
        canonical = str(transcripts[0]["id"]).split(".")[0]
    transcript_url = f"{ENSEMBL}/lookup/id/{quote(canonical)}?expand=0"
    return gene_info, api_json(transcript_url), gene_url, transcript_url


def promoter_context(species: str, transcript: dict, upstream: int, downstream: int) -> tuple[dict, str, str]:
    chrom = str(transcript["seq_region_name"])
    strand = int(transcript["strand"])
    tss = int(transcript["start"] if strand == 1 else transcript["end"])
    region_start = max(1, tss - upstream if strand == 1 else tss - downstream)
    region_end = tss + downstream if strand == 1 else tss + upstream
    region = f"{chrom}:{region_start}..{region_end}:{strand}"
    sequence_url = f"{ENSEMBL}/sequence/region/{quote(species)}/{region}"
    sequence = api_text(sequence_url)
    expected = region_end - region_start + 1
    if len(sequence) != expected:
        raise RuntimeError(f"Promoter length mismatch: expected {expected}, received {len(sequence)}")
    return {
        "chrom": chrom,
        "strand": strand,
        "tss_1based": tss,
        "region_start_1based": region_start,
        "region_end_1based": region_end,
        "upstream": upstream,
        "downstream": downstream,
        "assembly": transcript.get("assembly_name", "unknown"),
    }, sequence, sequence_url


def pssm_from_pfm(pfm: dict, pseudocount: float = 1e-4) -> list[list[float]]:
    widths = {len(pfm[base]) for base in BASES}
    if len(widths) != 1:
        raise ValueError("JASPAR PFM has inconsistent base-array lengths")
    rows = []
    for index in range(widths.pop()):
        counts = [float(pfm[base][index]) for base in BASES]
        total = sum(counts) + 4 * pseudocount
        rows.append([math.log2(((count + pseudocount) / total) / 0.25) for count in counts])
    return rows


def score_word(word: str, pssm: list[list[float]]) -> float | None:
    score = 0.0
    for base, row in zip(word, pssm):
        index = BASES.find(base)
        if index < 0:
            return None
        score += row[index]
    return score


def scan_hits(sequence: str, matrix: dict, context: dict, threshold: float) -> list[dict]:
    pssm = pssm_from_pfm(matrix["pfm"])
    width = len(pssm)
    minimum = sum(min(row) for row in pssm)
    maximum = sum(max(row) for row in pssm)
    hits = []
    for local_start in range(len(sequence) - width + 1):
        genomic_oriented = sequence[local_start : local_start + width]
        for motif_strand, motif_sequence in (("+", genomic_oriented), ("-", reverse_complement(genomic_oriented))):
            raw_score = score_word(motif_sequence, pssm)
            if raw_score is None:
                continue
            relative = (raw_score - minimum) / (maximum - minimum) if maximum != minimum else 0.0
            if relative < threshold:
                continue
            if context["strand"] == 1:
                genomic_start = context["region_start_1based"] - 1 + local_start
                genomic_end = genomic_start + width
                genomic_strand = motif_strand
            else:
                genomic_end = context["region_end_1based"] - local_start
                genomic_start = genomic_end - width
                genomic_strand = "-" if motif_strand == "+" else "+"
            hits.append(
                {
                    "chrom": "chr" + context["chrom"] if not context["chrom"].startswith("chr") else context["chrom"],
                    "bed_start": genomic_start,
                    "bed_end": genomic_end,
                    "tss_relative_start": local_start - context["upstream"],
                    "tss_relative_end": local_start + width - context["upstream"],
                    "motif_strand_transcript": motif_strand,
                    "genomic_strand": genomic_strand,
                    "relative_score": relative,
                    "raw_score": raw_score,
                    "motif_sequence_5to3": motif_sequence,
                    "genomic_sequence_transcript_oriented": genomic_oriented,
                    "matrix_id": matrix["matrix_id"],
                    "matrix_name": matrix["name"],
                    "matrix_version": matrix["version"],
                }
            )
    return sorted(hits, key=lambda item: (-item["relative_score"], item["bed_start"], item["genomic_strand"]))


def top_nonoverlapping(hits: list[dict], count: int = 2) -> list[dict]:
    chosen = []
    for hit in hits:
        if all(hit["bed_end"] <= prior["bed_start"] or hit["bed_start"] >= prior["bed_end"] for prior in chosen):
            chosen.append(hit)
        if len(chosen) == count:
            break
    return chosen


FIELDS = [
    "rank", "chrom", "bed_start", "bed_end", "tss_relative_start", "tss_relative_end",
    "genomic_strand", "motif_strand_transcript", "relative_score", "raw_score",
    "motif_sequence_5to3", "genomic_sequence_transcript_oriented", "matrix_id", "matrix_name", "matrix_version",
]


def write_tsv(path: Path, hits: list[dict], rank: bool = False) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for index, hit in enumerate(hits, 1):
            row = dict(hit)
            row["rank"] = index if rank else ""
            writer.writerow(row)


def write_bed(path: Path, hits: list[dict]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for index, hit in enumerate(hits, 1):
            name = f"rank{index}|{hit['matrix_id']}|{hit['relative_score']:.3f}|TSS{hit['tss_relative_start']:+d}"
            score = max(0, min(1000, round(hit["relative_score"] * 1000)))
            handle.write(f"{hit['chrom']}\t{hit['bed_start']}\t{hit['bed_end']}\t{name}\t{score}\t{hit['genomic_strand']}\n")


def svg_schematic(path: Path, context: dict, hits: list[dict], top: list[dict], gene: str, protein: str) -> None:
    width, height, left, right = 1600, 300, 130, 80
    plot = width - left - right
    total = context["upstream"] + context["downstream"] + 1
    def x(relative: int) -> float:
        return left + (relative + context["upstream"]) / total * plot
    axis_y = 160
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1f2937}.small{font-size:16px}.label{font-size:18px;font-weight:bold}</style>',
        f'<text x="40" y="45" class="label">{escape(protein)} motif candidates at {escape(gene)} promoter</text>',
        f'<text x="40" y="72" class="small">Canonical-transcript promoter; {context["assembly"]}; TSS-oriented 5′→3′</text>',
        f'<line x1="{left}" y1="{axis_y}" x2="{width-right}" y2="{axis_y}" stroke="#475569" stroke-width="4"/>',
    ]
    for relative, label in [(-context["upstream"], f"-{context['upstream']:,}"), (0, "TSS"), (context["downstream"], f"+{context['downstream']:,}")]:
        xx = x(relative)
        parts.append(f'<line x1="{xx:.1f}" y1="{axis_y-18}" x2="{xx:.1f}" y2="{axis_y+28}" stroke="#334155" stroke-width="2"/>')
        parts.append(f'<text x="{xx:.1f}" y="{axis_y+55}" text-anchor="middle" class="small">{label}</text>')
    for hit in hits[:80]:
        xx = x(hit["tss_relative_start"])
        color = "#d946ef" if hit in top else "#94a3b8"
        h = 42 if hit in top else 18
        parts.append(f'<rect x="{xx:.1f}" y="{axis_y-h/2:.1f}" width="{max(3, x(hit["tss_relative_end"])-xx):.1f}" height="{h}" rx="2" fill="{color}"/>')
    for index, hit in enumerate(top, 1):
        xx = x(hit["tss_relative_start"])
        label_y = 112 if index == 1 else 220
        arrow_y = axis_y - 24 if index == 1 else axis_y + 28
        label_x = xx - 14 if abs(hit["tss_relative_start"]) < 250 else xx
        anchor = "end" if label_x != xx else "middle"
        parts.append(f'<line x1="{xx:.1f}" y1="{arrow_y}" x2="{xx:.1f}" y2="{label_y + (-8 if index == 1 else 8)}" stroke="#a21caf" stroke-width="2"/>')
        parts.append(f'<text x="{label_x:.1f}" y="{label_y}" text-anchor="{anchor}" class="label">Site {index}: {hit["tss_relative_start"]:+d} ({hit["relative_score"]:.3f})</text>')
    parts.append('<text x="40" y="280" class="small">Grey: all threshold-passing motif matches (up to 80 shown). Magenta: top two non-overlapping candidates.</text>')
    parts.append('</svg>\n')
    path.write_text("\n".join(parts), encoding="utf-8")


def write_report(path: Path, protein: str, gene: str, matrix: dict, context: dict, top: list[dict], sources: dict, threshold: float) -> None:
    lines = [
        f"# {protein} → {gene} ChIP-qPCR motif candidates", "",
        "## Method", "",
        f"- JASPAR matrix: `{matrix['matrix_id']}` ({matrix['name']}, version {matrix['version']}).",
        f"- Assembly: `{context['assembly']}`; canonical transcript TSS: `{context['chrom']}:{context['tss_1based']}` ({'+' if context['strand'] == 1 else '-'} strand).",
        f"- Promoter: `{context['upstream']}` bp upstream to `{context['downstream']}` bp downstream of the TSS, oriented 5′→3′ in transcription direction.",
        f"- Scan threshold: relative PWM score ≥ `{threshold:.3f}` on both strands.",
        "- Relative score is a sequence-match metric, not binding probability, affinity, or evidence of occupancy.", "",
        "## Top two non-overlapping motif candidates", "",
        "| Rank | Position vs TSS | Genomic coordinate (BED) | Strand | Relative score | Motif-oriented sequence |",
        "|---:|---:|---|:---:|---:|---|",
    ]
    if top:
        for index, hit in enumerate(top, 1):
            lines.append(f"| {index} | {hit['tss_relative_start']:+d} to {hit['tss_relative_end']:+d} | {hit['chrom']}:{hit['bed_start']}-{hit['bed_end']} | {hit['genomic_strand']} | {hit['relative_score']:.3f} | `{hit['motif_sequence_5to3']}` |")
    else:
        lines.append("| — | No hit passed the selected threshold | — | — | — | — |")
    lines += ["", "## Source URLs", ""]
    lines += [f"- {label}: {url}" for label, url in sources.items()]
    lines += ["", f"Retrieved: {dt.datetime.now(dt.timezone.utc).isoformat()}", "", "## Validation note", "", "Use these as motif-based ChIP-qPCR candidates. Validate primer specificity, amplicon placement, biological replication, IgG/input controls, and chromatin/antibody performance before inferring occupancy."]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protein", help="Exact JASPAR TF name, unless --matrix-id is supplied")
    parser.add_argument("gene", help="Ensembl gene symbol")
    parser.add_argument("--species", default="human", help="Ensembl species alias (default: human)")
    parser.add_argument("--matrix-id", help="Versioned JASPAR matrix ID, e.g. MA1099.3")
    parser.add_argument("--transcript-id", help="Ensembl transcript ID; otherwise use canonical transcript")
    parser.add_argument("--upstream", type=int, default=2000)
    parser.add_argument("--downstream", type=int, default=500)
    parser.add_argument("--threshold", type=float, default=0.85, help="Relative PWM score threshold in [0,1]")
    parser.add_argument("--output-dir", type=Path, default=Path("chip_qpcr_candidates"))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.upstream < 0 or args.downstream < 0:
        raise ValueError("--upstream and --downstream must be non-negative")
    if not 0 <= args.threshold <= 1:
        raise ValueError("--threshold must be between 0 and 1")
    matrix, candidates, matrix_url = select_matrix(args.protein, args.matrix_id)
    gene_info, transcript, gene_url, transcript_url = resolve_transcript(args.species, args.gene, args.transcript_id)
    context, sequence, sequence_url = promoter_context(args.species, transcript, args.upstream, args.downstream)
    hits = scan_hits(sequence, matrix, context, args.threshold)
    top = top_nonoverlapping(hits)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    fasta_header = (
        f">{args.gene}|{args.protein}|{context['assembly']}|{context['chrom']}:{context['region_start_1based']}-{context['region_end_1based']}"
        f":{'+' if context['strand'] == 1 else '-'}|TSS={context['tss_1based']}|relative=-{args.upstream}/+{args.downstream}\n"
    )
    (output / "promoter.fasta").write_text(fasta_header + sequence + "\n", encoding="utf-8")
    selection = {"selected_matrix": matrix, "exact_name_candidates": candidates, "query_protein": args.protein}
    (output / "motif_selection.json").write_text(json.dumps(selection, indent=2) + "\n", encoding="utf-8")
    write_tsv(output / "all_hits.tsv", hits)
    write_bed(output / "all_hits.bed", hits)
    write_tsv(output / "top_two_sites.tsv", top, rank=True)
    write_bed(output / "top_two_sites.bed", top)
    svg_schematic(output / "promoter_schematic.svg", context, hits, top, args.gene, args.protein)
    sources = {"JASPAR matrix": matrix_url, "Ensembl gene": gene_url, "Ensembl transcript": transcript_url, "Ensembl promoter sequence": sequence_url}
    write_report(output / "candidate_report.md", args.protein, args.gene, matrix, context, top, sources, args.threshold)
    print(f"matrix={matrix['matrix_id']} ({matrix['name']}, v{matrix['version']})")
    print(f"promoter={context['assembly']} {context['chrom']}:{context['region_start_1based']}-{context['region_end_1based']} strand={context['strand']}")
    print(f"hits={len(hits)} top_two={len(top)} output={output.resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
