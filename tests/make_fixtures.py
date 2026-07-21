#!/usr/bin/env python3
"""Generate deterministic fixture data for the R analysis tests.

Uses a fixed seed so expected results are stable. The counts are synthetic and
exist only to exercise the tools; they are not biological data.
"""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path

FIXTURES = Path(__file__).resolve().parent / "fixtures"
# Real symbols so the ENTREZ mapping in enrichment tests has something to match.
GENES = [
    "TP53", "EGFR", "MYC", "AKT1", "PTEN", "KRAS", "BRAF", "PIK3CA", "CDKN2A", "RB1",
    "IL6", "TNF", "IL1B", "CXCL8", "CCL2", "NFKB1", "STAT3", "JAK2", "SOCS3", "IL10",
    "COL1A1", "COL3A1", "FN1", "MMP9", "MMP2", "TIMP1", "TGFB1", "SMAD3", "ACTA2", "VIM",
    "HIF1A", "VEGFA", "SLC2A1", "LDHA", "PKM", "HK2", "PGK1", "ENO1", "ALDOA", "GAPDH",
    "CDK1", "CCNB1", "MKI67", "PCNA", "TOP2A", "AURKA", "PLK1", "BUB1", "CDC20", "CCNA2",
]
UP = set(GENES[10:20])
DOWN = set(GENES[20:30])


def write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    FIXTURES.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)

    n_per_group = 6
    samples = [f"CTRL{i + 1}" for i in range(n_per_group)] + [f"DIS{i + 1}" for i in range(n_per_group)]
    groups = ["control"] * n_per_group + ["disease"] * n_per_group

    counts_rows = []
    for gene in GENES:
        row = [gene]
        base = rng.uniform(200, 2000)
        for group in groups:
            effect = 1.0
            if group == "disease":
                if gene in UP:
                    effect = 4.0
                elif gene in DOWN:
                    effect = 0.25
            mean = base * effect
            # Negative-binomial-ish noise so dispersion estimation has something to do.
            value = max(0, int(rng.gauss(mean, mean * 0.25)))
            row.append(value)
        counts_rows.append(row)
    write_csv(FIXTURES / "counts.csv", ["gene", *samples], counts_rows)

    # Same data log2-transformed: what a user pastes in by mistake. DESeq2/edgeR
    # must refuse this rather than silently returning wrong results.
    logged_rows = [[row[0], *[round(math.log2(value + 1), 4) for value in row[1:]]] for row in counts_rows]
    write_csv(FIXTURES / "logged.csv", ["gene", *samples], logged_rows)

    coldata_rows = [
        [sample, group, "M" if index % 2 else "F", 40 + (index * 3) % 30, f"batch{index % 2 + 1}"]
        for index, (sample, group) in enumerate(zip(samples, groups))
    ]
    write_csv(FIXTURES / "coldata.csv", ["sample", "group", "sex", "age", "batch"], coldata_rows)

    # Survival fixture: risk driven by the first three genes. Effects are strong
    # enough that LASSO retains them under the conservative lambda.1se rule.
    surv_rows = []
    n_subjects = 200
    for i in range(n_subjects):
        expr = [rng.gauss(0, 1) for _ in range(8)]
        linear = 1.8 * expr[0] + 1.4 * expr[1] - 1.2 * expr[2]
        scale = math.exp(-linear / 2)
        time = rng.expovariate(1 / (24 * scale))
        censor = rng.expovariate(1 / 40)
        observed = min(time, censor)
        event = 1 if time <= censor else 0
        surv_rows.append([f"S{i + 1}", round(observed, 3), event, *[round(v, 4) for v in expr]])
    write_csv(
        FIXTURES / "survival.csv",
        ["sample", "time", "event", *[f"gene{i + 1}" for i in range(8)]],
        surv_rows,
    )

    events = sum(row[2] for row in surv_rows)
    print(f"Wrote fixtures to {FIXTURES}")
    print(f"  counts.csv   {len(GENES)} genes x {len(samples)} samples")
    print(f"  coldata.csv  {len(samples)} samples, 2 groups")
    print(f"  survival.csv {n_subjects} subjects, {events} events")


if __name__ == "__main__":
    main()
