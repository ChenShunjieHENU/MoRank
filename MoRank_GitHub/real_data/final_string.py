"""Induce candidate and reference-neighbour networks from the pinned STRING files."""
from pathlib import Path
import csv
import gzip
import os
import sys

project = Path(sys.argv[1]).resolve()
parent = Path(os.environ["MORANK_INPUT_DIR"])
out = project / "results/final_model"
out.mkdir(parents=True, exist_ok=True)
source = parent / "round20_network/sources"


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write(name, columns, records):
    with (out / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(records)


ranks = rows(parent / "round33_model_optimization/results/TCGA_exploratory_ranks.csv")
final = {r["gene"] for r in ranks if r["final_candidate"] == "TRUE"}
candidates = {r["canonical_symbol"] for r in rows(parent / "round9/results/gene_mapping_audit.csv")
              if r["input_gene"] in final and r["mapping_status"] in {"exact_approved", "unique_previous"}}
anchors = ("ABCC3", "NQO1", "IRS2", "KIT")
assert set(anchors) <= candidates
with gzip.open(source / "9606.protein.info.v12.5.txt.gz", "rt", encoding="utf-8") as stream:
    names = {r["#string_protein_id"]: r["preferred_name"] for r in csv.DictReader(stream, delimiter="\t")}
mapped = sorted(candidates & set(names.values()))
write("string_mapped.csv", ("gene",), ((gene,) for gene in mapped))

neighbors = {g: {} for g in anchors}
for kind in ("functional", "physical"):
    stem = "protein.links" if kind == "functional" else "protein.physical.links"
    induced = {}
    with gzip.open(source / f"9606.{stem}.v12.5.txt.gz", "rt", encoding="utf-8") as stream:
        next(stream)
        for line in stream:
            a, b, score = line.split()
            strength = int(score)
            if a >= b or strength < 400:
                continue
            g, h = names[a], names[b]
            if g == h:
                continue
            if g in candidates and h in candidates:
                key = tuple(sorted((g, h)))
                induced[key] = max(strength, induced.get(key, 0))
            if kind == "functional" and strength >= 700:
                for anchor, neighbor in ((g, h), (h, g)):
                    if anchor in neighbors:
                        d = neighbors[anchor]
                        d[neighbor] = max(strength, d.get(neighbor, 0))
    write(f"string_{kind}.csv", ("from", "to", "score"),
          ((*edge, strength / 1000) for edge, strength in sorted(induced.items())))
    print(kind, len(mapped), "mapped nodes", sum(s >= 700 for s in induced.values()), "edges at 0.7")

write("string_reference_neighbors.csv", ("anchor", "neighbor", "score", "candidate"),
      ((anchor, gene, strength / 1000, gene in candidates) for anchor in anchors
       for gene, strength in sorted(neighbors[anchor].items(), key=lambda p: (-p[1], p[0]))[:5]))
