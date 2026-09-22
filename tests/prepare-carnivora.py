#!/usr/bin/env python3
"""Inspect the downloaded Carnivora archive and sample without extracting it.

Developer helper, specific to this archive's unquoted, one-tree-per-file Newick.
Run from the parent repository after downloading validation-data/carnivora-full/trees.tar.xz.
No tree distance matrices are constructed, even during the full-data census.
"""
import csv
import hashlib
import json
from pathlib import Path
import random
import re
import tarfile

root = Path(__file__).resolve().parents[1] / "validation-data/carnivora-full"
archive = root / "trees.tar.xz"
seed = 20260921
records = {}
with tarfile.open(archive, "r:xz") as handle:
    for member in handle:
        if member.isdir():
            continue
        if not member.isfile() or not member.name.endswith(".treefile"):
            raise ValueError(f"Unexpected archive member: {member.name}")
        if member.size > 1_000_000 or member.name in records:
            raise ValueError(f"Oversized or duplicate tree: {member.name}")
        tree = handle.extractfile(member).read().decode("ascii").strip()
        # Fail on unsupported syntax rather than silently miscounting taxa.
        assert not any(c in tree for c in "'\"[]"), member.name
        assert tree.count(";") == 1 and tree.endswith(";"), member.name
        taxa = re.findall(r"(?<=[(,])\s*([^():,;\s]+)\s*:", tree)
        assert len(taxa) == tree.count(",") + 1 == len(set(taxa)), member.name
        assert tree.count("(") == tree.count(")"), member.name
        records[member.name] = (tree, set(taxa))

names = sorted(records)
random.Random(seed).shuffle(names)
all_taxa = sorted(set.union(*(item[1] for item in records.values())))
metadata = {
    "source_url": "https://raw.githubusercontent.com/damiendevienne/phylter-data/main/Carnivora/data/trees.tar.xz",
    "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
    "genes": len(records), "taxa": len(all_taxa), "taxon_names": all_taxa,
    "taxa_per_gene_min": min(len(item[1]) for item in records.values()),
    "taxa_per_gene_max": max(len(item[1]) for item in records.values()),
    "seed": seed,
    "selection": "Python random.Random(seed).shuffle(sorted(archive tree paths)); nested prefixes, sorted within each sample",
    "samples": [],
}
for count in (100, 250, 1000, 3000):
    selected = sorted(names[:count])
    stem = f"sample-{count}"
    path = root / (stem + ".nwk")
    path.write_text("\n".join(records[name][0] for name in selected) + "\n")
    with (root / (stem + ".manifest.tsv")).open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "archive_member", "taxa", "tree_sha256"])
        for i, name in enumerate(selected, 1):
            writer.writerow([f"{stem}:{i}", name, len(records[name][1]),
                             hashlib.sha256(records[name][0].encode()).hexdigest()])
    metadata["samples"].append({
        "file": path.name, "genes": count,
        "taxa": len(set.union(*(records[name][1] for name in selected))),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    })
with (root / "census.tsv").open("w") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["archive_member", "taxa"])
    writer.writerows((name, len(records[name][1])) for name in sorted(records))
(root / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
print(json.dumps(metadata, indent=2))
