#!/usr/bin/env python3
"""Generate deterministic synthetic tip expansions; never alter the source data.

Each observed tip becomes a balanced clade containing that original tip and
factor-1 synthetic relatives. Original-to-original patristic distances remain
unchanged (within floating-point rounding). Missing taxa stay missing as clades.
This helper accepts the downloaded archive's simple unquoted Newick syntax.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
import re

root = Path(__file__).resolve().parents[1]
source = root / "validation-data/carnivora-full"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--genes",type=int,nargs="+",default=[100,250,1000],help="Counts with existing sample-K.nwk files")
parser.add_argument("--factors",type=int,nargs="+",default=[1,2,4])
parser.add_argument("--output",type=Path,default=root / "validation-data/carnivora-tip-scaling")
args = parser.parse_args()
if any(n < 1 for n in args.genes+args.factors) or len(set(args.genes))!=len(args.genes) or len(set(args.factors))!=len(args.factors):
    parser.error("counts and factors must be positive and unique")
for count in args.genes:
    if not (source / f"sample-{count}.nwk").is_file():
        parser.error(f"Prepare source sample-{count}.nwk first")
destination = args.output.resolve()
destination.relative_to(root.parent)
destination.mkdir(parents=True, exist_ok=False)
seed = 20260922
leaf = re.compile(r"(?<=[(,])([^():,;\s]+):([0-9.eE+\-]+)")
provenance = {"seed": seed, "source": json.loads((source / "provenance.json").read_text()),
              "method": "balanced terminal grafts; original distances preserved; gene-specific branch jitter; synthetic clade inherits missingness",
              "datasets": []}
for genes in args.genes:
    trees = (source / f"sample-{genes}.nwk").read_text().splitlines()
    with (source / f"sample-{genes}.manifest.tsv").open() as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    assert len(trees) == len(manifest) == genes
    for factor in args.factors:
        folder = destination / f"tips-{53*factor}"
        folder.mkdir(exist_ok=True)
        output = folder / f"sample-{genes}.nwk"
        generated = []
        union = set()
        for tree, entry in zip(trees, manifest):
            assert not any(x in tree for x in "'\"[]")
            rng = random.Random(int.from_bytes(hashlib.sha256(
                f"{seed}|{entry['archive_member']}|{factor}".encode()).digest(), "big"))
            def replace(match):
                name, length = match[1], float(match[2])
                assert length >= 0 and "__synthetic" not in name
                names = [name] + [f"{name}__synthetic{i:03}" for i in range(1, factor)]
                def clade(labels, budget):
                    if len(labels) == 1:
                        # The original path length stays exact; other pendant
                        # lengths vary so these are not zero-distance duplicates.
                        pendant = budget if labels[0] == name else max(1e-10, budget * rng.uniform(.75, 1.25))
                        return f"{labels[0]}:{pendant:.17g}"
                    stem = budget * rng.uniform(.35, .55)
                    half = len(labels) // 2
                    return f"({clade(labels[:half], budget-stem)},{clade(labels[half:], budget-stem)}):{stem:.17g}"
                return clade(names, length)
            expanded = tree if factor == 1 else leaf.sub(replace, tree)
            labels = leaf.findall(expanded)
            assert len(labels) == int(entry["taxa"]) * factor
            assert len(set(x[0] for x in labels)) == len(labels)
            union.update(x[0] for x in labels)
            generated.append(expanded)
        assert len(union) == 53 * factor
        output.write_text("\n".join(generated) + "\n")
        # Identifiers sample-K:i remain aligned with the real-data manifests.
        with (folder / f"sample-{genes}.manifest.tsv").open("w") as handle:
            writer = csv.writer(handle,delimiter="\t")
            writer.writerow(["gene","archive_member","source_taxa","expanded_taxa",
                             "source_tree_sha256","expanded_tree_sha256"])
            for entry, tree in zip(manifest,generated):
                writer.writerow([entry["gene"],entry["archive_member"],entry["taxa"],
                    int(entry["taxa"])*factor,entry["tree_sha256"],hashlib.sha256(tree.encode()).hexdigest()])
        provenance["datasets"].append({"genes": genes, "taxa": 53*factor, "factor": factor,
            "file": str(output.relative_to(destination)),
            "sha256": hashlib.sha256(output.read_bytes()).hexdigest()})
(destination / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
with (destination / "datasets.tsv").open("w") as handle:
    writer = csv.DictWriter(handle,fieldnames=list(provenance["datasets"][0]),delimiter="\t")
    writer.writeheader()
    writer.writerows(provenance["datasets"])
print(json.dumps(provenance["datasets"], indent=2))
