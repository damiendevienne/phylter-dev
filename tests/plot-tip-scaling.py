#!/usr/bin/env python3
"""Render standalone figures from measured scaling results (developer tool)."""
import argparse
import json
import os
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("results",type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
results = args.results.resolve()
results.relative_to(root.parent)
os.environ["MPLCONFIGDIR"] = str(results / "plot-cache")
os.environ["XDG_CACHE_HOME"] = str(results / "plot-cache")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

data = json.loads((results / "summary.json").read_text())
metadata = json.loads((results / "environment.json").read_text())
initial = metadata.get("initial_only",False)
repeats = metadata["cpp_repeats"]
genes = sorted(set(x["genes"] for x in data))
taxa = sorted(set(x["taxa"] for x in data))
lookup = {(x["genes"],x["taxa"]):x for x in data}
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"svg.fonttype":"none"})
fig, axes = plt.subplots(1,2,figsize=(11,5.4),layout="constrained")
for ax,field,title,unit in zip(axes,["elapsed_seconds","peak_rss_mib"],
                              ["Elapsed time","Peak resident memory"],["seconds","MiB"]):
    values = np.array([[lookup.get((g,n),{}).get("C++",{}).get(field,float("nan")) for n in taxa] for g in genes])
    norm = LogNorm(vmin=np.nanmin(values),vmax=np.nanmax(values))
    im = ax.imshow(values,norm=norm,cmap="viridis",aspect="auto")
    ax.set(xticks=range(len(taxa)),xticklabels=taxa,yticks=range(len(genes)),yticklabels=[f"{g:,}" for g in genes],
           xlabel="Taxa in dataset (union)",ylabel="Genes",title=title)
    for i,g in enumerate(genes):
        for j,n in enumerate(taxa):
            v = values[i,j]
            if not np.isfinite(v):
                label,color = "not measured","black"
            else:
                label = f"{v:.2f} s" if unit=="seconds" else f"{v:,.1f} MiB"
                if unit=="seconds":
                    states = lookup[g,n]["C++"]["accepted_states"]
                    label += f"\n{states:g} states"
                color = "white" if norm(v)<.6 else "black"
            ax.text(j,i,label,ha="center",va="center",color=color,fontsize=11)
    fig.colorbar(im,ax=ax,label=f"{unit} (log color scale)",shrink=.85)
phase = "Initial pass" if initial else "Full filtering"
fig.suptitle(f"Phylter C++: genes × taxa\n{phase} · 1 thread · medians of {repeats} runs",fontsize=15)
fig.supxlabel("53 taxa: real Carnivora; 106/212: synthetic terminal clades. States include the initial state.",fontsize=10)
folder = root / "figures"
folder.mkdir(exist_ok=True)
for extension in ("png","svg","pdf"):
    stem = "tip-scaling-initial" if initial else "tip-scaling"
    path = folder / f"{stem}.{extension}"
    fig.savefig(path,dpi=170)
    print(path)
