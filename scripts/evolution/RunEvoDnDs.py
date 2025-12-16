# scripts/EvoDnDsModule/plot.py
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scripts.utils.EvoScoring import scoring
from scripts.evolution.whole_dnds import run_dnds_parallel
from pathlib import Path
import random
import glob
import os


def load_evo_scores(file_list: str) -> dict[str, list[float]]:
    """
    Load Evo scores from output CSV files (e.g., geneX_NLL_score.csv).
    File format: first column = "name", other columns = scores (usually one).
    Returns {scripts: [scores...]}.
    """
    out = {}
    for csv in file_list:
        df = pd.read_csv(csv)
        gene = Path(csv).stem.replace("_NLL_score", "")
        gene_cols = [c for c in df.columns if c != "name"]
        if not gene_cols:
            continue
        vals = pd.to_numeric(df[gene_cols[0]], errors="coerce").dropna().to_list()
        out[gene] = vals
    return out

def load_dnds_m0_vertical(results):
    m0_map, sig_map = {}, {}
    for item in results:
        if len(item) == 2:
            gene, csv_path = item
        else:
            gene, _, csv_path = item
        p = Path(csv_path)
        if not p.exists():
            print(f"Result file not found: {p}")
            continue

        m0_val = None
        sig = ""
        # Robust line-by-line parsing (no dependency on pandas)
        for raw in p.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = raw.strip()
            if not line or "," not in line:
                continue
            key, val = line.split(",", 1)
            key = key.strip().lower()
            val = val.strip()
            if key == "m0_omega":
                try:
                    m0_val = float(val)
                except ValueError:
                    m0_val = None
            elif key == "lrt_sig":
                sig = val 

        if m0_val is not None:
            m0_map[gene] = m0_val
        sig_map[gene] = sig
    return m0_map, sig_map

def plot_evo_vs_m0_bars(
    evo_map: dict[str, list[float]],
    m0_map: dict[str, float],
    sig_map: dict[str, str],
    out_png: str,
    title: str = "Evo and dN/dS Scores",
    dpi: int = 300,
    order: str = "alpha",   
):


    all_genes = sorted(set(evo_map) & set(m0_map))
    genes = []
    for g in all_genes:
        xs = [v for v in (evo_map.get(g) or []) if v is not None and np.isfinite(v)]
        m0 = m0_map.get(g, None)
        if xs and (m0 is not None) and np.isfinite(m0):
            genes.append(g)
    if not genes:
        raise ValueError("No genes with both valid Evo and M0_omega values.")


    if order == "evo":
        genes.sort(key=lambda g: np.mean([v for v in evo_map[g] if v is not None and np.isfinite(v)]))
    elif order == "m0":
        genes.sort(key=lambda g: float(m0_map[g]))
    else:
        genes.sort()


    evo_mean, evo_std, m0_vals = [], [], []
    for g in genes:
        xs = [v for v in (evo_map.get(g) or []) if v is not None and np.isfinite(v)]
        mu = np.mean(xs)
        sd = np.std(xs, ddof=1) if len(xs) > 1 else 0.0
        evo_mean.append(mu)
        evo_std.append(sd)
        m0_vals.append(float(m0_map[g]))

    evo_mean = np.array(evo_mean, dtype=float)
    evo_std  = np.array(evo_std, dtype=float)
    m0_vals  = np.array(m0_vals, dtype=float)


    n = len(genes)
    fig_w = max(8, 0.55 * n + 3)
    fig_h = 4.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    centers = np.arange(n, dtype=float)


    width_outer = 0.4
    width_inner = 0.3


    bars_evo = ax.bar(
        centers, evo_mean, width=width_outer,
        # facecolor="none", edgecolor="#F0868C", linewidth=1.5,
        facecolor="none", edgecolor="#799ad5", linewidth=1.5,
        label="Evo (mean±sd)", zorder=2
    )
    ax.errorbar(
        centers, evo_mean, yerr=evo_std,
        fmt='none', ecolor="#799ad5", elinewidth=1.1, capsize=3, zorder=3
    )

                                                                                                                                                                      

    bars_m0 = ax.bar(
        centers, m0_vals, width=width_inner,
        # color="#92cbea", edgecolor="#92cbea", linewidth=1.0,
        color="#F3A2A2", linewidth=1.0,
        label="dN/dS (M0 ω)", zorder=1
    )

    for i, g in enumerate(genes):
        star = (sig_map.get(g) or "").strip()
        if not star or star.lower() == "ns":
            continue
        top = max(evo_mean[i] + evo_std[i], m0_vals[i])
        y_text = top * 1.01 if top > 0 else 0.02
        ax.text(centers[i], y_text, star, ha="center", va="bottom",
                fontsize=10, color="black", zorder=4)


    ax.set_xlim(-0.7, n - 0.3)
    ax.set_xticks(centers)
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_ylabel("Score / ω (M0)")
    if title:
        ax.set_title(title)
    ax.axhline(y=1, color="#F3A2A2", linestyle="--", linewidth=1.2)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)



    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] Saved: {out_png}")
    return out_png
    
def RunEvoDnDs(fasta_input, output_dir, fasta_tree_map=None):
    # outpout = "/home/shiyi/Output_dir/evo_dnds_compare"
    mapping = {}
    if fasta_tree_map is not None:
        df = pd.read_csv(fasta_tree_map, sep="\t")
        for _, row in df.iterrows():
            mapping[str(Path(row["fasta"]))] = None if pd.isna(row["tree"]) else str(row["tree"])
    # evo_output = Path(output_dir) / 'evo_output'
    out_png = Path(output_dir) / 'Evo_dNdS.png'

    dnds_record,evo_results = run_dnds_parallel(fasta_input, output_dir, mapping=mapping)
    print(evo_results)
    evo_map = load_evo_scores(evo_results)
    m0_map, sig_map = load_dnds_m0_vertical(dnds_record)
    return plot_evo_vs_m0_bars(evo_map, m0_map,sig_map, out_png=out_png)

if __name__ == '__main__':
    fasta_input = '/Users/sy/Gene2Struct-main/scripts/EvoDnDsModule/test_data'
    output_dir = ('/Users/sy/Gene2Struct-main/scripts/EvoDnDsModule/test_out')
    RunEvoDnDs(fasta_input, output_dir)