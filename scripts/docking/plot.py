import re, os
from pathlib import Path
from matplotlib import font_manager as fm
import pandas as pd
import numpy as np
from scripts.utils.font import times_new_roman
import Bio.Phylo as Phylo
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scripts.evolution.plot import root_and_force_bifurcation, set_topological_ultrametric, sort_subtrees_by_size, get_leaf_order_topdown

fm.fontManager.addfont(times_new_roman)
font_prop = fm.FontProperties(fname=times_new_roman)
plt.rcParams['font.family'] = font_prop.get_name()
plt.rcParams["mathtext.fontset"] = "stix"

def best_affinity(log_file):
    _aff_pat = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)")
    with open(log_file) as f:
        for line in f:
            m = _aff_pat.match(line)
            if m:
                return float(m.group(1))
            
    print("Could not parse log file:", log_file)
    return None
    

def build_table(docking_dir, gene_ligand_map, out_csv=None):
    docking_dir = Path(docking_dir)
    mapping = {}
    for gene, info in gene_ligand_map.items():
        roles = {}
        for role in ["substrate", "product", "cofactor"]:
            ligs = info.get(role) or []
            if ligs:
                roles[role] = ligs
        if roles:
            mapping[gene] = roles
    columns = [(gene, role, ligand)
               for gene, rmap in mapping.items()
               for role, ligands in rmap.items()
               for ligand in ligands]
    col_map = {col: idx for idx, col in enumerate(columns)}
    rows_data = {}
    species_seen = set()
    for gene_dir in docking_dir.iterdir():
        if not gene_dir.is_dir():
            continue
        gene = gene_dir.name
        if gene not in mapping:
            continue
        for run_dir in gene_dir.iterdir():
            if not run_dir.is_dir():
                continue
            parts = run_dir.name.split("__", 2)
            if len(parts) != 3:
                print(f"Directory {run_dir} does not follow naming rule: gene__species__ligand, skipped.")
                continue
            g, species, ligand = parts
            if g != gene:
                continue
            species_seen.add(species)

            role = None
            for r, lst in mapping[gene].items():
                if ligand.lower() in (l.lower() for l in lst):
                    role = r
                    break
            if role is None:
                print(f"Ligand {ligand} of {gene} not found in substrate/product mapping, skipped.")
                continue
            log_file = next(run_dir.glob("*.log"), None)
            if not log_file:
                continue
            aff = best_affinity(log_file)
            if aff is None:
                # print(f"Failed to parse affinity from {log_file}")
                continue
            col_key = (gene, role, ligand)
            if col_key not in col_map:
                print(f"{col_key} not in predefined columns, skipped.")
                continue
            col_idx = col_map[col_key]
            if species not in rows_data:
                rows_data[species] = [None] * len(columns)
            rows_data[species][col_idx] = aff

    df = pd.DataFrame.from_dict(rows_data, orient="index",columns=pd.MultiIndex.from_tuples(columns, names=["Gene", "Role", "Ligand"])).sort_index(axis=1).sort_index(axis=0)
    if out_csv is not None and out_csv.strip() != "":
        df.to_csv(out_csv, float_format="%.6g")
        print(f"Docking summary table saved to {out_csv}")
    return df


def standardize_z_score(df):
    df_standardized = df.copy()
    index_cols = [col for col in df.columns if col != 'Species' and col != 'Gene']
    for col in index_cols:
        data = df_standardized[col].astype(float)
        if data.empty or data.isnull().all():
            continue
        mean = data.mean()
        std = data.std()
        if std != 0:
            z_score = (data - mean) / std
            if col == "Cofactor Binding Energy":
                z_score = -z_score
            df_standardized[f'{col}_ZScore'] = z_score.apply(lambda x: round(x,2) if not pd.isna(x) else np.nan)
        else:
            df_standardized[f'{col}_ZScore'] = 0.0
    return df_standardized


def min_max_normalize(df):
    df_normalized = df.copy()
    index_cols = [col for col in df.columns if col != 'Species' and col != 'Gene']
    for col in index_cols:
        data = df_normalized[col].astype(float)
        if data.empty or data.isnull().all():
            continue
        min_val = data.min()
        max_val = data.max()
        if max_val != min_val:
            normalized = (data - min_val) / (max_val - min_val)
            if col == "Cofactor Binding Energy":
                normalized = 1 - normalized
            df_normalized[f'{col}_Normalized'] = normalized.apply(lambda x: round(x,2) if not pd.isna(x) else np.nan)
        else:
            df_normalized[f'{col}_Normalized'] = 0.0
    return df_normalized


def calculate_indices(raw_df):
    results = []
    species_list = raw_df.index.tolist()
    is_cofactor = 'cofactor' in raw_df.columns.get_level_values('Role')
    for gene in raw_df.columns.get_level_values(0).unique():
        for species in species_list:
            substrate_energy = raw_df.loc[species,(gene,'substrate')].iloc[0]
            product_energy = raw_df.loc[species,(gene,'product')].iloc[0]
            cofactor_energy = None
            if is_cofactor:
                cofactor_energy = raw_df.loc[species,(gene,'cofactor')].iloc[0]

            inhibition_diff = substrate_energy - product_energy if substrate_energy is not None and product_energy is not None else None
            inhibition_diff = round(inhibition_diff, 2) if inhibition_diff is not None else None
            results.append({
                'Gene': gene,
                'Species': species,
                "Cofactor Binding Energy": cofactor_energy,
                "Inhibition Difference (Substrate - Product)": inhibition_diff,
            })
    results_df = pd.DataFrame(results)
    standardize_df = min_max_normalize(results_df)
    cofactor_df = pd.DataFrame()
    col = 'Cofactor Binding Energy_Normalized'
    if is_cofactor and col in standardize_df.columns:
        cofactor_df = standardize_df.pivot(index='Species', columns='Gene', values='Cofactor Binding Energy_Normalized')
    else:
        cofactor_df = pd.DataFrame()
    inhibition_diff_df = standardize_df.pivot(index='Species', columns='Gene', values='Inhibition Difference (Substrate - Product)_Normalized')

    return cofactor_df, inhibition_diff_df, is_cofactor


def plot_single_figure(tree_path, df, title, out_path):
    if df is None or df.empty or df.dropna().empty:
        return
    fig = plt.figure(figsize=(16, 10), dpi=150)
    gs = GridSpec(1, 2, width_ratios=[1, 5], wspace=0.2)
    ax_tree = fig.add_subplot(gs[0])
    ax_heatmap = fig.add_subplot(gs[1])
    # plot tree
    tree = Phylo.read(tree_path, "newick")
    max_depth = max(c.topo_depth if hasattr(c, 'topo_depth') else 0 for c in tree.find_clades(order='postorder'))
    unit_length = 1.0 / (max_depth or 1.0)
    aligned_tree = set_topological_ultrametric(tree, unit_length)
    sort_subtrees_by_size(aligned_tree.root)
    leaf_labels = get_leaf_order_topdown(aligned_tree)
    Phylo.draw(aligned_tree, axes=ax_tree, do_show=False,
            branch_labels=None,
            show_confidence=False,
            label_func=lambda x: x.name if x.name else "")
    ax_tree.axis('off')
    # adjust tree position
    fig.canvas.draw()
    ys = np.array([t.get_position()[1] for t in ax_tree.texts if t.get_text()])
    ys.sort()
    if ys.size > 0:
        ymin, ymax = float(ys[0]), float(ys[-1])
        ax_tree.set_ylim(ymax + 0.5, ymin - 0.5)
    # 
    leaf_labels_lower = [label.lower() for label in leaf_labels]
    df.index = [idx.lower() for idx in df.index]
    df = df.reindex(leaf_labels_lower)
    # plot heatmap
    sns.heatmap(df, annot=True, cmap="vlag",  fmt=".2f", linewidths=0.1, ax=ax_heatmap)
    ax_heatmap.set_title(title, fontsize=12)
    ax_heatmap.set_ylabel("")
    ax_heatmap.set_yticklabels([])
    ax_heatmap.set_yticks([])
    # adjust heatmap position
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    x1_fig = [t.get_window_extent(renderer=renderer).transformed(fig.transFigure.inverted()).x1 for t in ax_tree.texts]
    tree_right_edge = max(x1_fig) if x1_fig else ax_tree.get_position().x1
    idea_hm_x0 = tree_right_edge + 0.005
    hm_pos = ax_heatmap.get_position()
    ax_heatmap.set_position([idea_hm_x0, hm_pos.y0, hm_pos.width, hm_pos.height])

    plt.savefig(out_path, dpi=300)
    print(f"[INFO] Plot saved to {out_path}")
    plt.close(fig)


def fig(tree_path, cofactor_df, inhibition_diff_df, is_cofactor,out_dir):
    out_path = os.path.join(out_dir, "InhibitionDiff.png")
    plot_single_figure(tree_path=tree_path,
                       df=inhibition_diff_df,
                       title="Inhibition Difference (Substrate - Product) Z-Score",
                       out_path=out_path
    )
    if is_cofactor:
        out_path = os.path.join(out_dir, "CofactorBindingEnergy.png")
        plot_single_figure(tree_path=tree_path,
                           df=cofactor_df,
                           title="Cofactor Binding Energy Z-Score",
                           out_path=out_path
        )