# import pandas as pd
# import numpy as np
# from scripts.utils.TreeFunction import load_tree, compute_leaf_positions
# import matplotlib.pyplot as plt
# import seaborn as sns
# from matplotlib.gridspec import GridSpec
# from Bio import Phylo
# import pandas as pd
# import numpy as np
# from typing import Optional

# def compute_universal_activity_matrix(
#     csv_path: str, 
#     out_csv_path: Optional[str] = None,
#     activity_method: str = "ratio",
#     standardize: bool = "False"
# ) -> pd.DataFrame:


#     R_kcal = 1.987e-3      # kcal·mol⁻¹·K⁻¹
#     T = 298                # K
#     RT = R_kcal * T        # ≈ 0.592 kcal/mol
    
#     df_raw = pd.read_csv(csv_path, header=None)
#     gene_row = df_raw.iloc[0, 1:]
#     role_row = df_raw.iloc[1, 1:].str.lower()
#     ligand_row = df_raw.iloc[2, 1:]
    
#     col_names = [f"{g}_{r}_{l}" for g, r, l in zip(gene_row, role_row, ligand_row)]
    
#     df_energy = df_raw.iloc[3:].copy()
#     df_energy.columns = ["Sample"] + col_names
#     df_energy.set_index("Sample", inplace=True)
    
#     tidy = (
#         df_energy.stack()
#                  .reset_index()
#                  .rename(columns={0: "ΔG", "level_1": "GeneRole"})
#     )
#     tidy[["Gene", "Role", "Ligand"]] = tidy["GeneRole"].str.split("_", n=2, expand=True)
#     tidy["ΔG"] = pd.to_numeric(tidy["ΔG"], errors="coerce") 
#     tidy = tidy[(tidy["ΔG"].notna()) & (tidy["ΔG"] < 0)]
    
#     results = []
#     for (sample, scripts), sub_df in tidy.groupby(["Sample", "Gene"]):
#         prod_energy = sub_df[sub_df["Role"] == "product"]["ΔG"].values
#         subs_energy = sub_df[sub_df["Role"] == "substrate"]["ΔG"].values

#         if (
#             len(prod_energy) == 0 or len(subs_energy) == 0 or
#             np.any(prod_energy >= 0) or np.any(subs_energy >= 0)
#         ):
#             continue

#         mean_prod = np.mean(prod_energy)
#         mean_subs = np.mean(subs_energy)

#         if activity_method == "ratio":
#             activity = abs(mean_prod / mean_subs)
#         elif activity_method == "delta":
#             activity = mean_subs - mean_prod
#         elif activity_method == "exp_delta":
#             delta_g = mean_subs - mean_prod
#             activity = np.exp(delta_g / RT)
#         else:
#             raise ValueError(f"Unknown activity_method: {activity_method}")

#         results.append({
#             "Sample": sample,
#             "Gene": scripts,
#             "Activity": activity,
#             "Mean_Substrate": mean_subs,
#             "Mean_Product": mean_prod
#         })
    
#     df_activity = pd.DataFrame(results)
    
#     if standardize and not df_activity.empty:
#         df_activity["Activity_Z"] = df_activity.groupby("Gene")["Activity"].transform(
#             lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
#         )
    
#     value_col = "Activity_Z" if standardize else "Activity"
#     df_matrix = df_activity.pivot(index="Sample", columns="Gene", values=value_col)
    

#     df_matrix = df_matrix.reset_index()
    

#     if out_csv_path:
#         df_matrix.to_csv(out_csv_path, index=False, float_format="%.6g")
    
#     return df_matrix

# def plot_clade(clade, ax, depths, max_depth, leaf_positions, tick=0.25):
#     x0 = depths[clade]
#     if clade.is_terminal():
#         y = leaf_positions[clade.name]
#         ax.hlines(y, x0, max_depth, colors='gray', linestyles='--', linewidth=0.8)
#         ax.vlines(max_depth, y - tick, y + tick, colors='black', linewidth=0.8)
#         return y

#     ys = []
#     for child in clade.clades:
#         y_child = plot_clade(child, ax, depths, max_depth, leaf_positions, tick)
#         ys.append(y_child)
#         ax.hlines(y_child, x0, depths[child], colors='black', linestyles='-', linewidth=1)

#     y_min, y_max = min(ys), max(ys)
#     ax.vlines(x0, y_min, y_max, colors='black', linestyles='-', linewidth=1)
#     return 0.5 * (y_min + y_max)

# def draw_tree_and_heatmap(
#     tree, depths, max_depth,
#     leaf_positions, heatmap_df,
#     picture_path: str,
#     name_limit: int = 20,
#     cmap: str = "RdBu_r",
#     norm_label: str = "kd-based (norm)",
# ):
#     fig = plt.figure(figsize=(25, 0.4 * len(leaf_positions) + 2), dpi=150)
#     gs = GridSpec(1, 2, width_ratios=[1, 6], wspace=0.15)

#     ax_tree = fig.add_subplot(gs[0])
#     ax_heat = fig.add_subplot(gs[1])

#     plot_clade(tree.root, ax_tree, depths, max_depth, leaf_positions)
#     ax_tree.set_ylim(-0.5, len(leaf_positions) - 0.5)
#     ax_tree.set_xlim(0, max_depth * 1.02)
#     ax_tree.invert_yaxis()
#     ax_tree.axis("off")

#     for name, y in leaf_positions.items():
#         ax_tree.text(max_depth * 1.005, y, name[:name_limit], va='center', ha='left', fontsize=8)

#     sns.heatmap(
#         heatmap_df.values,
#         ax=ax_heat,
#         cmap=cmap,
#         cbar_kws={"label": norm_label},
#         linewidths=0.5,
#         linecolor="white",
#     )
#     ax_heat.set_yticklabels([])
#     ax_heat.tick_params(axis="y", which="both", labelleft=False, left=True)
#     ax_heat.yaxis.set_ticks_position("left")
#     ax_heat.set_yticks(np.arange(len(heatmap_df)) + 0.5)

#     ax_heat.set_xticks(np.arange(len(heatmap_df.columns)) + 0.5)
#     ax_heat.set_xticklabels(heatmap_df.columns, rotation=45, ha='right', fontsize=8)
#     ax_heat.set_xlabel("Gene", fontsize=10)

#     plt.tight_layout()
#     plt.savefig(picture_path)
#     plt.close(fig)



# def plot(csv_path, out_csv_path, tree_path, pic_path):
#     tree, depths, max_depths = load_tree(tree_path)
#     leaf_position = compute_leaf_positions(tree)
#     leaf_position_lower = {k.lower(): v for k, v in leaf_position.items()}
#     kd_martrix = compute_universal_activity_matrix(csv_path, out_csv_path)
#     df_heatmap = kd_martrix.set_index("Sample")
#     df_heatmap = df_heatmap.loc[leaf_position_lower.keys()]

#     draw_tree_and_heatmap(tree, depths, max_depths, leaf_position, df_heatmap, pic_path)


import re, os
from pathlib import Path
from matplotlib import font_manager as fm
import pandas as pd
# from pyparsing import pythonStyleComment
# from scripts.DockingModule.AutoDocking import parse_mapping_csv
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
    results_df.to_csv('/home/shiyi/Gene2Struct-main/准备删除/binding_energy_results.csv', index=False)
    standardize_df = min_max_normalize(results_df)
    cofactor_df = pd.DataFrame()
    if is_cofactor:
        cofactor_df = standardize_df.pivot(index='Species', columns='Gene', values='Cofactor Binding Energy_Normalized')
    inhibition_diff_df = standardize_df.pivot(index='Species', columns='Gene', values='Inhibition Difference (Substrate - Product)_Normalized')

    return cofactor_df, inhibition_diff_df, is_cofactor


def plot_single_figure(tree_path, df, title, out_path):
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