from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from Bio.Phylo.BaseTree import Clade
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import FixedLocator, NullFormatter
import matplotlib.font_manager as fm
from pathlib import Path
import re
import pandas as pd

def _collapse_mean_by_k_cols(arr, k=3, strict=True):

    nrows, ncols = arr.shape
    if strict and (ncols % k != 0):
        raise ValueError(f"Number of columns {ncols} is not divisible by {k}, cannot collapse by codons.")
    pad = (-ncols) % k
    if pad:
        arr = np.concatenate([arr, np.full((nrows, pad), np.nan)], axis=1)
    g = arr.reshape(nrows, -1, k)          # (rows, groups, k)
    return np.nanmean(g, axis=2)           # (rows, groups)

def root_and_force_bifurcation(tree, outgroups):
    """
    Roots the tree using the outgroup's MRCA and surgically forces a binary root split
    by merging all scattered outgroup and ingroup branches.
    """
    if outgroups is None:
        outgroups = []
    # 1. Prepare nodes and names
    valid_og_names = [og for og in outgroups if og in [c.name for c in tree.get_terminals()]]
    if not valid_og_names:
        return tree

    valid_og_nodes = [c for c in tree.get_terminals() if c.name in valid_og_names]

    # 2. Initial rooting using the MRCA
    try:
        monophyletic_og_root = tree.common_ancestor(valid_og_nodes)
        tree.root_with_outgroup(monophyletic_og_root)
    except Exception:
        return tree

    # 3. Classify root's current children
    outgroup_branches = [
        branch for branch in tree.clade.clades
        if any(t.name in valid_og_names for t in branch.get_terminals())
    ]
    innergroup_branches = [
        branch for branch in tree.clade.clades
        if not any(t.name in valid_og_names for t in branch.get_terminals())
    ]

    if not outgroup_branches:
        return tree

    # 4. Surgical Reassembly (Creating new Clades and rebuilding the root)
    outgroup_clade = Clade()
    outgroup_clade.clades = outgroup_branches

    innergroup_clade = Clade()
    innergroup_clade.clades = innergroup_branches

    tree.clade.clades = [outgroup_clade, innergroup_clade]

    # 5. Set branch lengths to avoid None-value errors
    if outgroup_clade.branch_length is None:
        outgroup_clade.branch_length = 0.0
    innergroup_clade.branch_length = outgroup_clade.branch_length

    return tree


def set_topological_ultrametric(tree, unit_length=1.0):
    """
    根据最深子簇的拓扑层级数来设置分支长度，实现您要求的对齐和延伸。
    """

    # 1. 后序遍历：计算每个节点到其最深叶子的拓扑层级数 (Clade Depth)
    # 这个层级数将作为新的距离度量
    for clade in tree.find_clades(order='postorder'):
        if clade.is_terminal():
            # 叶子节点，拓扑深度为 0
            clade.topo_depth = 0
        else:
            # 内部节点到最深叶子的拓扑距离 = max(子节点拓扑深度) + 1
            max_child_depth = max(child.topo_depth for child in clade.clades)
            clade.topo_depth = max_child_depth + 1

    # 2. 找出最大的拓扑深度（即总距离）
    max_topo_depth = tree.root.topo_depth if hasattr(tree.root, 'topo_depth') else 1

    # 3. 前序遍历：设置新的分支长度
    # 新分支长度 = (父节点到最深叶子的层级数) - (子节点到最深叶子的层级数)

    for clade in tree.find_clades(order='preorder'):

        # 找出该节点到其最深叶子的总距离
        parent_depth_to_deepest_leaf = clade.topo_depth

        for child in clade.clades:
            # 子节点到其最深叶子的总距离
            child_depth_to_deepest_leaf = child.topo_depth

            # 新分支长度 = (父节点拓扑深度 - 子节点拓扑深度) * 单位长度
            # 这个差值就是该分支需要“延伸”的层级数
            new_branch_length = (parent_depth_to_deepest_leaf - child_depth_to_deepest_leaf) * unit_length

            # 设置新的分支长度
            child.branch_length = max(0.000001, new_branch_length)

            # 根节点分支长度设为 0
    tree.root.branch_length = 0
    return tree


def sort_subtrees_by_size(clade):
    """
    Recursively sort the subtrees of a clade based on the size (number of terminals)
    of the sub-clades. Used to order branches for visualization.
    """
    if clade.clades:
        # 递归调用自身，先对子树进行排序
        for subclade in clade.clades:
            sort_subtrees_by_size(subclade)

        # 根据子分支的叶子数进行排序 (从小到大)
        clade.clades.sort(key=lambda c: -len(c.get_terminals()))


def get_leaf_order_topdown(tree):
    """
    Returns a list of terminal names in the order they appear when traversing
    the tree from top to bottom (based on the current sorted layout).
    """
    # tree.get_terminals() 默认按 postorder 排序，我们需要 pre-order 或 top-down 顺序
    leaf_names = []

    def traverse_clade(clade):
        if clade.is_terminal():
            leaf_names.append(clade.name)
        else:
            # 遍历子分支，保持 sort_subtrees_by_size 设定的顺序（从小到大）
            for subclade in clade.clades:
                traverse_clade(subclade)

    traverse_clade(tree.clade)
    return leaf_names

def parse_rst_beb_p_positive(rst_path):
    with open(rst_path, "r", encoding="UTF-8", errors="ignore") as f:
        text = f.read()

        # 定位 BEB 段落的标题行
    hdr_pat = re.compile(
        r"Bayes\s+Empirical\s+Bayes\s+\(BEB\)\s+probabilities\s+for\s+(\d+)\s+classes.*",
        re.IGNORECASE
    )
    hdr_m = hdr_pat.search(text)
    if not hdr_m:
        raise ValueError("未在 rst 中找到 BEB 段落标题。")

    n_classes_decl = int(hdr_m.group(1))
    # 目标正选择类：默认取“最后一类”
    target_idx_1based = n_classes_decl
    if not (1 <= target_idx_1based <= n_classes_decl):
        raise ValueError(f"positive_class_index 超出范围 1..{n_classes_decl}")

    # 从标题行后开始逐行解析
    lines = text[hdr_m.end():].splitlines()

    rows = []
    started = False
    for line in lines:
        s = line.strip()
        if not s:
            # 如果已经开始采集并遇到空行，通常说明 BEB 表结束
            if started:
                break
            else:
                continue

        # 匹配一行：开头是位点编号 + 氨基酸
        # 例： "  10 P   0.00005 0.00027 ... 0.90234 (11)  1.486 +-  0.345"
        m_head = re.match(r"^\s*(\d+)\s+([A-Za-z\-\?])\s+(.*)$", line)
        if not m_head:
            # 若已经开始采集，遇到不再匹配的行，就视为表结束
            if started:
                break
            else:
                continue

        started = True
        pos = int(m_head.group(1))
        aa = m_head.group(2)
        rest = m_head.group(3)

        # 抓取该行中的所有浮点数（包含小数 / 科学计数）
        floats = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?", rest)
        # 必须至少包含：n_classes 个概率 + postmean + sd  → 共 n_classes + 2
        if len(floats) < n_classes_decl + 2:
            # 这行可能不是一个完整的 BEB 行；跳过
            continue

        probs = list(map(float, floats[:n_classes_decl]))
        postmean_w = float(floats[n_classes_decl])
        postmean_sd = float(floats[n_classes_decl + 1])

        # 提取括号里的最可能类别（可选）
        m_cls = re.search(r"\(\s*(\d+)\s*\)", rest)
        class_ml = int(m_cls.group(1)) if m_cls else None

        prob_pos = probs[target_idx_1based - 1]  # 1-based → 0-based

        rows.append({
            "pos": pos,
            "aa": aa,
            "prob": prob_pos,  # 正选择类（通常第11类）后验概率
            "class_ml": class_ml,
            "postmean_w": postmean_w,
            "postmean_sd": postmean_sd,
            "probs": probs,
            "n_classes": n_classes_decl
        })

    if not rows:
        raise ValueError("未能解析到任何 BEB 行。请检查 rst 是否为 codeml 的 BEB 输出格式。")

    df = pd.DataFrame(rows).sort_values("pos").reset_index(drop=True)
    return df


def add_beb_track_all(ax_heatmap, beb_df, ncols):
    gap = 0.01
    fig = ax_heatmap.figure
    box = ax_heatmap.get_position()
    h = box.height * 0.1
    ax_b = fig.add_axes([box.x0,
                         box.y0 - h - gap,
                         box.width,
                         h],
                        sharex=ax_heatmap)
    pr = np.zeros(int(ncols), dtype=float)
    if beb_df is not None and not beb_df.empty:
        for pos, prob in zip(beb_df['pos'].astype(int),
                             beb_df['prob'].astype(float)):
            i = pos - 1  # 1-based → 0-based
            if 0 <= i < ncols:
                pr[i] = prob
    x = np.arange(0, ncols)
    ax_b.bar(x, pr, width=0.9, color="black", edgecolor='none',
             align='center', zorder=1)
    ax_b.set_ylim(0, 1.10)
    ax_b.axhline(0.95, ls='--', lw=0.8, color='black')
    ax_b.set_ylabel("Pr(ω>1)", fontsize=8)
    ax_b.grid(axis='y', ls=':', lw=0.5, color='#e2e8f0')
    ax_b.tick_params(axis='x', bottom=True, labelbottom=True, labelrotation=60, pad=2)
    plt.setp(ax_b.get_xticklabels(), ha='right', rotation_mode='anchor')
    for s in ("right", "top"):
        ax_b.spines[s].set_visible(False)


def Visualization(paml_treefile, heatmap_df, outpng_path, rst_path=None, name_limit=10, thr = 1.5, site_model=False, outgroups=None):
    # 2. 设置全局字体为 Times New Roman
    plt.rcParams["font.family"] = "Times New Roman"
    # 3. 确保数学符号字体也匹配 Times
    plt.rcParams["mathtext.fontset"] = "stix"
    # 创建画布
    fig = plt.figure(figsize=(20, 10), dpi=150)
    gs = GridSpec(1, 3, width_ratios=[1, 5, 0.5], wspace=0)
    ax_tree = fig.add_subplot(gs[0])
    ax_heatmap = fig.add_subplot(gs[1])
    ax_legend = fig.add_subplot(gs[2])
    # 绘制进化树
    tree = Phylo.read(paml_treefile, "newick")
    tree = root_and_force_bifurcation(tree, outgroups)
    max_depth = max(c.topo_depth if hasattr(c, 'topo_depth') else 0 for c in tree.find_clades(order='postorder'))
    unit_length = 1.0 / (max_depth or 1.0)
    aligned_tree = set_topological_ultrametric(tree, unit_length)
    sort_subtrees_by_size(aligned_tree.root)
    leaf_labels = get_leaf_order_topdown(aligned_tree)
    Phylo.draw(aligned_tree, axes=ax_tree, do_show=False,
               branch_labels=None,
               show_confidence=False,
               label_func=lambda x: x.name[:name_limit] if x.name else "")

    for text_obj in ax_tree.texts:
        text_obj.set_fontname("Times New Roman")
        text_obj.set_fontsize(8)
        text_obj.set_fontweight('normal')
    ax_tree.axis('off')
    fig.canvas.draw()  # 让文本排版完成，拿到真实坐标

    # 叶标签的 y 坐标（数据坐标）
    ys = np.array([t.get_position()[1] for t in ax_tree.texts if t.get_text()])
    ys.sort()
    ymin, ymax = float(ys[0]), float(ys[-1])
    # 只调树：把树的可视范围设置成“上下各半格”的长度，和热图等长
    ax_tree.set_ylim(ymax + 0.5 , ymin -0.5)

    # retu
    raw = heatmap_df.reindex(leaf_labels)
    raw = raw.values.astype(float)
    agg = _collapse_mean_by_k_cols(raw, k=3, strict=True)
    cat = np.full_like(agg, np.nan, dtype=float)
    mask = ~np.isnan(agg)
    cat[mask] = (agg[mask] >= thr).astype(int)
    cmap = ListedColormap(['#d24735', "#a9be7b"])
    cmap.set_bad('#e6e6e6')
    norm = BoundaryNorm([0, 0.5, 1.5], cmap.N)
    ax_heatmap.set_facecolor('white')
    ax_heatmap.imshow(cat, aspect='auto', cmap=cmap, norm=norm,
                           interpolation='nearest', origin='upper')

    legend_elements = [
        Patch(facecolor="#d24735", edgecolor='k', label=f'Highly variable (<{thr})'),
        Patch(facecolor="#a9be7b", edgecolor='k', label=f'Highly conserved (≥{thr})')
    ]
    leg = ax_legend.legend(handles=legend_elements,
        title="  Amino acid site  \nconservation (bits)",
        loc="upper center", frameon=True, fancybox=True,
        framealpha=0.9, edgecolor="0.3",
        fontsize=10, title_fontsize=11,
        handlelength=1.4, handleheight=1.0, borderaxespad=0.0
    )
    ax_legend.axis('off')


    if site_model:
        ncols = cat.shape[1]
        ticks = np.arange(0, ncols , 10)
        ax_heatmap.set_xticks(ticks)
        ax_heatmap.set_xticklabels(ticks+1, rotation=60, fontsize=6)
        ax_heatmap.tick_params(
            axis='x',
            top=False, labeltop=False,  # 顶部显示刻度与标签
            bottom=False, labelbottom=False,  # 关闭底部刻度与标签
            pad=2  # 标签与轴的间距，按需调
        )
    else:
        ncols = cat.shape[1]
        ticks = np.arange(0, ncols , 10)
        ax_heatmap.set_xticks((ticks))
        ax_heatmap.set_xticklabels(ticks + 1, fontsize=6)
        ax_heatmap.tick_params(
            axis='x',
            top=False, labeltop=False,
            bottom=True, labelbottom=True,
        )
    nrows = cat.shape[0]
    ax_heatmap.yaxis.set_major_locator(FixedLocator(np.arange(nrows)))  # 0,1,2,...
    ax_heatmap.yaxis.set_major_formatter(NullFormatter())  # 不显示数值
    ax_heatmap.tick_params(axis='y', which='major',
                           left=True, right=False, length=2.5, width=0.8,
                           labelleft=False)

    # --- 5. AXES ALIGNMENT: Dynamic Layout (Heatmap First) ---

    # 定义布局常量
    EDGE_MARGIN = 0.01  # 画布左右边缘保留的距离
    CONTENT_PAD = 0.005  # 元素之间的间距 (树标签/热图)

    MIN_HM_WIDTH = 0.10  # 热图最小宽度（防止被挤压得太小）

    # ----------------------------------------------------
    # 步骤 1: 设置树的左边界，并获取标签的精确位置
    # ----------------------------------------------------
    tree_pos = ax_tree.get_position()
    ax_tree.set_position([
        EDGE_MARGIN,  # 左边缘固定
        tree_pos.y0+0.05,
        tree_pos.width,
        tree_pos.height
    ])

    legend_pos = ax_legend.get_position()
    ax_legend.set_position([
        1-EDGE_MARGIN-0.1,  # 左边缘固定
        legend_pos.y0+0.05,
        0.1,
        legend_pos.height
    ])

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # A. 树标签的右边界 (T_R)
    x1_fig = [t.get_window_extent(renderer=renderer).transformed(fig.transFigure.inverted()).x1 for t in ax_tree.texts]
    tree_right_labels = max(x1_fig) if x1_fig else ax_tree.get_position().x1

    # ----------------------------------------------------
    # 步骤 2: 计算热图的最终绝对位置 (左侧依赖树，右侧依赖画布边缘)
    # ----------------------------------------------------
    hm_pos = ax_heatmap.get_position()
    legend_pos = ax_legend.get_position()
    # 目标热图左边界 (HM_L): 紧贴树标签 + pad
    ideal_hm_x0 = tree_right_labels + CONTENT_PAD

    # 目标热图右边界 (HM_R): 紧贴画布右边缘 - margin
    ideal_hm_width = legend_pos.x0-ideal_hm_x0 -0.0001
    ax_heatmap.set_position([
        ideal_hm_x0,
        hm_pos.y0+0.05,
        ideal_hm_width,
        hm_pos.height
    ])

    if rst_path:
        beb_df = parse_rst_beb_p_positive(rst_path)  # 不筛阈值，拿全体 BEB 位点
        ncols = cat.shape[1]
        add_beb_track_all(ax_heatmap, beb_df, ncols=ncols)


    plt.savefig(outpng_path)
    plt.close(fig)

    return str(outpng_path)


if __name__ == "__main__":
    import ast
    import pandas as pd
    treefile1 = '/Users/sy/Gene2Struct-main/scripts/evolution/OUTPUT/4HPAAS/paml_input/4HPAAS.paml.tree'
    heatmap_df1 = pd.read_csv('//evolution/OUTPUT/4HPAAS/evo_dir/Conservation_score.csv', index_col=0, header=0)
    Visualization(treefile1, heatmap_df1, '//evolution', '4HPAAS', thr=1.0, rst_path="//evolution/OUTPUT/4HPAAS/M7M8/__work_M8/rst", outgroups=['Rhodiola_himalensis'], name_limit=20)
    # Visualization(treefile1, heatmap_df1,'/Users/sy/Gene2Struct-main/scripts/evolution','4HPAAS',thr=1.0,outgroups=[ 'Rhodiola_himalensis'],name_limit=20)

    # treefile2 = '/Users/sy/Gene2Struct-main/scripts/evolution/out/rbcl.treefile'
    # heatmap_df2 = pd.read_csv('/Users/sy/Gene2Struct-main/scripts/evolution/out/rbcl_entropy.csv', index_col=0,header=0)
    # heatmap_df2 = heatmap_df2.map(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    # Visualization(treefile2, heatmap_df2,'/Users/sy/Gene2Struct-main/scripts/evolution','rbcl1',thr=1.6,outgroups=['Dacrycarpus_imbricatus', 'Dacrycarpus_dacrydioides', 'Dacrydium_elatum', 'Ginkgo_biloba'])

    # print(heatmap_df.shape)

# if __name__ == "__main__":
#     # 创建画布
#     BASE_COLS = 420
#     BASE_ROWS = 70
#     BASE_W_IN, BASE_H_IN = 20.0, 10.0
#     DPI_FIXED = 150
#     PX_PER_COL = BASE_W_IN * DPI_FIXED / BASE_COLS
#     PX_PER_ROW = BASE_H_IN * DPI_FIXED / BASE_ROWS
#     R = 50
#     C = 400
#     # 可选画布尺寸上下限（英寸）
#     MIN_W, MAX_W = 8.0, 48.0
#     MIN_H, MAX_H = 6.0, 36.0
#     w_in_A = (PX_PER_COL * C) / DPI_FIXED
#     h_in_A = (PX_PER_ROW * R) / DPI_FIXED
#     w_in_A = float(np.clip(w_in_A, MIN_W, MAX_W))
#     h_in_A = float(np.clip(h_in_A, MIN_H, MAX_H))
#     print(w_in_A, h_in_A)