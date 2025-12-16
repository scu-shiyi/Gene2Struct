from Bio import Phylo
import subprocess
import os
from pathlib import Path
def build_tree(fasta_file, temp_tree_dir):

    gene_name = Path(fasta_file).stem
    gene_outdir = os.path.join(temp_tree_dir, gene_name)
    os.makedirs(gene_outdir, exist_ok=True)

    tree_prefix = os.path.join(gene_outdir, gene_name)
    final_tree_path = f"{tree_prefix}.treefile"

    if os.path.exists(final_tree_path):
        return final_tree_path

    iqtree_cmd = ["iqtree2", "-s", fasta_file, "-pre", tree_prefix]
    subprocess.run(iqtree_cmd, check=True)

    return final_tree_path


def load_tree(tree_path: str):

    tree = Phylo.read(tree_path, "newick")
    tree.ladderize()
    depths = tree.depths()
    max_depth = max(depths.values())
    return tree, depths, max_depth

def compute_leaf_positions(tree) -> dict:

    terminals = tree.get_terminals()
    leaf_positions = {term.name: idx for idx, term in enumerate(terminals)}
    return leaf_positions
