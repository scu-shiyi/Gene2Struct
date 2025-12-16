from pathlib import Path
from Bio import SeqIO
from scripts.utils.Phylip_Prepare import  prepare_paml_input
from scripts.utils.EvoScoring import position_entropy
import pandas as pd
import ast
from scripts.utils.site_model import run_pair_model
from scripts.evolution.plot import Visualization

def Run_Site(fasta_path, output_dir, tree_path,is_site, name_limit, thr, outgroups, heatmap_path=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_name = Path(fasta_path).stem

    file_input = output_dir / gene_name / "file_input"
    file_input.mkdir(parents=True, exist_ok=True)
    phylip_file, paml_treefile, species = prepare_paml_input(fasta_path, file_input, tree_path=tree_path)

    if heatmap_path is None:
        evo_output = output_dir / gene_name / "evo_output"
        evo_output.mkdir(parents=True, exist_ok=True)
        heatmap_path = position_entropy(phylip_file, str(evo_output))

    heatmap_df = pd.read_csv(heatmap_path, index_col=0, header=0)
    heatmap_df = heatmap_df.map(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    if is_site:
        paml_output = output_dir / gene_name / "paml_output"
        paml_output.mkdir(parents=True, exist_ok=True)
        m0m3_base_dir = paml_output / "M0M3"
        M0M3_result = run_pair_model(phylip_file, paml_treefile, m0m3_base_dir, "M0M3",gene_name=gene_name)
        if M0M3_result["p"] < 0.05:
            print("LRT (M0 vs M3) is significant at the 0.05/0.01 level: evidence of site-specific ω heterogeneity.")
            m7m8_base_dir = paml_output / "M7M8"
            m7m8_result = run_pair_model(phylip_file, paml_treefile, m7m8_base_dir, "M7M8", gene_name=gene_name)
            mlc_path = m7m8_base_dir / 'result' / 'M8_mlc'
            rst_path = m7m8_base_dir / '__work_M8' / 'rst'
            if m7m8_result["p"] < 0.05:
                print("LRT (M7 vs M8) is significant at the 0.05/0.01 level: positive selection sites detected.")
            else:
                print("LRT (M7 vs M8) is not significant: no positive selection sites detected.")
        else:
            mlc_path = None
            rst_path = None
    else:
        rst_path = None
    outpng_path = output_dir / gene_name / f"{gene_name}.png"
    Visualization(paml_treefile, heatmap_df, outpng_path, rst_path=rst_path, name_limit=name_limit, thr = thr, site_model=is_site, outgroups=outgroups)

if __name__ == "__main__":
    Run_Site('//evolution/FPS.fasta',
            '//evolution/OUTPUT',
            '//evolution/FPS.paml.tree',
            is_site=False, name_limit=10, thr=1.4, outgroups=['XM_001693116'])