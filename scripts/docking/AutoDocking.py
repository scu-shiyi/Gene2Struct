from pathlib import Path
import os
import shutil
import pandas as pd
from scripts.docking.DockingExecutor import collect_tasks,run_task
import re
import pandas as pd
from scripts.utils.PreparePDBQT import ligand2pdbqt, receptor2pdbqt, downloading,convert_cif_to_pdb,clean_filename,download_reference_pdb
from functools import partial
from multiprocessing import Pool, cpu_count
from scripts.docking.plot import build_table, calculate_indices, fig

def split_cell(raw):
    if raw is None or pd.isna(raw):
        return []
    text = str(raw).strip()
    if not text:
        return []
    parts = re.split(r"[;,|/]+", text)
    return [p.strip().replace(" ", "-") for p in parts if p.strip()]

def parse_mapping_csv(mapping_csv):
    df = pd.read_csv(mapping_csv)
    df.columns = [col.strip().capitalize() for col in df.columns]
    required_cols =['Gene', 'Receptordir', 'Substrate', 'Product']
    if not {'Gene','Receptordir','Substrate', 'Product'}.issubset(df.columns):
        raise ValueError("Missing required columns: Gene, Receptordir,Substrate, Product")

    gene_ligand_map = {}
    ligand_set = set()
    reference_set = set()

    for idx, row in df.iterrows():
        is_missing = False
        for col in required_cols:
            if pd.isna(row.get(col)) or str(row.get(col)).strip()=="":
                is_missing = True
                break
        if is_missing:
            print(f"Skipping row {idx + 2} due to missing required data in {required_cols} columns.")
            continue

        gene = str(row['Gene']).strip()
        receptor_dir = str(row['Receptordir']).strip()
        substrates = split_cell(row['Substrate'])
        products = split_cell(row['Product'])
        cofactors = split_cell(row.get('Cofactor', ''))
        reference = split_cell(row.get('Reference', ''))
        if reference:
            ref = reference[0]
            reference_set.add(ref)
        else:
            ref= None
        gene_ligand_map[gene] = {
            "receptordir": receptor_dir,
            "substrate": substrates,
            "product": products,
            "cofactor": cofactors,
            "reference": ref
        }
        ligands = list(set(substrates + products+cofactors))
        ligand_set.update(ligands)

    return gene_ligand_map, ligand_set, reference_set

def cif2pdb(gene, input_receptor_dir, receptor_pdb_dir):
    out_dir = os.path.join(receptor_pdb_dir, gene)
    os.makedirs(out_dir, exist_ok=True)

    for file in os.listdir(input_receptor_dir):
        src_path = os.path.join(input_receptor_dir, file)
        if os.path.isdir(src_path):
            continue
        base,_ = os.path.splitext(file)
        basename_cleaned = clean_filename(base, gene)
        dst_path_base = os.path.join(out_dir, basename_cleaned)
        dst_path = f"{dst_path_base}.pdb"
        if os.path.exists(dst_path):
            continue
        if file.lower().endswith('.pdb'):
            shutil.copy2(src_path, dst_path)
        elif file.lower().endswith(".cif"):
            convert_cif_to_pdb(src_path, dst_path)

def pdb2pdbqt(receptor_dir):
    gene_receptor_map = {}
    for root, _, files in os.walk(receptor_dir):
        if os.path.normpath(root) == os.path.normpath(receptor_dir):
            continue
        gene_name = os.path.basename(root)
        current_receptor_info = []
        for file in files:
            if file.lower().endswith('.pdb'):
                base_name = Path(file).stem
                species_name = base_name
                src_pdb_path = os.path.join(root, file)
                dst_pdbqt_path = os.path.join(root, f"{base_name}.pdbqt")
                pdbqt_file_path = None

                if os.path.exists(dst_pdbqt_path):
                    pdbqt_file_path = dst_pdbqt_path
                else:
                    try:
                        pdbqt_file_path = receptor2pdbqt(src_pdb_path, dst_pdbqt_path)
                    except Exception as e:
                        print(f"[ERROR] PDBQT conversion failed for {file} in {gene_name}: {e}")
                if pdbqt_file_path:
                    current_receptor_info.append({'species':species_name,
                                                  'pdb_file':src_pdb_path,
                                                  'pdbqt_file':pdbqt_file_path})
        
        if current_receptor_info:
            gene_receptor_map[gene_name] = current_receptor_info
    return gene_receptor_map


def AutoDocking(mapping_csv, tree_path, output_dir):
    receptor_dir = os.path.join(output_dir, "receptors")
    ligand_dir = os.path.join(output_dir, "ligands")
    docking_dir = os.path.join(output_dir, "docking_results")
    pocket_dir = os.path.join(output_dir, "pockets")
    os.makedirs(docking_dir, exist_ok=True)
    os.makedirs(pocket_dir, exist_ok=True)
    os.makedirs(receptor_dir, exist_ok=True)
    os.makedirs(ligand_dir, exist_ok=True)

    FPOCKET_BIN = shutil.which("fpocket")
    if FPOCKET_BIN is None:
        raise FileNotFoundError("fpocket not found in PATH. Please check your environment setup.")
    VINA_BIN = shutil.which('vina')
    if VINA_BIN is None:
        raise FileNotFoundError("vina not found in PATH. Please check your Conda environment.")
    
    gene_ligand_map, ligand_set, reference_set = parse_mapping_csv(mapping_csv)
    ref_pdb_dir = download_reference_pdb(list(reference_set), receptor_dir)
    sdf_file_list = downloading(list(ligand_set), ligand_dir)
    ligand_path_map = ligand2pdbqt(sdf_file_list, ligand_dir)
    for gene in gene_ligand_map:
        cif2pdb(gene, gene_ligand_map[gene]['receptordir'], receptor_dir)
    gene_receptor_map = pdb2pdbqt(receptor_dir)
    tasks = collect_tasks(gene_ligand_map, gene_receptor_map, ligand_path_map,ref_pdb_dir)
    n_tasks = len(tasks)
    total_cpus = cpu_count()
    if total_cpus >= 2:
        cpu_per_task = 8
    else:
        cpu_per_task = 1
    n_parallel = max(1, total_cpus // cpu_per_task)
    n_parallel = min(n_parallel, n_tasks)

    run_task_partial = partial(
        run_task,
        docking_dir=docking_dir,
        pocket_dir=pocket_dir,
        VINA_BIN=VINA_BIN,
        FPOCKET_BIN=FPOCKET_BIN,
        cpu_per_task=cpu_per_task)
        
    print(f"[INFO] Collected {len(tasks)} tasks.")
    print(f"[INFO] Starting parallel docking using {total_cpus} CPU cores...")
    print(f"[INFO] Running {n_parallel} tasks in parallel, total {n_tasks / n_parallel} tasks.")
    # from concurrent.futures import ProcessPoolExecutor
    with Pool(processes=n_parallel) as pool:
        pool.map(run_task_partial, tasks)
    print("[INFO] All docking tasks finished.") 
    print("[INFO] All docking tasks finished.")


    docking_summary_csv = os.path.join(output_dir, 'docking_summary.csv')
    orign_df = build_table(docking_dir, gene_ligand_map, docking_summary_csv)
    cofactor_df, inhibition_diff_df, is_cofactor = calculate_indices(orign_df)
    fig(tree_path=tree_path,
        cofactor_df=cofactor_df, 
        inhibition_diff_df=inhibition_diff_df,
        is_cofactor=is_cofactor,
        out_dir = output_dir)


# def AutoDocking(data_mapping_csv, output_dir, tree_path):
    # pdb_dir = os.path.join(output_dir, 'pdb')
    # receptor_pdb = os.path.join(pdb_dir, 'receptor_pdb')
    # ligand_pdb = os.path.join(pdb_dir, 'ligand_pdb')
    # os.makedirs(receptor_pdb, exist_ok=True)
    # os.makedirs(ligand_pdb, exist_ok=True)

    # temp_dir = os.path.join(output_dir, "temp")
    # temp_ligand = os.path.join(temp_dir, 'ligand')
    # temp_center = os.path.join(temp_dir, 'center')
    # os.makedirs(temp_ligand, exist_ok=True)
    # os.makedirs(temp_center, exist_ok=True)
    
    # pdbqt_dir = os.path.join(output_dir, 'pdbqt')
    # receptor_pdbqt = os.path.join(pdbqt_dir, "receptor_pdbqt")
    # ligand_pdbqt = os.path.join(pdbqt_dir, "ligand_pdbqt")
    # os.makedirs(receptor_pdbqt, exist_ok=True)
    # os.makedirs(ligand_pdbqt, exist_ok=True)

    # docking_dir = os.path.join(output_dir, 'docking')

    # gene_list = [] 
    # for scripts in sorted(os.listdir(receptor_pdb)):
    #     gene_path = os.path.join(receptor_pdb, scripts)
    #     if os.path.isdir(gene_path):
    #         valid = any(f.endswith(('.pdb', '.cif')) for f in os.listdir(gene_path))
    #         if valid:
    #             gene_list.append(scripts)

    # gene_ligand_map, ligand_set, genes = parse_mapping_csv(data_mapping_csv)

    # if set(gene_list) == set(genes):
    #     process_ligand(ligand_set, temp_ligand, ligand_pdb, ligand_pdbqt)
    # else:
    #     print(set(gene_list))
    #     print(set(genes))
    #     raise ValueError("Gene names in protein folder and mapping CSV must match exactly.")


    # for root,_,files in os.walk(receptor_pdb):
    #     for file in files:
    #         if file.endswith(".pdb"):
    #             src_file = os.path.join(root, file)
    #             relative_path = os.path.relpath(root, receptor_pdb)
    #             out_dir = os.path.join(receptor_pdbqt, relative_path)
    #             os.makedirs(out_dir, exist_ok=True)
    #             process_receptors(src_file, out_dir)


    # executor = DockingExecutor(
    #     receptor_pdb_dir = receptor_pdb,
    #     receptor_pdbqt_dir = receptor_pdbqt,
    #     ligand_pdbqt_dir = ligand_pdbqt,
    #     docking_dir= docking_dir,
    #     temp_center=temp_center,
    #     gene_ligand_map = gene_ligand_map,
    # )
    # executor.run_all()

    # row_binding_energy = os.path.join(output_dir, 'raw_binding_energy.csv')
    # catalytic_activity_matrix = os.path.join(output_dir,'catalytic_activity_matrix.csv')
    # phylo_activity_heatmap = os.path.join(output_dir, 'phylo_activity_heatmap.png')
    # build_table(docking_dir, mapping_csv, row_binding_energy)

    # plot(row_binding_energy, catalytic_activity_matrix, tree_path, phylo_activity_heatmap)

# if __name__ == '__main__':
#     data_mapping_csv = '/home/shiyi/Gene2Struct-main/scripts/DockingModule/data.csv'
#     out = '/home/shiyi/Gene2Struct-main/scripts/DockingModule/final_out'
#     tree_path = "/home/shiyi/Gene2Struct-main/scripts/evolution/OUTPUT/4HPAAS/paml_input/4HPAAS.paml.tree"
#     AutoDocking(data_mapping_csv, out,tree_path)


    # gene_ligand_map, ligand_set,x = parse_mapping_csv(data_mapping_csv)
    # # print(type(gene_ligand_map))
    # for scripts, info in gene_ligand_map.items():
    #     print(scripts, info)
    # print(x)
    # sdf_file_list = downloading(list(ligand_set), '/home/shiyi/Gene2Struct-main/scripts/DockingModule/output/ligands')
    # ligand_path_map = ligand2pdbqt(sdf_file_list,'/home/shiyi/Gene2Struct-main/scripts/DockingModule/output/ligands')
    # for scripts in gene_ligand_map:
    #     gene_ligand_map[scripts]['receptordir']
    #     cif2pdb(scripts, gene_ligand_map[scripts]['receptordir'], '/home/shiyi/Gene2Struct-main/scripts/DockingModule/output/receptors')
    # gene_receptor_map = pdb2pdbqt('/home/shiyi/Gene2Struct-main/scripts/DockingModule/output/receptors')
    # tasks = collect_tasks(gene_ligand_map, gene_receptor_map,ligand_path_map)
