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
    os.makedirs(receptor_dir, exist_ok=True)
    os.makedirs(ligand_dir, exist_ok=True)


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


