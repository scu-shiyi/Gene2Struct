import os
import subprocess
from pathlib import Path
import shutil
import tempfile
from weakref import ref
import numpy as np
from Bio.PDB import PDBParser, Superimposer, Selection
from Bio.PDB.Polypeptide import standard_aa_names
import re
def best_affinity(log_file):
    _aff_pat = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)")
    with open(log_file) as f:
        for line in f:
            m = _aff_pat.match(line)
            if m:
                return float(m.group(1))
            
    print("Could not parse log file:", log_file)
    return None
def run_fpocket(receptor_pdb, dest,receptor_name, FPOCKET_BIN):
    # FPOCKET_BIN = shutil.which("fpocket")
    # if FPOCKET_BIN is None:
    #     raise FileNotFoundError("fpocket not found in PATH. Please check your environment setup.")
    receptor_pdb = Path(receptor_pdb)
    if Path(dest).exists():
        return dest
    os.makedirs(dest, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        tmp_receptor = tmpdir / f"{receptor_name}.pdb"
        shutil.copy2(receptor_pdb, tmp_receptor)

        command = [FPOCKET_BIN, "-f", str(tmp_receptor)]
        try:
            subprocess.run(command, check=True, cwd=tmpdir)
        except subprocess.CalledProcessError as e:
            print(f"fpocket failed for {receptor_pdb.name}: {e}")
            shutil.rmtree(dest)
            return None
        
        out_dir = tmpdir / f"{receptor_name}_out"
        if out_dir.exists():
            try:
                shutil.copytree(str(out_dir), str(dest), dirs_exist_ok=True)
            except shutil.Error as e:
                print(f"[REEOR] Copying fpocket output failed: {e}")
                shutil.rmtree(dest, ignore_errors=True)
                return None
            return dest
        else:
            print(f"fpocket completed but no output directory generated: {out_dir}")
            shutil.rmtree(dest, ignore_errors=True)
            return None


def parse_fpocket_pqr(output_file):

    pqr_file = os.path.join(output_file, 'pockets/pocket2_vert.pqr')
    if not os.path.exists(pqr_file):
        print("No pockets.pqr file found. Please check fpocket results.")
        return None, None
    coords = []
    with open(pqr_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
                coords.append([x, y, z])
    if len(coords) == 0:
        print("PQR file contains no atom data. Cannot calculate pocket center/size.")
        return None, None

    coords = np.array(coords)
    center = np.mean(coords, axis=0)  
    size = np.ptp(coords, axis=0) + 10  

    return center, size



def align_and_get_pocket_params(
    target_pdb_path,
    reference_pdb_path,
    chain_id: str = "A"):
    """
    使用 BioPython 执行结构比对，将参考蛋白移动到目标蛋白的空间位置，
    然后根据移动后的参考蛋白上的 HET 配体计算 Vina 盒子中心和尺寸。

    Args:
        target_pdb_path: 目标受体 (UGT3) 的 PDB 文件路径。
        reference_pdb_path: 参考共晶结构 (8ITA) 的 PDB 文件路径。
        chain_id: 用于比对的链 ID，默认为 'A'。
        min_ligand_atoms: 视为配体的最小原子数。
        box_margin: 盒子尺寸的额外边距 (Å)。

    Returns:
        (center, size) tuple of numpy arrays, or (None, None) if alignment fails.
    """
    min_ligand_atoms = 5
    box_margin = 6.0 
    parser = PDBParser(QUIET=True)
    try:
        # 1. 加载结构
        target_structure = parser.get_structure('target', target_pdb_path)
        reference_structure = parser.get_structure('reference', reference_pdb_path)

        target_model = target_structure[0]
        ref_model = reference_structure[0]

        if chain_id not in target_model or chain_id not in ref_model:
            print(f"[ERROR] Chain {chain_id} not found in target or reference.")
            return None, None

        target_chain = target_model[chain_id]
        ref_chain = ref_model[chain_id]

        # 2. 提取 Cα 列表（按残基 ID 匹配）
        
        # 辅助函数：提取 (res_id, CA_atom) 列表
        def get_ca_map(chain):
            ca_map = {}
            for residue in chain:
                # 仅考虑标准氨基酸（残基ID的插入码部分为 ' '）
                if residue.id[0] == ' ' and residue.get_resname() in standard_aa_names:
                    if 'CA' in residue:
                        ca_map[residue.id] = residue['CA']
            return ca_map

        target_ca_map = get_ca_map(target_chain)
        ref_ca_map = get_ca_map(ref_chain)

        # 找出两个结构中都存在的残基 ID (res_id)
        common_res_ids = sorted(list(set(target_ca_map.keys()) & set(ref_ca_map.keys())))

        if not common_res_ids:
            print("[ERROR] No common Cα atoms found for standard amino acids in the specified chain.")
            return None, None

        # 构建最终比对列表
        target_ca = [target_ca_map[res_id] for res_id in common_res_ids]
        ref_ca = [ref_ca_map[res_id] for res_id in common_res_ids]

        # 3. 对齐：把 reference -> target
        #    set_atoms(fixed, moving) ； fixed = target, moving = reference
        super_imposer = Superimposer()
        super_imposer.set_atoms(target_ca, ref_ca)
        rot, tran = super_imposer.rotran  # 旋转矩阵 & 平移向量
        
        # 必须对参考结构的所有原子应用变换，包括配体
        reference_structure.transform(rot, tran)

        print(f"[INFO] Structural alignment based on {len(common_res_ids)} Cα pairs. RMSD: {super_imposer.rms:.3f} Å")

        # 4. 在参考结构中自动识别“配体”
        #    规则：非水 (HOH/WAT)、非标准氨基酸、原子数 >= min_ligand_atoms
        ligand_coords_tgt = [] # 注意：现在配体坐标已经在目标坐标系中了

        # 一些常见核酸残基，防止被当成配体
        nucleic_resnames = {"DA", "DT", "DG", "DC", "A", "U", "G", "C"}
        found_ligands = set()

        # 遍历已变换的参考结构
        for model in reference_structure:
            for chain in model:
                for residue in chain:
                    resn = residue.get_resname().strip()
                    # PDB Residue ID结构: ('H_', residue_number, insertion_code)
                    # 'H_' 代表 HETATM
                    het_flag = residue.id[0]

                    # 过滤掉水 (HOH/WAT)
                    if resn in ("HOH", "WAT"):
                        continue
                    # 过滤掉标准氨基酸 (标准残基 flag 是 ' ')
                    if het_flag == ' ' and resn in standard_aa_names:
                        continue
                    # 过滤掉核酸残基
                    if resn in nucleic_resnames:
                        continue
                    # 过滤掉太小的东西
                    if len(residue) < min_ligand_atoms:
                        continue

                    # 剩下的就视为“配体的一部分”
                    found_ligands.add(resn)
                    for atom in residue:
                        # 坐标已经是变换后的
                        ligand_coords_tgt.append(atom.get_coord())

        if not ligand_coords_tgt:
            print("[ERROR] No ligand-like residues found (non-water, non-AA, size>=min_ligand_atoms).")
            return None, None

        print(f"[INFO] Dynamically identified reference ligands: {', '.join(found_ligands)}")
        ligand_coords_tgt = np.array(ligand_coords_tgt, dtype=float)

        # 5. 计算盒子中心和尺寸
        min_coords = ligand_coords_tgt.min(axis=0)
        max_coords = ligand_coords_tgt.max(axis=0)

        center = (min_coords + max_coords) / 2.0
        size = (max_coords - min_coords) + box_margin # 使用 box_margin 参数

        return center, size

    except Exception as e:
        print(f"[FATAL ERROR] Alignment/Pocket calculation failed: {e}")
        return None, None
    
    
def collect_tasks(gene_ligand_map, gene_receptor_map,ligand_path_map, ref_pdb_dir):
    tasks = []
    for gene_name, receptors in gene_receptor_map.items():

        if gene_name not in gene_ligand_map:
            print(f"[WARN] Gene {gene_name} found in receptor map but missing in ligand map. Skipping.")
            continue
        ligand_names = (gene_ligand_map[gene_name].get('substrate', []) +  gene_ligand_map[gene_name].get('product', []))
        reference = gene_ligand_map[gene_name].get('reference', None)
        reference_path = os.path.join(ref_pdb_dir, f"{reference}.pdb") if reference else None
        if not ligand_names:
            print(f"[WARN] No substrates or products defined for scripts {gene_name}. Skipping.")
            continue

        for receptor_info in receptors:
            receptor_pdb_file = receptor_info['pdb_file']
            receptor_pdbqt_file = receptor_info['pdbqt_file']
            receptor_name = receptor_info['species']
            for ligand_name in ligand_names:
                if ligand_name in ligand_path_map:
                    ligand_pdbqt_file = ligand_path_map[ligand_name]
                    tasks.append((receptor_pdb_file, receptor_pdbqt_file, ligand_pdbqt_file, gene_name, ligand_name, receptor_name,reference_path))
                else:
                    print(f"[WARN] Missing ligand PDBQT: {ligand_pdbqt_file}")
    return tasks

def run_task(task, docking_dir, pocket_dir, VINA_BIN, FPOCKET_BIN,cpu_per_task):
    # VINA_BIN = shutil.which('vina')
    # if VINA_BIN is None:
    #     raise FileNotFoundError("vina not found in PATH. Please check your Conda environment.")

    receptor_pdb_file, receptor_pdbqt_file, ligand_pdbqt_file, gene_name, ligand_name, receptor_name, reference_path= task
    docking_result = os.path.join(docking_dir, gene_name, f"{gene_name}__{receptor_name}__{ligand_name}")
    os.makedirs(docking_result, exist_ok=True)
    log_file = os.path.join(docking_result, "result.log")
    if os.path.isfile(log_file):
        aff = best_affinity(Path(log_file))
        if aff is not None:
            return 
    os.makedirs(docking_result, exist_ok=True)

    if reference_path and os.path.exists(reference_path):
        center, size = align_and_get_pocket_params(receptor_pdb_file, reference_path)
    if center is None or size is None:
        fpocket_out = os.path.join(pocket_dir, gene_name, f"{receptor_name}")
        fpocket_out = run_fpocket(receptor_pdb_file, fpocket_out, receptor_name,FPOCKET_BIN) 
        if fpocket_out is None:
            print(f"[ERROR] fpocket failed: {receptor_name}")
            return
        else:
            center, size = parse_fpocket_pqr(fpocket_out)
            if center is None:
                print(f"[ERROR] fpocket parsing failed: {receptor_name}")
                return

        


    vina_config = os.path.join(docking_result, "vina_config.txt")
    with open(vina_config, 'w') as f:
        f.write(f"receptor = {receptor_pdbqt_file}\n")
        f.write(f"ligand = {ligand_pdbqt_file}\n")
        f.write(f"center_x = {center[0]:.3f}\ncenter_y = {center[1]:.3f}\ncenter_z = {center[2]:.3f}\n")
        f.write(f"size_x = {size[0]:.3f}\nsize_y = {size[1]:.3f}\nsize_z = {size[2]:.3f}\n")
        f.write(f"out = {docking_result}/result.pdbqt\nlog = {docking_result}/result.log\n")
        f.write(f"exhaustiveness = 32\nnum_modes = 20\n")

    try:
        subprocess.run([VINA_BIN, "--config", vina_config, "--cpu", str(cpu_per_task)], check=True)
        print(f"[OK] Docking completed: {gene_name}_{receptor_name} x {ligand_name}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] vina failed: {e}")




if __name__ == "__main__":
    # output_file = '/home/shiyi/Gene2Struct-main/scripts/DockingModule/output/pockets/UGT3/rhodiola_himalensis'
    # center, size = parse_fpocket_pqr(output_file)
    # print(center)
    # print(size)
    target_pdb_path = '/home/shiyi/Gene2Struct-main/scripts/DockingModule/output/receptors/UGT3/rhodiola_himalensis.pdb'
    reference_pdb_path = '/home/shiyi/Gene2Struct-main/scripts/DockingModule/cif_data/8ITA.pdb'
    center,size = align_and_get_pocket_params(target_pdb_path, reference_pdb_path)
    print(center)
    print(size)
# class DockingExecutor:
#     def __init__(self, receptor_pdb_dir, receptor_pdbqt_dir, ligand_pdbqt_dir, docking_dir,temp_center, gene_ligand_map):
#         self.receptor_pdbqt_dir = receptor_pdbqt_dir
#         self.ligand_pdbqt_dir = ligand_pdbqt_dir
#         self.docking_dir = docking_dir
#         self.gene_ligand_map = gene_ligand_map
#         self.receptor_pdb_dir =receptor_pdb_dir
#         self.temp_center = temp_center

#         self.VINA_BIN = shutil.which('vina')
#         if self.VINA_BIN is None:
#             raise FileNotFoundError("vina not found in PATH. Please check your Conda environment.")
        
#     def collect_tasks(self):
#         tasks = []
#         for root, _, files in os.walk(self.receptor_pdb_dir):
#             for file in files:
#                 if file.endswith('.pdb'):
#                     receptor_pdb_file = os.path.join(root, file)
#                     gene_name = Path(root).name  

#                     if gene_name not in self.gene_ligand_map:
#                         continue


#                     relative_subfolder = os.path.relpath(root, self.receptor_pdb_dir)
#                     receptor_pdbqt_file = os.path.join(self.receptor_pdbqt_dir, relative_subfolder, file.replace(".pdb", ".pdbqt"))
#                     if not Path(receptor_pdbqt_file).exists():
#                         print(f"[ERROR] Missing receptor PDBQT: {receptor_pdbqt_file}")
#                         continue

#                     ligands = self.gene_ligand_map[gene_name]
#                     for ligand_name in ligands:
#                         ligand_pdbqt_file = os.path.join(self.ligand_pdbqt_dir, f"{ligand_name}.pdbqt")
#                         if Path(ligand_pdbqt_file).exists():
#                             tasks.append((receptor_pdb_file, receptor_pdbqt_file, ligand_pdbqt_file, gene_name, ligand_name))
#                         else:
#                             print(f"[WARN] Missing ligand PDBQT: {ligand_pdbqt_file}")
#         return tasks

#     def run_task(self, args):
#         receptor_pdb_file, receptor_pdbqt_file, ligand_pdbqt_file, gene_name, ligand_name = args
#         receptor_name = Path(receptor_pdb_file).stem

#         fpocket_dir = os.path.join(self.temp_center, gene_name, f"{receptor_name}_out")
#         if Path(fpocket_dir).exists():
#             fpocket_out = fpocket_dir
#             center, size = parse_fpocket_pqr(fpocket_out)
#             if center is None:
#                 fpocket_out = run_fpocket(str(receptor_pdb_file), str(self.temp_center))
#                 if not fpocket_out:
#                     print(f"[ERROR] fpocket failed: {receptor_name}")
#                     return
#                 center, size = parse_fpocket_pqr(fpocket_out)
#                 if center is None:
#                     print(f"[ERROR] fpocket parsing failed: {receptor_name}")
#                     return
#         else:
#             fpocket_out = run_fpocket(str(receptor_pdb_file), str(self.temp_center), gene_name)
#             if not fpocket_out:
#                 print(f"[ERROR] fpocket failed: {receptor_name}")
#                 return
#             center, size = parse_fpocket_pqr(fpocket_out)
#             if center is None:
#                 print(f"[ERROR] fpocket parsing failed: {receptor_name}")
#                 return

#         out_dir = os.path.join(self.docking_dir, gene_name, f"{gene_name}__{receptor_name}__{ligand_name}")
#         os.makedirs(out_dir, exist_ok=True)

#         vina_config = os.path.join(out_dir, "vina_config.txt")
#         with open(vina_config, 'w') as f:
#             f.write(f"receptor = {receptor_pdbqt_file}\n")
#             f.write(f"ligand = {ligand_pdbqt_file}\n")
#             f.write(f"center_x = {center[0]:.3f}\ncenter_y = {center[1]:.3f}\ncenter_z = {center[2]:.3f}\n")
#             f.write(f"size_x = {size[0]:.3f}\nsize_y = {size[1]:.3f}\nsize_z = {size[2]:.3f}\n")
#             f.write(f"out = {out_dir}/result.pdbqt\nlog = {out_dir}/result.log\n")
#             f.write(f"exhaustiveness = 32\nnum_modes = 20\n")

#         try:
#             subprocess.run([self.VINA_BIN, "--config", vina_config], check=True)
#             print(f"[OK] Docking completed: {gene_name} x {ligand_name}")
#         except subprocess.CalledProcessError as e:
#             print(f"[ERROR] vina failed: {e}")

#     def run_all(self):
#         for root, _, files in os.walk(self.receptor_pdb_dir):
#             for file in files:
#                 if file.endswith('.pdb'):
#                     receptor_pdb_file = os.path.join(root, file)
#                     receptor_name = Path(file).stem
#                     gene_name = os.path.basename(root)
#                     fpocket_dir = os.path.join(self.temp_center, f"{receptor_name}_out")
#                     if not Path(fpocket_dir).exists():
#                         print(f"[INFO] Running fpocket: {receptor_name}")
#                         run_fpocket(receptor_pdb_file, self.temp_center, gene_name)

#         tasks = self.collect_tasks()
#         with Pool(cpu_count()) as pool:
#             pool.map(self.run_task, tasks)

