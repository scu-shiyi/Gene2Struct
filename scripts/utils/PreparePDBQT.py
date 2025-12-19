import subprocess
import os
from pathlib import Path
import requests
import os
from openbabel import pybel
import re
from urllib.parse import quote
from pathlib import Path
from scripts.mgltools import prepare_ligand4, prepare_receptor4
import sys
PYTHON = sys.executable



def run_subprocess(command):
	"""Run a subprocess command with error handling."""
	try:
		subprocess.run(command, check=True)
	except subprocess.CalledProcessError as e:
		raise ValueError(f"Conversion failed: {e}")


def ligand2pdbqt(ligand_input, ligand_dir):
	prepare_ligand_script = str(Path(prepare_ligand4.__file__).resolve())
	ligand_pdb_list = []
	if isinstance(ligand_input, (list, tuple, set)):
		for item in ligand_input:
			p = os.path.abspath(str(item))
			if not os.path.isfile(p):
				raise FileNotFoundError(f"ligand_input item is not a file: {p}")
			if not p.lower().endswith(".pdb"):
				raise ValueError(f"ligand_input item is not a .pdb file: {p}")
			ligand_pdb_list.append(p)
	else:
		ligand_input = os.path.abspath(ligand_input)
		if os.path.isfile(ligand_input):
			if not ligand_input.lower().endswith('.pdb'):
				raise ValueError("format converting from sdf to pdb has failed")
			ligand_pdb_list = [ligand_input]
		elif os.path.isdir(ligand_input):
			ligand_pdb_list = [os.path.join(ligand_input, f) for f in os.listdir(ligand_input) if
							   f.lower().endswith('.pdb')]
			if not ligand_pdb_list:
				raise FileNotFoundError(f"No .pdb files found in directory: {ligand_input}")
		else:
			raise FileNotFoundError(f"ligand_input does not exist: {ligand_input}")

	ligand_pdbqt_list = []
	for ligand_file in ligand_pdb_list:
		file_name = Path(ligand_file).stem + ".pdbqt"
		output_file = os.path.join(ligand_dir, file_name)
		command = [PYTHON, prepare_ligand_script,
				   "-l", ligand_file,
				   "-o", output_file,
				   "-U", "waters",
				   "-A", "hydrogens"
				   ]
		run_subprocess(command)
		ligand_pdbqt_list.append(output_file)
		print(f"{file_name} was successfully converted and saved to {output_file}")

	ligand_path_map = {}
	for path in ligand_pdbqt_list:
		ligand_name = Path(path).stem
		ligand_path_map[ligand_name] = path
	return ligand_path_map


def receptor2pdbqt(receptor_file, receptor_pdbqt_out):
	prepare_receptor_script = str(Path(prepare_receptor4.__file__).resolve())
	if os.path.exists(receptor_pdbqt_out):
		return receptor_pdbqt_out

	command = [PYTHON, prepare_receptor_script,
			   "-r", receptor_file,
			   "-o", receptor_pdbqt_out,
			   "-U", "waters",
			   "-A", "hydrogens"]
	run_subprocess(command)
	print(f"Receptor was successfully converted and saved to {receptor_pdbqt_out}")
	return receptor_pdbqt_out


def downloading(ligand_list, ligand_dir):
	NAME_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
	CID_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
	sdf_file_list = []
	for ligand in ligand_list:
		ligand = str(ligand).strip()
		out_file = os.path.join(ligand_dir, f"{ligand}.sdf")
		if os.path.exists(out_file):
			sdf_file_list.append(out_file)
			print(f"{ligand}.sdf already exists, skipping download.")
			continue
		if ligand.isdigit():
			cid = ligand
			file_name = f"{cid}.sdf"
		else:
			url = NAME_URL.format(name=quote(ligand))
			res = requests.get(url, timeout=10)
			res.raise_for_status()
			data = res.json()
			try:
				cid = str(data["IdentifierList"]["CID"][0])
			except (KeyError, IndexError):
				raise ValueError(f"cannot find CID for {ligand} from PubChem")
			safe_name = ligand.replace(" ", "_")
			file_name = f"{safe_name}.sdf"

		sdf_url = CID_URL.format(cid=cid)
		params = {"record_type": "3d"}
		res = requests.get(sdf_url, params=params, timeout=20)
		res.raise_for_status()

		out_path = os.path.join(ligand_dir, file_name)
		with open(out_path, "wb") as f:
			f.write(res.content)
		sdf_file_list.append(out_path)
		print(f"{file_name} has been downloaded!")

	ligand_pdb_list = []
	for sdf_file in sdf_file_list:
		name = Path(sdf_file).stem
		try:
			mol = next(pybel.readfile("sdf", sdf_file))
		except StopIteration:
			raise ValueError(f"Invalid sdf: {sdf_file}")
		pdb_file = os.path.join(ligand_dir, f"{name}.pdb")
		mol.write('pdb', pdb_file, overwrite=True)
		ligand_pdb_list.append(pdb_file)

	return ligand_pdb_list


def download_reference_pdb(reference_list, receptor_dir):
	ref_pdb_dir = os.path.join(receptor_dir, "reference_pdb")
	os.makedirs(ref_pdb_dir, exist_ok=True)
	for reference in reference_list:
		reference = str(reference).strip()
		if not reference:
			continue
		pdb_id = reference.upper()
		pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
		pdb_path = os.path.join(ref_pdb_dir, f"{pdb_id}.pdb")
		if os.path.exists(pdb_path):
			continue
		try:
			response = requests.get(pdb_url, timeout=10)
		except requests.RequestException as e:
			print(f"[ERROR] Failed to download {pdb_id}: {e}")
			continue
		if response.status_code != 200:
			print(f"[ERROR] Failed to download {pdb_id}: HTTP {response.status_code}")
			continue
		with open(pdb_path, "wb") as f:
			f.write(response.content)
	return ref_pdb_dir


def clean_filename(basename, gene_name):
	if gene_name:
		pattern = re.compile(rf"[_-]?{re.escape(gene_name)}[_-]?", flags=re.I)
		basename = pattern.sub("", basename)
	for tag in ['fold', 'model', 'structure', 'protein', 'scripts']:
		pattern = re.compile(rf"[_-]?{tag}(?:[_-]?\d+)?", flags=re.I)
		basename = pattern.sub("", basename)

	basename = re.sub(r"[_-]{2,}", "_", basename).strip("_-")
	if not basename:
		raise ValueError("Invalid receptor naming inside the protein directory.")
	return basename


def convert_cif_to_pdb(src_path, dst_path):
	try:
		mol = next(pybel.readfile("cif", src_path))
	except StopIteration:
		print(f"Invalid CIF file: {src_path}")
		return None

	mol.write("pdb", dst_path, overwrite=True)
	print(f"Converted {src_path} to {dst_path}")
	return dst_path


if __name__ == "__main__":
	name = "8ITA"
	download_reference_pdb([name], "/home/shiyi/Gene2Struct-main/scripts/DockingModule/final_out/receptors")