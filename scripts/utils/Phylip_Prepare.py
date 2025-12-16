import os
from pdb import run
import sys
import shutil
import warnings
from pathlib import Path
from Bio import SeqIO
from Bio import BiopythonDeprecationWarning
from scripts.utils.TreeLoad import _sanitize_id_strict, parse_tree
# from scripts.utils.TreeLoad import TreeLoad  # 你已有的模块
import subprocess
from Bio.SeqRecord import SeqRecord
import re
from Bio import SeqIO
import unicodedata




def check_cds(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        raise ValueError(f"{fasta_file} contains no valid sequences")

    cleaned_records = []
    stops = {"TAA", "TAG", "TGA"}
    id_counter = {}
    for rec in records:
        seq = str(rec.seq).upper().replace("-", "")
        if len(seq) % 3 != 0:
            raise ValueError(f"{rec.id}: sequence length is not a multiple of 3")

        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]

        # Check premature stops
        if any(codon in stops for codon in codons[:-1]):
            premature_stop_index = next((i for i, codon in enumerate(codons[:-1]) if codon in stops), -1)
            raise ValueError(f"{rec.id}: contains a premature stop codon at postion {premature_stop_index * 3 + 1}")

        # Remove terminal stop codon if present
        if codons[-1] in stops:
            print(f"[Warning] {rec.id}: terminal stop codon removed")
            seq = seq[:-3]

        rec.seq = rec.seq.__class__(seq)  # keep Seq type

        header = rec.description if rec.description else rec.id
        parts = header.split('|')
        base_name = parts[0].strip() if parts else header
        cleaned_base_id = _sanitize_id_strict(base_name)
        final_id = cleaned_base_id
        if cleaned_base_id in id_counter:
            id_counter[cleaned_base_id] += 1
            final_id = f"{cleaned_base_id}_{id_counter[cleaned_base_id] - 1}"
        else:
            id_counter[cleaned_base_id] = 1

        rec.id = rec.name = final_id
        rec.description = ""  # 关键：清空描述，避免空格再次混入
        cleaned_records.append(rec)

    return cleaned_records

def run_pal2nal_revised(cds_fasta, output_prefix):
    cds_path = Path(cds_fasta)
    out_prefix = Path(output_prefix)
    aa_fasta_path = out_prefix.with_suffix(".aa.fasta")
    aa_aln_path = out_prefix.with_suffix(".aa.aln.fasta")
    codon_alignment_path = out_prefix.with_suffix(".codon.aln.fasta")

    aa_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        # 1. Translate CDS → Protein
        new_records = []
        for rec in SeqIO.parse(cds_path, "fasta"):
            aa_seq = rec.seq.translate(to_stop=True)
            new_rec = SeqRecord(aa_seq, id=rec.id, description="")
            new_records.append(new_rec)
        SeqIO.write(new_records, aa_fasta_path, "fasta")


        # 2. Muscle alignment (protein)
        muscle_path = shutil.which("muscle")
        if not muscle_path:
            raise FileNotFoundError(
            "MUSCLE not found. Please install MUSCLE and make sure it is on your PATH, "
            "or set the MUSCLE environment variable to the absolute path of the executable.\n"
            "Examples:\n"
            "  conda install -c bioconda muscle\n"
            "  export MUSCLE=/path/to/muscle"
    )
        cmd_mafft = [muscle_path, "-align", str(aa_fasta_path), '-output', aa_aln_path]
        try:
            subprocess.run(cmd_mafft, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"MUSCLE failed with exit code {e.returncode}.\n"
                f"Command: {' '.join(cmd_mafft)}\n"
                f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
            ) from e
        
        if not aa_aln_path.exists() or aa_aln_path.stat().st_size == 0:
            raise RuntimeError("MUSCLE produced no output alignment file.")


        # 3. PAL2NAL: map back to codon alignment
        pal2nal_path = shutil.which("pal2nal.pl")
        if not pal2nal_path:
            raise FileNotFoundError(
                "PAL2NAL not found. Please install PAL2NAL and make sure it is on your PATH, "
                "or set the PAL2NAL environment variable to the absolute path of the executable.\n"
                "Examples:\n"
                "  conda install -c bioconda pal2nal\n"
            )
        with open(codon_alignment_path, "w") as fout:
            cmd_p2n = [pal2nal_path, str(aa_aln_path), str(cds_path), "-output", "fasta"]
            try:
                subprocess.run(cmd_p2n,stdout=fout,check=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"PAL2NAL failed. Exit code {e.returncode}.\n"
                    f"Command: {' '.join(cmd_p2n)}\n"
                    f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
                ) from e
        if not codon_alignment_path.exists() or codon_alignment_path.stat().st_size == 0:
            raise RuntimeError("PAL2NAL produced no output codon alignment file.")
        return str(codon_alignment_path)
    finally:
        # Clean up temporary files
        if aa_fasta_path.exists():
            aa_fasta_path.unlink()


def run_trimal(input_fasta, output_phy, mode='automated1'):
    trimal_exe = shutil.which('trimal')
    if not trimal_exe:
        raise FileNotFoundError(
            "trimalnot found. Please install trimAl and ensure it is on your PATH, "
            "or set the TRIMAL environment variable to the absolute path of the executable.\n"
            "Examples:\n"
            "  conda install -c bioconda trimal\n"
            "  export TRIMAL=/path/to/trimal"
        )
    if mode == "automated1":
        cmd = [trimal_exe, "-in", str(input_fasta), "-out", str(output_phy), "-phylip_paml", "-gt", "0.8"]
    else:
        cmd = [trimal_exe, "-in", str(input_fasta), "-out", str(output_phy), "-phylip_paml", "-gt", str(mode)]
    subprocess.run(cmd, check=True)
    return Path(output_phy)

def convert(fasta_file, phylip_file):
    trimal_exe = shutil.which('trimal')
    if not trimal_exe:
        raise FileNotFoundError(
            "trimalnot found. Please install trimAl and ensure it is on your PATH, "
            "or set the TRIMAL environment variable to the absolute path of the executable.\n"
            "Examples:\n"
            "  conda install -c bioconda trimal\n"
            "  export TRIMAL=/path/to/trimal"
        )
    cmd = [trimal_exe, "-in", str(fasta_file), "-out", str(phylip_file), "-phylip"]
    subprocess.run(cmd, check=True)
    return Path(phylip_file)

def check_species_match_before_alignment(fasta_file, species_from_tree):
    records = list(SeqIO.parse(fasta_file, "phylip-relaxed"))
    if not records:
        raise ValueError(f"{fasta_file} contains no sequences")

    seq_ids = {rec.id for rec in records}
    tree_species = set(species_from_tree)

    if seq_ids != tree_species:
        missing_in_seq = tree_species - seq_ids
        missing_in_tree = seq_ids - tree_species
        msg = []
        if missing_in_seq:
            msg.append(f"Missing in fasta: {missing_in_seq}")
        if missing_in_tree:
            msg.append(f"Missing in tree: {missing_in_tree}")
        raise ValueError("Tree and fasta species do not match!\n" + "\n".join(msg))
    else:
        print("[Check] Sequence IDs match tree species.")


def prepare_paml_input(fasta_path, output_dir, tree_path=None):
    # output_dir = "../file_input"
    fasta_name = Path(fasta_path).stem
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Checking CDS and sanitizing FASTA IDs...")
    clean_id_fasta = output_dir / "cleaned.fasta"
    cleaned_records = check_cds(fasta_path)
    SeqIO.write(cleaned_records, clean_id_fasta, "fasta")
    print("Aligning and trimming...")
    codon_alignment_path= run_pal2nal_revised(clean_id_fasta, output_dir / fasta_name)
    if Path(output_dir / f"{fasta_name}.phy").exists():
        fasta_phy_path = output_dir / f"{fasta_name}.phy"
    else:
        fasta_phy_path = run_trimal(codon_alignment_path,output_dir / f"{fasta_name}.phy")
    # fasta_phy_path = convert(codon_alignment_path, output_dir / f"{fasta_name}.phy")
    # if tmp_clean_fasta.exists():
    #     tmp_clean_fasta.unlink()

    # --- tree ---
    print("Locating or generating phylogeny...")
    paml_treefile = output_dir / f"{fasta_name}.paml.tree"
    parsed_tree_data = None
    tree_candidates = []
    if tree_path:
        tree_candidates.append(Path(tree_path))
    tree_candidates.extend([paml_treefile, output_dir / f"{fasta_name}.treefile"])

    for cand_path in tree_candidates:
        if Path(cand_path).exists() and Path(cand_path).stat().st_size > 0:
            print(f"[Status] Found existing tree: {cand_path.name}")
            species, newick_tree = parse_tree(cand_path)
            parsed_tree_data = (species, newick_tree)
            break  # 找到第一个有效的就跳出

    if not parsed_tree_data:
        print("[Status] Tree file not found. Running IQ-TREE 2 to generate phylogeny...")
        cmd = ["iqtree", "-s", str(codon_alignment_path), "-m", "MFP", "-fast", '-pre', str(output_dir / fasta_name)]
        try:
            subprocess.run(cmd, cwd=output_dir, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"IQ-TREE 2 failed: {e}")

        raw_tree = output_dir / f"{fasta_name}.treefile"
        species, newick_tree = parse_tree(raw_tree)
        parsed_tree_data = (species, newick_tree)
        for suffix in ["model.gz", "log", "iqtree", "ckp.gz", "mldist", "bionj"]:
            f = output_dir / f"{fasta_name}.{suffix}"
            if f.exists():
                f.unlink()

    species, newick_tree = parsed_tree_data
    n_species = len(species)
    header = f"   {n_species}  1\n"
    with open(paml_treefile, 'w') as f:
        f.write(header)
        f.write(newick_tree)
    print("Checking species consistency...")
    check_species_match_before_alignment(fasta_phy_path, species)

    print(f"[Success] PAML input prepared: {fasta_phy_path.name} and {paml_treefile.name}")
    return fasta_phy_path, paml_treefile, species



#
# if __name__ == "__main__":
#     file = '/Users/sy/Library/Containers/com.tencent.xinWeChat/Data/Documents/xwechat_files/wxid_567rvuytb70a22_22c0/msg/file/2025-11/MIKC1.fa'
#     fasta_phy_path, paml_treefile, species = prepare_paml_input(file, '/Users/sy/Gene2Struct-main/scripts/evolution/test1')
#     # file = '/Users/sy/Library/Containers/com.tencent.xinWeChat/Data/Documents/xwechat_files/wxid_567rvuytb70a22_22c0/msg/file/2025-11/MIKC(1).fa'
#     # fasta_phy_path, paml_treefile, species = prepare_paml_input(file, '/Users/sy/Gene2Struct-main/scripts/evolution/test2')
#     # file = '/Users/sy/Library/Containers/com.tencent.xinWeChat/Data/Documents/xwechat_files/wxid_567rvuytb70a22_22c0/msg/file/2025-11/Type_1.1.fa'
#     # fasta_phy_path, paml_treefile, species = prepare_paml_input(file, '/Users/sy/Gene2Struct-main/scripts/evolution/test')
if __name__ == '__main__':
    new_records = []
    for rec in SeqIO.parse('/Users/sy/Library/Containers/com.tencent.xinWeChat/Data/Documents/xwechat_files/wxid_567rvuytb70a22_22c0/msg/file/2025-11/A_E.fa', "fasta"):
        aa_seq = rec.seq.translate(to_stop=True)
        new_rec = SeqRecord(aa_seq, id=rec.id, description="")
        new_records.append(new_rec)
    SeqIO.write(new_records, "/Users/sy/Library/Containers/com.tencent.xinWeChat/Data/Documents/xwechat_files/wxid_567rvuytb70a22_22c0/msg/file/2025-11/A_E.aa.fa", "fasta")