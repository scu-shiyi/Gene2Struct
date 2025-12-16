
import os
import shutil
import re
from pathlib import Path
from typing import List, Tuple
from Bio import SeqIO

def clean_sample_id(file_name: str) -> str:
    basename = Path(file_name).stem
    cleaned = re.sub(r'[\s\-]+', '_', basename)
    cleaned = re.sub(r'[^A-Za-z0-9_]+', '_', cleaned)
    cleaned = re.sub(r'_+', '_', cleaned).strip('_')
    return cleaned


def _filter_fasta(in_path: str, out_path: str, outgroups: List[str]) -> bool:
    records = list(SeqIO.parse(in_path, "fasta"))
    if outgroups:
        records = [r for r in records if all(og not in r.id for og in outgroups)]
    if not records:                 return False
    SeqIO.write(records, out_path, "fasta")
    return True


def preprocess_fasta_dir(
        fasta_files: List[str],
        output_dir: str,
        outgroups: List[str] | None = None,
) -> Tuple[List[str], str]:

    outgroups = outgroups or []

    processed_dir = os.path.join(output_dir, "processed_fasta")
    os.makedirs(processed_dir, exist_ok=True)


    processed_files = []
    for fasta_path in fasta_files:
        clean_gene_name = clean_sample_id(fasta_path)
        dst_path = os.path.join(processed_dir, f"{clean_gene_name}.fa")

        ok = _filter_fasta(str(fasta_path), str(dst_path), outgroups)
        if ok:
            processed_files.append(str(dst_path))

    processed_files = list(dict.fromkeys(processed_files)) 
    return processed_files, str(processed_dir)

