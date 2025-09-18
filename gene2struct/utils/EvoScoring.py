import requests
import pandas as pd
import numpy as np
from Bio import SeqIO
import os
from pathlib import Path
import re


API_SCORE_URL   = "https://681030def768.ngrok-free.app/score_only"
API_ENTROPY_URL = "https://681030def768.ngrok-free.app/entropy_only"

def collect_fasta_files(fasta_input):

    if isinstance(fasta_input, str):
        p = Path(fasta_input)
        if p.is_dir():
            return [str(p / f) for f in os.listdir(p) if f.endswith(('.fa', '.fasta', '.fas', '.txt', '.phy'))]
        else:
            return fasta_input.split(';')
    elif isinstance(fasta_input, list):
        return fasta_input
    else:
        raise ValueError("Invalid fasta input")
def clean_header_to_base_id(header: str) -> str:
    header = header or ""
    parts = header.split('|')
    species = parts[0].strip().replace(' ', '_')
    pid = None
    for p in parts[1:]:
        if p.startswith('protein_id:'):
            pid = p.split(':', 1)[1].strip()
            break
    base = f"{species}_{pid}" if pid else species

    return re.sub(r'[^A-Za-z0-9_.-]', '_', base)

def remove_outgroups(fasta_file, outgroups=None):

    outgroups = outgroups or []
    records, used_fmt = [], None

    if fasta_file.endswith('.phy'):
        for fmt in ('phylip-relaxed', 'phylip'):
            try:
                records = list(SeqIO.parse(fasta_file, fmt))
                if records:
                    used_fmt = fmt
                    break
            except Exception:
                pass
    if not records:
        try:
            records = list(SeqIO.parse(fasta_file, 'fasta'))
            used_fmt = 'fasta'
        except Exception:
            pass
    if not records: 
        raise ValueError(f"Failed to parse {fasta_file} as phylip or fasta.")


    cleaned = []
    for r in records:
        base = clean_header_to_base_id(r.description or r.id)
        r.id = r.name = base
        r.description = ""
        cleaned.append(r)
    records = cleaned

    if used_fmt == "phylip":
        og_set = {og[:10] for og in (outgroups or [])}
    else:
        og_set = set(outgroups or [])

    filtered_records = [r for r in records if r.id not in og_set]
    if not filtered_records:
        raise ValueError("No sequences left after removing outgroups. Please check IDs.")

    names = [r.id for r in filtered_records]

    seqs = [str(r.seq).replace('-', '').replace('?', '').upper() for r in filtered_records]

    return names, seqs, filtered_records


def scoring(fasta_input, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = collect_fasta_files(fasta_input)
    generated = []
    for fasta_path in fasta_files:
        fasta_path = str(fasta_path)
        gene_name = Path(fasta_path).stem
        per_gene_csv = Path(output_dir) / f"{gene_name}_NLL_score.csv"
        if per_gene_csv.exists() and per_gene_csv.stat().st_size > 0:
            generated.append(str(per_gene_csv))
            continue
        records, used_fmt = [], None
        if fasta_path.endswith('.phy'):
            for fmt in ('phylip-relaxed', 'phylip'):
                try:
                    records = list(SeqIO.parse(fasta_path, fmt))
                    if records:
                        used_fmt = fmt
                        break
                except Exception:
                     pass
        if not records:
            try:
                records = list(SeqIO.parse(fasta_path, 'fasta'))
                used_fmt = 'fasta'
            except Exception:
                pass
        if not records: 
            raise ValueError(f"Failed to parse {fasta_path} as phylip or fasta.")
        cleaned = []
        for r in records:
            base = clean_header_to_base_id(r.description or r.id)
            r.id = r.name = base
            r.description = ""   
            cleaned.append(r)
        records = cleaned

        names = [r.id for r in records]
        seqs = [str(r.seq).replace('-', '').replace('?', '').upper() for r in records]

        request_data = {'seq_contents': seqs}
        response = requests.post(API_SCORE_URL, json=request_data)

        if response.status_code == 200:
            scores = response.json()['results']
            per_gene_df = pd.DataFrame({'name': names, gene_name: [-round(score, 2) for score in scores]})
            per_gene_df.to_csv(per_gene_csv, index=False)
            generated.append(str(per_gene_csv))
        else:
            raise ValueError(f"API Error {response.status_code}: {response.text}")
    return generated


def align_scores_to_msa(msa_seq, score_list):
    result = []
    idx = 0
    for aa in msa_seq:
        if aa == "-":
            result.append(np.nan)
        else:
            result.append(score_list[idx])
            idx += 1
    return result
    
def position_entropy(fasta_file, output_dir, outgroups=None):
    gene_name = Path(fasta_file).stem
    names, seqs, records = remove_outgroups(fasta_file, outgroups=outgroups)
    request_data = {'seq_contents': seqs}
    response = requests.post(API_ENTROPY_URL, json=request_data)
    if response.status_code == 200:
        entropies = response.json()['entropies']
        conservation = [[round(2 - entropy, 4) for entropy in entropy_sequence] for entropy_sequence in entropies]

        seqs = [str(record.seq).upper() for record in records]
        aligned = {
            name: align_scores_to_msa(seq, score)
            for name, seq, score in zip(names, seqs, conservation)
        }

        df_heatmap = pd.DataFrame.from_dict(aligned, orient="index")
        os.makedirs(output_dir, exist_ok=True)
        out_file = os.path.join(output_dir, f"{gene_name}_entropy.csv")
        df_heatmap.to_csv(out_file, index=True)
        return out_file
    else:
        raise ValueError(f"API Error {response.status_code}: {response.text}")
