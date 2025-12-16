import requests
import pandas as pd
import numpy as np
from Bio import SeqIO
import os
from pathlib import Path
from scripts.utils.TreeLoad import _sanitize_id_strict


API_SCORE_URL   = "http://life-bioinfo.cn:3423/score_only"
API_ENTROPY_URL = "http://life-bioinfo.cn:3423/entropy_only"

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
def scoring(fasta_file, output_dir):
    gene_name = Path(fasta_file).stem
    output_file = Path(output_dir) / f"{gene_name}_NLL_score.csv"
    if output_file.exists() and output_file.stat().st_size > 0:
        return output_file
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

    names = [r.id for r in records]
    seqs = [str(r.seq).replace('-', '').replace('?', '').upper() for r in records]

    request_data = {'seq_contents': seqs}
    response = requests.post(API_SCORE_URL, json=request_data)

    if response.status_code == 200:
        scores = response.json()['results']
        score_df = pd.DataFrame({'name': names, gene_name: [-round(score, 2) for score in scores]})
        score_df.to_csv(output_file, index=False)
        return output_file
    else:
        raise ValueError(f"API Error {response.status_code}: {response.text}")



# def scoring(fasta_input, output_dir):
#     output_dir = Path(output_dir)
#     output_dir.mkdir(parents=True, exist_ok=True)
#
#     fasta_files = collect_fasta_files(fasta_input)
#     generated = []
#     for fasta_path in fasta_files:
#         fasta_path = str(fasta_path)
#         gene_name = Path(fasta_path).stem
#         per_gene_csv = Path(output_dir) / f"{gene_name}_NLL_score.csv"
#         if per_gene_csv.exists() and per_gene_csv.stat().st_size > 0:
#             generated.append(str(per_gene_csv))
#             continue
#         records, used_fmt = [], None
#         if fasta_path.endswith('.phy'):
#             for fmt in ('phylip-relaxed', 'phylip'):
#                 try:
#                     records = list(SeqIO.parse(fasta_path, fmt))
#                     if records:
#                         used_fmt = fmt
#                         break
#                 except Exception:
#                      pass
#         if not records:
#             try:
#                 records = list(SeqIO.parse(fasta_path, 'fasta'))
#                 used_fmt = 'fasta'
#             except Exception:
#                 pass
#         if not records:
#             raise ValueError(f"Failed to parse {fasta_path} as phylip or fasta.")
#         cleaned = []
#         for r in records:
#             base = clean_header_to_base_id(r.description or r.id)
#             r.id = r.name = base
#             r.description = ""
#             cleaned.append(r)
#         records = cleaned
#
#         names = [r.id for r in records]
#         seqs = [str(r.seq).replace('-', '').replace('?', '').upper() for r in records]
#
#         request_data = {'seq_contents': seqs}
#         response = requests.post(API_SCORE_URL, json=request_data)
#
#         if response.status_code == 200:
#             scores = response.json()['results']
#             per_gene_df = pd.DataFrame({'name': names, gene_name: [-round(score, 2) for score in scores]})
#             per_gene_df.to_csv(per_gene_csv, index=False)
#             generated.append(str(per_gene_csv))
#         else:
#             raise ValueError(f"API Error {response.status_code}: {response.text}")
#     return generated


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
    
def position_entropy(fasta_file, output_dir):
    # otuput_dir = ../evo2output
    gene_name = Path(fasta_file).stem
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records =[]

    fmt = ["phylip-relaxed", "phylip", "fasta"]
    for fmt in fmt:
        try:
            with open(fasta_file, 'r') as handle:
                records = list(SeqIO.parse(handle, fmt))
            if records and len(records) > 0 and len(records[0].seq) > 0:
                break
        except Exception:
            continue
    if not records:
        raise ValueError(f"Could not parse the multiple sequence alignment file {fasta_file}. Please ensure it is a valid FASTA or Phylip format.")

    names = [rec.id for rec in records]
    gapped_seqs = [str(rec.seq).upper() for rec in records]

    ungapped_seqs = [s.replace('-', '').replace('.', '').replace('*', '') for s in gapped_seqs]
    request_data = {'seq_contents': ungapped_seqs}
    response = requests.post(API_ENTROPY_URL, json=request_data)
    if response.status_code == 200:
        entropies = response.json()['entropies']
        conservation = [[round(2 - entropy, 4) for entropy in entropy_sequence] for entropy_sequence in entropies]
        aligned = {}
        for i in range(len(names)):
            name = names[i]
            seq = gapped_seqs[i]
            score = conservation[i]
            aligned[name] = align_scores_to_msa(seq, score)

        df_heatmap = pd.DataFrame.from_dict(aligned, orient="index")
        df_heatmap.columns = [f'{i+1}' for i in df_heatmap.columns]
        out_file = os.path.join(output_dir, f"{gene_name}_entropy.csv")
        df_heatmap.to_csv(out_file, index=True)
        return out_file
    else:
        raise ValueError(f"API Error {response.status_code}: {response.text}")


if __name__ == "__main__":
    file = '/Users/sy/Gene2Struct-main/scripts/EvoDnDsModule/out/psbB/file_input/psbB.phy'
    print(scoring(file, '/Users/sy/Gene2Struct-main/scripts/EvoDnDsModule/out/psbB/result'))