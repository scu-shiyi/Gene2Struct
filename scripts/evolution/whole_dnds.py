from scripts.utils.site_model import run_pair_model
from pathlib import Path
import os
from scripts.utils.Phylip_Prepare import prepare_paml_input
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scripts.utils.EvoScoring import scoring

def parse_tree_branch_values(tree_text: str) -> dict[str, float]:
    """
    Extract branch values for leaves from a dS tree or dN tree.
    Returns {species: value}.
    """
    out = {}
    # Match pattern: SpeciesName: number
    pat = re.compile(r'([\w\.\'\-]+):\s*([0-9.eE+-]+)')
    for sp, val in pat.findall(tree_text):
        try:
            out[sp] = float(val)
        except ValueError:
            continue
    return out


def get_omega_from_freeratio_mlc(freeratio_mlc: Path, species: list[str],min_ds: float = 5e-4,min_dn: float = 1e-3) -> dict[str, float | None]:
    """
    Extract per-species ω (dN/dS) from a FreeRatio .mlc file,
    and filter unreliable results using dN/dS tree values.
    """
    text = Path(freeratio_mlc).read_text(errors="ignore")

    # 1) Capture "w ratios as node labels"
    omega_map = {}
    for sp in species:
        pat = re.compile(rf'{re.escape(sp)}\s*#\s*([0-9.eE+-]+)')
        m = pat.search(text)
        if m:
            try:
                omega_map[sp] = float(m.group(1))
            except ValueError:
                omega_map[sp] = None
        else:
            omega_map[sp] = None

    # 2) Capture dS tree
    ds_match = re.search(r'dS tree:\s*(\(.+?\));', text, re.S)
    ds_map = parse_tree_branch_values(ds_match.group(1)) if ds_match else {}

    # 3) Capture dN tree
    dn_match = re.search(r'dN tree:\s*(\(.+?\));', text, re.S)
    dn_map = parse_tree_branch_values(dn_match.group(1)) if dn_match else {}

    # 4) Filtering
    out = {}
    for sp in species:
        omega = omega_map.get(sp)
        ds = ds_map.get(sp, None)
        dn = dn_map.get(sp, None)

        if omega is None or not np.isfinite(omega):
            out[sp] = None
            continue

        # Discard if dS or dN is too small
        if ds is None or dn is None:
            out[sp] = None
            continue
        if ds < 0.0009 or ds > 3:
            out[sp] = None
            continue
        if dn < 0.0009 or dn > 1:
            out[sp] = None
            continue
        if omega <= 0 or omega > 3:
            out[sp] = None
        else:
            out[sp] = omega

    return out

def parse_m0_omega(m0_mlc: str | Path) -> float | None:
    """
    Parse overall omega from an M0 .mlc file.
    Supports both 'omega (dN/dS) =' and 'w (dN/dS) =' formats.
    """
    text = Path(m0_mlc).read_text(errors="ignore")
    m = re.search(r'\bomega\s*\(dN/dS\)\s*=\s*([0-9.eE+-]+)', text)
    if not m:
        m = re.search(r'\bw\s*\(dN/dS\)\s*=\s*([0-9.eE+-]+)', text)
    if not m:
        return None
    try:
        return float(m.group(1))
    except ValueError:
        return None
    
def save_omega_row(tsv_path: Path, gene: str, species: list[str], omega_map,m0_omega,lrt_result):
    """
    Write one row to the summary table: first column = Gene,
    then species-specific ω values in the given order.
    If the file does not exist, create it with header. Missing values = blank.
    """
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    # Significance annotation
    if lrt_result and "p" in lrt_result:
        p_val = lrt_result["p"]
        if lrt_result["sig"].get("0.01"):
            sig = "**"
        elif lrt_result["sig"].get("0.05"):
            sig = "*"
        else:
            sig = "ns"
    else:
        p_val, sig = "", "ns"

    with tsv_path.open("w", encoding="utf-8") as f:
        f.write(f"name,{gene}\n")
        for sp in species:
            val = omega_map.get(sp)
            f.write(f"{sp},{'' if val is None else val}\n")
        f.write("\n")
        f.write(f"M0_omega,{'' if m0_omega is None else m0_omega}\n")
        f.write("\n")
        f.write(f"LRT_p,{'' if p_val == '' else p_val}\n")
        f.write(f"LRT_sig,{sig}\n")

def run_one_gene_freeratio(fasta_file: str, base_output: str, tree_file: str | None = None,omp_threads: int = 1):
    """
    Run full FreeRatio analysis for a single scripts:
    prepare input → run FreeRatio → parse ω.
    Returns (gene_name, mlc_path, summary_row_path).
    """
    gene = Path(fasta_file).stem
    root = Path(base_output) / gene
    paml_input = root / 'file_input'
    paml_output_base = root / 'paml_output' / 'FreeRatio'
    result_dir = root / 'result'
    env = os.environ.copy()
    # Prevent codeml from overusing threads
    env["OMP_NUM_THREADS"] = str(omp_threads)

    # 1) Prepare input (remove outgroups and normalize IDs)
    phylip_path, newick_path, species = prepare_paml_input(
        fasta_file, str(paml_input), tree_path=tree_file
    )

    # 2) Run FreeRatio (model name must be uppercase "FREERATIO")
    lrt_result ,mo_mlc, free_mlc= run_pair_model(phylip_path, newick_path, str(paml_output_base), model="FREERATIO", gene_name=gene)

    # 3) Parse ω and save one row for this scripts
    omega_map = get_omega_from_freeratio_mlc(free_mlc, species)
    mo_omega = parse_m0_omega(mo_mlc)

    omega_csv = Path(result_dir) / f"{gene}_omega.csv"
    save_omega_row(omega_csv, gene, species, omega_map,mo_omega,lrt_result)
    evo_csv = scoring(fasta_file, result_dir)
    return gene, omega_csv, evo_csv

def collect_fasta_files(fasta_input):
    """Accept directory/semicolon-separated string/list, return fasta file list."""
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
    
def run_dnds_parallel(fasta_files, output_dir: str,max_workers: int | None = None,omp_threads_each_codeml: int = 4, mapping: dict[str, str | None] | None = None,):

    fasta_input = collect_fasta_files(fasta_files)
    output_dir = str(Path(output_dir).resolve())
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    norm_map: dict[str, str | None] = {}
    if mapping:
        for k, v in mapping.items():
            fk = str(Path(k).resolve())
            if v is None or (isinstance(v, str) and v.strip() == ""):
                norm_map[fk] = None
            else:
                tv = str(Path(v).resolve())
                if not Path(tv).exists():
                    print(f"Tree in mapping not found, fallback to auto-build: {fk} -> {tv}")
                    norm_map[fk] = None
                else:
                    norm_map[fk] = tv

    results = []
    evo_results = []
    with ProcessPoolExecutor(max_workers=max_workers or os.cpu_count()) as ex:
        futures = []
        for fa in fasta_input:
            fa_norm = str(Path(fa).resolve())
            tree_for_fa = norm_map.get(fa_norm) if norm_map else None
            futures.append(
                ex.submit(run_one_gene_freeratio, fa_norm, output_dir,tree_for_fa, omp_threads_each_codeml)
            )
        for fut in as_completed(futures):
            try:
                gene, per_gene_tsv, per_evo_csv = fut.result()
                print(gene, per_gene_tsv, per_evo_csv)
                results.append((gene, per_gene_tsv))
                evo_results.append(per_evo_csv)
            except Exception as e:
                print(f"Task failed: {e}")
    return results, evo_results

if __name__ == "__main__":
    fasta_files = ['/Users/sy/Gene2Struct-main/scripts/EvoDnDsModule/ndhB.fa','/Users/sy/Gene2Struct-main/scripts/EvoDnDsModule/psbB.fa']
    out = '/Users/sy/Gene2Struct-main/scripts/EvoDnDsModule/out'
    run_dnds_parallel(fasta_files, out)
    # run_one_gene_freeratio(fasta_files, out)