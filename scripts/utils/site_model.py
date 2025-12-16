from string import Template
from pathlib import Path
import subprocess
import shutil
import re
from scipy.stats import chi2
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

def _run_codeml_async(codeml_bin, ctl_path, work_dir):
    work_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "4")
    return subprocess.Popen([codeml_bin, str(ctl_path)], cwd=str(work_dir), env=env)


TEMPLATES = {
    "M0": Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 0
NSsites = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 3
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
fix_blength = 0        
method = 0
cleandata = 0'''),
    "M3": Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 0
NSsites = 3
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 3
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
fix_blength = 0        
method = 0
cleandata = 0'''),
    'M7': Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 0
NSsites = 7
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
fix_blength = 0
method = 0
cleandata = 0'''),
    'M8': Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 0
NSsites = 8
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
fix_blength = 0        
method = 0
cleandata = 0'''),
    'FREERATIO': Template('''seqfile = $seq_path
treefile = $tree_path
outfile = $output_mlc
noisy = 9
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
ndata = 1
clock = 0
icode = 0
Mgene = 0
model = 1
NSsites = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8
fix_rho = 1
rho = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 0
method = 0
fix_blength = 0''')
}


def _mlc_complete(p: Path) -> bool:

    try:
        if not p.exists() or p.stat().st_size < 200:
            return False
        txt = p.read_text(errors="ignore")
        return ("lnL(" in txt) and ("Time used" in txt)
    except Exception:
        return False
def write_ctl(model, seq_path, tree_path, mlc_path, ctl_path):
    model = model.upper()
    if model not in TEMPLATES:
        raise ValueError(f"Unknown model: {model}")
    text = TEMPLATES[model].substitute(
        seq_path=str(seq_path),
        tree_path=str(tree_path),
        output_mlc=str(mlc_path)
    )
    ctl_path.write_text(text)
    return ctl_path, mlc_path


def parse_mlc_np_lnl(path: Path):

    NUM = r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'
    PAT = re.compile(
    rf'lnL\s*\(\s*ntime:\s*\d+\s+np:\s*(\d+)\s*\)\s*:\s*({NUM})')

    text = Path(path).read_text(errors="ignore")
    matches = list(PAT.finditer(text))
    if not matches:
        raise ValueError(f"无法在 {path} 里找到 lnL(ntime: .. np: ..): 行。")
    m = matches[-1]
    np_ = int(m.group(1))
    lnL = float(m.group(2))
    return lnL,np_

def lrt(null_lnl, null_np, alt_lnl, alt_np):
    """
    Likelihood Ratio Test: H0 = simpler model (null), H1 = more complex model (alt)
    Inputs: lnL and np from .mlc files
    Returns: test statistic, degrees of freedom, p-value, significance decision
    """
    df = alt_np - null_np
    if df <= 0:
        raise ValueError(f"Degrees of freedom (df)={df} <= 0, please check input np values.")
    stat = 2.0 * (alt_lnl - null_lnl)
    if stat < 0:
        stat = 0.0
    p = chi2.sf(stat, df)
    crit_005 = chi2.ppf(0.95, df)  
    crit_001 = chi2.ppf(0.99, df)  

    return {"stat": stat, "df": df, "p": p,
        "critical": {"0.05": crit_005, "0.01": crit_001},
        "sig": {"0.05": stat > crit_005, "0.01": stat > crit_001}}



def run_pair_model(seq_path,tree_path,out_dir, model, only_branch=True, gene_name=None):
    # output_dir = 传入base，手动创建res

    codeml_bin = shutil.which('codeml')
    if not codeml_bin:
        raise FileNotFoundError(
            "PAML not found. Please install PAML (codeml) and ensure it is on your PATH, "
            "or set the CODEML environment variable to the absolute path of the executable.\n"
            "Examples:\n"
            "  conda install -c bioconda paml\n"
            "  export CODEML=/path/to/codeml"
        )
    seq_path  = Path(seq_path)
    tree_path = Path(tree_path)
    if gene_name is None:
        gene_name = seq_path.stem


    base = Path(out_dir)
    res = base / "result"
    base.mkdir(parents=True, exist_ok=True)
    res.mkdir(parents=True, exist_ok=True)
    if model.upper() == "M0M3":
        null_ctl = base / "M0.ctl"
        alt_ctl = base / "M3.ctl"
        null_mlc = res / "M0.mlc"
        alt_mlc = res / "M3.mlc"
        # 写 ctl
        write_ctl("M0", seq_path, tree_path, null_mlc, null_ctl)
        write_ctl("M3", seq_path, tree_path, alt_mlc, alt_ctl)
        # run codeml
        procs =[]
        if not _mlc_complete(null_mlc):
            # work_m0 = base / "__work_M0"
            # work_m0.mkdir(parents=True, exist_ok=True)
            print(f"找不到{null_mlc}")
            procs.append(_run_codeml_async(codeml_bin, null_ctl, base / "__work_M0"))
            # subprocess.run([codeml_bin, str(null_ctl)], cwd=work_m0, check=True)
        if not _mlc_complete(alt_mlc):
            # work_m3 = base / "__work_M3"
            # work_m3.mkdir(parents=True, exist_ok=True)
            procs.append(_run_codeml_async(codeml_bin, alt_ctl, base / "__work_M3"))
            # subprocess.run([codeml_bin, str(alt_ctl)], cwd=work_m3, check=True)
        for p in procs:
            ret = p.wait()
            if ret != 0:
                raise RuntimeError(f"codeml failed with exit code {ret}")
        # parse lnL/np
        null_lnl, null_np = parse_mlc_np_lnl(null_mlc)
        alt_lnl, alt_np = parse_mlc_np_lnl(alt_mlc)
        # scipy test
        return lrt(null_lnl, null_np, alt_lnl, alt_np)

    elif model.upper() == "M7M8":
        null_ctl = base / "M7.ctl"
        alt_ctl = base / "M8.ctl"
        null_mlc = res / "M7.mlc"
        alt_mlc = res / "M8.mlc"
        write_ctl("M7", seq_path, tree_path, null_mlc, null_ctl)
        write_ctl("M8", seq_path, tree_path, alt_mlc, alt_ctl)

        # run codeml
        procs = []
        if not _mlc_complete(null_mlc):
            # work_m7 = base / "__work_M7"
            # work_m7.mkdir(parents=True, exist_ok=True)
            # subprocess.run([codeml_bin, str(null_ctl)], cwd=work_m7, check=True)
            procs.append(_run_codeml_async(codeml_bin, null_ctl, base / "__work_M7"))
        if not _mlc_complete(alt_mlc):
            # work_m8 = base / "__work_M8"
            # work_m8.mkdir(parents=True, exist_ok=True)
            # subprocess.run([codeml_bin, str(alt_ctl)], cwd=work_m8, check=True)
            procs.append(_run_codeml_async(codeml_bin, alt_ctl,  base / "__work_M8"))
        for p in procs:
            ret = p.wait()
            if ret != 0:
                raise RuntimeError(f"codeml failed with exit code {ret}")
        # parse lnL/np
        null_lnl, null_np = parse_mlc_np_lnl(null_mlc)
        alt_lnl, alt_np = parse_mlc_np_lnl(alt_mlc)
        # scipy test
        return lrt(null_lnl, null_np, alt_lnl, alt_np)

    elif model.upper() == "FREERATIO":
        alt_ctl = base / "FREERATIO.ctl"
        alt_mlc = res / "FREERATIO.mlc"
        null_ctl = base / "M0.ctl"
        null_mlc = res / "M0.mlc"

        write_ctl("M0", seq_path, tree_path, null_mlc, null_ctl)
        write_ctl("FREERATIO", seq_path, tree_path, alt_mlc, alt_ctl)
        procs = []
        if not _mlc_complete(null_mlc):
            procs.append(_run_codeml_async(codeml_bin, null_ctl, base / "__work_M0"))
        if not _mlc_complete(alt_mlc):
            procs.append(_run_codeml_async(codeml_bin, alt_ctl,  base / "__work_FreeRatio"))
        for p in procs:
            ret = p.wait()
            if ret != 0:
                raise RuntimeError(f"codeml failed with exit code {ret}")
        # parse lnL/np
        null_lnl, null_np = parse_mlc_np_lnl(null_mlc)
        alt_lnl, alt_np = parse_mlc_np_lnl(alt_mlc)
        # scipy test
        return lrt(null_lnl, null_np, alt_lnl, alt_np),str(null_mlc), str(alt_mlc)

    return None


