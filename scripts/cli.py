import argparse
import os
import tarfile
from scripts.evolution.core import RunSite, RunBranch
from scripts.seq_mining.core import run as run_geneminer
from scripts.docking.core import run as run_docking
import shutil
import subprocess
import shlex
import zlib
from urllib.parse import urljoin

# ===== MGLTools Constants and Configuration =====
MGL_URL_BASE = os.environ.get("MGL_URL_BASE", "https://ccsb.scripps.edu/download/532/")
if not MGL_URL_BASE.endswith('/'):
    MGL_URL_BASE += '/'
MGL_ARCHIVE = os.environ.get("MGL_ARCHIVE_NAME", "mgltools_x86_64Linux2.1.5.7.tar.gz")
MGL_DIR_NAME = "mgltools_x86_64Linux2_1.5.7"
MGL_SENTINEL_FILE = ".install_complete"

# ===== Helper Functions (as refined in our discussions) =====

def _is_gzip_magic(path: str) -> bool:
    """Checks if a file has a valid gzip header."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False

def _archive_is_complete(archive_path: str) -> bool:
    """Checks if the archive is a complete and valid tar.gz file."""
    if not os.path.isfile(archive_path) or not _is_gzip_magic(archive_path):
        return False
    try:
        with tarfile.open(archive_path, "r:gz") as tar:
            names = set(tar.getnames())
            req_dir = f"{MGL_DIR_NAME}/"
            # Now, we only check for the MGLToolsPckgs.tar.gz sub-archive.
            req_file = f"{MGL_DIR_NAME}/MGLToolsPckgs.tar.gz"
            return any(n.startswith(req_dir) for n in names) and req_file in names
    except (tarfile.ReadError, OSError, EOFError, zlib.error):
        return False

def _safe_extract(tar: tarfile.TarFile, path: str):
    """Safely extracts a tar archive, preventing path traversal."""
    base = os.path.realpath(path)
    for member in tar.getmembers():
        target = os.path.realpath(os.path.join(path, member.name))
        if not target.startswith(base + os.sep):
            raise RuntimeError(f"Unsafe path in tar: {member.name}")
    tar.extractall(path=path)

def _extract_archive(archive_path: str, utils_dir: str):
    print("Extracting MGLTools …")
    with tarfile.open(archive_path, "r:gz") as tar:
        _safe_extract(tar, utils_dir)
    
    print("Extracting MGLTools sub-archives...")
    # Get the extracted MGLTools root directory
    extracted_root = os.path.join(utils_dir, MGL_DIR_NAME)
    
    # Extract the MGLToolsPckgs.tar.gz sub-archive
    pkg_tarball_path = os.path.join(extracted_root, "MGLToolsPckgs.tar.gz")
    if os.path.isfile(pkg_tarball_path):
        with tarfile.open(pkg_tarball_path) as tar:
            _safe_extract(tar, extracted_root)
        os.remove(pkg_tarball_path) # Clean up the sub-archive
    else:
        raise FileNotFoundError(f"Sub-archive {pkg_tarball_path} not found.")
def _patch_prepare_ligand(target_dir: str):
    """Applies a patch to the prepare_ligand4.py script."""
    patch_path = os.path.join(target_dir, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_ligand4.py")
    if not os.path.isfile(patch_path):
        # 如果文件不存在，引发异常让调用者处理
        raise FileNotFoundError(f"Patch target not found: {patch_path}")
    print(f"Patching {patch_path} …")
    with open(patch_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    patched_lines = []
    for line in lines:
        s = line.strip()
        if s == "ligand_filename = os.path.basename(a)":
            patched_lines.append("# " + line if not line.lstrip().startswith("#") else line)
        elif s == "#ligand_filename = a":
            indentation = len(line) - len(line.lstrip())
            patched_lines.append(f"{' ' * indentation}ligand_filename = a\n")
        else:
            patched_lines.append(line)
    with open(patch_path, "w", encoding="utf-8") as f:
        f.writelines(patched_lines)

def _download_fast(url: str, dest: str) -> bool:
    """Uses an external downloader to fetch the file."""
    prefer = os.environ.get("MGL_DOWNLOADER", "").lower()
    order = [prefer] if prefer in {"aria2c", "wget", "curl"} else ["wget", "aria2c", "curl"]
    extra = shlex.split(os.environ.get("MGL_DOWNLOADER_ARGS", ""))
    ipv4 = os.environ.get("MGL_IPV4") == "1"

    for tool in order:
        if tool == "wget" and shutil.which("wget"):
            cmd = ["wget", "-c", "--tries=10", "--timeout=30", "-O", dest]
            if ipv4: cmd.insert(1, "-4")
            cmd += extra + [url]
        elif tool == "aria2c" and shutil.which("aria2c"):
            cmd = ["aria2c", "-x16", "-s16", "-k", "1M", "--continue=true", "-o", os.path.basename(dest), "--dir", os.path.dirname(dest), url]
            cmd += extra
        elif tool == "curl" and shutil.which("curl"):
            cmd = ["curl", "-L", "--retry", "10", "--retry-delay", "2", "-C", "-", "-o", dest, url]
            if ipv4: cmd.insert(1, "-4")
            cmd += extra
        else:
            continue
        try:
            print(f"Trying external downloader: {tool} …")
            subprocess.run(cmd, check=True)
            return True
        except Exception as e:
            print(f"[WARN] {tool} failed: {e}; fallback…")
            continue
    return False

# ===== Main MGLTools Installation Logic =====

def ensure_mgltools():
    """
    Checks for MGLTools directory, then archive, and downloads only if necessary.
    """
    base_dir = os.path.abspath(os.path.dirname(__file__))
    utils_dir = os.path.join(base_dir, "utils")
    target_dir = os.path.join(utils_dir, MGL_DIR_NAME)
    archive_path = os.path.join(utils_dir, MGL_ARCHIVE)
    sentinel_path = os.path.join(target_dir, MGL_SENTINEL_FILE)

    # Make sure the base directory for utils exists
    os.makedirs(utils_dir, exist_ok=True)
    # 0 
    if os.path.isfile(sentinel_path):
        return target_dir
    print("MGLTools installation check initiated.")

    is_install = False
    # 1. Check if the directory is already installed and patched
    patch_target_path = os.path.join(target_dir, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_ligand4.py")
    if os.path.isdir(target_dir) and os.path.isfile(patch_target_path):
        print("MGLTools already present, skipping download and extract.")
        try:
            _patch_prepare_ligand(target_dir)
            is_install = True
        except FileNotFoundError:
            print(f"[WARN] Incomplete MGLTools installation found during patching. Reinstalling.")
            shutil.rmtree(target_dir, ignore_errors=True)
        except Exception as e:
            print(f"[ERROR] An unexpected error occurred during patching: {e}. Reinstalling.")
            shutil.rmtree(target_dir, ignore_errors=True)
    if is_install:
        with open(sentinel_path, "w") as f:
            f.write("Installation complete.\n")
        return target_dir

    # 2. Check if a local, complete archive exists and extract it
    if os.path.isfile(archive_path) and _archive_is_complete(archive_path):
        print("Found local, complete archive. Extracting...")
        try:
            _extract_archive(archive_path, utils_dir)
            print("MGLTools installed from local archive.")
            _patch_prepare_ligand(target_dir)
            with open(sentinel_path, "w") as f:
                f.write("Installation complete.\n")
            return target_dir
        except (tarfile.ReadError, RuntimeError, FileNotFoundError) as e:
            print(f"[ERROR] Failed to extract local archive: {e}. Removing...")
            os.remove(archive_path)

    # 3. If no directory and no complete archive, download the file
    print("MGLTools not found. Downloading...")
    # Prioritize the direct file link for downloaders like wget and curl
    url = urljoin(MGL_URL_BASE, MGL_ARCHIVE)
    if not _download_fast(url, archive_path):
        raise RuntimeError(
            "Download failed. Please check your network connection or try again later."
        )

    # 4. After download, perform final checks and installation
    if not _archive_is_complete(archive_path):
        raise RuntimeError("Downloaded archive appears invalid or incomplete.")

    _safe_extract(tarfile.open(archive_path, "r:gz"), utils_dir)
    print("MGLTools installed under utils.")
    _patch_prepare_ligand(target_dir)
    with open(sentinel_path, "w") as f:
        f.write("Installation complete.\n")
    return target_dir



def main():
    ensure_mgltools()
    parser = argparse.ArgumentParser(
        prog="phyloselect",
        description="PhyloSelect: A comprehensive toolkit for phylogenetic analysis and evolutionary selection.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # GeneMiner
    COMMAND_HELP = '''
    filter    Reference-based filtering of raw reads
    refilter  Refinement of filtered reads
    assemble  Gene assembly using wDBG
    consensus Consensus generation on heterozygous sites
    trim      Flank sequence removal
    combine   Gene alignment, concatenation and cleanup
    tree      Phylogenetic tree reconstruction
    '''
    parser_geneminer = subparsers.add_parser("geneminer", help="Extract phylogenetic marker genes(CDS) from genomic data for evolutionary studies.",
                                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_geneminer.add_argument('command',
                        choices=('filter', 'assemble', 'consensus', 'trim', 'combine', 'tree',[]),
                        help='One or several of the following actions, separated by space:' + COMMAND_HELP,
                        metavar='command',
                        nargs='*')

    parser_geneminer.add_argument('-f', help='Sample list file', metavar='FILE', required=True)
    parser_geneminer.add_argument('-r', help='Reference directory', metavar='DIR', required=True)
    parser_geneminer.add_argument('-o', help='Output directory', metavar='DIR', required=True)
    parser_geneminer.add_argument('-p', default=1, help='Number of parallel processes', metavar='INT', type=int)

    parser_geneminer.add_argument('-kf', default=31, help='Filter k-mer size', metavar='INT', type=int)
    parser_geneminer.add_argument('-ka', default=0, help='Assembly k-mer size (default = auto)', metavar='INT', type=int)
    parser_geneminer.add_argument('-s', '--step-size', default=4, help='Filter step size', metavar='INT', type=int)
    parser_geneminer.add_argument('-e', '--error-threshold', default=2, help='Error threshold', metavar='INT', type=int)
    parser_geneminer.add_argument('-sb', '--soft-boundary', choices=('0', 'auto', 'unlimited'), default='auto', help='Soft boundary (default = auto)', type=str)
    parser_geneminer.add_argument('-i', '--iteration', default=4096, help='Search depth', metavar='INT', type=int)

    parser_geneminer.add_argument('-c', '--consensus-threshold', default='0.75', help='Consensus threshold (default = 0.75)', metavar='FLOAT', type=float)

    parser_geneminer.add_argument('-ts', '--trim-source', choices=('assembly', 'consensus'), default=None, help='Whether to trim the primary assembly or the consensus sequence (default = output of last step, assembly if no other command given)')
    parser_geneminer.add_argument('-tm', '--trim-mode', choices=('all', 'longest', 'terminal', 'isoform'), default='terminal', help='Trim mode (default = terminal)', type=str)
    parser_geneminer.add_argument('-tr', '--trim-retention', default=0, help='Retention length threshold (default = 0.0)', metavar='FLOAT', type=float)

    parser_geneminer.add_argument('-cs', '--combine-source', choices=('assembly', 'consensus', 'trimmed'), default=None, help='Whether to combine the primary assembly, the consensus sequences or the trimmed sequences (default = output of last step, assembly if no other command given)')
    parser_geneminer.add_argument('-cd', '--clean-difference', default=1, help='Maximum acceptable pairwise difference in an alignment (default = 1.0)', metavar='FLOAT', type=float)
    parser_geneminer.add_argument('-cn', '--clean-sequences', default=0, help='Number of sequences required in an alignment (default = 0)', metavar='INT', type=int)

    parser_geneminer.add_argument('-m', '--tree-method', choices=('coalescent', 'concatenation'), default='coalescent', help='Multi-scripts tree reconstruction method (default = coalescent)')
    parser_geneminer.add_argument('-b', '--bootstrap', default=1000, help='Number of bootstrap replicates', metavar='INT', type=int)

    parser_geneminer.add_argument('--max-reads', default=0, help='Maximum reads per file', metavar='INT', type=int)
    parser_geneminer.add_argument('--min-depth', default=50, help='Minimum acceptable depth during re-filtering', metavar='INT', type=int)
    parser_geneminer.add_argument('--max-depth', default=768, help='Maximum acceptable depth during re-filtering', metavar='INT', type=int)
    parser_geneminer.add_argument('--max-size', default=6, help='Maximum file size during re-filtering', metavar='INT', type=int)
    parser_geneminer.add_argument('--min-ka', default=21, help='Minimum auto-estimated assembly k-mer size', metavar='INT', type=int)
    parser_geneminer.add_argument('--max-ka', default=51, help='Maximum auto-estimated assembly k-mer size', metavar='INT', type=int)
    parser_geneminer.add_argument('--msa-program', choices=('clustalo', 'mafft', 'muscle'), default='mafft', help='Program for multiple sequence alignment', type=str)
    parser_geneminer.add_argument('--no-alignment', action='store_true', default=True, help='Do not perform multiple sequence alignment')
    parser_geneminer.add_argument('--no-trimal', action='store_true', default=False, help='Do not run trimAl on alignments')
    parser_geneminer.add_argument('--phylo-program', choices=('raxmlng', 'iqtree', 'fasttree', 'veryfasttree'), default='fasttree', help='Program for phylogenetic tree reconstruction', type=str)

    # siteview
    parser_site = subparsers.add_parser("siteview", 
                                        help="Perform phylogenetic analysis to study site-wise heterogeneity, positive selection, and entropy in evolutionary trees.",
                                        description = ("This command performs phylogenetic analysis to study the site-wise heterogeneity, "
                                                        "positive selection, and entropy within evolutionary trees."),
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_site.add_argument("-s", required=True, metavar="FILE",
                             help="Path to a single aligned FASTA file (e.g., -f gene1.aln.fasta)")
    parser_site.add_argument("-o", required=True, metavar="DIR",
                            help="Directory to save results")
    parser_site.add_argument("-t", default=None, metavar="FILE",
                             help="Optional: Newick-format phylogenetic tree file (.treefile). If not provided, a tree will be inferred automatically.")
    parser_site.add_argument("-n", "--name-limit", type=int, default=20, metavar="CHAR_LIMIT",
                             help="Maximum number of characters shown per leaf label on the tree.")
    parser_site.add_argument("-thr", type=float, default=1.4, metavar="THRESHOLD",
                             help="Threshold for distinguishing conserved vs. poorly conserved sites (default: 1.4).")
    parser_site.add_argument("-g", "--og", nargs="+", default=None, metavar="SPECIES",
                             help="Optional: One or more outgroup species names for tree rooting (e.g., -g Sp1 Sp2 Sp3).")
    parser_site.add_argument("--site", action="store_true",
                             help="Enable site-level conservation calculation (default: disabled).")
    parser_site.add_argument("--heatmap_path", default=None, metavar="ENTROPY_CSV",
                             help="Optional: Precomputed entropy matrix CSV (for debugging only).")

    # selection
    """可用线程数"""
    parser_dnds = subparsers.add_parser("selection", help="Analyze evolutionary selection on genetic sequences to identify positive selection and evolutionary trends.",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_dnds.add_argument("-i", "--fasta_input", required=True, metavar="FASTA_INPUT",
                             help="FASTA input: a directory of FASTA files or a text file listing fasta paths.")
    parser_dnds.add_argument("-o", "--output_dir", required=True,metavar="OUT_DIR",
                             help="Directory to save results ")
    parser_dnds.add_argument("-m", "--tree_map",metavar="TREE_MAP",default=None,
                             help = ("Text file mapping each FASTA to its phylogenetic tree. "
                                     "Each line should contain: <fasta_file> <tree_file>."))


    # Docking 
    parser_dock = subparsers.add_parser("docking", help="Simulate molecular interactions and evaluate ligand binding affinities in enzyme-substrate systems.",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_dock.add_argument("-in", required=True, metavar="FILE",
                             help="CSV file mapping scripts to substrate/product CIDs.")
    parser_dock.add_argument("-t", required=True, metavar="FILE",
                             help="Phylogenetic tree in Newick format.")
    parser_dock.add_argument("-o", required=True, metavar="DIR",
                             help="Output directory")



    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    if args.command == "siteview":
        RunSite(args)
    elif args.command == "selection":
        RunBranch(args)
    elif args.command == "docking":
        run_docking(args)
    elif args.command == "geneminer":
        run_geneminer(args)

if __name__ == "__main__":
    main()
