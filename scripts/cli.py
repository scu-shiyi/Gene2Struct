import argparse
import os
import tarfile
import shutil
import subprocess
import shlex
import zlib
from urllib.parse import urljoin


def main():
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
                        choices=('filter', 'assemble', 'consensus', 'trim', 'combine', 'tree'),
                        help='One or several of the following actions, separated by space:' + COMMAND_HELP,
                        metavar='command',
                        nargs='*')

    parser_geneminer.add_argument('-f', help='Sample list file', metavar='FILE', required=True)
    parser_geneminer.add_argument('-r', help='Reference directory', metavar='DIR', required=True)
    parser_geneminer.add_argument('-o', dest = "output_dir", help='Output directory', metavar='DIR', required=True)
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
    parser_site.add_argument("-o", dest = "output_dir", required=True, metavar="DIR",
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
    parser_dnds.add_argument("-o", dest = "output_dir", required=True,metavar="OUT_DIR",
                             help="Directory to save results ")
    parser_dnds.add_argument("-m", "--tree_map",metavar="TREE_MAP",default=None,
                             help = ("Text file mapping each FASTA to its phylogenetic tree. "
                                     "Each line should contain: <fasta_file> <tree_file>."))


    # Docking 
    parser_dock = subparsers.add_parser("docking", help="Simulate molecular interactions and evaluate ligand binding affinities in enzyme-substrate systems.",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_dock.add_argument("-in", dest = 'mapping_csv', required=True, metavar="FILE",
                             help="CSV file mapping scripts to substrate/product CIDs.")
    parser_dock.add_argument("-t", dest='tree_path', required=True, metavar="FILE",
                             help="Phylogenetic tree in Newick format.")
    parser_dock.add_argument("-o", dest = 'output_dir', required=True, metavar="DIR",
                             help="Output directory")



    args = parser.parse_args()

    if args.command == "siteview":
        os.makedirs(args.output_dir, exist_ok=True)
        from scripts.evolution.core import RunSite
        RunSite(args)
    elif args.command == "selection":
        os.makedirs(args.output_dir, exist_ok=True)
        from scripts.evolution.core import RunBranch
        RunBranch(args)
    elif args.command == "docking":
        os.makedirs(args.output_dir, exist_ok=True)
        from scripts.docking.core import run as run_docking
        run_docking(args)
    elif args.command == "geneminer":
        os.makedirs(args.output_dir, exist_ok=True)
        from scripts.seq_mining.core import run as run_geneminer
        run_geneminer(args)

if __name__ == "__main__":
    main()
