import argparse
from scripts.seq_mining.unix_command import prepare_workdir, execute_tasks, COMMAND_HELP

def run(argv=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='seq_mining is a tool for extracting phylogenetic marker genes.')
    parser.add_argument('command',
                        choices=('filter', 'assemble', 'consensus', 'trim', 'combine', 'tree'),
                        help='One or several of the following actions, separated by space:' + COMMAND_HELP,
                        metavar='command',
                        nargs='*')

    parser.add_argument('-f', help='Sample list file', metavar='FILE', required=True)
    parser.add_argument('-r', help='Reference directory', metavar='DIR', required=True)
    parser.add_argument('-o', help='Output directory', metavar='DIR', required=True)
    parser.add_argument('-p', default=1, help='Number of parallel processes', metavar='INT', type=int)

    parser.add_argument('-kf', default=31, help='Filter k-mer size', metavar='INT', type=int)
    parser.add_argument('-ka', default=0, help='Assembly k-mer size (default = auto)', metavar='INT', type=int)
    parser.add_argument('-s', '--step-size', default=4, help='Filter step size', metavar='INT', type=int)
    parser.add_argument('-e', '--error-threshold', default=2, help='Error threshold', metavar='INT', type=int)
    parser.add_argument('-sb', '--soft-boundary', choices=('0', 'auto', 'unlimited'), default='auto', help='Soft boundary (default = auto)', type=str)
    parser.add_argument('-i', '--iteration', default=4096, help='Search depth', metavar='INT', type=int)

    parser.add_argument('-c', '--consensus-threshold', default='0.75', help='Consensus threshold (default = 0.75)', metavar='FLOAT', type=float)

    parser.add_argument('-ts', '--trim-source', choices=('assembly', 'consensus'), default=None, help='Whether to trim the primary assembly or the consensus sequence (default = output of last step, assembly if no other command given)')
    parser.add_argument('-tm', '--trim-mode', choices=('all', 'longest', 'terminal', 'isoform'), default='terminal', help='Trim mode (default = terminal)', type=str)
    parser.add_argument('-tr', '--trim-retention', default=0, help='Retention length threshold (default = 0.0)', metavar='FLOAT', type=float)

    parser.add_argument('-cs', '--combine-source', choices=('assembly', 'consensus', 'trimmed'), default=None, help='Whether to combine the primary assembly, the consensus sequences or the trimmed sequences (default = output of last step, assembly if no other command given)')
    parser.add_argument('-cd', '--clean-difference', default=1, help='Maximum acceptable pairwise difference in an alignment (default = 1.0)', metavar='FLOAT', type=float)
    parser.add_argument('-cn', '--clean-sequences', default=0, help='Number of sequences required in an alignment (default = 0)', metavar='INT', type=int)

    parser.add_argument('-m', '--tree-method', choices=('coalescent', 'concatenation'), default='coalescent', help='Multi-scripts tree reconstruction method (default = coalescent)')
    parser.add_argument('-b', '--bootstrap', default=1000, help='Number of bootstrap replicates', metavar='INT', type=int)

    parser.add_argument('--max-reads', default=0, help='Maximum reads per file', metavar='INT', type=int)
    parser.add_argument('--min-depth', default=50, help='Minimum acceptable depth during re-filtering', metavar='INT', type=int)
    parser.add_argument('--max-depth', default=768, help='Maximum acceptable depth during re-filtering', metavar='INT', type=int)
    parser.add_argument('--max-size', default=6, help='Maximum file size during re-filtering', metavar='INT', type=int)
    parser.add_argument('--min-ka', default=21, help='Minimum auto-estimated assembly k-mer size', metavar='INT', type=int)
    parser.add_argument('--max-ka', default=51, help='Maximum auto-estimated assembly k-mer size', metavar='INT', type=int)
    parser.add_argument('--msa-program', choices=('clustalo', 'mafft', 'muscle'), default='mafft', help='Program for multiple sequence alignment', type=str)
    parser.add_argument('--no-alignment', action='store_true', default=True, help='Do not perform multiple sequence alignment')
    parser.add_argument('--no-trimal', action='store_true', default=False, help='Do not run trimAl on alignments')
    parser.add_argument('--phylo-program', choices=('raxmlng', 'iqtree', 'fasttree', 'veryfasttree'), default='fasttree', help='Program for phylogenetic tree reconstruction', type=str)


    args = parser.parse_args(argv)

    # 保持上面的默认和 --concat 逻辑
    args.command = args.command or ('filter','assemble','trim','combine')
    samples = prepare_workdir(args)
    if samples:
        print(f'Running tasks: {", ".join(args.command)}')
        print()
        execute_tasks(args, samples)
    else:
        print('Sample list is empty or invalid, exiting')


if __name__ == '__main__':
    run()