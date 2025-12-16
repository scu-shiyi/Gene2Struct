
import argparse
from scripts.docking.core import run
import os
        
def main():
    parser = argparse.ArgumentParser(
        prog="scripts docking",
        usage = "%(prog)s [-m mapping_csv] [-t tree_path] [-o <output_dir>] ...",
        description= ("DockingModule: A molecular docking pipeline that calculates binding affinities, "
            "assesses substrate/product inhibition potential, and visualizes "
            "results on a phylogenetic tree."),
        formatter_class=argparse.MetavarTypeHelpFormatter)


    parser.add_argument("-m", required=True, type= str,
                        metavar="FILE",
                        help="configuration file for scripts information.")

    parser.add_argument("-t", required=True,type = str,
                        metavar="FILE",
                        help="Phylogenetic tree file in Newick format (used for evolutionary visualization).")
    
    parser.add_argument("-o", required=True, type=str,
                        metavar="DIR",
                        help="Directory to save results")

    args = parser.parse_args()

    run(args)

if __name__ == "__main__":
    main()
