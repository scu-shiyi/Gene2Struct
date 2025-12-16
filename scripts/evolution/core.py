from scripts.evolution.runsite import Run_Site
from scripts.evolution.RunEvoDnDs import RunEvoDnDs
from pathlib import Path


def RunSite(args):
    fasta_path = str(Path(args.fasta_path).resolve())
    output_dir = str(Path(args.output_dir).resolve())
    tree_path = str(Path(args.tree_path).resolve()) if args.tree_path else None
    is_calculate_site = args.site
    name_limit = args.name_limit if args.name_limit is not None else 20
    thr = args.thr if args.thr is not None else 1.4
    outgroups = args.og if args.og else []
    heatmap_path = str(Path(args.heatmap_path).resolve()) if args.heatmap_path else None

    output_path = Run_Site(
        fasta_path=fasta_path,
        output_dir=output_dir,
        tree_path=tree_path,
        is_site=is_calculate_site,
        name_limit=name_limit,
        thr=thr,
        outgroups=outgroups,
        heatmap_path=heatmap_path,)

    return output_path

def RunBranch(args):
    fasta_input = str(Path(args.fasta_input).resolve())
    output_dir = str(Path(args.output_dir).resolve())
    fasta_tree_map = str(Path(args.tree_map).resolve()) if args.tree_map else None
    output_png = RunEvoDnDs(
        fasta_input=fasta_input,
        output_dir=output_dir,
        fasta_tree_map=fasta_tree_map
    )
    return output_png