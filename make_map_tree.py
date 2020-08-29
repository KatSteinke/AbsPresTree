from argparse import ArgumentParser
import sys

from map_tree import draw_cluster_tree, make_abs_pres, make_abs_pres_mibig

# TODO: find good phrasing for how it works or change it
parser = ArgumentParser(description="""Generate a tree with absence/presence matrices from BiG-SCAPE data and a tree.""")
parser.add_argument("cluster_dir", action="store", help="Directory with BiG-SCAPE network files")
parser.add_argument("tree_file", action="store", help="Tree to add absence/presence data to")
parser.add_argument("out_tree", action="store", help="Output file (image)")
parser.add_argument("--mibig", action="store_true", help="Process data that includes MIBiG clusters")
parser.add_argument("-c", "--cutoff", action="store", default=0.3, help="Cutoff of BiG-SCAPE run of interest (default 0.3)")
args = parser.parse_args()

cluster_dir = args.cluster_dir
tree_file = args.tree_file
out_tree = args.out_tree
mibig = args.mibig
cutoff = float(args.cutoff)

try:
    if mibig:
        absence_presence_matrix = make_abs_pres_mibig.bigscape_to_matrix(cluster_dir, cutoff)
    else:
        absence_presence_matrix = make_abs_pres.bigscape_to_matrix(cluster_dir, cutoff)
except ValueError as value_err:
    sys.stderr.write(str(value_err))

draw_cluster_tree.make_cluster_tree(tree_file, absence_presence_matrix, out_tree)
