"""This module visualizes trees with cluster absence/presence information.
It is intended as part of make_map_tree for creating a tree using absence/presence data.
However, it can be used on its own from the commandline.

Example:
    python3 draw_cluster_tree.py treefile matrix outfile

Note that the matrix parameter can either be a matrix file or a tab-separated matrix
as a string.
"""
from argparse import ArgumentParser
from typing import Optional, List

from ete3 import AttrFace, ClusterTree, TreeStyle, NodeStyle, TextFace, CircleFace
from ete3.treeview.faces import add_face_to_node


def make_cluster_tree(tree_file: str, matrix: str, out_file: str, outgroup: Optional[List[str]] = None) -> None:
    """Draw a tree with cluster absence/presence information from an existing
    tree file and absence/presence matrix, and save it as an image under the
    supplied file name.

    Arguments:
        tree_file: the name of the file containing the tree to annotate
        matrix:    either a tab-separated absence/presence matrix, or the name
                   of a file containing such a matrix.
        out_file:  the name under which to save the resulting image
        outgroup:  the organism(s) to use as an outgroup, if any
    """
    # ClusterTree needs tab-separated, but that can't be exported cleanly
    matrix = matrix.replace(",", "\t")
    # tree with clustering analysis
    tree = ClusterTree(tree_file, text_array=matrix)

    # rerooting the tree
    if outgroup:
        ancestor = tree.get_common_ancestor(outgroup)
        tree.set_outgroup(ancestor)
        tree.ladderize(direction=1)

    # set drawing line width to 2
    my_node_style = NodeStyle()
    my_node_style["vt_line_width"] = 2
    my_node_style["hz_line_width"] = 2
    my_node_style["size"] = 5

    # layout function
    def sel_mylayout(node):
        node.set_style(my_node_style)

        if node.is_leaf():
            # add names in larger font + italics
            species_name = AttrFace("name", fsize=12, fstyle="italic")
            add_face_to_node(species_name, node, column=0, position="branch-right")
            # add absence/presence matrix
            for i, value in enumerate(getattr(node, "profile", [])):
                if value > 0:
                    color = "#FF0000"
                else:
                    color = "#EEEEEE"
                my_face = CircleFace(8, color, style="circle")
                my_face.margin_right = 3
                my_face.margin_bottom = 3
                add_face_to_node(my_face, node, position="aligned", column=i)

    # Use my layout to visualize the tree
    my_tree_style = TreeStyle()

    # Add header
    for j, name in enumerate(tree.arraytable.colNames):
        name_face = TextFace(name, fsize=11)
        name_face.rotation = -90
        name_face.hz_align = 1
        name_face.vt_align = 1
        name_face.margin_bottom = 10
        my_tree_style.aligned_header.add_face(name_face, column=j)

    my_tree_style.scale_length = 0.1
    # myTreeStyle.show_branch_support = True
    # don't auto-show leaf names, since we dealt with that above
    my_tree_style.show_leaf_name = False

    # set layout function for my_tree_style
    my_tree_style.layout_fn = sel_mylayout

    #tree.render(out_file, w=183, units="mm", dpi=600, tree_style=my_tree_style)
    tree.render(out_file, dpi=600, tree_style=my_tree_style)


if __name__ == "__main__":
    parser = ArgumentParser(description="Generate a tree with absence/presence matrix with optional outgroup.")
    parser.add_argument('in_tree', action="store", help="Treefile to visualize")
    parser.add_argument('matrix_file', action="store", help="Absence/presence matrix")
    parser.add_argument('out_tree', action="store", help="Filename of image to output.")
    parser.add_argument("--outgroup", action="store", nargs="*", help="Outgroup(s) for rerooting the tree, optional.")
    args = parser.parse_args()
    IN_TREE = args.in_tree
    MATRIX_FILE = args.matrix_file
    OUT_TREE = args.out_tree
    OUTGROUP = args.outgroup
    make_cluster_tree(IN_TREE, MATRIX_FILE, OUT_TREE, OUTGROUP)
