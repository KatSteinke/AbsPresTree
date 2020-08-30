"""This module generates an absence/presence matrix of cluster families from BiG-SCAPE results.

Example:
    python3 make_abs_pres_networkx.py result_dir/

NOTE: module currently designed to work with autoMLST trees, so organism names are sanitized.
"""
import pathlib
import re
import sys

from argparse import ArgumentParser
from collections import OrderedDict
from typing import Dict, List, Set, Tuple

import pandas as pd
import networkx as nx
import scipy.cluster


class CompoundFamily:
    """Records type of compound with a set of clusters."""
    def __init__(self, compound_class: str, clusters: Set[str]):
        self.compound_class = compound_class
        self.clusters = clusters

    def __len__(self):
        return len(self.clusters)


def get_connected_from_df(network_df: pd.DataFrame, cluster_types: Dict[str, str]) -> List[CompoundFamily]:
    """"From a network as a dataframe, extract connected components.
    Arguments:
        network_df: dataframe representing the network
        cluster_types: a dictionary mapping cluster IDs to Bigscape secondary
                        metabolite types
        Returns:
            The connected components in their network with their compounds
    """
    nx_network = nx.from_pandas_edgelist(network_df,
                                         source="Clustername 1",
                                         target="Clustername 2",
                                         edge_attr="Raw distance")
    trimmed_components = []
    for components in nx.connected_components(nx_network):
        if not all([component.startswith("BGC") for component in components]):
            compound_type = "_".join({cluster_types[component] for component in components})
            trimmed_components.append(CompoundFamily(compound_type, components))
    return trimmed_components


def map_clusters_to_names(summary_file: pathlib.Path) -> Tuple[Dict[str, Set[str]], Dict[str, str], Dict[str, str]]:
    """Read the Bigscape summary file and map clusters to genomes.
    Arguments:
        summary_file: pathlib.Path to the summary file
    Returns:
        A dict of genomes with the clusters they contain;
        another dict mapping BGC IDs to cluster names;
        a third mapping clusters to Bigscape cluster types.
    """
    genomes_with_clusters = {}
    bgcs_with_names = {}
    clusters_with_types = {}
    wrong_line_format = False

    with summary_file.open(mode="r") as infile:
        line = infile.readline()
        if not line.startswith("BGC"):
            raise ValueError("Header line of summary file missing.")
        # move to the next line
        line = infile.readline()
        while line:
            line = line.rstrip().split("\t")
            if len(line) != 7:  # can happen in functional files  but also be a sign of wrongness
                wrong_line_format = True

            cluster = line[0]
            # quick fix for MIBiG families
            if cluster.startswith("BGC"):
                cluster_name = line[2].replace(" biosynthetic gene cluster", "")
                # sanitize cluster name
                cluster_name = re.sub(r"[^A-Za-z0-9_-]", "", cluster_name.replace(" ", "_"))
                bgcs_with_names[cluster] = cluster_name
            else:
                genome = line[5]
                # sanitize to match autoMLST rules: remove anything not an alphanumeric character, - or _
                genome = re.sub(r"[^A-Za-z0-9_-]", "", genome.replace(" ", "_"))
                if genome not in genomes_with_clusters:
                    genomes_with_clusters[genome] = {cluster}
                else:
                    genomes_with_clusters[genome].add(cluster)
            metabolite_type = line[4]
            clusters_with_types[cluster] = metabolite_type

            line = infile.readline()
    if wrong_line_format:
        format_warn = "Line(s) in the summary file contain a nonstandard number of columns (standard is 7). \n" \
                      "In case of errors/strange results, please check the format of your summary file.\n"
        sys.stderr.write(format_warn)
    return genomes_with_clusters, bgcs_with_names, clusters_with_types


def merge_all_bigscape_networks(results: List[pathlib.Path], bigscape_cutoff: float = 0.3) -> pd.DataFrame:
    """Merge all .network files for a single Bigscape run into a DataFrame.
    Arguments:
        results: a list of the paths of all folders from which to extract
        network files.
        bigscape_cutoff: the cutoff of the Bigscape run to use. Default 0.3.
    Returns:
        A DataFrame with all raw cluster-to-cluster distances. Where a pairing
        appears in multiple folders, raw distance is the estimated average of
        all raw distance values.
    """
    network_values = pd.DataFrame()
    network_name_pattern = "*_c{:.2f}.network".format(bigscape_cutoff)
    for result_dir in results:
        # catch duplicates and missing files
        network_found = False
        for dir_entry in result_dir.iterdir():
            if dir_entry.is_file() and dir_entry.match(network_name_pattern):
                if not network_found:
                    network_values = network_values.append(pd.read_csv(dir_entry, sep="\t")[['Clustername 1',
                                                                                             'Clustername 2',
                                                                                             'Raw distance']])
                    network_found = True
                # if there's already a cluster file found, complain
                else:
                    raise ValueError("Multiple network files in directory {}.".format(result_dir.name))
    # get all combinations of cluster matches
    network_values['Match_name'] = network_values["Clustername 1"] + network_values["Clustername 2"]
    # extract average raw distance for matches appearing more than once (= in more than one file)
    network_matches = network_values.groupby('Match_name')['Raw distance'].mean()
    # replace original raw distances with new averages via inner join; get rid of duplicate lines
    networks_combined = network_values[["Match_name",
                                        "Clustername 1",
                                        "Clustername 2"]].join(network_matches.to_frame(name="Raw distance"),
                                                               on="Match_name",
                                                               how="inner").drop_duplicates()
    # match name is no longer needed in the final output
    networks_combined = networks_combined.drop("Match_name", axis=1)

    return networks_combined


def get_families(connected_families: List[CompoundFamily], bgc_names: Dict[str, str]) -> Dict[str, Set[str]]:
    """Name compound families according to their type and any Mibig clusters
    they contain.
    Arguments:
        connected_families: the compound families to parse
        bgc_names: a dict mapping Mibig accession numbers to known compounds.
    Returns:
        A dict mapping informative family names to the clusters comprising the
        family.
    """
    families_names = {}
    for index, components in enumerate(sorted(connected_families,
                                              key=len,
                                              reverse=True),
                                       start=1):
        family_name = str(index) + "_" + components.compound_class
        for component in components.clusters:
            # don't add the same name over and over
            if component.startswith("BGC") \
            and not bgc_names[component] in family_name:
                family_name += "_" + bgc_names[component]
        families_names[family_name] = components.clusters
    return families_names


def get_absence_presence_matrix(families_with_names: Dict[str, Set[str]], clusters_in_genomes: Dict[str, Set[str]]) \
        -> str:
    """Create an absence/presence matrix for a given set of families in a given
    set of genomes.
    Arguments:
        families_with_names: a dict sorting clusters into named families
        clusters_in_genomes: a dict giving the clusters contained in each genome
        to examine.
    Returns:
        An absence/presence matrix as a comma-separated string.
    """
    families_clusters = OrderedDict()
    # for each family, get how many members of the family are in each genome
    for family in sorted(families_with_names.keys()):
        families_clusters[family] = [len(clusters.intersection(families_with_names[family]))
                                     for genomes, clusters in sorted(clusters_in_genomes.items())]
    # use scipy's linkage, like seaborn does, to cluster the genes by their absence/presence in the species,
    # then rearrange columns in dataframe accordingly
    absence_presence_tmp = pd.DataFrame(data=families_clusters, index=sorted(clusters_in_genomes.keys()))
    absence_presence_by_clusters = absence_presence_tmp.T
    linkage_transposed = scipy.cluster.hierarchy.linkage(absence_presence_by_clusters, optimal_ordering=True)
    transposed_leaves = scipy.cluster.hierarchy.leaves_list(linkage_transposed)
    absence_presence_reordered = absence_presence_tmp[[absence_presence_by_clusters.index.values[leaf_ind]
                                                       for leaf_ind in transposed_leaves]]
    absence_presence_matrix = "#NAMES" + absence_presence_reordered.to_csv()
    return absence_presence_matrix


def parse_family_order(order_file: pathlib.Path) -> List[str]:
    """Extract the desired order of families from a file listing family numbers
    separated by commas on a single line.
    Arguments:
        order_file: the path to the file

    Returns:
        The order of families in a list.
    """
    with open(order_file, "r") as family_order:
        family_numbers = family_order.readline().rstrip().split(",")
    # sanity check: are there duplicates in this? - TODO: remove if it bloats code too much
    numbers_seen = set()
    duplicates = []
    for number in family_numbers:
        if number in numbers_seen:
            duplicates.append(number)
        numbers_seen.add(number)
    if duplicates:
        raise ValueError("Family order list contains duplicates: {}".format(", ".join(duplicates)))
    return family_numbers


def get_preordered_matrix(families_with_names: Dict[str, Set[str]], clusters_in_genomes: Dict[str, Set[str]],
                          family_order: List[str]) -> str:
    """Create an absence/presence matrix for a given set of families in a given
        set of genomes, sorted by a previously specified order.
        Arguments:
            families_with_names: a dict sorting clusters into named families
            clusters_in_genomes: a dict giving the clusters contained in each genome
            to examine.
            family_order:        the order in which the clusters are to be sorted
        Returns:
            An absence/presence matrix as a comma-separated string.
        """
    # get matrix as above, but sort by list rather than family names
    families_clusters = OrderedDict()
    families_by_order = []
    if not len(families_with_names) == len(family_order):
        raise ValueError("""Length of family order list must match amount of families.
        Amount of families: {} Length of family order list: {}""".format(len(families_with_names), len(family_order)))
    # for each number, in the order given in the file
    for family_number in family_order:
        # add all families from the families_with_names list that start with the number + "_"
        # this should only be one
        families_by_order.extend(filter(lambda x: x.startswith(str(family_number)+"_"), families_with_names.keys()))
    for family in families_by_order:
        families_clusters[family] = [len(clusters.intersection(families_with_names[family]))
                                     for genomes, clusters in sorted(clusters_in_genomes.items())]
    # TODO: make sure genomes still match lists
    absence_presence_df = pd.DataFrame(data=families_clusters, index=sorted(clusters_in_genomes.keys()))
    absence_presence_matrix = "#NAMES" + absence_presence_df.to_csv()
    return absence_presence_matrix



if __name__ == "__main__":
    parser = ArgumentParser(description="""Generate an absence/presence matrix and the families represented in it from
    BiG-SCAPE data.""")
    parser.add_argument("cluster_dir", action="store",
                        help="Directory with BiG-SCAPE network files")
    parser.add_argument("--out_families", action="store", default="families.txt",
                        help="Name of file to store connected component families (default families.txt)")
    parser.add_argument("-c", "--cutoff", action="store", default=0.3,
                        help="Cutoff of BiG-SCAPE run of interest (default 0.3)")
    parser.add_argument("--family_order", action="store",
                        help="Optional file listing desired order of families in matrix as comma-separated numbers.")
    args = parser.parse_args()

    cluster_dir = args.cluster_dir
    out_families = args.out_families
    cutoff = float(args.cutoff)
    family_order = args.family_order

    # parse dirs
    base_path = pathlib.Path(cluster_dir).resolve()

    # find network file and result dirs
    all_networks = None
    result_dirs = []
    for entry in base_path.iterdir():
        if entry.is_file() and entry.match("Network_Annotations_Full*"):
            # catch duplicates
            if not all_networks:
                all_networks = entry
            else:
                raise ValueError("Multiple summary files found.")
        if entry.is_dir():
            result_dirs.append(entry)

    # sanity checks
    if not all_networks:
        raise IOError("No summary file found.")
    if not result_dirs:
        raise IOError("No results found.")
    # extract names
    genomes_to_clusters, bgcs_to_names, clusters_to_types = map_clusters_to_names(all_networks)
    # extract networks
    networks_merged = merge_all_bigscape_networks(result_dirs)
    # get connected components
    connected_components = get_connected_from_df(networks_merged, clusters_to_types)
    # get named families
    named_families = get_families(connected_components, bgcs_to_names)
    with open(out_families, "w") as family_file:
        family_file.write("Cluster\tFamily")
        for named_family, family_clusters in sorted(named_families.items()):
            for fam_cluster in family_clusters:
                family_file.write("\n{}\t{}".format(fam_cluster, named_family))
    # generate matrix
    # if families are to be in specific order, find and apply
    if family_order:
        order_of_families = parse_family_order(pathlib.Path(family_order).resolve())
        abs_pres_matrix = get_preordered_matrix(named_families, genomes_to_clusters, order_of_families)
    # otherwise cluster as networkx does
    else:
        abs_pres_matrix = get_absence_presence_matrix(named_families, genomes_to_clusters)
    print(abs_pres_matrix)
