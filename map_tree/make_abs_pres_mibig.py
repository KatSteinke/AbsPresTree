"""This module generates an absence/presence matrix of cluster families from BiG-SCAPE results.
It is intended as part of make_map_tree for creating a tree using absence/presence data.
However, it can be used on its own from the commandline.

Example:
    python3 make_abs_pres.py result_dir/

NOTE: module currently designed to work with autoMLST trees, so organism names are sanitized.
"""

import pathlib
import re
import sys

from typing import Dict, Set, Tuple


def get_families(cluster_file: pathlib.Path, families: Dict[str, Set]) -> None:
    """Read a Bigscape clustering file and map families to clusters, and updates
    an existing dictionary with the families extracted.
    Arguments:
        cluster_file: pathlib.Path to the clustering file
        families: dictionary to update with the clusters found
    """
    with cluster_file.open(mode="r") as infile:
        line = infile.readline()
        # check presence of the comment line
        if not line.startswith("#"):
            raise ValueError("Header line of file {} missing.".format(cluster_file.name))
        # move to the next line
        line = infile.readline()
        while line:
            line = line.rstrip().split("\t")
            if len(line) != 2:
                raise ValueError("Cluster file must contain two columns. File {} contains {}.".format(
                                                                                                      cluster_file.name,
                                                                                                      len(line)))
            cluster = line[0]
            family = line[1]
            if family not in families:
                families[family] = {cluster}
            else:
                families[family].add(cluster)
            line = infile.readline()


def map_clusters_to_names(summary_file: pathlib.Path) -> Tuple[Dict[str, Set[str]], Dict[str, str]]:
    """Read the Bigscape summary file and map clusters to genomes.
    Arguments:
        summary_file: pathlib.Path to the summary file
    Returns:
        A dict of genomes with the clusters they contain;
        another dict mapping BGC IDs to cluster names.
    """
    genomes_with_clusters = {}
    bgcs_with_names = {}
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

            line = infile.readline()
    if wrong_line_format:
        format_warn = "Line(s) in the summary file contain a nonstandard number of columns (standard is 7). \n" \
                      "In case of errors/strange results, please check the format of your summary file.\n"
        sys.stderr.write(format_warn)
    return genomes_with_clusters, bgcs_with_names


def bigscape_to_matrix(network_result_dir: str, cutoff: float = 0.3) -> str:
    """Read a set of Bigscape results and create an absence/presence matrix
    for the cluster families in each organism.
    Arguments:
        network_result_dir: directory containing the Bigscape networking results
        cutoff: the cutoff of the Bigscape run of interest, if the folder 
                contains multiple ones. Default: 0.3 (Bigscape default)
    Returns:
        An absence/presence matrix for the cluster families, as a tab-separated
        string for processing in ETE3.
    """
    base_path = pathlib.Path(network_result_dir).resolve()

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

    # parsing individual results
    families_with_clusters = {}
    cluster_name_pattern = "*_clustering_c{:.2f}.tsv".format(cutoff) # TODO: something happens here
    for result_dir in result_dirs:
        # catch duplicates and missing files - TODO: make prettier?
        clusters_found = False
        for entry in result_dir.iterdir():
            if entry.is_file() and entry.match(cluster_name_pattern):
                if not clusters_found:
                    get_families(entry, families_with_clusters)
                    clusters_found = True
                # if there's already a cluster file found, complain
                else:
                    raise ValueError("Multiple clustering files in directory {}.".format(result_dir.name))
        # if there's no clustering file found in the entire dir, complain - but this can happen? TODO: what even is happening
        #if not clusters_found:
        #    raise ValueError("No clustering files in directory {}".format(result_dir.name))


    # parsing networking file
    genomes_to_clusters, bgcs_to_names = map_clusters_to_names(all_networks)

    # create name-to-family mapping - TODO: efficiency tweaks?
    named_families = {}
    for family, cluster_list in families_with_clusters.items():
        family_name = family
        # only add if it's not just a MiBIG family
        if not all([found_cluster.startswith("BGC") for found_cluster in cluster_list]):
            for found_cluster in cluster_list:
            # don't add the same name over and over
                if found_cluster.startswith("BGC") \
                and not bgcs_to_names[found_cluster] in family_name:
                    family_name += "_" + bgcs_to_names[found_cluster]
            named_families[family_name] = cluster_list



    # get a list of the keys for fixed order - TODO: find a way to have it properly sorted?
    all_families = sorted(named_families.keys())

    # construct matrix string
    absence_presence_matrix = "#NAMES\t" + "\t".join(all_families)
    for genome, clusters in sorted(genomes_to_clusters.items()):
        absence_presence_matrix += "\n" + genome + "\t"  # avoids trailing newlines
        absence_presence_row = [str(len(clusters.intersection(named_families[family])))
                                for family in all_families]
        absence_presence_matrix += "\t".join(absence_presence_row)

    return absence_presence_matrix


if __name__ == "__main__":
    if len(sys.argv) == 2:
        CURRENT_PATH = sys.argv[1]
        CUTOFF = 0.3
    elif len(sys.argv) == 3:
        CURRENT_PATH = sys.argv[1]
        CUTOFF = sys.argv[2]
    else:
        print("Please supply Bigscape result directory. Usage: make_abs_pres_mibig.py result_dir [cutoff]")
        sys.exit(1)
    try:
        print(bigscape_to_matrix(CURRENT_PATH, CUTOFF))
    except ValueError as err:  # TODO - raise more specifically!
        sys.stderr.write(err)
