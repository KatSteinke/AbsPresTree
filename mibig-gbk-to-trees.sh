#!/bin/bash

set -e

# virtualenvs - replace placeholders with own virtualenvs
antismash_env=/path/to/.virtualenvs/antismash-dmz_markers/bin
automlst_env=/path/to/.virtualenvs/automlst_env/bin
ete_env=/path/to/.virtualenvs/ete3_env/bin
general_env=/path/to/.virtualenvs/base_env/bin
networkx_env=/path/to/.virtualenvs/networkx_env/bin/

# scripts - replace placeholders with own script locations
bigscape_base=/path/to/BiG-SCAPE-master
pfam_base=/path/to/Pfam-A
automlst_base=/path/to/autoMLST/ziemertlab-automlst-7b2b5a9a8961
antismash_base=/path/to/antismash-dmz_markers
abs_pres=/path/to/absprestree


# arguments: base_dir, accessions, out_tree, outgroups

base_dir=/path/to/default/dir
accessions="curated-tiny"
out_tree="matrix_tree.png"
outgroups=""

usage="$(basename "$0") [-h] [-b -a -t -g] -- download genomes from a list of accessions, run antiSMASH and BiG-SCAPE to detect BGCs, and visualize absence/presence of these in an autoMLST-generated tree

where:
	-h show this help text
	-b base directory to generate files in
	-a file listing accessions to download (one accession per line)
	-t name of the tree file to output
	-g outgroup(s) for the tree (if more than one, separate by spaces and wrap the entire list in double quotes)

"

while getopts hb:a:t:g: option; do
	case "$option" in
		h) echo "$usage"
		   exit
		   ;;
		b) base_dir=$(realpath ${OPTARG});;
		a) accessions=$(realpath ${OPTARG});;
		t) out_tree=${OPTARG};;
		g) outgroups="--outgroup ${OPTARG}";;
	esac
done

if [[ ! -d $base_dir ]]; then
	echo "$base_dir is not a directory."
	exit 1
fi

if [[ ! -f $accessions ]]; then
	echo "$accessions is not a file or does not exist."
	exit 1
fi

cd $base_dir
# download all the files without stopping for broken ones
echo "Downloading files"
for accession in $(cat $accessions); do
    $general_env/ncbi-acc-download -e all $accession
done

# create results folder if it doesn't exist
if [ ! -d "antismash_results" ]; then
    mkdir "antismash_results"
fi
cd "antismash_results"
# run antismash
for gbk_file in $(ls ..); do
    if [ "$(file -b "../$gbk_file")" = "ASCII text" ]; then
        # check for regex match: does this look like a Genbank file?
        # must start with LOCUS and end with a DD-MMM-YYYY date
        if [[ "$(head -n 1 ../$gbk_file)" =~ ^LOCUS.*[0-9]{2}-[A-Z]{3}-[0-9]{4}$ ]]; then
	    echo "Fixing strain naming for $gbk_file"
            $general_env/python3 $abs_pres/utility/rename_strainless_organisms.py "../$gbk_file" --overwrite
            echo "Running antiSMASH on $gbk_file"
            $antismash_env/python3 $antismash_base/run_antismash.py --data /opt/antismash/data --cpus 4 --minimal "../$gbk_file" 
        fi 
    fi
done
# run Bigscape
echo "Running Bigscape"
# generate unique prefix for run
bigscape_id=$(uuidgen)
$general_env/python3 $bigscape_base/bigscape.py -i . -o ../bigscape-results --pfam_dir $pfam_base --mibig --label $bigscape_id
# run autoMLST
echo "Creating tree"
$automlst_env/python $automlst_base/simplified_wrapper.py ".."
# absence/presence matrix
echo "Creating absence/presence matrix"
cd ../bigscape-results/network_files/
# find results directory by unique identifier
cd $(find -name "*_$bigscape_id*" -type d)

network_matrix=$($networkx_env/python $abs_pres/map_tree/make_abs_pres_networkx.py .)
$ete_env/python $abs_pres/map_tree/draw_cluster_tree.py $base_dir/raxmlpart.txt.treefile "$network_matrix" $base_dir/$out_tree $outgroups
