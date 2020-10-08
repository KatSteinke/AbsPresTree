#!/bin/bash

set -e

# virtualenv juggling
antismash_env=/home/kat/.virtualenvs/antismash-dmz_markers/bin
automlst_env=/home/kat/.virtualenvs/automlst_env/bin
ete_env=/home/kat/.virtualenvs/ete3_test/bin
general_env=/home/kat/.virtualenvs/bioeng_env/bin
networkx_env=/home/kat/.virtualenvs/networkx-env/bin/

# script juggling
general_base=~/Documents/work/DTU/bioengineering10
automlst_base=~/Documents/autoMLST/ziemertlab-automlst-7b2b5a9a8961
antismash_base=~/Documents/antismash-dmz_markers
abs_pres=~/Documents/PycharmProjects/bacillus_job


# arguments: base_dir, accessions, out_tree, outgroups

base_dir=~/Documents/work/DTU/bioengineering10/test_cleaned_script
accessions="curated-tiny"
out_tree="matrix_tree.png"
outgroups=""
while getopts b:a:t:g: option; do
	case "$option" in
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
# download all the files - now hopefully without stopping for broken ones
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
            $general_env/python3 $general_base/rename_strainless_organisms.py "../$gbk_file" --overwrite
            echo "Running antiSMASH on $gbk_file"
            $antismash_env/python3 $antismash_base/run_antismash.py --data /opt/antismash/data --cpus 4 --minimal "../$gbk_file" 
        fi 
    fi
done
# run Bigscape
echo "Running Bigscape"
# generate unique prefix for run
bigscape_id=$(uuidgen)
$general_env/python3 $general_base/BiG-SCAPE-master/bigscape.py -i . -o ../bigscape-results --pfam_dir $general_base/Pfam-A --mibig --label $bigscape_id
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
