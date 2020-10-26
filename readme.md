#Absence/Presence Tree script
-----------------------------

This script generates a phylogenetic tree, as well as a matrix showing absence/presence of biosynthetic
gene clusters alongside this tree, for a set of genomes. 
It takes a list of accession numbers of genomes, predicts biosynthetic gene clusters for them using
antiSMASH, groups these BGCs into families with BiG-SCAPE and NetworkX, and visualizes absence/presence
of these families for each organism alongside a phylogenetic tree with autoMLST.

##Prerequisites
The following programs need to be installed:
* [ncbi-acc-download](https://github.com/kblin/ncbi-acc-download/tree/master/ncbi_acc_download)
* [BiG-SCAPE](https://git.wageningenur.nl/medema-group/BiG-SCAPE/)
* the [cluster-splitting fork of antiSMASH](https://github.com/KatSteinke/dmz-antismash)
* the [simplified wrapper fork of autoMLST](https://github.com/KatSteinke/automlst-simplified-wrapper)

Additionally, several Python virtual environments are required. Requirements for three of these are
given in `requirements-ete-env.txt`, `requirements-general-env.txt`, and `requirements-networkx-env.txt`.
Beyond this, virtual environments for autoMLST and antiSMASH are also required; requirements for these 
can be found in their respective repositories.

##Setup
Before running mibig-gbk-to-trees.sh for the first time, several placeholders need to be replaced. These are:
###Virtual environments
* `/path/to/.virtualenvs/antismash-dmz_markers/bin`: replace with path to antiSMASH virtualenv
* `/path/to/.virtualenvs/automlst_env/bin`: replace with path to autoMLST virtualenv
* `/path/to/.virtualenvs/ete3_env/bin`: replace with path to virtualenv specified in `requirements-ete-env.txt`
* `/path/to/.virtualenvs/base_env/bin`: replace with path to virtualenv specified in `requirements-general-env.txt`
* `/path/to/.virtualenvs/networkx_env/bin/`: replace with path to virtualenv specified in `requirements-networkx-env.txt`
###Programs
* `/path/to/BiG-SCAPE-master`: replace with path to dir containing BiG-SCAPE
* `/path/to/Pfam-A`: replace with path to PFAM directory used for BiG-SCAPE (normally, but not always, in BiG-SCAPE dir)
* `/path/to/autoMLST/ziemertlab-automlst-7b2b5a9a8961`: replace with path to dir containing autoMLST fork
* `/path/to/antismash-dmz_markers`: replace with path to dir containing antiSMASH fork
* `/path/to/absprestree`: replace with path to dir containing this script
###Directories
* `/path/to/default/dir`: replace with path to the directory you wish to run the script in by default
If desired, you can also change `curated-accessions-tiny` to a list of accessions you wish to run by default.

##Running the script
To run the script with default settings, simply run `bash mibig-gbk-to-trees.sh`. All that is required is a plain
text file listing accession numbers of the genomes you wish to generate an absence/presence matrix and tree for, one per 
line, under the name set as a default in the default base directory.
Additional arguments are:
* `-b`: Base directory to run script in; all files and directories are normally generated here
* `-a`: plain text file containing accession numbers of your desired genomes, one per line
* `-t`: name (not path!) of the final tree image. Formats supported by ete3 are .png, .pdf and .svg.
* `-g`: name (not accession!) of outgroup(s) to be used. Replace spaces with \_ and remove all non-alphanumeric characters 
  except - and \_. If you wish to use more than one outgroup, separate by spaces and wrap the entire list in double quotes
  (e.g. `"Bacillus_cereus_ATCC_14579 Bacillus_megaterium_NBRC_15308__ATCC_14581"`)
