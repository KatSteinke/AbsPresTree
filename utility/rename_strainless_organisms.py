#!/usr/bin/env python
import os
import sys
from argparse import ArgumentParser

from Bio import SeqIO

parser = ArgumentParser(description="Add strain names to organism fields.")
parser.add_argument("gbk_file", action="store", help="Genbank file to edit")
parser.add_argument("--overwrite", action="store_true", help="Overwrite files (default: create new ones)")
args = parser.parse_args()
gbk_file = args.gbk_file
overwrite = args.overwrite

def rename_organism(genbank_file: str, overwrite=False) -> None:
    """Extracts strain name from a genbank file and adds it to the file's
    organism attribute if it's not alread present.
    """
    out_name = genbank_file
    if not overwrite:
        file_name, extension = os.path.splitext(genbank_file)
        out_name = file_name + "_renamed" + extension
    strain_name = ""
    # find strain attribute
    with open(genbank_file, "r") as infile:
        line = infile.readline()
        while line and not strain_name:
            if "/strain" in line:
                strain_parts = line.strip().split("=")
                strain_name = strain_parts[1].replace('"', "")
            line = infile.readline()
    # change organism attribute
    records = list(SeqIO.parse(genbank_file, "genbank"))
    for record in records:
        if strain_name not in record.annotations["organism"]:
            record.annotations["organism"] += (" " + strain_name)
        if strain_name not in record.annotations["source"]:
            record.annotations["source"] += (" " + strain_name)
    SeqIO.write(records, out_name, "genbank")

if __name__ == "__main__":
    rename_organism(gbk_file, overwrite)
