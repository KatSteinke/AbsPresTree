import pathlib
import re
import sys

if len(sys.argv) == 2:
    antismash_dir = pathlib.Path(sys.argv[1]).resolve()
else:
    print("Supply the folder with antiSMASH results to rename. Usage: rename_clusters.py folder_name")
    sys.exit(1)

for entry in antismash_dir.iterdir():
    if entry.is_dir():
        for result in entry.iterdir():  # TODO: this looks so ugly
            if result.is_file():
                # if the name has only digits and "region"
                plain_name = re.search(r"^\d{1,7}\.region\d{3}\.gbk", result.name)
                if plain_name:
                    new_name = entry / (entry.name + "_" + result.name)
                    result.rename(new_name)
