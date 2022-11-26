import sys, os, time
import biorun.libs.placlib as plac
from tqdm import tqdm

URL = "http://data.biostarhandbook.com/make/code.tar.gz"

@plac.flg('update', "update existing files", abbrev="u")
def run(update=False):
    """
    Downloads the data from the URL.
    """
    import requests
    import tarfile
    import io

    # Download the data.
    resp = requests.get(URL)
    resp.raise_for_status()

    # Extract the data.
    tar = tarfile.open(fileobj=io.BytesIO(resp.content))

    # Flags each tarinfo for extraction
    def check_file(tarinfo):
        return os.path.isfile(tarinfo.name), tarinfo

    pairs = list(map(check_file, tar))
    pairs = list(filter(lambda x: x[1].isfile(), pairs))

    found = filter(lambda x: x[0], pairs)
    found = list(map(lambda x: x[1], found))

    not_found = filter(lambda x: not x[0], pairs)
    not_found = list(map(lambda x: x[1], not_found))

    for tf in found:
        if update:
            print(f"# Updating: {tf.name}")
            tar.extract(member=tf)

    for tf in not_found:
        print(f"# Writing: {tf.name}")
        tar.extract(member=tf)

    if not_found:
        count = len(not_found)
        label = "file" if count == 1 else "files"
        print(f"# Created {count} {label}.")

    if found:
        if not update:
            print(f"# Skipped {len(found)} files that already exist.")
            print(f"# Use -u to update existing files.")
        else:
            print(f"# Updated {len(found)} files.")
    print("# Biostar Workflows: https://www.biostarhandbook.com/")

def main():
    """
    Entry point for the script.
    """
    plac.call(run)

if __name__ == '__main__':
    main()
