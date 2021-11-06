# PolyGFS

Analysis of GFS metagenome assembled genomes (MAGs) in a coassembly to assess polymorphism in multiple populations, across genomic regions.

## Downlaoding the repository

```bash
git clone --recurse-submodules https://github.com/Mass23/PolyGFS.git
cd PolyGFS
```

## Note(s):
- After the repo is downloaded, complete the setup as follows:
    - Adjust the paths in the `mantis.config` file to where the databases need to be downloaded
    - The databases for running `MANTIS` need to by running the below:
```bash
# Adjust the paths in the `mantis.config` file to where the databases need to be downloaded
python submodules/mantis/ setup_databases --mantis_config mantis.config
```
