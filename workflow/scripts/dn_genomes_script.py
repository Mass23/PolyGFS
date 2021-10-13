import os
import glob
import pandas as pd

for line in open(snakemake.input[0], 'r').readlines():
    if line.startswith('GCA'):
        args = ['ncbi-genome-download','-o',snakemake.output[0],'--formats','fasta','--flat-output','-s','genbank','-A',line,'all']
        subprocess.call(' '.join(args), shell = True)
    elif input.startswith('GCF'):
        args = ['ncbi-genome-download','-o',snakemake.output[0],'--formats','fasta','--flat-output','-s','refseq','-A',line,'all']
        subprocess.call(' '.join(args), shell = True)
    else:
        print(str(line) + ' is not a refseq or genbank accession!')
os.system('touch ' + str(snakemake.output[1]))
