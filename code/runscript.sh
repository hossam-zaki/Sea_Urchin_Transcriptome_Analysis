#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -t 2:00:00

python3 miRNAsearch.py --seqs ../data/motifSearchData/depleted_sequences.txt --mirna ../data/miRNAData/miRNASeqs.txt --enriched
