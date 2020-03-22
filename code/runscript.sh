#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -t 1:00:00

python3 miRNAsearch.py --seqs ../data/depleted_sequences_cds.txt --mirna ../data/miRNASeqs.txt
