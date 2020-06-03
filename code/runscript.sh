#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -t 2:00:00

python3 miRNAsearch.py --mirna ../data/miRNAData/miRNASeqs.txt
