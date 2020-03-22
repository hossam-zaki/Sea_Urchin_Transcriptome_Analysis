import sys
import re
import argparse
import operator 
from argparse import ArgumentParser
import KnuthMorrisPratt
import PyPDF2
import local_align

class miRNASearch:
    def __init__(self, seqsFile, miRNA):  
        self.pathToSeq = seqsFile
        self.pathToMiRNA = miRNA
        self.seqsDict = {}
        self.miRnaDict = {}
        self.scoringMatrix = {
                ('C', 'C'): 1,
                ('C', 'A'): -1,
                ('C', 'G'): -1,
                ('C', 'T'): -1,
                ('A', 'A'): 1,
                ('A', 'C'): -1,
                ('A', 'G'): -1,
                ('A', 'T'): -1,
                ('G', 'G'): 1,
                ('G', 'C'): -1,
                ('G', 'A'): -1,
                ('G', 'T'): -1,
                ('T', 'T'): 1,
                ('T', 'C'): -1,
                ('T', 'A'): -1,
                ('T', 'G'): -1,
                ('N', 'G'): 0,
                ('N', 'C'): 0,
                ('N', 'A'): 0,
                ('N', 'T'): 0,
                ('T', 'N'): 0,
                ('A', 'N'): 0,
                ('G', 'N'): 0,
                ('C', 'N'): 0   
        }
    def retrieve_mirna(self):
        di = {}
        with open(self.pathToMiRNA) as file:
            for line in file:
                split_line = line.strip().split(" ")
                di[split_line[0]] = self.getComplement(split_line[1])
        return di
    def retrieve_seqs(self, seqs, dict, boolean1):
        with open(seqs) as f:
            boolean = False
            label = None
            string = ""
            for line in f:
                matches = re.match(r">(.*)-cds", line)
                if matches != None:
                    boolean = True
                    label = matches[1]
                    continue
                if not line.strip():
                    boolean = False
                    dict[label] = string
                    label = None
                    string = ""
                    continue
                if boolean == True:
                    string += line.strip()
                    continue
        return dict 
    def miRNAQuery(self):
        print(self.miRnaDict)
        with(open('../results/miRNAdata.txt', 'w+')) as f:
            for mirna in self.miRnaDict:
                for cds in self.seqsDict:
                    align = local_align.local_aligning(self.seqsDict[cds], self.miRnaDict[mirna], self.scoringMatrix, -1)
                    score = align.align_locally()
                    if score > 45:
                        f.write(f"{mirna} was found  in {cds} with the score {score} and the following alignment \n")
                        align.trace_back()
                        f.write(f"{align.final_top_seq} \n")
                        f.write(f"{align.final_bottom_seq} \n")
                        f.write("\n")
                    
    def searchForMiRNA(self):
        self.miRnaDict = self.retrieve_mirna()
        self.seqsDict = self.retrieve_seqs(self.pathToSeq, self.seqsDict, True)
        self.miRNAQuery()
    def getComplement(self, sequence):
        returnedSeq = ""
        print(sequence)
        for char in sequence:
            if char == 'u':
                returnedSeq += "T"
            if char == 'a':
                returnedSeq += "A"
            if char == 'g':
                returnedSeq += "G"
            if char == 'c':
                returnedSeq += "C"
        print(len(returnedSeq))
        return returnedSeq
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--seqs',type=str,required=True)
    parser.add_argument('--mirna',type=str,required=True)
    config = parser.parse_args()
    kmp = miRNASearch(config.seqs, config.mirna)
    kmp.searchForMiRNA()
