import sys
import re
import argparse
from argparse import ArgumentParser

class kmp:
    def __init__(self, seqs, pattern, cds):  
        self.pathToSeq = seqs
        self.pathToPattern = pattern
        self.cds = cds
        self.pattern = None
        self.failure_fun = {}
        self.seq1 = None
        self.successArray = []
        self.seqsDict = {}
        self.cdsDict = {}
    def failure_function(self):
        self.failure_fun[0] = 0
        i = 0
        for j in range (1, len(self.pattern)):
            i = self.failure_fun[j-1]
            while self.pattern[j] != self.pattern[i] and i > 0:
                i = self.failure_fun[i-1]
            if self.pattern[j] != self.pattern[i] and i==0:
                self.failure_fun[j] = 0
            else:
                self.failure_fun[j] = i+1
    def retrieve_pattern(self):
        with open(self.pathToPattern) as file:
            for line in file:
                self.pattern = line.strip()
    def retrieve_seqs(self, seqs, dict, boolean1):
        with open(seqs) as f:
            boolean = False
            label = None
            string = ""
            for line in f:
                if boolean1:
                    matches = re.match(r">(.*)-mrna", line)
                else:
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
    def kmpsearching(self, seq):
        patInd = 0
        textInd = 0
        lenSeq = len(seq)
        lenPattern = len(self.pattern)
        self.successArray = []
        while lenSeq > textInd:
            if self.pattern[patInd] == seq[textInd]:
                patInd = patInd + 1
                textInd = textInd + 1
                if patInd == lenPattern:
                    self.successArray.append(textInd - lenPattern)
                    patInd = self.failure_fun[patInd-1]
            if textInd < lenSeq and self.pattern[patInd] is not seq[textInd]:
                if patInd != 0:
                    patInd = self.failure_fun[patInd-1]
                else:
                    textInd +=1
        if len(self.successArray) == 1:
            return (self.successArray[0])
        else: 
            if len(self.successArray) == 0:
                self.successArray = -1
            return (self.successArray)
        
    def kmp(self):
        self.retrieve_pattern()
        self.seqsDict = self.retrieve_seqs(self.pathToSeq, self.seqsDict, True)
        self.cdsDict = self.retrieve_seqs(self.cds, self.cdsDict, False) 
        self.failure_function()
        for i in self.seqsDict:
            ind = self.kmpsearching(self.seqsDict[i])
            if (ind) != -1:
                    if self.kmpsearching(self.cdsDict[i]) != -1:
                        with open("depleted.txt", "a+") as f:
                            f.write(f"{i}, {ind} \n")  
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--seqs',type=str,required=True)
    parser.add_argument('--pattern',type=str,required=True)
    parser.add_argument('--cds', type=str, required=True)
    config = parser.parse_args()
    kmp = kmp(config.seqs, config.pattern, config.cds)
    kmp.kmp()