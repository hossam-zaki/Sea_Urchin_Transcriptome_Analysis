import sys
import re
import argparse
import operator
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
        self.foundSeqs = []
        self.seqs5 = []
        self.seqs3 = []
    def failure_function(self, pattern):
        self.failure_fun={}
        self.failure_fun[0] = 0
        i = 0
        for j in range (1, len(pattern)):
            i = self.failure_fun[j-1]
            while pattern[j] != pattern[i] and i > 0:
                i = self.failure_fun[i-1]
            if pattern[j] != pattern[i] and i==0:
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
    def kmpsearching(self, seq, pattern):
        patInd = 0
        textInd = 0
        lenSeq = len(seq)
        lenPattern = len(pattern)
        self.successArray = []
        while lenSeq > textInd:
            if pattern[patInd] == seq[textInd]:
                patInd = patInd + 1
                textInd = textInd + 1
                if patInd == lenPattern:
                    self.successArray.append(textInd - lenPattern)
                    patInd = self.failure_fun[patInd-1]
            if textInd < lenSeq and pattern[patInd] is not seq[textInd]:
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
        self.failure_function(self.pattern)
        for i in self.seqsDict:
            ind = self.kmpsearching(self.seqsDict[i], self.pattern)
            if (ind) != -1:
                if self.kmpsearching(self.cdsDict[i], self.pattern) == -1:
                    self.foundSeqs.append((i, ind))
                    with open("enriched.txt", "a+") as f:
                        f.write(f"{i}, {ind} \n")  
        for j in self.foundSeqs:
            self.failure_function(self.cdsDict[j[0]])
            start = self.kmpsearching(self.seqsDict[j[0]], self.cdsDict[j[0]])
            end = start + len(self.cdsDict[j[0]])
            counter_3 = 0
            counter_5 = 0
            if type(j[1]) == list:
                for ind in j[1]:
                    if ind < start:
                        counter_5 +=1
                    if ind > end:
                        counter_3 +=1
            else:
                if j[1] < start:
                    counter_5 +=1
                if j[1] > end:
                    counter_3 +=1
            if counter_5 > 0:
                self.seqs5.append((j[0], counter_5))
            if counter_3 > 0:
                self.seqs3.append((j[0], counter_3))
        self.seqs5.sort(key=lambda tup: tup[1], reverse=True)
        self.seqs3.sort(key=lambda tup: tup[1], reverse=True)
        for k in self.seqs5:
            with open("enriched_utr_5.txt", "a+") as file:
                file.write(f"{k[0]} has {k[1]} PRE element's in 5' utr \n") 
        for l in self.seqs3:
            with open("enriched_utr_3.txt", "a+") as file:
                file.write(f"{l[0]} has {l[1]} PRE element's in 3' utr \n") 
    def find_motif(self):
        for i in range(4, 16):
            with open("enriched_most_common_kmers.txt", "a+") as file:
                file.write(f"Length of k-mer = {i} \n")
            kmers = {}
            seqs_vals_list = []
            for j in self.seqs3:
                self.failure_function(self.cdsDict[j[0]])
                start = self.kmpsearching(self.seqsDict[j[0]], self.cdsDict[j[0]])
                end = start + len(self.cdsDict[j[0]])
                string = self.seqsDict[j[0]][end:len(self.seqsDict[j[0]])]
                for l in range (0, len(string)):
                    if string[l: (l + i)] in kmers:
                        kmers[string[l: (l + i)]] += 1
                        continue
                    elif len(string[l: (l + i)]) != i:
                        continue
                    else: 
                        kmers[string[l: (l + i)]] = 1
            for m in kmers:
                seqs_vals_list.append((m, kmers[m]))
            seqs_vals_list.sort(key=lambda tup: tup[1], reverse=True)
            lis = []
            for val in range(0, 20):
                seq = seqs_vals_list[val][0]
                self.failure_function(seq)
                counter = 0
                for h in self.seqs3:
                    if self.kmpsearching(self.seqsDict[h[0]], seq) != -1:
                        counter += 1
                lis.append((seq, counter))
                print(seq, counter)
                
            lis.sort(key=lambda tup: tup[1], reverse=True)
            for q in range(0, len(lis)):
                with open("enriched_most_common_kmers.txt", "a+") as fi:
                    fi.write(f"kmer: {lis[q][0]} appears: {lis[q][1]} out of {len(self.foundSeqs)} \n")



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--seqs',type=str,required=True)
    parser.add_argument('--pattern',type=str,required=True)
    parser.add_argument('--cds', type=str, required=True)
    config = parser.parse_args()
    kmp = kmp(config.seqs, config.pattern, config.cds)
    kmp.kmp()
    kmp.find_motif()
