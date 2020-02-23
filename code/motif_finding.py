import sys
import re
import argparse
import operator
from argparse import ArgumentParser
from KnuthMorrisPratt import *

class motifFinding:
    def __init__(self, seqs, pattern, cds, enriched):  
        self.pathToSeq = seqs
        self.pathToPattern = pattern
        self.cds = cds
        self.pattern = None
        self.failure_fun = {}
        self.seqsDict = {}
        self.cdsDict = {}
        self.foundSeqs = []
        self.seqs5 = []
        self.seqs3 = []
        if enriched == "True":
            self.classification = "enriched"
        else:
            self.classification = "depleted"
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
    def searchForMotif(self):
        self.retrieve_pattern()
        self.seqsDict = self.retrieve_seqs(self.pathToSeq, self.seqsDict, True)
        self.cdsDict = self.retrieve_seqs(self.cds, self.cdsDict, False) 
        for i in self.seqsDict:
            ind = kmpMatching(self.seqsDict[i], self.pattern)
            if ind != -1:
                if kmpMatching(self.cdsDict[i], self.pattern) == -1:
                    self.foundSeqs.append((i, ind))
        for j in self.foundSeqs:
            start = kmpMatching(self.seqsDict[j[0]], self.cdsDict[j[0]])
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
            with open(f"../results/{self.classification}_utr_5.txt", "a+") as file:
                file.write(f"{k[0]} has {k[1]} PRE element's in 5' utr \n") 
        for l in self.seqs3:
            with open(f"../results/{self.classification}_utr_3.txt", "a+") as file:
                file.write(f"{l[0]} has {l[1]} PRE element's in 3' utr \n") 
    def find_motif(self):
        for i in range(4, 16):
            with open(f"../results/{self.classification}_most_common_kmers.txt", "a+") as file:
                file.write(f"Length of k-mer = {i} \n")
            kmers = {}
            seqs_vals_list = []
            for j in self.seqs3:
                start = kmpMatching(self.seqsDict[j[0]], self.cdsDict[j[0]])
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
                    #if (string[l: (l + i)]) == "TGTAAAT":
                    #    quit()
            for m in kmers:
                seqs_vals_list.append((m, kmers[m]))
            seqs_vals_list.sort(key=lambda tup: tup[1], reverse=True)
            lis = []
            for val in range(0, 50):
                seq = seqs_vals_list[val][0]
                counter = 0
                for h in self.seqs3:
                    start = kmpMatching(self.seqsDict[h[0]], self.cdsDict[h[0]])
                    end = start + len(self.cdsDict[h[0]])
                    string = self.seqsDict[h[0]][end:len(self.seqsDict[h[0]])]
                    if kmpMatching(string, seq) != -1:
                        counter += 1
                lis.append((seq, counter))
            lis.sort(key=lambda tup: tup[1], reverse=True)
            for q in range(0, len(lis)):
                with open(f"../results/{self.classification}_most_common_kmers.txt", "a+") as fi:
                    fi.write(f"kmer: {lis[q][0]} appears: {lis[q][1]} out of {len(self.seqs3)} \n")



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--seqs',type=str,required=True)
    parser.add_argument('--pattern',type=str,required=True)
    parser.add_argument('--cds', type=str, required=True)
    parser.add_argument('--enriched', type=str, required=True)
    config = parser.parse_args()
    kmp = motifFinding(config.seqs, config.pattern, config.cds, config.enriched)
    kmp.searchForMotif()
    kmp.find_motif()
