import argparse
import operator
import re
import sys
from argparse import ArgumentParser

import KnuthMorrisPratt
import local_align


class miRNASearch:
    def __init__(self, seqsFile, miRNA, enriched):
        self.enrichedBool = enriched
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
                # used to get all of the mrna that the mirna will downregulate
                matches = re.match(r">(.*)-mrna", line)
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
        file = ''
        miRNAtoTranscript = {}
        transcripttoMiRNA = {}
        miRNACount = {}
        transcriptMiRNACount = {}
        if self.enrichedBool:
            mainFile = '../results/miRNAResults/miRNAEnriched.txt'
            miRNAtoTranscriptFile = '../results/miRNAResults/EnrichedmiRNAtoTranscript.txt'
            transcripttoMiRNAFile = '../results/miRNAResults/EnrichedTranscripttoMiRNA.txt'
        else:
            mainFile = '../results/miRNAResults/miRNADepleted.txt'
            miRNAtoTranscriptFile = '../results/miRNAResults/DepletedmiRNAtoTranscript.txt'
            transcripttoMiRNAFile = '../results/miRNAResults/DepletedTranscripttoMiRNA.txt'
        with(open(mainFile, 'w+')) as f:
            for mirna in self.miRnaDict:
                if mirna not in miRNAtoTranscript:  # which transcripts are bound to the mirna
                    miRNAtoTranscript[mirna] = []
                if mirna not in miRNACount:  # number of transcripts the mirna bind to
                    miRNACount[mirna] = 0
                for seq in self.seqsDict:
                    if seq not in transcripttoMiRNA:  # which mirna bind to which transcripts
                        transcripttoMiRNA[seq] = []
                    if seq not in transcriptMiRNACount:  # number of mirna bound to the transcript
                        transcriptMiRNACount[seq] = 0
                    align = local_align.local_aligning(
                        self.seqsDict[seq], self.miRnaDict[mirna], self.scoringMatrix, -1)
                    score = align.align_locally()
                    if score > 45:
                        f.write(
                            f"{mirna} was found  in {seq} with the score {score} and the following alignment \n")
                        align.trace_back()
                        f.write(f"{align.final_top_seq} \n")
                        f.write(f"{align.final_bottom_seq} \n")
                        f.write("\n")
                        miRNAtoTranscript[mirna].append(seq)
                        transcripttoMiRNA[seq].append(mirna)
                        transcriptMiRNACount[seq] += 1
                        miRNACount[mirna] += 1
        miRNACount.sort(key=lambda tup: tup[1], reverse=True)
        with(open(miRNAtoTranscriptFile, "w+")) as f:
            f.write("-------------------miRNA counts------------------- \n")
            f.write("\n")
            for miRNA in miRNACount:
                f.write(f"{miRNA} binds to {miRNACount[miRNA]} transcripts \n")
            f.write("\n")
            f.write(
                "-------------------miRNA found in transcripts------------------- \n")
            f.write("\n")
            for miRNAinTranscript in miRNAtoTranscript:
                f.write(
                    f"{miRNAinTranscript} binds to these sequences {miRNAtoTranscript[miRNAinTranscript]} \n")
        transcriptMiRNACount.sort(key=lambda tup: tup[1], reverse=True)
        with(open(transcripttoMiRNAFile, "w+")) as f:
            f.write("-------------------Transcript counts------------------- \n")
            f.write("\n")
            for transcript in transcriptMiRNACount:
                f.write(
                    f"{transcript} is bound by {transcriptMiRNACount[transcript]} miRNA's \n")
            f.write("\n")
            f.write(
                "-------------------Transcript bound by miRNA------------------- \n")
            f.write("\n")
            for transcript in transcripttoMiRNA:
                f.write(
                    f"{transcript} was found to be bound by these miRNA's {transcripttoMiRNA[transcript]} \n")

    def searchForMiRNA(self):
        self.miRnaDict = self.retrieve_mirna()
        self.seqsDict = self.retrieve_seqs(self.pathToSeq, self.seqsDict, True)
        self.miRNAQuery()

    def getComplement(self, sequence):
        returnedSeq = ""
        for char in sequence:
            if char == 'u':
                returnedSeq += "T"
            if char == 'a':
                returnedSeq += "A"
            if char == 'g':
                returnedSeq += "G"
            if char == 'c':
                returnedSeq += "C"
        return returnedSeq


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--seqs', type=str, required=True)
    parser.add_argument('--mirna', type=str, required=True)
    parser.add_argument('--enriched', action="store_true")
    config = parser.parse_args()
    kmp = miRNASearch(config.seqs, config.mirna, config.enriched)
    kmp.searchForMiRNA()
