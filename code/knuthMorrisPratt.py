import sys
import re
import argparse
import operator
from argparse import ArgumentParser

class kmpSearch:
    def __init__(self):
        self.failure_fun = {}

    def computeFailureFunction(self, pattern):
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
    
    def kmpMatching(self, seq, pattern):
        self.computeFailureFunction(pattern)
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
            return self.successArray[0]
        else: 
            if len(self.successArray) == 0:
                self.successArray = -1
            return self.successArray