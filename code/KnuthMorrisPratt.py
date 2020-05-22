import argparse
import operator
import re
import sys
from argparse import ArgumentParser


def computeFailureFunction(pattern):
    failure_fun = {}
    failure_fun[0] = 0
    i = 0
    for j in range(1, len(pattern)):
        i = failure_fun[j-1]
        while pattern[j] != pattern[i] and i > 0:
            i = failure_fun[i-1]
        if pattern[j] != pattern[i] and i == 0:
            failure_fun[j] = 0
        else:
            failure_fun[j] = i+1
    return failure_fun


def kmpMatching(seq, pattern):
    failure_fun = computeFailureFunction(pattern)
    patInd = 0
    textInd = 0
    lenSeq = len(seq)
    lenPattern = len(pattern)
    successArray = []
    while lenSeq > textInd:
        if pattern[patInd] == seq[textInd]:
            patInd = patInd + 1
            textInd = textInd + 1
            if patInd == lenPattern:
                successArray.append(textInd - lenPattern)
                patInd = failure_fun[patInd-1]
        if textInd < lenSeq and pattern[patInd] is not seq[textInd]:
            if patInd != 0:
                patInd = failure_fun[patInd-1]
            else:
                textInd += 1
    if len(successArray) == 1:
        return successArray[0]
    else:
        if len(successArray) == 0:
            successArray = -1
        return successArray
