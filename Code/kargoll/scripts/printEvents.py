#!/usr/bin/env python

import ROOT
import argparse
import numpy

files = []
files.append(ROOT.TFile("/net/scratch_cms/institut_3b/kargoll/TestPUBinning/workdirAnalysis_Apr_04_2014/MuTauSyncTree_VBFHiggs_TestPUBinning.root","READ"))
files.append(ROOT.TFile("/net/scratch_cms/institut_3b/kargoll/SyncExercise/VBFHToTauTauM125/workdirAnalysis_Apr_02_2014/MuTauSyncTree_VBFHToTauTauM125_2014_Mar_20.root","READ"))

names =["fineBinning","AachenDefault"]

trees = [files[0].Get("syncTree"), files[1].Get("syncTree")]

NEvents = [trees[0].GetEntries(), trees[1].GetEntries()]
if NEvents[0] != NEvents[1]:
    print "WARNING: Different number of events: ", NEvents[0], NEvents[1]

# read from tree 1
vars1 = []
for event in trees[0]:
    vars1.append(event.pt_1)

vars2 = []
for event in trees[1]:
    vars2.append(event.pt_1)

vars1.sort()
vars2.sort()

fmt = "{0:10} : {1:20}, {2:20}"
print fmt.format("Event", names[0], names[1])

unmatched1 = []
unmatched2 = []
for i in range(len(vars1)):
    #print fmt.format(i, vars1[i], vars2[i])
    if vars1[i] not in vars2:
        unmatched1.append(vars1[i])
    if vars2[i] not in vars1:
        unmatched2.append(vars2[i])
        
print "Unmatched values in", names[0], ":"
print unmatched1    
print "Unmatched values in", names[1], ":"
print unmatched2

