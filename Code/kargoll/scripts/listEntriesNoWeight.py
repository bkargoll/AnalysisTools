#!/usr/bin/env python

import ROOT

files = []
files.append(ROOT.TFile("/net/scratch_cms/institut_3b/kargoll/SyncExercise/VBFHToTauTauM125/workdirAnalysis_Apr_02_2014/LOCAL_COMBINED_mutausync_default.root","READ"))
files.append(ROOT.TFile("/net/scratch_cms/institut_3b/kargoll/TestPUBinning/workdirAnalysis_Apr_04_2014/LOCAL_COMBINED_mutausync_default.root","READ"))
files.append(ROOT.TFile("/net/scratch_cms/institut_3b/kargoll/TestPUBinning/ClaudiasSample/workdirAnalysis_Apr_07_2014/LOCAL_COMBINED_mutausync_default.root","READ"))
files.append(ROOT.TFile("/net/scratch_cms/institut_3b/kargoll/TestPUBinning/ClaudiasSample/forcePUweightToOne/workdirAnalysis_Apr_07_2014/LOCAL_COMBINED_mutausync_default.root","READ"))
files.append(ROOT.TFile("/net/scratch_cms/institut_3b/kargoll/TestPUBinning/switchToFloatNumInter/workdirAnalysis_Apr_10_2014/LOCAL_COMBINED_mutausync_default.root"))
files.append(ROOT.TFile("/net/scratch_cms/institut_3b/kargoll/TestPUBinning/switchToFloatNumInter/Test2ndRun/workdirAnalysis_Apr_10_2014/LOCAL_COMBINED_mutausync_default.root"))

names =["Vladimir","fineBinning","Claudia","forcePUWeightTo1","fixedPUBug_1","fixedPUBug_2"]

hists = []
for file in files:
    hists.append(file.Get("mutausync_default_NPassMC_VBFHTauTauM125_noweight"))


for i in range(len(names)):
    print names[i], " : "
    for bin in range(16):
        print "\t", bin, " : " , hists[i].GetBinContent(bin) 
