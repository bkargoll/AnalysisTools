#!/usr/bin/env python

import argparse

from ROOT import TH1F, TCanvas, gROOT, gStyle

parser = argparse.ArgumentParser(description='Analyze output of printTreeSize.C')
parser.add_argument('inputLogs', nargs='+', help="One or more files containing the output of printTreeSize.C")
args = parser.parse_args()

# create a dictionary to hold the lists of variable names
variables = {"total":[], "rest":[],'Vtx':[], 'Muon':[], 'PFTau':[], 'Electron':[], 'PFJet':[], 'MET':[], 'Track':[], 'MC':[]}

# loop over and analyze input files
for logName in args.inputLogs:
    log = open(logName, 'r')
    nEvts = 0
    lBranches = []
    sumOfSorted = 0
    # parse lines
    for line in log:
        if line.strip().startswith("Its size is"):
            nEvts = int(line.strip().split()[3].strip(","))
            variables["total"].append( int(line.strip().split()[5]) )
        if line.strip().startswith("Its branch"):
            name = line.split('\"')[1]
            bytes = float(line.strip().split()[4])
            lBranches.append([name, bytes])
    # sort numbers
    for key in variables.keys():
        if key in ["total","rest"]: continue
        nBytesVar = 0
        for b in lBranches:
            if b[0].startswith(key):
                nBytesVar += b[1]
        # add one entry to list for each logfile
        variables[key].append(nBytesVar / nEvts)
        sumOfSorted += nBytesVar / nEvts
    variables["rest"].append(variables["total"][-1] - sumOfSorted)
        
for key in variables.keys():
        print "\n++++++++++", key
        print variables[key]
        
# create histograms
for key in variables.keys():
    gROOT.Reset()
    gStyle.SetOptStat(0)
    c1 = TCanvas("c"+key)
    hist = TH1F(key, key, len(variables[key]), -0.5, len(variables[key])-0.5)
    for i_sam in range(len(variables[key])):
        hist.SetBinContent(i_sam+1, variables[key][i_sam])
        hist.GetXaxis().SetBinLabel(i_sam+1, args.inputLogs[i_sam].split("/")[-1] )
        hist.SetMinimum(0)
    hist.Draw()
    c1.SaveAs("sizeOf" + key + ".pdf")
