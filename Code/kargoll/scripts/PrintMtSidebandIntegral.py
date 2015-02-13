#!/usr/bin/env python

execfile('loadPyRoot.py')
import argparse

parser = argparse.ArgumentParser(description='Obtain integral of mT sideband region for individual categories, to be used for WJet BG estimation.')
parser.add_argument('inputFile', help="location of the input root file")
args = parser.parse_args()

file = ROOT.TFile(args.inputFile,"READ")

categories =["0JetLow","0JetHigh","1JetHigh","1JetLow","1JetBoost","VBFLoose","VBFTight"]

for cat in categories:
    hist = file.Get("htotaumutauh_default_Cat" + cat + "MtSidebandData")
    print "Category", cat, ": Events in data mT sideband:", hist.Integral()

print "\n"
lowBin = -1
upBin = -1
for cat in categories:
    histW1 = file.Get("htotaumutauh_default_Cat" + cat + "MtMC_W_lnu")
    histW2 = file.Get("htotaumutauh_default_Cat" + cat + "MtMC_W_taunu")
    
    if (lowBin == -1):
        lowBin = histW1.FindFixBin(0.0)
    if (upBin == -1):
        upBin = histW1.FindFixBin(30.0) 
    
    Nwlnu = histW1.Integral(lowBin,upBin)
    Nwtaunu = histW2.Integral(lowBin,upBin)
    
    print 