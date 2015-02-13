#!/usr/bin/env python

execfile('loadPyRoot.py')
from glob import glob
import argparse

parser = argparse.ArgumentParser(description='Get WJet Yield from Sets. To be cross-checked against yield obtained from data-driven method.')
parser.add_argument('filePattern', help="Pattern of the root input files. CATEGORY will be replaced in script by the category names.")
args = parser.parse_args()

categories = ["ZeroJetLow","ZeroJetHigh","OneJetHigh","OneJetLow","OneJetBoost","VBFLoose","VBFTight","Inclusive"]
catWJetYield = []

for cat in categories:
    #print "*********", cat, "*********"
    filePattern = args.filePattern.replace("CATEGORY",cat.lower())
    fileList = glob(filePattern)
    nEvts = 0
    for file in fileList:
        #print file
        f = ROOT.TFile(file)
        hist = f.Get(cat.lower()+"_default_NPassMC_W_lnu")
        nEvts += hist.GetBinContent(18)
        #print "WLnu: ", nEvts
        hist = f.Get(cat.lower()+"_default_NPassMC_W_taunu")
        nEvts += hist.GetBinContent(18)
        #print "WTaunu: ", nEvts
    print "Category", cat, "has", nEvts 
        
print "Done"



