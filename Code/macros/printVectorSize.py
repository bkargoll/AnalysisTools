#!/usr/bin/env python

import argparse
import ROOT
import sys

# allow accessing additional data types
ROOT.gROOT.ProcessLine('.L Loader.C+')
#ROOT.gInterpreter.GenerateDictionary("std::vector<std::vector<float>>", "vector;vector;float")

parser = argparse.ArgumentParser(description='Compare size of 2 vectors in Ntuple.')
parser.add_argument('NtupleFile', help="Location of the input root file containing the Ntuple.")
parser.add_argument('Vector1', help="Name of the first vector to check")
parser.add_argument('Vector2', help="Name of the second vector to check")
parser.add_argument('--treeName', dest="treeName", default="t",help="Name of the tree in the Ntuple file.")
parser.add_argument('--verbose', action='store_true',help="Print more stuff.")
args = parser.parse_args()

## load the ROOT file
ntupleFile = ROOT.TFile(args.NtupleFile, "READ")
tree = ntupleFile.Get(args.treeName)
nevts = tree.GetEntries()

wrongSizes = []

for event in tree:
    size1 = 0
    size2 = 0
    size1 = eval("event."+args.Vector1+".size()")   
    size2 = eval("event."+args.Vector2+".size()")
    if (size1 != size2):
        #print "Warning! The 2 vectors have different size!\n\tsize(",args.Vector1,") =", size1, ", size(",args.Vector2,") =", size2
        wrongSizes.append([size1,size2])
    
print nevts, "have been investigated"
if (wrongSizes == []):
    print "The vectors have always the same length"
else:
    print "The vectors are of DIFFERENT length"
    if (args.verbose):
        print wrongSizes
