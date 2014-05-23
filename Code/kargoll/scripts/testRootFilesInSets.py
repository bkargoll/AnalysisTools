#!/usr/bin/env python

import argparse
import ROOT

parser = argparse.ArgumentParser(description='Test if all the Sets have a valid root file containing the skimming tree.')
parser.add_argument('firstSet', help="Offset value, first set to start looking.", type=int)
parser.add_argument('nSets', help="Amount of Sets you have in your directory.", type=int)
args = parser.parse_args()

for i in range(args.firstSet-1,args.nSets):
    set = i+1
    print "Analyze Set", set
    
    file = ROOT.TFile("Set_"+str(set)+"/SKIMMED_NTUP.root","read")
    tree = file.Get("t")
    print "   entries:" , tree.GetEntries()
    file.Close()
