#!/usr/bin/env python

execfile('loadPyRoot.py')
import argparse
import os
from sys import stdout
from array import array

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout


# parse arguments from command line
parser = argparse.ArgumentParser(description='Create shape datacards (ROOT files) from HToTauhTaumu analysis results.')
parser.add_argument('inputFolder', help="location of the input root files")
parser.add_argument('outputFolder', help="location to save the output datacards")
parser.add_argument('--categories', nargs="+",
                    default=["zerojetlow", "zerojethigh", "onejetlow", "onejethigh", "onejetboost", "vbfloose", "vbftight"],
                    help="list of categories to be processed (default: all of them)")
parser.add_argument('--fitVariable', default ='visibleMass') #default="shape_SVfitM")
parser.add_argument('--verbose', action='store_true', help="Print more stuff")
args = parser.parse_args()

# define template structure of files, histograms, ...
inFileTemplate = 'LOCAL_COMBINED_<CAT>_default_LumiScaled.root'
histoTemplate = '<CAT>_default_<VAR><PROC>'
outFileTemplate = 'htt_<CHANNEL>.inputs-sm-<DATASET>.root'
# specify some fields in templates
histoTemplate = histoTemplate.replace('<VAR>', args.fitVariable)
outFileTemplate = outFileTemplate.replace('<CHANNEL>', 'mt').replace('<DATASET>', '8TeV')

# define binning to be used in mass plots
binning = {'nonVBF' :   array('d', [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]),#, 225, 250, 275, 300, 325, 350]),
           'VBF'    :   array('d', [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200])}#, 250, 300, 350])}

# create dictionary to map output histo names against list of input histo names (processes)
procs = {'ZTT'      :   ['MC_DY_embedded'],
         'W'        :   ['MC_W_lnu', 'MC_W_taunu'],
         'QCD'      :   ['QCD'],
         'ZL'       :   ['MC_DY_ll'],
         'ZJ'       :   [],  # todo: split ZL and ZJ
         'TT'       :   ['MC_ttbar', 'MC_tw', 'MC_tbarw'],
         'VV'       :   ['MC_WW_2l2nu', 'MC_WZ_2l2q', 'MC_WZ_3l1nu', 'MC_ZZ_4l', 'MC_ZZ_2l2nu', 'MC_ZZ_2l2q'],
         'ggH125'   :   ['MC_ggFHTauTauM125'], # todo: allow for different masses
         'qqH125'   :   ['MC_VBFHTauTauM125'],
         'WH125'    :   ['MC_WHZHTTHTauTauM125'],
         'ZH125'    :   [],  # todo: split WH from ZH?
         'data_obs' :   ['Data']  
         }

# create dictionary to map Ntuple category names to official category names
translate = {"inclusive"    :   "inclusive",
             "zerojetlow"   :   "0jet_medium",
             "zerojethigh"  :   "0jet_high",
             "onejetlow"    :   "1jet_medium",
             "onejethigh"   :   "1jet_high_lowhiggs",
             "onejetboost"  :   "1jet_high_mediumhiggs",
             "vbfloose"     :   "vbf_loose",
             "vbftight"     :   "vbf_tight"
             }

# create output root file
if not os.path.exists(args.outputFolder):
    os.makedirs(args.outputFolder)
outFileName = args.outputFolder + '/' + outFileTemplate
outFile = ROOT.TFile(outFileName, "RECREATE")

# loop over categories
for cat in args.categories:
    # read in output from analysis
    inFile = ROOT.TFile(args.inputFolder + '/' + inFileTemplate.replace('<CAT>', cat), "READ")
    
    # update directory structure in output file
    outFile.cd()
    newDir = outFile.mkdir('muTau_' + translate[cat])
    newDir.cd()
    
    # choose correct binning
    if 'vbf' in cat:
        binKey = 'VBF'
    else:
        binKey = 'nonVBF'
    
    # loop over processes to be put in the datacard
    for procKey in procs.keys():
        outHist = ROOT.TH1D(procKey, procKey, len(binning[binKey])-1, binning[binKey])
        outHist.Sumw2()
        
        # loop over samples to be combined into this process
        for sample in procs[procKey]:
            histName = histoTemplate.replace('<CAT>', cat).replace('<PROC>', sample)
            if args.verbose : print 'Loading histogram', histName
            inHist = inFile.Get(histName)
            inHistReb = inHist.Rebin(len(binning[binKey])-1, 'inHistReb', binning[binKey])
            outHist.Add(inHistReb)
        
        # save histogram into output file
        outHist.Write()
        if args.verbose:
            print '----> Histogram ', outHist.GetName(), " written to ", newDir.GetPath()


print 'Results have been stored in', outFileName



