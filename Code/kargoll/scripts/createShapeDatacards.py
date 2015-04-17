#!/usr/bin/env python

execfile('loadPyRoot.py')
import argparse
import os
from sys import stdout
from array import array
from itertools import takewhile

# context manager to collect output from stdout
# printouts to stdout will be collected and appended to a list
# Usage:
# with Capturing() as output:
#     do_something(my_object)
# output is now a list containing the output
# http://stackoverflow.com/questions/16571150/how-to-capture-stdout-output-from-a-python-function-call
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

# save histogram to file
def SaveHist(hist):
    "save histogram into output file"
    global args, newDir
    hist.Write()
    if args.verbose:
        print '----> Histogram ', hist.GetName(), " written to ", newDir.GetPath()
        
# Scale all bins up to value upTo of a histogram hist by scaleBy
def ScaleHist(hist, scaleBy, upTo):
    "Scale a TH1 hist by a given amount scaleBy up to value upTo"
    global binning, binKey
    for binEdge in takewhile(lambda b: b <= upTo, binning[binKey]): # loop only over bins < 'upTo' GeV
        if binEdge == 0: continue
        binNum = hist.FindFixBin(binEdge - 0.1)
        hist.SetBinContent(binNum, scaleBy * hist.GetBinContent(binNum))

# QCD scaling
# put in a "function" (accessing global variables) for readability
def QCDScaling():
    "Scale QCD bins below 50GeV up by 10% for some categories"
    global cat, args, binning, binKey, inHistReb
    if cat in ["onejetlow", "onejethigh"]:
        if args.verbose: print '  -> Rescale by 1.1 for mtt < 50 GeV'
        ScaleHist(inHistReb, 1.1, 50.)

###### 

# parse arguments from command line
parser = argparse.ArgumentParser(description='Create shape datacards (ROOT files) from HToTauhTaumu analysis results.')
parser.add_argument('inputFolder', help="location of the input root files")
parser.add_argument('outputFolder', help="location to save the output datacards")
parser.add_argument('--categories', nargs="+",
                    default=["zerojetlow", "zerojethigh", "onejetlow", "onejethigh", "onejetboost", "vbfloose", "vbftight"],
                    help="list of categories to be processed (default: all of them)")
parser.add_argument('--fitVariable', default ='visibleMass') #default="shape_SVfitM")
parser.add_argument('--noQCDScaling', action='store_true', help="Disable scaling QCD up by 10%% in 1-jet low and high categories for mtt<50GeV.")
parser.add_argument('--verbose', action='store_true', help="Print more stuff")
args = parser.parse_args()

# define template structure of files, histograms, ...
inFileTemplate = 'LOCAL_COMBINED_<CAT>_default_LumiScaled.root'
histoTemplate = '<CAT>_default_<VAR><PROC>'
outFileTemplate = 'htt_<CHANNEL>.inputs-sm-<DATASET>.root'
qcdShapeUncTemplate = 'CMS_htt_QCDShape_<CHANNEL>_<CAT>_<DATASET><UPDOWN>'
# specify some fields in templates
histoTemplate = histoTemplate.replace('<VAR>', args.fitVariable)
outFileTemplate = outFileTemplate.replace('<CHANNEL>', 'mt').replace('<DATASET>', '8TeV')
qcdShapeUncTemplate = qcdShapeUncTemplate.replace('<CHANNEL>', 'mutau').replace('<DATASET>', '8TeV')

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
        
    # store scaled QCD histogram for uncertainties
    qcdHist = 0
    
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
            
            
            if sample == "QCD":
                if not args.noQCDScaling:
                    # scale QCD histogram up by 10% for mtt<50
                    QCDScaling()
                    
                qcdHist = inHistReb
            
            outHist.Add(inHistReb)
        
        # save histogram into output file
        SaveHist(outHist)
            
    # create QCD shape uncertainties
    histName = qcdShapeUncTemplate.replace('<CAT>', translate[cat])
    
    qcdScaleUp   = qcdHist.Clone( histName.replace('<UPDOWN>', 'Up'))
    ScaleHist(qcdScaleUp, 1.1, 50.)
    SaveHist(qcdScaleUp)
    
    qcdScaleDown = qcdHist.Clone( histName.replace('<UPDOWN>', 'Down'))
    ScaleHist(qcdScaleDown, 0.9, 50.)
    SaveHist(qcdScaleDown)


print 'Results have been stored in', outFileName



