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
    norm = hist.Integral()
    global binning, binKey
    for binEdge in takewhile(lambda b: b <= upTo, binning[binKey]): # loop only over bins < 'upTo' GeV
        if binEdge == 0: continue
        binNum = hist.FindFixBin(binEdge - 0.1)
        hist.SetBinContent(binNum, scaleBy * hist.GetBinContent(binNum))
    # scale to keep normalization
    hist.Scale(norm / hist.Integral())

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
parser.add_argument('outputFolder', help="location to save the output datacards")
parser.add_argument('--categories', nargs="+",
                    default=["inclusive", "zerojetlow", "zerojethigh", "onejetlow", "onejethigh", "onejetboost", "vbfloose", "vbftight"],
                    help="list of categories to be processed (default: all of them)")
parser.add_argument('--fitVariable', default ='SVfitM')
parser.add_argument('--noQCDScaling', action='store_true', help="Disable scaling QCD up by 10%% in 1-jet low and high categories for mtt<50GeV.")
parser.add_argument('--verbose', action='store_true', help="Print more stuff")
args = parser.parse_args()

# define template structure of files, histograms, ...
inFileTemplate = 'LOCAL_COMBINED_<CAT>_default_LumiScaled.root'
inHistoTemplate = '<CAT>_default_shape_<VAR><SAMPLE>'
outFileTemplate = 'htt_<CHANNEL>.inputs-sm-<DATASET>.root'
qcdShapeUncTemplate = 'QCD_CMS_htt_QCDShape_<CHANNEL>_<CAT>_<DATASET><UPDOWN>'
zlShapeUncTemplate = 'ZL_CMS_htt_ZLScale_<CHANNEL>_<DATASET><UPDOWN>'
tauESShapeUncTemplate = '<PROC>_CMS_scale_t_mutau_<DATASET><UPDOWN>'
# specify some fields in templates
histoTemplate = inHistoTemplate.replace('<VAR>', args.fitVariable)
inputZlShapeTemplate = inHistoTemplate.replace('<VAR>', args.fitVariable + '_ZLScale<UPDOWN>')
outFileTemplate = outFileTemplate.replace('<CHANNEL>', 'mt').replace('<DATASET>', '8TeV')
qcdShapeUncTemplate = qcdShapeUncTemplate.replace('<CHANNEL>', 'mutau').replace('<DATASET>', '8TeV')
zlShapeUncTemplate = zlShapeUncTemplate.replace('<CHANNEL>', 'mutau').replace('<DATASET>', '8TeV')
tauESShapeUncTemplate = tauESShapeUncTemplate.replace('<DATASET>', '8TeV')

# create map of categories to input folders

baseDir = "/net/scratch_cms/institut_3b/kargoll/"
inputDirs = {"inclusive"    :   baseDir + "workdirInclusive3ProngAug_05_2016",
             "zerojetlow"   :   baseDir + "workdirZeroJetAug_05_2016",
             "zerojethigh"  :   baseDir + "workdirZeroJetAug_05_2016",
             "onejetlow"    :   baseDir + "workdirOneJetAug_05_2016",
             "onejethigh"   :   baseDir + "workdirOneJetAug_05_2016",
             "onejetboost"  :   baseDir + "workdirOneJetAug_05_2016",
             "vbfloose"     :   baseDir + "workdirVBFAug_05_2016",
             "vbftight"     :   baseDir + "workdirVBFAug_05_2016"
             }
inputDirsTauESUp = {"inclusive"    :   baseDir + "workdirInclusive3Prong_TauESUpAug_06_2016",
             "zerojetlow"   :   baseDir + "workdirZeroJet_TauESUpAug_06_2016",
             "zerojethigh"  :   baseDir + "workdirZeroJet_TauESUpAug_06_2016",
             "onejetlow"    :   baseDir + "workdirOneJet_TauESUpAug_06_2016",
             "onejethigh"   :   baseDir + "workdirOneJet_TauESUpAug_06_2016",
             "onejetboost"  :   baseDir + "workdirOneJet_TauESUpAug_06_2016",
             "vbfloose"     :   baseDir + "workdirVBF_TauESUpAug_06_2016",
             "vbftight"     :   baseDir + "workdirVBF_TauESUpAug_06_2016"
             }
inputDirsTauESDown = {"inclusive"    :   baseDir + "workdirInclusive3Prong_TauESDownAug_06_2016",
             "zerojetlow"   :   baseDir + "workdirZeroJet_TauESDownAug_06_2016",
             "zerojethigh"  :   baseDir + "workdirZeroJet_TauESDownAug_06_2016",
             "onejetlow"    :   baseDir + "workdirOneJet_TauESDownAug_07_2016",
             "onejethigh"   :   baseDir + "workdirOneJet_TauESDownAug_07_2016",
             "onejetboost"  :   baseDir + "workdirOneJet_TauESDownAug_07_2016",
             "vbfloose"     :   baseDir + "workdirVBF_TauESDownAug_07_2016",
             "vbftight"     :   baseDir + "workdirVBF_TauESDownAug_07_2016"
             }

# define binning to be used in mass plots
binning = {'nonVBF' :   array('d', [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 225, 250, 275, 300, 325, 350]),
           'VBF'    :   array('d', [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 250, 300, 350])}

# list of Higgs masses
masses = [110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160]

# create dictionary to map output histo names against list of input histo names (processes)
procs = {'ZTT'      :   ['MC_DY_embedded'],
         'W'        :   ['MC_W_lnu', 'MC_W_taunu'],
         'QCD'      :   ['QCD'],
         'ZL'       :   ['MC_DY_ll'],
         'ZJ'       :   [],  # todo: split ZL and ZJ
         'TT'       :   ['MC_ttbar', 'MC_tw', 'MC_tbarw'],
         'VV'       :   ['MC_WW_2l2nu', 'MC_WZ_2l2q', 'MC_WZ_3l1nu', 'MC_ZZ_4l', 'MC_ZZ_2l2nu', 'MC_ZZ_2l2q'],
         'ggH<M>'   :   ['MC_ggFHTauTauM<M>'],
         'qqH<M>'   :   ['MC_VBFHTauTauM<M>'],
         'WH<M>'    :   ['MC_WHZHTTHTauTauM<M>'],
         'ZH<M>'    :   [],  # todo: split WH from ZH?
         'data_obs' :   ['Data']  
         }

# list of Higgs masses
masses = [110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160]
massProcesses = []
for procKey in procs:
    if '<M>' in procKey:
        massProcesses.append(procKey)
                
for procKey in massProcesses:
    entry = procs.pop(procKey)
    for m in masses:
        procs[procKey.replace('<M>', str(m))] = [e.replace('<M>', str(m)) for e in entry]

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
    inFileName = inFileTemplate.replace('<CAT>', cat)
    inFiles = { "main"      : ROOT.TFile(inputDirs[cat] + '/' + inFileName, "READ"),
                "TauESUp"   : ROOT.TFile(inputDirsTauESUp[cat] + '/' + inFileName, "READ"),
                "TauESDown" : ROOT.TFile(inputDirsTauESDown[cat] + '/' + inFileName, "READ")
                }
    
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
        
        tauESName = tauESShapeUncTemplate.replace('<PROC>', procKey)
        tauESUpHist = ROOT.TH1D(tauESName.replace('<UPDOWN>', 'Up'), tauESName.replace('<UPDOWN>', 'Up'), len(binning[binKey])-1, binning[binKey])
        tauESUpHist.Sumw2()
        tauESDownHist = ROOT.TH1D(tauESName.replace('<UPDOWN>', 'Down'), tauESName.replace('<UPDOWN>', 'Down'), len(binning[binKey])-1, binning[binKey])
        tauESDownHist.Sumw2()
        
        # loop over samples to be combined into this process
        for sample in procs[procKey]:
            histName = histoTemplate.replace('<CAT>', cat).replace('<SAMPLE>', sample)
            if args.verbose : print 'Loading histogram', histName
            inHist = inFiles["main"].Get(histName)
            inHistReb = inHist.Rebin(len(binning[binKey])-1, 'inHistReb', binning[binKey])
                        
            if sample == "QCD":
                if not args.noQCDScaling:
                    # scale QCD histogram up by 10% for mtt<50
                    QCDScaling()
                    
                qcdHist = inHistReb
            
            outHist.Add(inHistReb)
            
            # tau energy scale
            inHistTauESUp = inFiles["TauESUp"].Get(histName)
            inHistRebTauESUp = inHistTauESUp.Rebin(len(binning[binKey])-1, 'inHistRebTauESUp', binning[binKey])
            tauESUpHist.Add(inHistRebTauESUp)
            inHistTauESDown = inFiles["TauESDown"].Get(histName)
            inHistRebTauESDown = inHistTauESDown.Rebin(len(binning[binKey])-1, 'inHistRebTauESDown', binning[binKey])
            tauESDownHist.Add(inHistRebTauESDown)
        
        # save histogram into output file
        SaveHist(outHist)
        SaveHist(tauESUpHist)
        SaveHist(tauESDownHist)
        
        # ZL shape uncertainty
        if procKey == "ZL":
            for direction in ['Up', 'Down']:
                if len(procs[procKey]) > 1: print "WARNING: ZL should be only 1 sample"
                sample = procs[procKey][0]
                histName = inputZlShapeTemplate.replace('<CAT>', cat).replace('<SAMPLE>', sample).replace('<UPDOWN>', direction)
                if args.verbose : print 'Loading histogram', histName
                inHist = inFiles["main"].Get(histName)
                inHistReb = inHist.Rebin(len(binning[binKey])-1, 'inHistReb', binning[binKey])
                
                outName = zlShapeUncTemplate.replace('<UPDOWN>', direction)
                outHist = ROOT.TH1D(outName, outName, len(binning[binKey])-1, binning[binKey])
                outHist.Sumw2()
                outHist.Add(inHistReb)
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



