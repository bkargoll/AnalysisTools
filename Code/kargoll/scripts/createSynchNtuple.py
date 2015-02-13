#!/usr/bin/env python

execfile('loadPyRoot.py')
import argparse
import numpy

# allow accessing additional data types
ROOT.gROOT.ProcessLine('.L ../../macros/Loader.C+')

parser = argparse.ArgumentParser(description='Create a HTauTau Mu+Tau synchronization Ntuple from an Aachen3B TauNtuple.')
parser.add_argument('inputFile', help="location of the input root file")
parser.add_argument('outputFile', help="location to save the output root file")
parser.add_argument('--verbose', action='store_true',help="Print more stuff.")
args = parser.parse_args()

# create dictionary to hold mapping between sync Ntuple and TauNtuple
# list taken from https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/HiggsTauTau/src/HTTSync.cc#L48-L201
names = {"run":"inTree.Event_RunNumber",
         "lumi":"inTree.Event_luminosityBlock",
         "evt":"inTree.Event_EventNumber",
         # Event Variables
         "npv":"inTree.Vtx_ndof.size()",
         "npu":"inTree.PileupInfo_NumInteractions_n0",
         "rho":"inTree.RhoIsolationAllInputTags",
         # Event Weights
         "mcweight":"",
         "puweight":"inTree.EvtWeight3D",
         "trigweight_1":"",
         "trigweight_2":"",
         "idweight_1":"",
         "idweight_2":"",
         "isoweight_1":"",
         "isoweight_2":"",
         "fakeweight":"",
         "effweight":"",
         "weight":"",
         "embeddedWeight":"",
         "signalWeight":"",
         # SV Fit variables
         "mvis":"(p4_1+p4_2).M()",
         "m_sv":"",
         "pt_sv":"",
         "eta_sv":"",
         "phi_sv":"",
         "m_sv_Up":"",
         "m_sv_Down":"",
         # First lepton : muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
         "pt_1":"p4_1.Pt()",
         "phi_1":"p4_1.Phi()",
         "eta_1":"p4_1.Eta()",
         "m_1":"p4_1.M()",
         "q_1":"inTree.Muon_charge.at(0)",
         "iso_1":"",
         "mva_1":"",
         "d0_1":"",
         "dZ_1":"",
         "passid_1":"",
         "passiso_1":"",
         "mt_1":"",
         # Second lepton : hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT) Tau for Tau Tau
         "pt_2":"p4_2.Pt()",
         "phi_2":"p4_2.Phi()",
         "eta_2":"p4_2.Eta()",
         "m_2":"p4_2.M()",
         "q_2":"inTree.PFTau_Charge.at(0)",
         "iso_2":"",
         "d0_2":"",
         "dZ_2":"",
         "pt_tt":"(p4_1+p4_2).Pt()",
         "byCombinedIsolationDeltaBetaCorrRaw3Hits_2":"inTree.PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits.at(0)",
         "againstElectronMVA3raw_2":"",
         "byIsolationMVA2raw_2":"",
         "againstMuonLoose2_2":"inTree.PFTau_isHPSAgainstMuonLoose2.at(0)",
         "againstMuonMedium2_2":"inTree.PFTau_isHPSAgainstMuonMedium2.at(0)",
         "againstMuonTight2_2":"inTree.PFTau_isHPSAgainstMuonTight2.at(0)",
         "mva_2":"",
         "passid_2":"",
         "passiso_2":"",
         "mt_2":"",
         # Met related variables
         "met":"inTree.MET_Uncorr_et",
         "metphi":"inTree.MET_Uncorr_phi",
         "l1met":"",
         "l1metphi":"",
         "l1metcorr":"",
         "calomet":"",
         "calometphi":"",
         "calometcorr":"",
         "calometphicorr":"",
         "mvamet":"inTree.MET_CorrMVA_et",
         "mvametphi":"inTree.MET_CorrMVA_phi",
         "pzetavis":"",
         "pzetamiss":"",
         # met covariance matrices
         "metcov00":"inTree.MET_Uncorr_significance_xx",
         "metcov01":"inTree.MET_Uncorr_significance_xy",
         "metcov10":"inTree.MET_Uncorr_significance_xy",
         "metcov11":"inTree.MET_Uncorr_significance_yy",
         # mva met covariance matrices
         "mvacov00":"inTree.MET_CorrMVA_significance_xx",
         "mvacov01":"inTree.MET_CorrMVA_significance_xy",
         "mvacov10":"inTree.MET_CorrMVA_significance_xy",
         "mvacov11":"inTree.MET_CorrMVA_significance_yy",
         # First jet: leading jet after applying Jet energy corrections (excluding hadronic Tau)
         "jpt_1":"jetP4_1.Pt()",
         "jeta_1":"jetP4_1.Eta()",
         "jphi_1":"jetP4_1.Phi()",
         "jptraw_1":"",
         "jptunc_1":"",
         "jmva_1":"inTree.PFJet_PUJetID_discr.at(0)",
         "jlrm_1":"",
         "jctm_1":"",
         "jpass_1":"inTree.PFJet_PUJetID_looseWP.at(0)",
         # Second Jet : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
         "jpt_2":"jetP4_2.Pt()",
         "jeta_2":"jetP4_2.Eta()",
         "jphi_2":"jetP4_2.Phi()",
         "jptraw_2":"",
         "jptunc_2":"",
         "jmva_2":"inTree.PFJet_PUJetID_discr.at(1)",
         "jlrm_2":"",
         "jctm_2":"",
         "jpass_2":"inTree.PFJet_PUJetID_looseWP.at(1)",
         # B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
         "bpt":"",
         "beta":"",
         "bphi":"",
         # Di Jet kinematic variables for VBF selection ==> Two leading pT Jets
         "mjj":"",
         "jdeta":"",
         "njetingap":"",
         "mva":"",
         # variables that go into the VBF MVA
         "jdphi":"",
         "dijetpt":"",
         "dijetphi":"",
         "hdijetphi":"",
         "visjeteta":"",
         "ptvis":"",
         # number of btags passing btag id (pt > 20)
         "nbtag":"",
         # number of jets passing jet id ( pt > 30 )
         "njets":"",
         "njetspt20":"",
         # mva output for e+mu channel
         "mva_gf":"",
         "mva_vbf":""
         }

# create empty dictionary to hold mapping between sync Ntuple name and variable
vars = {}

if (args.verbose):
    print "Input file: ", args.inputFile
    print "Output file: ", args.outputFile

inFile = ROOT.TFile(args.inputFile, "READ")
outFile = ROOT.TFile(args.outputFile, "RECREATE")

inTree = inFile.Get("t")
outTree = ROOT.TTree("syncTree","syncTree")

# create variables to copy tree content to
for key in names.keys():
    if names[key] == "":
        if (args.verbose):
            print "Variable", key, "will be ignored"
        continue
    addVar = {key:numpy.zeros(1, dtype='float')}
    vars.update(addVar)
    outTree.Branch(key,vars[key],key+'/D')
    if (args.verbose):
        print "added variable", key, "to hold contents of", names[key]

print "Start Event Loop"
for event in inTree:
    # calculate 4vectors to be used for individual objects
    p4_1 = ROOT.TLorentzVector(inTree.Muon_p4.at(0).at(1),inTree.Muon_p4.at(0).at(2),inTree.Muon_p4.at(0).at(3),inTree.Muon_p4.at(0).at(0))
    p4_2 = ROOT.TLorentzVector(inTree.PFTau_p4.at(0).at(1),inTree.PFTau_p4.at(0).at(2),inTree.PFTau_p4.at(0).at(3),inTree.PFTau_p4.at(0).at(0))
    if inTree.PFJet_p4.size() > 0:
        jetP4_1  = ROOT.TLorentzVector(inTree.PFJet_p4.at(0).at(1),inTree.PFJet_p4.at(0).at(2),inTree.PFJet_p4.at(0).at(3),inTree.PFJet_p4.at(0).at(0))
        jet1 = True
    else:
        jetP4_1 = ROOT.TLorentzVector()
        jet2 = False
        
    if inTree.PFJet_p4.size() > 1:
        jetP4_2  = ROOT.TLorentzVector(inTree.PFJet_p4.at(1).at(1),inTree.PFJet_p4.at(1).at(2),inTree.PFJet_p4.at(1).at(3),inTree.PFJet_p4.at(1).at(0))
        jet2 = True
    else:
        jetP4_2 = ROOT.TLorentzVector()
        jet2 = False
    
    # fill branches
    for key in names.keys():
        if names[key] == "":
            continue
        if ("PFJet" in names[key]):
            if ("at(0)" in names[key]) and not jet1:
                continue
            if ("at(1)" in names[key]) and not jet2:
                continue
        #print "vars[", key, "][0] = eval(", names[key], ")"
        vars[key][0] = eval(names[key])
    outTree.Fill()

outFile.Write()
outFile.Close()

print "The file \"", args.outputFile, "\" was created to hold your sync Ntuple."
