#!/usr/bin/env python

#execfile('loadPyRoot.py')

from ROOT import TFile, TCanvas, TGraphErrors, TH1D, TF1, gStyle, TAxis, TH1F
from tdrstyle import setTDRStyle
setTDRStyle(False)

gStyle.SetOptFit(True)

file_official = TFile("/home/home2/institut_3b/kargoll/TauAnalysis/H2TauLimits/CMSSW_7_1_5/src/auxiliaries/shapes/CERN/htt_mt.inputs-sm-8TeV.root","READ")
channels = ["0jet_high", "0jet_medium",
            "1jet_high_mediumhiggs", "1jet_high_lowhiggs", "1jet_medium",
            "vbf_loose", "vbf_tight"]

mLow   = 90
mHigh  = 160
mStep  = 5

for channel in channels:
    can = TCanvas()
    
    NPoints = (mHigh-mLow)/mStep + 1
    graph = TGraphErrors(NPoints)
    graph.SetTitle(channel)
    
    m = mLow
    while m <= mHigh:
        i = (m-mLow)/5
        x = m
        
        hist = file_official.Get("muTau_" + channel + "/ggH" + str(m));
        graph.SetPoint(i, x, hist.GetMean())
        graph.SetPointError(i, 0, hist.GetMeanError())
        
        m +=mStep
    
    graph.GetHistogram().GetXaxis().SetRangeUser(75, mHigh+mStep)
    graph.GetHistogram().GetYaxis().SetRangeUser(75, mHigh+mStep)
    graph.Draw("APL")
    graph.Fit("pol1", "", "", mLow, mHigh)
    
    hist = file_official.Get("muTau_" + channel + "/ZTT");
    dy_graph = TGraphErrors(1)
    dy_graph.SetPoint(0, 91, hist.GetMean())
    dy_graph.SetPointError(0, 0, hist.GetMeanError())
    dy_graph.SetMarkerColor(3)
    dy_graph.Draw("P")
    
    line = TF1("line", "x", 0, mHigh+40)
    line.SetLineColor(16)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.Draw("same")
    
    graph.Draw("PL")

    can.Print("svfMReg_" + channel + ".pdf")