#!/usr/bin/env python

#execfile('loadPyRoot.py')

import sys
from ROOT import TFile, TCanvas, TH1D, TColor, TLegend, gStyle, TGraph
from tdrstyle import setTDRStyle
setTDRStyle(False)

gStyle.SetLegendBorderSize(0)
folder = "/net/scratch_cms/institut_3b/kargoll/workdirSVFitSmallMassSampleRangeAug_26_2015/Set_1/SVFitLikelihoods/"
colors = [TColor.GetColor(0,84,159), # dark blue
          TColor.GetColor(204,7,30), # red
          TColor.GetColor(87,171,39), # green
          TColor.GetColor(0,97,101), # petrol
          TColor.GetColor(246,168,0)] # orange

def drawLikelihoods(files, names, outputName, genMass = None, xMax = 300):
    can = TCanvas()
    legend = TLegend(.54,.72,.94,.85)
    hists = []
    
    if len(files) != len(names):
        sys.exit("Files and Names do not have equal length. Abort.")
    
    for i in range(len(files)):
        hist = files[i].Get("SVfitStandaloneAlgorithmMarkovChain_histogramMass")
        hist.SetName("Hist"+str(i))
        if hist.GetMaximum(): hist.Scale(1./hist.GetMaximum())
        hist.SetLineColor( colors[i] )

        if i == 0:
            hist.GetXaxis().SetRangeUser(0.,xMax)
            hist.GetXaxis().SetTitle("m_{#tau#tau} / GeV")
            hist.GetXaxis().SetTitleOffset(0.98)
            hist.GetYaxis().SetTitle("Likelihood / a.u.")
            hist.GetYaxis().SetTitleOffset(1.2)
            hist.SetTitle("")
            hist.Draw()
        else:
            hist.Draw("same")
            
        legend.AddEntry(hist.GetName(), names[i], "l")
            
    if genMass:
        graph = TGraph(1)
        graph.SetPoint(0, genMass, 1.03)
        graph.SetMarkerStyle(23)
        graph.SetMarkerSize(2)
        graph.SetMarkerColor(TColor.GetColor(97,33,88))
        graph.Draw("P")
        legend.AddEntry(graph, "gen. Mass", "p")
    
    legend.Draw()
    can.Print(outputName + ".pdf")
    
######     
Event = "TruePlus_Rotated"
genM =  90.4558
Files = [ TFile(folder+"SVFLikelihood_Mass82.963897705078125_8kP07.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass64.963993072509766_hMFSv.root","READ")] 
Names = [ "3-prong SVFit Rotated",
          "Default SVFit"]

drawLikelihoods(Files, Names, Event, genM)

######     
Event = "TrueMinus_Rotated"
genM =  90.4302
Files = [ TFile(folder+"SVFLikelihood_Mass170.00048065185547_2tK0L.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass94.475265502929688_KWrO7.root","READ")] 
Names = [ "3-prong SVFit Rotated",
          "Default SVFit"]

drawLikelihoods(Files, Names, Event, genM)

######     
Event = "TrueMinus_Rotated_2"
genM =  89.5857
Files = [ TFile(folder+"SVFLikelihood_Mass146.70893096923828_tjH22.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass97.278984069824219_e1IX5.root","READ")] 
Names = [ "3-prong SVFit Rotated",
          "Default SVFit"]

drawLikelihoods(Files, Names, Event, genM)

######     
Event = "TrueMinus_Ambiguity"
genM =  89.097
Files = [ TFile(folder+"SVFLikelihood_Mass341.85054016113281_4ogHF.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass121.29331588745117_ZQxwQ.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass91.262699127197266_Iwo7J.root","READ")] 
Names = [ "3-prong SVFit Plus Sol.",
          "3-prong SVFit Minus Sol.",
          "Default SVFit"]

drawLikelihoods(Files, Names, Event, genM, 400)

######     
Event = "TrueMinus_Ambiguity_2"
genM =  91.603
Files = [ TFile(folder+"SVFLikelihood_Mass205.82180023193359_8b7mt.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass99.432666778564453_KLHIG.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass94.247104644775391_lFrXL.root","READ")] 
Names = [ "3-prong SVFit Plus Sol.",
          "3-prong SVFit Minus Sol.",
          "Default SVFit"]

drawLikelihoods(Files, Names, Event, genM)

######     
Event = "TrueMinus_Ambiguity_3"
genM =  111.691
Files = [ TFile(folder+"SVFLikelihood_Mass133.78108978271484_YL8rq.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass0_DiZSk.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass112.15602111816406_iNJj8.root","READ")] 
Names = [ "3-prong SVFit Plus Sol.",
          "3-prong SVFit Minus Sol.",
          "Default SVFit"]

drawLikelihoods(Files, Names, Event, genM)

######     
Event = "TrueMinus_Ambiguity_4"
genM =   90.9831
Files = [ TFile(folder+"SVFLikelihood_Mass94.164077758789062_E2fMt.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass89.527969360351562_VXa2z.root","READ"), 
          TFile(folder+"SVFLikelihood_Mass76.925491333007812_Cd5WK.root","READ")] 
Names = [ "3-prong SVFit Plus Sol.",
          "3-prong SVFit Minus Sol.",
          "Default SVFit"]

drawLikelihoods(Files, Names, Event, genM)

sys.exit(0)
