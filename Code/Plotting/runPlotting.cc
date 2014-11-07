/////////////////////////////////////////////////
//
// plotting.cc: main macro to run in ROOT
// call like "root runPlotting.cc+" in your shell
//
/////////////////////////////////////////////////

#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TColor.h"

#include <vector>

#include "tdrstyle.C"
#include "CMS_lumi.C"

#include "userConfig.h"
#include "defineSamplesAndPlots.h"
#include "plottingHelpers.h"

void runPlotting(){
	if(verbose) std::cout << "--> plotting()" << std::endl;
	gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle();
	gROOT->LoadMacro("CMS_lumi.C");
	writeExtraText = true;
	gStyle->SetOptStat(0);

	// define which samples to use
	std::vector<sample> samples = defineSamples();

	// define which plots to draw
	std::vector<plotInfo> plots = definePlots();

	// todo: implement syst

	// test validity of structs
	testInputs(conf, samples, plots);
	
	// create plots
	if(testPlotting){
		plotInfo testPlot = plots.at(0);
		TH1D* datahist = getHisto(conf, testPlot.identifier+"Data", 1, 1);
		manipulateHisto(datahist, testPlot);
		drawPlot(conf, testPlot, datahist, samples);
	}else{
		std::vector<TH1D*> datahists;
		for(unsigned p = 0; p < plots.size(); p++){
			TH1D* datahist = getHisto(conf, plots.at(p).identifier+"Data", 1, 1);
			manipulateHisto(datahist, plots.at(p));
			datahists.push_back(datahist);
		}
		for(unsigned p = 0; p < plots.size(); p++){
			drawPlot(conf, plots.at(p), datahists.at(p), samples);
		}
	}

	// convert all created eps files into pdf
	std::cout << "Plots created. Convert all eps into pdf..." << std::endl;
	system("for i in `ls -1 *.eps`; do epstopdf $i; done");
}

