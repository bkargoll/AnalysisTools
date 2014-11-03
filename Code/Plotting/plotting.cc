#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"

#include <vector>

#include "tdrstyle.C"
#include "CMS_lumi.C"

#include "userConfig.h"
#include "plottingHelpers.h"

void plotting(){
	if(verbose) std::cout << "--> plotting()" << std::endl;
	gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle();
	gROOT->LoadMacro("CMS_lumi.C");
	writeExtraText = true;
	gStyle->SetOptStat(0);
		
	// define samples to use
	sample s_ggHiggs	("H_{gg}"				, 632, "MC_ggFHTauTauM125");
	sample s_vbfHiggs	("H_{VBF}"				, 633, "MC_VBFHTauTauM125");
	sample s_apHiggs	("H_{AP}"				, 634, "MC_WHZHTTHTauTauM125");
	sample s_qcd		("QCD"					, 606, "MC_QCD");
	sample s_ztt		("Z#rightarrow#tau#tau"	, 400, "MC_DY_tautau");
	sample s_zll		("Z#rightarrowll"		, 800, "MC_DY_ll");
	sample s_ww			("WW"					, 417, "MC_WW_2l2nu");
	sample s_wz2l2q		("WZ#rightarrowllqq"	, 418, "MC_WZ_2l2q");
	sample s_wz3l1nu	("WZ#rightarrowlll#nu"	, 419, "MC_WZ_3l1nu");
	sample s_zz4l		("ZZ#rightarrowllll"	, 420, "MC_ZZ_4l");
	sample s_zz2l2nu	("ZZ#rightarrowll#nu#nu", 416, "MC_ZZ_2l2nu");
	sample s_zz2l2q		("ZZ#rightarrowllqq"	, 415, "MC_ZZ_2l2q");
	sample s_ttbar		("t#bar{t}"				, 600, "MC_ttbar");
	sample s_tw			("tW"					, 603, "MC_tw");
	sample s_tbarW		("#bar{t}W"				, 591, "MC_tbarw");
	sample s_Wlnu		("W#rightarrowl#nu"		, 876, "MC_W_lnu");
	sample s_Wtaunu		("W#rightarrow#tau#nu"	, 874, "MC_W_taunu");
	
	// *** examples for sample creation when scaling manually ***
	// std::map<int, int> MnEvents = readSkimSummary("mySkimSummary.log");
	// sample s_ggHiggs(11,"H_{gg}", 1.233664, lumi, 632, "MC_ggFHTauTauM125", MnEvents );
	// sample s_vbfHigg(12,"H_{VBF}", 0.0997296, lumi, 633, "MC_VBFHTauTauM125", MnEvents );
	// sample s_apHiggs(13,"H_{AP}", 0.0771792, lumi, 634, "MC_WHZHTTHTauTauM125", MnEvents );

	// combine samples for plotting
	sample s_higgs(s_ggHiggs, "H", 630);
	s_higgs += s_vbfHiggs;
	s_higgs += s_apHiggs;
	sample s_diboson(s_ww, "diboson", 416);
	s_diboson += s_wz2l2q;
	s_diboson += s_wz3l1nu;
	s_diboson += s_zz4l;
	s_diboson += s_zz2l2nu;
	s_diboson += s_zz2l2q;
	sample s_top(s_ttbar, "top", 600);
	s_top += s_tw;
	s_top += s_tbarW;

	// define which samples to plot and in which order
	std::vector<sample> samples;
	samples.push_back(s_higgs);
	// samples.push_back(s_qcd); // no qcd in background plots
	samples.push_back(s_ztt);
	samples.push_back(s_zll);
	samples.push_back(s_diboson);
	samples.push_back(s_top);
	samples.push_back(s_Wlnu);
	samples.push_back(s_Wtaunu);

	if (verbose){
		std::cout << "###>-- The following samples will be plotted:" << std::endl;
		for(unsigned i = 0; i<samples.size(); i++){
			printSample(samples.at(i), conf);
		}
	}

	/* THIS STUFF STILL NEEDS TO BE ADAPTED
	// todo: implement syst
	std::vector<double> syst;
	// lumi + xsec + eid + muid + pileup + trigger
	const int nsyst = 6;
	double qcd[nsyst] = {0.026,0.387,0.083,0.083,0.,0.083};
	double zz4l[nsyst] = {0.026,0.15,0.002,0.006,0.008,0.007};
	double zz2l2q[nsyst] = {0.026,0.15,0.0044,0.0042,0.004,0.007};
	double zz2l2nu[nsyst] = {0.026,0.15,0.0044,0.007,0.011,0.009};
	double wz3lnu[nsyst] = {0.026,0.15,0.0018,0.006,0.042,0.035};
	double wz2l2q[nsyst] = {0.026,0.15,0.007,0.006,0.015,0.014};
	double ww2l2nu[nsyst] = {0.026,0.15,0.0016,0.0055,0.013,0.015};
	double ttbar[nsyst] = {0.026,0.047,0.0014,0.0053,0.021,0.010};
	double tw[nsyst] = {0.026,0.09,0.0011,0.0055,0.294,0.325};
	double tbarw[nsyst] = {0.026,0.09,0.0011,0.0055,0.029,0.054};
	double dyll[nsyst] = {0.026,0.033,0.0015,0.0063,0.060,0.075};
	double dytt[nsyst] = {0.026,0.033,0.0016,0.0055,0.067,0.066};
	double dyemu[nsyst+1] = {0.026,0.055,0.0016,0.0055,0.011,0.011};
	syst.push_back(QuadraticSum(nsyst,qcd));
	syst.push_back(QuadraticSum(nsyst,zz4l));
	syst.push_back(QuadraticSum(nsyst,zz2l2q));
	syst.push_back(QuadraticSum(nsyst,zz2l2nu));
	syst.push_back(QuadraticSum(nsyst,wz3lnu));
	syst.push_back(QuadraticSum(nsyst,wz2l2q));
	syst.push_back(QuadraticSum(nsyst,ww2l2nu));
	syst.push_back(QuadraticSum(nsyst,ttbar));
	syst.push_back(QuadraticSum(nsyst,tw));
	syst.push_back(QuadraticSum(nsyst,tbarw));
	syst.push_back(QuadraticSum(nsyst,dyll));
	syst.push_back(QuadraticSum(nsyst,dytt));
	syst.push_back(QuadraticSum(nsyst+1,dyemu));
	*/
	
	// define which plots to draw
	std::vector<plotInfo> plots;
	plots.push_back( plotInfo("CatFired", "") );
	plots.push_back( plotInfo("MuPt", "GeV") );

	// test validity of structs
	testInputs(conf, samples, plots);
	
	// create plots
	if(testPlotting){
		plotInfo testPlot = plots.at(0);
		TH1D* datahist = getHisto(conf, testPlot.identifier+"Data", 1, 1);
		drawPlot(conf, testPlot, datahist, samples);
	}else{
		std::vector<TH1D*> datahists;
		for(unsigned p = 0; p < plots.size(); p++){
			datahists.push_back(getHisto(conf, plots.at(p).identifier+"Data", 1, 1));
		}
		for(unsigned p = 0; p < plots.size(); p++){
			drawPlot(conf, plots.at(p), datahists.at(p), samples);
		}
	}

	/*
	//TH1D* onejet_signal = getHisto("onejetMC_emu_DY",mcscale.at(mcscale.size()-1),0,infile,syst.at(syst.size()-1));
	//TH1D* onejet_dy = getHisto("onejetMC_tautau_DY",mcscale.at(mcscale.size()-2),2345,infile,syst.at(syst.size()-2));
	//onejet_signal->Scale(1./onejet_signal->Integral());
	//onejet_dy->Scale(1./onejet_dy->Integral());
	//TH1D* onejet_ratio = getDataMC(onejet_signal,onejet_dy);
	//drawPlot(onejet_signal, onejet_dy, onejet_ratio, "Custom MC", "Official MC", "Single Jet Pt", "GeV");

	TH1D* zmass_signal = getHisto("zmassMC_emu_DY",1,0,infile,0);
	TH1D* zmass_dy = getHisto("zmassMC_DY",1,2345,infile,0);
	zmass_signal->Scale(1./zmass_signal->Integral());
	zmass_dy->Scale(1./zmass_dy->Integral());
	TH1D* zmass_ratio = getDataMC(zmass_signal,zmass_dy);
	//drawPlot(zmass_signal,zmass_dy,zmass_ratio,"Custom MC","Official MC","Z mass","GeV");

	TH1D* zpt_signal = getHisto("zptMC_emu_DY",1,0,infile,0);
	TH1D* zpt_dy = getHisto("zptMC_tautau_DY",1,2345,infile,0);
	zpt_signal->Scale(1./zpt_signal->Integral());
	zpt_dy->Scale(1./zpt_dy->Integral());
	TH1D* zpt_ratio = getDataMC(zpt_signal,zpt_dy);
	//for(unsigned i=1;i<=zpt_ratio->GetNbinsX();i++){
	//	std::cout << "bin " << i << ": " << zpt_ratio->GetBinContent(i) << std::endl;
	//}
	//drawPlot(zpt_signal,zpt_dy,zpt_ratio,"Custom MC","Official MC","Z pt","GeV");

	TH1D* zeta_signal = getHisto("zetaMC_emu_DY",1,0,infile,0);
	TH1D* zeta_dy = getHisto("zetaMC_tautau_DY",1,2345,infile,0);
	zeta_signal->Scale(1./zeta_signal->Integral());
	zeta_dy->Scale(1./zeta_dy->Integral());
	TH1D* zeta_ratio = getDataMC(zeta_signal,zeta_dy);
	//drawPlot(zeta_signal,zeta_dy,zeta_ratio,"Custom MC","Official MC","Z eta","GeV");

	TH1D* njets_signal = getHisto("NJetsMC_emu_DY",mcscale.at(mcscale.size()-1),0,infile,syst.at(syst.size()-1));
	TH1D* njets_dy = getHisto("NJetsMC_tautau_DY",mcscale.at(mcscale.size()-2),2345,infile,syst.at(syst.size()-2));
	njets_signal->Scale(1./njets_signal->Integral());
	njets_dy->Scale(1./njets_dy->Integral());
	TH1D* njets_ratio = getDataMC(njets_signal,njets_dy);
	//drawPlot(njets_signal,njets_dy,njets_ratio,"Custom MC","Official MC","Number of jets","");
*/
}

