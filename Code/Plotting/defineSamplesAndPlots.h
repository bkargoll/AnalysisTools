#ifndef defineSamplesAndPlots_h
#define defineSamplesAndPlots_h

/////////////////////////////////////////////////
//
// defineSamplesAndPlots:
// declare which plots you want to draw using which samples
// this cannot go into userConfig.h, as it would cause circular dependencies
//
/////////////////////////////////////////////////

#include "userConfig.h"
#include "sample.h"
#include "plotInfo.h"
#include "plottingHelpers.h"

// define samples to use
std::vector<sample> defineSamples(){
	sample s_ggHiggs	("H_{gg}"				, 632, "MC_ggFHTauTauM125");
	sample s_vbfHiggs	("H_{VBF}"				, 633, "MC_VBFHTauTauM125");
	sample s_apHiggs	("H_{AP}"				, 634, "MC_WHZHTTHTauTauM125");
	sample s_qcd		("QCD"					, col_htt_qcd, "QCD");
	sample s_ztt		("Z#rightarrow#tau#tau"	, col_htt_Ztt, "MC_DY_tautau");
	sample s_zll		("Z#rightarrowll"		, col_rwth_lightblue, "MC_DY_ll");
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
	sample s_ewk(s_ww, "electroweak", col_htt_W);
	s_ewk += s_wz2l2q;
	s_ewk += s_wz3l1nu;
	s_ewk += s_zz4l;
	s_ewk += s_zz2l2nu;
	s_ewk += s_zz2l2q;
	s_ewk += s_Wlnu;
	s_ewk += s_Wtaunu;
	sample s_top(s_ttbar, "top", col_htt_tt);
	s_top += s_tw;
	s_top += s_tbarW;

	// define which samples to plot and in which order
	std::vector<sample> samples;
	samples.push_back(s_qcd);
	samples.push_back(s_top);
	samples.push_back(s_ewk);
	samples.push_back(s_zll);
	samples.push_back(s_ztt);
	samples.push_back(s_higgs);

	if (verbose){
		std::cout << "###>-- The following samples will be plotted:" << std::endl;
		for(unsigned i = 0; i<samples.size(); i++){
			printSample(samples.at(i), conf);
		}
	}

	return samples;
}


// define which plots to draw
std::vector<plotInfo> definePlots(){
	std::vector<plotInfo> plots;
	plots.push_back( plotInfo("CatFired", "") );
	plots.push_back( plotInfo("MuPt", "GeV") );
	plots.push_back( plotInfo("Cut_10_Nminus1_MT_", "GeV") );
	plots.push_back( plotInfo("MtAfterOppCharge", "GeV", 2) );
	plots.push_back( plotInfo("MetPt", "GeV", 1, 0., 80.));

	return plots;
}


#endif
