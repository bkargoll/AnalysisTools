/*
 * Category.cxx
 *
 *  Created on: Apr 20, 2015
 *      Author: kargoll
 */

#include "Category.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

Category::Category(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	// run Categories using embedding
	useEmbedding = true;

	// run Categories using data-driven WJets BG
	wJetsBGSource = "Data";

	// run Categories using data-driven QCD BG
	qcdShapeFromData = true;
}

Category::~Category() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

void Category::Configure() {
	Logger(Logger::Verbose) << "Start." << std::endl;

	// configure base class
	HToTaumuTauh::Setup();
	// configure categorization
	categoryConfiguration();

	// make everything known to the framework
	Selection::ConfigureHistograms();
	HConfig.GetHistoInfo(types, CrossSectionandAcceptance, legend, colour);
}

void Category::doEvent() {
	Logger(Logger::Verbose) << ">>>>>>>>>>>>>>>>" << std::endl;
	Logger(Logger::Verbose) << "\tCategory: " << categoryFlag << std::endl;

	// setup event
	HToTaumuTauh::doEventSetup();

	// run selection part of analysis
	HToTaumuTauh::doSelection(false);
	catPassed = categorySelection();

	// execute AnalysisCuts
	status = AnalysisCuts(t,w,wobs);

	// do plotting part of analysis
	HToTaumuTauh::doPlotting();
	categoryPlotting();
}

void Category::Store_ExtraDist() {
	Logger(Logger::Verbose) << std::endl;

	// store histograms from base class
	HToTaumuTauh::Store_ExtraDist();

	// store histograms from this category
	categoryExtradist();
}
