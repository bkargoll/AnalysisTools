/*
 * OneJetHigh.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "OneJetHigh.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

OneJetHigh::OneJetHigh(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class OneJetHigh" << std::endl;
	// run OneJetHigh category
	categoryFlag = "OneJetHigh";

	// run Categories using embedding
	useEmbedding = true;

	// run Categories using data-driven WJets BG
	wJetsBGSource = "Data";

	// run Categories using data-driven QCD BG
	qcdShapeFromData = true;
}

OneJetHigh::~OneJetHigh() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

