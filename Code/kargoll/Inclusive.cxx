/*
 * Inclusive.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "Inclusive.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

Inclusive::Inclusive(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class Inclusive" << std::endl;
	// run Inclusive category
	categoryFlag = "Inclusive";

	// run Categories using embedding
	useEmbedding = true;

	// run Categories using data-driven WJets BG
	wJetsBGSource = "Data";

	// run Categories using data-driven QCD BG
	qcdShapeFromData = true;
}

Inclusive::~Inclusive() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

