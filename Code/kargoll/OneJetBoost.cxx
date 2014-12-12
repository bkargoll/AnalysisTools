/*
 * OneJetBoost.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "OneJetBoost.h"

OneJetBoost::OneJetBoost(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class OneJetBoost" << std::endl;
	// run OneJetBoost category
	categoryFlag = "OneJetBoost";

	// run Categories using embedding
	useEmbedding = true;

	// run Categories using data-driven WJets BG
	wJetsBGSource = "Data";

	// run Categories using data-driven QCD BG
	qcdShapeFromData = true;
}

OneJetBoost::~OneJetBoost() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "OneJetBoost::~OneJetBoost Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "OneJetBoost::~OneJetBoost()" << std::endl;
}

