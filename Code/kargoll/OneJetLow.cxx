/*
 * OneJetLow.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "OneJetLow.h"

OneJetLow::OneJetLow(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class OneJetLow" << std::endl;
	// run OneJetLow category
	categoryFlag = "OneJetLow";

	// run Categories using embedding
	useEmbedding = true;

	// run Categories using data-driven WJets BG
	wJetsBGSource = "Data";

	// run Categories using data-driven QCD BG
	qcdShapeFromData = true;
}

OneJetLow::~OneJetLow() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "OneJetLow::~OneJetLow Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "OneJetLow::~OneJetLow()" << std::endl;
}

