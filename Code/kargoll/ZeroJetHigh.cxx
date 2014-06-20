/*
 * ZeroJetHigh.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "ZeroJetHigh.h"

ZeroJetHigh::ZeroJetHigh(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class ZeroJetHigh" << std::endl;
	// run ZeroJetHigh category
	categoryFlag = "ZeroJetHigh";

	// run Categories using data-driven WJets BG
	wJetsBGSource = "Data";
}

ZeroJetHigh::~ZeroJetHigh() {
	  for(int j=0; j<Npassed.size(); j++){
	    std::cout << "ZeroJetHigh::~ZeroJetHigh Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "ZeroJetHigh::~ZeroJetHigh()" << std::endl;
}

