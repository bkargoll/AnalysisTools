/*
 * ZeroJetLow.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "ZeroJetLow.h"

ZeroJetLow::ZeroJetLow(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class ZeroJetLow" << std::endl;
	// run ZeroJetLow category
	categoryFlag = "ZeroJetLow";

}

ZeroJetLow::~ZeroJetLow() {
	  for(int j=0; j<Npassed.size(); j++){
	    std::cout << "ZeroJetLow::~ZeroJetLow Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "ZeroJetLow::~ZeroJetLow()" << std::endl;
}

