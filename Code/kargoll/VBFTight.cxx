/*
 * VBFTight.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "VBFTight.h"

VBFTight::VBFTight(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class VBFTight" << std::endl;
	// run VBFTight category
	categoryFlag = "VBFTight";

}

VBFTight::~VBFTight() {
	  for(int j=0; j<Npassed.size(); j++){
	    std::cout << "VBFTight::~VBFTight Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "VBFTight::~VBFTight()" << std::endl;
}

