/*
 * VBFLoose.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "VBFLoose.h"

VBFLoose::VBFLoose(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class VBFLoose" << std::endl;
	// run VBFLoose category
	categoryFlag = "VBFLoose";

}

VBFLoose::~VBFLoose() {
	  for(int j=0; j<Npassed.size(); j++){
	    std::cout << "VBFLoose::~VBFLoose Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "VBFLoose::~VBFLoose()" << std::endl;
}

