/*
 * Inclusive.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "Inclusive.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

Inclusive::Inclusive(TString Name_, TString id_):
	Category(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class Inclusive" << std::endl;
	// run Inclusive category
	categoryFlag = "Inclusive";
}

Inclusive::~Inclusive() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

bool Inclusive::categorySelection(){
	bool categoryPass = true;
	std::vector<float> value_Inclusive(HToTaumuTauh::NCuts,-10);
	std::vector<float> pass_Inclusive(HToTaumuTauh::NCuts,false);

	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("Inclusive",value_Inclusive, pass_Inclusive,HToTaumuTauh::CatCut1) && categoryPass;
	return categoryPass;
}
