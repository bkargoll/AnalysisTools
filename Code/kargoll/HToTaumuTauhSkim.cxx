/*
 * HToTaumuTauhSkim.cxx
 *
 *  Created on: May 14, 2014
 *      Author: kargoll
 */

#include "HToTaumuTauhSkim.h"

HToTaumuTauhSkim::HToTaumuTauhSkim(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_),
	cSkim_Mu_relIso(0.5),
	cSkim_Tau_rawIso(10.0)
{
	std::cout << "Setting up the class HToTaumuTauhSkim" << std::endl;
	// always run without category for skimming
	categoryFlag = "NoCategory";

	// run Skim always using MC for WJets BG
	wJetsBGSource = "MC";

	// check if skimming cuts are looser than main analysis cuts
	if (cSkim_Mu_relIso < cMu_relIso)
		std::cout << "WARNING: Mu_relIso cut in skim tighter than in main analysis! Using main analysis value." << std::endl;
	else
		cMu_relIso = cSkim_Mu_relIso;
	if (cSkim_Tau_rawIso < cTau_rawIso)
		std::cout << "WARNING: Tau_rawIso cut in skim tighter than in main analysis! Using main analysis value." << std::endl;
	else
		cTau_rawIso = cSkim_Tau_rawIso;
}

HToTaumuTauhSkim::~HToTaumuTauhSkim() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "HToTaumuTauhSkim::~HToTaumuTauhSkim Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "HToTaumuTauhSkim::~HToTaumuTauhSkim()" << std::endl;
}

void  HToTaumuTauhSkim::Configure(){
	HToTaumuTauh::Configure();

	// disable cuts which should not be applied on skimming level
	cut.at(OppCharge) = 999; // set to 999 to disable opp. charge cut
	title.at(OppCharge) += " DISABLED";
	cut.at(MT) = 999; // set to 999 to disable mt cut
	title.at(MT) += " DISABLED";

}

void HToTaumuTauhSkim::doEvent() {
}

void HToTaumuTauhSkim::Store_ExtraDist() {
}
