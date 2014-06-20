/*
 * HToTaumuTauhBackgrounds.cxx
 *
 *  Created on: Jun 20, 2014
 *      Author: kargoll
 */

#include "HToTaumuTauhBackgrounds.h"

HToTaumuTauhBackgrounds::HToTaumuTauhBackgrounds(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)

{
	std::cout << "Setting up the class HToTaumuTauhBackgrounds" << std::endl;
	// always run without category for background methods
	// the numbers will be produced for all categories individually
	categoryFlag = "NoCategory";

	// run Skim always using MC for WJets BG
	wJetsBGSource = "MC";
}

HToTaumuTauhBackgrounds::~HToTaumuTauhBackgrounds() {
	  for(int j=0; j<Npassed.size(); j++){
	    std::cout << "HToTaumuTauhBackgrounds::~HToTaumuTauhBackgrounds Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "HToTaumuTauhBackgrounds::~HToTaumuTauhBackgrounds()" << std::endl;
}

void HToTaumuTauhBackgrounds::Finish() {
	if (verbose) std::cout << "HToTaumuTauhBackgrounds::Finish()" << std::endl;

	Selection::Finish();

	// print MC scales
	unsigned histo;
	std::map<int, double> scales;
	for (unsigned id = 2; id < 100; id++){
			if (HConfig.GetHisto(false,id,histo)){
				double scale = Lumi * HConfig.GetCrossSection(id) / Npassed.at(histo).GetBinContent(0);
				scales.insert(std::pair<int, double>(id, scale));
				printf("ID = %2i has scale of %4f \n", id, scale);
				//printf("ID = %2i has scale of %4f (%4f * %4f / %4f) \n", id, scale, Lumi, HConfig.GetCrossSection(id), Npassed.at(histo).GetBinContent(0));
			}
	}
	if(wJetsBGSource != "MC") std::cout << "Please set wJetsBGSource = \"MC\" to obtain background yields.";
	else {
		// calculate W+Jet yield for categories
		const unsigned nCat = 8;
		TString n[nCat] = {"0-Jet Low","0-Jet High","1-Jet Low","1-Jet High","1-Jet Boost","VBF Loose","VBF Tight","Inclusive"};
		std::vector<TString> catNames(n,n+nCat);
		// extrapolation factor from MC
		std::vector<double> catEPSignal(nCat,0.0);
		std::vector<double> catEPSideband(nCat,0.0);
		std::vector<double> catWJetMCPrediction(nCat,0.0);
		int lowBin = Cat0JetLowMt.at(0).FindFixBin(0.0);
		int highBin = Cat0JetLowMt.at(0).FindFixBin(30.0) - 1;
		printf("lowBin = %i, highBin = %i \n", lowBin, highBin);
		for (unsigned id = 20; id < 24; id++){ //only for WJets processes
			if (HConfig.GetHisto(false,id,histo)){
				catEPSignal.at(0) += Cat0JetLowMtExtrapolation.at(histo).GetBinContent(1);
				catEPSideband.at(0) += Cat0JetLowMtExtrapolation.at(histo).GetBinContent(2);
				catWJetMCPrediction.at(0) += Cat0JetLowMt.at(histo).Integral(lowBin, highBin) * scales[id];

				catEPSignal.at(1) += Cat0JetHighMtExtrapolation.at(histo).GetBinContent(1);
				catEPSideband.at(1)+= Cat0JetHighMtExtrapolation.at(histo).GetBinContent(2);
				catWJetMCPrediction.at(1) += Cat0JetHighMt.at(histo).Integral(lowBin, highBin) * scales[id];

				catEPSignal.at(2) += Cat1JetLowMtExtrapolation.at(histo).GetBinContent(1);
				catEPSideband.at(2) += Cat1JetLowMtExtrapolation.at(histo).GetBinContent(2);
				catWJetMCPrediction.at(2) += Cat1JetLowMt.at(histo).Integral(lowBin, highBin) * scales[id];

				catEPSignal.at(3)+= Cat1JetHighMtExtrapolation.at(histo).GetBinContent(1);
				catEPSideband.at(3) += Cat1JetHighMtExtrapolation.at(histo).GetBinContent(2);
				catWJetMCPrediction.at(3) += Cat1JetHighMt.at(histo).Integral(lowBin, highBin) * scales[id];

				catEPSignal.at(4) += Cat1JetBoostMtExtrapolation.at(histo).GetBinContent(1);
				catEPSideband.at(4) += Cat1JetBoostMtExtrapolation.at(histo).GetBinContent(2);
				catWJetMCPrediction.at(4) += Cat1JetBoostMt.at(histo).Integral(lowBin, highBin) * scales[id];

				catEPSignal.at(5) += CatVBFLooseRelaxMtExtrapolation.at(histo).GetBinContent(1);
				catEPSideband.at(5) += CatVBFLooseRelaxMtExtrapolation.at(histo).GetBinContent(2);
				catWJetMCPrediction.at(5) += CatVBFLooseMt.at(histo).Integral(lowBin, highBin) * scales[id];

				catEPSignal.at(6) += CatVBFTightRelaxMtExtrapolation.at(histo).GetBinContent(1);
				catEPSideband.at(6) += CatVBFTightRelaxMtExtrapolation.at(histo).GetBinContent(2);
				catWJetMCPrediction.at(6) += CatVBFTightMt.at(histo).Integral(lowBin, highBin) * scales[id];

				catEPSignal.at(7)+= CatInclusiveMtExtrapolation.at(histo).GetBinContent(1);
				catEPSideband.at(7) += CatInclusiveMtExtrapolation.at(histo).GetBinContent(2);
				catWJetMCPrediction.at(7) += CatInclusiveMt.at(histo).Integral(lowBin, highBin) * scales[id];
			}
		}
		std::vector<double> catEPFactor(nCat,0.0);
		for (unsigned int icat = 0; icat < nCat; icat++){
			catEPFactor.at(icat) = catEPSignal.at(icat)/catEPSideband.at(icat);
		}
		// mT sideband events from data
		std::vector<double> catSBData(nCat,0.0);
		if (HConfig.GetHisto(true,1,histo)){
			catSBData.at(0) = Cat0JetLowMtSideband.at(histo).Integral();
			catSBData.at(1) = Cat0JetHighMtSideband.at(histo).Integral();
			catSBData.at(2) = Cat1JetLowMtSideband.at(histo).Integral();
			catSBData.at(3) = Cat1JetHighMtSideband.at(histo).Integral();
			catSBData.at(4) = Cat1JetBoostMtSideband.at(histo).Integral();
			catSBData.at(5) = CatVBFLooseMtSideband.at(histo).Integral();
			catSBData.at(6) = CatVBFTightMtSideband.at(histo).Integral();
			catSBData.at(7) = CatInclusiveMtSideband.at(histo).Integral();
		}
		std::vector<double> catSBBackgrounds(nCat,0.0);
		for (unsigned id = 30; id < 80; id++){ //remove DY, diboson, top from MC
			if (id == 60) continue; // 60 = QCD
			if (HConfig.GetHisto(false,id,histo)){
				catSBBackgrounds.at(0) += Cat0JetLowMtSideband.at(histo).Integral() * scales[id];
				catSBBackgrounds.at(1) += Cat0JetHighMtSideband.at(histo).Integral() * scales[id];
				catSBBackgrounds.at(2) += Cat1JetLowMtSideband.at(histo).Integral() * scales[id];
				catSBBackgrounds.at(3) += Cat1JetHighMtSideband.at(histo).Integral() * scales[id];
				catSBBackgrounds.at(4) += Cat1JetBoostMtSideband.at(histo).Integral() * scales[id];
				catSBBackgrounds.at(5) += CatVBFLooseMtSideband.at(histo).Integral() * scales[id];
				catSBBackgrounds.at(6) += CatVBFTightMtSideband.at(histo).Integral() * scales[id];
				catSBBackgrounds.at(7) += CatInclusiveMtSideband.at(histo).Integral() * scales[id];
			}
		}
		std::vector<double> catWJetsInSB(nCat,0.0);
		std::vector<double> catWJetsYield(nCat,0.0);
		std::vector<double> catWJetsMCRatio(nCat,0.0);
		for (unsigned int icat = 0; icat < nCat; icat++){
			catWJetsInSB.at(icat) = catSBData.at(icat) - catSBBackgrounds.at(icat);
			catWJetsYield.at(icat) = catWJetsInSB.at(icat) * catEPFactor.at(icat);
			catWJetsMCRatio.at(icat) = catWJetsYield.at(icat) / catWJetMCPrediction.at(icat);
		}


		// print results
		std::cout << "  ##########################################################" << std::endl;
		std::cout << "  ############# W+Jets MC extrapolation factor #############" << std::endl;
		printf("%12s  %13s : %13s = %12s \n","Category","Signal Region", "Sideband", "Extr. factor");
		const char* format = "%12s  %13.1f : %13.1f = %12f \n";
		for (unsigned int icat = 0; icat < nCat; icat++){
			printf(format,catNames.at(icat).Data(), catEPSignal.at(icat), catEPSideband.at(icat), catEPFactor.at(icat));
		}
		std::cout << "  ############# W+Jets Events in Sideband ##################" << std::endl;
		printf("%12s  %13s - %13s = %14s \n","Category","Nevts Data SB", "Nevts MC SB", "Nevts WJets SB");
		format = "%12s  %13.3f - %13.3f = %14.3f \n";
		for (unsigned int icat = 0; icat < nCat; icat++){
			printf(format,catNames.at(icat).Data(), catSBData.at(icat), catSBBackgrounds.at(icat), catWJetsInSB.at(icat));
		}
		std::cout << "  ############# W+Jets Yield ###############################" << std::endl;
		printf("%12s  %14s * %14s = %14s \n","Category","Nevts WJets SB", "Extr. factor", "WJets Yield");
		format = "%12s  %14.3f * %14f = %14.1f\n";
		for (unsigned int icat = 0; icat < nCat; icat++){
			printf(format,catNames.at(icat).Data(), catWJetsInSB.at(icat), catEPFactor.at(icat), catWJetsYield.at(icat));
		}
		std::cout << "  ############# W+Jets MC Comparison #######################" << std::endl;
		printf("%12s  %14s <-> %14s || %14s\n","Category","WJets Yield", "MC Pred.", "Data/MC");
		format = "%12s  %14.1f <-> %14.1f || %14.6f\n";
		for (unsigned int icat = 0; icat < nCat; icat++){
			printf(format,catNames.at(icat).Data(), catWJetsYield.at(icat), catWJetMCPrediction.at(icat), catWJetsMCRatio.at(icat));
		}
		std::cout << "  ##########################################################" << std::endl;
	}
}