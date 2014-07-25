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

	if(wJetsBGSource != "MC"){
		std::cout << "Please set wJetsBGSource = \"MC\" to obtain background yields. Abort...";
		return;
	}

	// calculate W+Jet yield for categories
	unsigned histo;
	const unsigned nCat = 8;
	TString n[nCat] = {"0-Jet Low","0-Jet High","1-Jet Low","1-Jet High","1-Jet Boost","VBF Loose","VBF Tight","Inclusive"};
	std::vector<TString> catNames(n,n+nCat);
	std::vector<double> catEPSignal(nCat,0.0);
	std::vector<double> catEPSideband(nCat,0.0);
	std::vector<double> catWJetMCPrediction(nCat,0.0);
	std::vector<double> catWJetRelaxedMCPrediction(nCat,0.0);
	std::vector<double> catEPFactor(nCat,-9.9);
	std::vector<double> catSBData(nCat,-9.9);
	std::vector<double> catSBBackgrounds(nCat,0.0);
	std::vector<double> catWJetsInSB(nCat,-9.9);
	std::vector<double> catWJetsYield(nCat,-9.9);
	std::vector<double> catWJetsMCRatio(nCat,-9.9);
	std::vector<double> catWJetsRelaxedMCRatio(nCat,-9.9);

	int lowBin = Cat0JetLowMt.at(0).FindFixBin(0.0);
	int highBin = Cat0JetLowMt.at(0).FindFixBin(30.0) - 1;
	//printf("lowBin = %i, highBin = %i \n", lowBin, highBin);

	// extrapolation factor from MC
	// use unscaled MC events for this, thus do it before Selection::Finish() is called
	for (unsigned id = 20; id < 24; id++){ //only for WJets processes
		if (HConfig.GetHisto(false,id,histo)){
			catEPSignal.at(0) += Cat0JetLowMtExtrapolation.at(histo).GetBinContent(1);
			catEPSideband.at(0) += Cat0JetLowMtExtrapolation.at(histo).GetBinContent(2);

			catEPSignal.at(1) += Cat0JetHighMtExtrapolation.at(histo).GetBinContent(1);
			catEPSideband.at(1)+= Cat0JetHighMtExtrapolation.at(histo).GetBinContent(2);

			catEPSignal.at(2) += Cat1JetLowMtExtrapolation.at(histo).GetBinContent(1);
			catEPSideband.at(2) += Cat1JetLowMtExtrapolation.at(histo).GetBinContent(2);

			catEPSignal.at(3)+= Cat1JetHighMtExtrapolation.at(histo).GetBinContent(1);
			catEPSideband.at(3) += Cat1JetHighMtExtrapolation.at(histo).GetBinContent(2);

			catEPSignal.at(4) += Cat1JetBoostMtExtrapolation.at(histo).GetBinContent(1);
			catEPSideband.at(4) += Cat1JetBoostMtExtrapolation.at(histo).GetBinContent(2);

			catEPSignal.at(5) += CatVBFLooseRelaxMtExtrapolation.at(histo).GetBinContent(1);
			catEPSideband.at(5) += CatVBFLooseRelaxMtExtrapolation.at(histo).GetBinContent(2);

			catEPSignal.at(6) += CatVBFTightRelaxMtExtrapolation.at(histo).GetBinContent(1);
			catEPSideband.at(6) += CatVBFTightRelaxMtExtrapolation.at(histo).GetBinContent(2);

			catEPSignal.at(7)+= CatInclusiveMtExtrapolation.at(histo).GetBinContent(1);
			catEPSideband.at(7) += CatInclusiveMtExtrapolation.at(histo).GetBinContent(2);
		}
	}
	for (unsigned int icat = 0; icat < nCat; icat++){
		catEPFactor.at(icat) = catEPSignal.at(icat)/catEPSideband.at(icat);
	}


	// do plotting and scale histograms
	Selection::Finish();
	// all histograms below this are scaled to luminosity!

	// MC prediction of WJet in signal region (for cross-checking)
	for (unsigned id = 20; id < 24; id++){ //only for WJets processes
		if (HConfig.GetHisto(false,id,histo)){
			catWJetMCPrediction.at(0) += Cat0JetLowMt.at(histo).Integral(lowBin, highBin);
			catWJetMCPrediction.at(1) += Cat0JetHighMt.at(histo).Integral(lowBin, highBin);
			catWJetMCPrediction.at(2) += Cat1JetLowMt.at(histo).Integral(lowBin, highBin);
			catWJetMCPrediction.at(3) += Cat1JetHighMt.at(histo).Integral(lowBin, highBin);
			catWJetMCPrediction.at(4) += Cat1JetBoostMt.at(histo).Integral(lowBin, highBin);
			catWJetMCPrediction.at(5) += CatVBFLooseMt.at(histo).Integral(lowBin, highBin);
			catWJetRelaxedMCPrediction.at(5) += CatVBFLooseRelaxMt.at(histo).Integral(lowBin, highBin);
			catWJetMCPrediction.at(6) += CatVBFTightMt.at(histo).Integral(lowBin, highBin);
			catWJetRelaxedMCPrediction.at(6) += CatVBFTightRelaxMt.at(histo).Integral(lowBin, highBin);
			catWJetMCPrediction.at(7) += CatInclusiveMt.at(histo).Integral(lowBin, highBin);
		}
	}
	// mT sideband events from data
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
	for (unsigned id = 30; id < 80; id++){ //remove DY, diboson, top from MC
		if (id == 60) continue; // 60 = QCD
		if (HConfig.GetHisto(false,id,histo)){
			catSBBackgrounds.at(0) += Cat0JetLowMtSideband.at(histo).Integral();
			catSBBackgrounds.at(1) += Cat0JetHighMtSideband.at(histo).Integral();
			catSBBackgrounds.at(2) += Cat1JetLowMtSideband.at(histo).Integral();
			catSBBackgrounds.at(3) += Cat1JetHighMtSideband.at(histo).Integral();
			catSBBackgrounds.at(4) += Cat1JetBoostMtSideband.at(histo).Integral();
			catSBBackgrounds.at(5) += CatVBFLooseMtSideband.at(histo).Integral();
			catSBBackgrounds.at(6) += CatVBFTightMtSideband.at(histo).Integral();
			catSBBackgrounds.at(7) += CatInclusiveMtSideband.at(histo).Integral();
		}
	}
	for (unsigned int icat = 0; icat < nCat; icat++){
		catWJetsInSB.at(icat) = catSBData.at(icat) - catSBBackgrounds.at(icat);
		catWJetsYield.at(icat) = catWJetsInSB.at(icat) * catEPFactor.at(icat);
		catWJetsMCRatio.at(icat) = catWJetsYield.at(icat) / catWJetMCPrediction.at(icat);
		catWJetsRelaxedMCRatio.at(icat) = catWJetsYield.at(icat) / catWJetRelaxedMCPrediction.at(icat);
	}

	// print MC scales
	std::map<int, double> scales;
	for (unsigned id = 2; id < 100; id++){
			if (HConfig.GetHisto(false,id,histo)){
				double scale = Lumi * HConfig.GetCrossSection(id) / Npassed.at(histo).GetBinContent(0);
				scales.insert(std::pair<int, double>(id, scale));
				printf("ID = %2i has scale of %4f \n", id, scale);
				//printf("ID = %2i has scale of %4f (%4f * %4f / %4f) \n", id, scale, Lumi, HConfig.GetCrossSection(id), Npassed.at(histo).GetBinContent(0));
			}
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
		if (icat == 5 || icat == 6) continue;
		printf(format,catNames.at(icat).Data(), catWJetsYield.at(icat), catWJetMCPrediction.at(icat), catWJetsMCRatio.at(icat));
	}
	printf("%12s  %14s <-> %14s || %14s    %14s %14s\n","Category","WJets Yield", "MC Pred.", "Data/MC", "Rel. MC Pred.", "Data/(Rel. MC)");
	format = "%12s  %14.1f <-> %14.1f || %14.6f    %14.1f %14.6f\n";
	for (unsigned int icat = 5; icat <= 6; icat++){
		printf(format,catNames.at(icat).Data(), catWJetsYield.at(icat), catWJetMCPrediction.at(icat), catWJetsMCRatio.at(icat), catWJetRelaxedMCPrediction.at(icat), catWJetsRelaxedMCRatio.at(icat));
	}
	std::cout << "  ##########################################################" << std::endl;
	printf("Please copy the following numbers in order to use the data driven WJets yield:\n");
	for (unsigned int icat = 0; icat < nCat; icat++){
		if (icat != 5 && icat != 6) printf("%12s : %14.8f\n", catNames.at(icat).Data(), catWJetsMCRatio.at(icat));
		else printf("%12s : %14.8f\n", catNames.at(icat).Data(), catWJetsYield.at(icat));
	}
}
