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

	//////////// WJets ////////////
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

	// do the stuff in SS region for QCD method
	std::vector<double> catEPSignalSS(nCat,0.0);
	std::vector<double> catEPSidebandSS(nCat,0.0);
	std::vector<double> catEPFactorSS(nCat,-9.9);
	std::vector<double> catSBDataSS(nCat,-9.9);
	std::vector<double> catSBBackgroundsSS(nCat,0.0);
	std::vector<double> catWJetsInSBSS(nCat,-9.9);
	std::vector<double> catWJetsYieldSS(nCat,-9.9);

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

			catEPSignalSS.at(0) += Cat0JetLowMtExtrapolationSS.at(histo).GetBinContent(1);
			catEPSidebandSS.at(0) += Cat0JetLowMtExtrapolationSS.at(histo).GetBinContent(2);
			catEPSignalSS.at(1) += Cat0JetHighMtExtrapolationSS.at(histo).GetBinContent(1);
			catEPSidebandSS.at(1)+= Cat0JetHighMtExtrapolationSS.at(histo).GetBinContent(2);
			catEPSignalSS.at(2) += Cat1JetLowMtExtrapolationSS.at(histo).GetBinContent(1);
			catEPSidebandSS.at(2) += Cat1JetLowMtExtrapolationSS.at(histo).GetBinContent(2);
			catEPSignalSS.at(3)+= Cat1JetHighMtExtrapolationSS.at(histo).GetBinContent(1);
			catEPSidebandSS.at(3) += Cat1JetHighMtExtrapolationSS.at(histo).GetBinContent(2);
			catEPSignalSS.at(4) += Cat1JetBoostMtExtrapolationSS.at(histo).GetBinContent(1);
			catEPSidebandSS.at(4) += Cat1JetBoostMtExtrapolationSS.at(histo).GetBinContent(2);
			catEPSignalSS.at(5) += CatVBFLooseRelaxMtExtrapolation.at(histo).GetBinContent(1); // no OS for VBT mT extrapolation
			catEPSidebandSS.at(5) += CatVBFLooseRelaxMtExtrapolation.at(histo).GetBinContent(2); // no OS for VBT mT extrapolation
			catEPSignalSS.at(6) += CatVBFTightRelaxMtExtrapolation.at(histo).GetBinContent(1); // no OS for VBT mT extrapolation
			catEPSidebandSS.at(6) += CatVBFTightRelaxMtExtrapolation.at(histo).GetBinContent(2); // no OS for VBT mT extrapolation
			catEPSignalSS.at(7)+= CatInclusiveMtExtrapolationSS.at(histo).GetBinContent(1);
			catEPSidebandSS.at(7) += CatInclusiveMtExtrapolationSS.at(histo).GetBinContent(2);
		}
	}
	for (unsigned int icat = 0; icat < nCat; icat++){
		catEPFactor.at(icat) = catEPSignal.at(icat)/catEPSideband.at(icat);
		catEPFactorSS.at(icat) = catEPSignalSS.at(icat)/catEPSidebandSS.at(icat);
	}

	// print MC scales before scaling
	for (unsigned id = 2; id < 100; id++){
			if (HConfig.GetHisto(false,id,histo)){
				double scale = Lumi * HConfig.GetCrossSection(id) / Npassed.at(histo).GetBinContent(0);
				printf("ID = %2i will be scaled by %4f \n", id, scale);
			}
	}

	// do plotting and scale histograms
	Selection::Finish();
	// all histograms below this are scaled to luminosity!

	// print MC scales after scaling
	for (unsigned id = 2; id < 100; id++){
			if (HConfig.GetHisto(false,id,histo)){
				double scale = Lumi * HConfig.GetCrossSection(id) / Npassed.at(histo).GetBinContent(0);
				printf("ID = %2i has scale of %4f after scaling\n", id, scale);
			}
	}

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

		catSBDataSS.at(0) = Cat0JetLowMtSidebandSS.at(histo).Integral();
		catSBDataSS.at(1) = Cat0JetHighMtSidebandSS.at(histo).Integral();
		catSBDataSS.at(2) = Cat1JetLowMtSidebandSS.at(histo).Integral();
		catSBDataSS.at(3) = Cat1JetHighMtSidebandSS.at(histo).Integral();
		catSBDataSS.at(4) = Cat1JetBoostMtSidebandSS.at(histo).Integral();
		catSBDataSS.at(5) = CatVBFLooseMtSidebandSS.at(histo).Integral();
		catSBDataSS.at(6) = CatVBFTightMtSidebandSS.at(histo).Integral();
		catSBDataSS.at(7) = CatInclusiveMtSidebandSS.at(histo).Integral();
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

			catSBBackgroundsSS.at(0) += Cat0JetLowMtSidebandSS.at(histo).Integral();
			catSBBackgroundsSS.at(1) += Cat0JetHighMtSidebandSS.at(histo).Integral();
			catSBBackgroundsSS.at(2) += Cat1JetLowMtSidebandSS.at(histo).Integral();
			catSBBackgroundsSS.at(3) += Cat1JetHighMtSidebandSS.at(histo).Integral();
			catSBBackgroundsSS.at(4) += Cat1JetBoostMtSidebandSS.at(histo).Integral();
			catSBBackgroundsSS.at(5) += CatVBFLooseMtSidebandSS.at(histo).Integral();
			catSBBackgroundsSS.at(6) += CatVBFTightMtSidebandSS.at(histo).Integral();
			catSBBackgroundsSS.at(7) += CatInclusiveMtSidebandSS.at(histo).Integral();
		}
	}
	for (unsigned int icat = 0; icat < nCat; icat++){
		catWJetsInSB.at(icat) = catSBData.at(icat) - catSBBackgrounds.at(icat);
		catWJetsYield.at(icat) = catWJetsInSB.at(icat) * catEPFactor.at(icat);
		catWJetsMCRatio.at(icat) = catWJetsYield.at(icat) / catWJetMCPrediction.at(icat);
		catWJetsRelaxedMCRatio.at(icat) = catWJetsYield.at(icat) / catWJetRelaxedMCPrediction.at(icat);

		catWJetsInSBSS.at(icat) = catSBDataSS.at(icat) - catSBBackgroundsSS.at(icat);
		catWJetsYieldSS.at(icat) = catWJetsInSBSS.at(icat) * catEPFactorSS.at(icat);
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

	//////////// QCD ////////////
	// calculate the OS/SS factor
	std::vector<TH2D> catQcdABCD(nCat,TH2D());
	std::vector<double> catOsSsRatio(nCat,0.0);
	if (HConfig.GetHisto(true,1,histo)){
		catQcdABCD.at(0) = Cat0JetLowQcdAbcd.at(histo);
		catQcdABCD.at(1) = Cat0JetHighQcdAbcd.at(histo);
		catQcdABCD.at(2) = Cat1JetLowQcdAbcd.at(histo);
		catQcdABCD.at(3) = Cat1JetHighQcdAbcd.at(histo);
		catQcdABCD.at(4) = Cat1JetBoostQcdAbcd.at(histo);
		catQcdABCD.at(5) = CatVBFLooseQcdAbcd.at(histo);
		catQcdABCD.at(6) = CatVBFTightQcdAbcd.at(histo);
		catQcdABCD.at(7) = CatInclusiveQcdAbcd.at(histo);
	}
	for (unsigned int icat = 0; icat < nCat; icat++){
		catOsSsRatio.at(icat) = catQcdABCD.at(icat).GetBinContent(2,1) / catQcdABCD.at(icat).GetBinContent(2,2);
	}

	// get events in SS region
	std::vector<double> catQcdSSYieldData(nCat,0.0);
	std::vector<double> catQcdSSYieldWJets(nCat,0.0);
	std::vector<double> catQcdSSYieldMCBG(nCat,0.0);
	std::vector<double> catQcdSSYieldBGCleaned(nCat,0.0);
	std::vector<double> catQcdOSYield(nCat,0.0);
	for (unsigned id = 30; id < 80; id++){ //remove DY, diboson, top from MC
		if (id == 60) continue; // 60 = QCD
		if (HConfig.GetHisto(false,id,histo)){
			catQcdSSYieldMCBG.at(0) += Cat0JetLowQcdAbcd.at(histo).GetBinContent(1,2);
			catQcdSSYieldMCBG.at(1) += Cat0JetHighQcdAbcd.at(histo).GetBinContent(1,2);
			catQcdSSYieldMCBG.at(2) += Cat1JetLowQcdAbcd.at(histo).GetBinContent(1,2);
			catQcdSSYieldMCBG.at(3) += Cat1JetHighQcdAbcd.at(histo).GetBinContent(1,2);
			catQcdSSYieldMCBG.at(4) += Cat1JetBoostQcdAbcd.at(histo).GetBinContent(1,2);
			catQcdSSYieldMCBG.at(5) += CatVBFLooseQcdAbcd.at(histo).GetBinContent(1,2);
			catQcdSSYieldMCBG.at(6) += CatVBFTightQcdAbcd.at(histo).GetBinContent(1,2);
			catQcdSSYieldMCBG.at(7) += CatInclusiveQcdAbcd.at(histo).GetBinContent(1,2);
		}
	}
	for (unsigned int icat = 0; icat < nCat; icat++){
		catQcdSSYieldData.at(icat) = catQcdABCD.at(icat).GetBinContent(1,2);
		catQcdSSYieldWJets.at(icat) = catWJetsYieldSS.at(icat);

		catQcdSSYieldBGCleaned.at(icat) = catQcdSSYieldData.at(icat) - catQcdSSYieldWJets.at(icat) - catQcdSSYieldMCBG.at(icat);
		catQcdOSYield.at(icat) = catQcdSSYieldBGCleaned.at(icat) * catOsSsRatio.at(icat);
	}

	std::cout << "  ############# QCD: OS/SS ratio #######################" << std::endl;
	printf("%12s  %12s / %12s = %12s\n", "Category", "N(OS)", "N(SS)", "OS/SS ratio");
	format = "%12s  %12.1f / %12.1f = %12f\n";
	for (unsigned int icat = 0; icat < nCat; icat++){
		double os = catQcdABCD.at(icat).GetBinContent(2,1);
		double ss = catQcdABCD.at(icat).GetBinContent(2,2);
		printf(format, catNames.at(icat).Data(), os, ss, catOsSsRatio.at(icat));
	}

	std::cout << "  ############# QCD: SS Yield #######################" << std::endl;
	printf("%12s  %12s - %12s - %12s = %12s\n", "Category", "N(Data)", "N(WJets)", "N(MC BG)", "QCD SS Yield");
	format = "%12s  %12.1f - %12.1f - %12.1f = %12f\n";
	for (unsigned int icat = 0; icat < nCat; icat++){
		printf(format, catNames.at(icat).Data(), catQcdSSYieldData.at(icat), catQcdSSYieldWJets.at(icat), catQcdSSYieldMCBG.at(icat), catQcdSSYieldBGCleaned.at(icat));
	}

	std::cout << "  ############# QCD: OS Yield #######################" << std::endl;
	printf("%12s  %12s * %12s = %12s\n", "Category", "SS Yield", "OS/SS ratio", "QCD OS Yield");
	format = "%12s  %12.1f * %12.5f = %12f\n";
	for (unsigned int icat = 0; icat < nCat; icat++){
		printf(format, catNames.at(icat).Data(), catQcdSSYieldBGCleaned.at(icat), catOsSsRatio.at(icat), catQcdOSYield.at(icat));
	}

	std::cout << "  ##########################################################\n" << std::endl;
	printf("Please copy the following numbers in order to use the data driven WJets yield:\n");
	for (unsigned int icat = 0; icat < nCat; icat++){
		printf("%12s : %14.8f\n", catNames.at(icat).Data(), catWJetsYield.at(icat));
	}
	printf("Please copy the following numbers in order to use the data driven QCD yield:\n");
	for (unsigned int icat = 0; icat < nCat; icat++){
		printf("%12s : %14.8f\n", catNames.at(icat).Data(), catQcdOSYield.at(icat));
	}
}
