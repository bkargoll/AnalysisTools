/*
 * HToTaumuTauhBackgrounds.cxx
 *
 *  Created on: Jun 20, 2014
 *      Author: kargoll
 */

#include "HToTaumuTauhBackgrounds.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

HToTaumuTauhBackgrounds::HToTaumuTauhBackgrounds(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)

{
	std::cout << "Setting up the class HToTaumuTauhBackgrounds" << std::endl;
	// always run without category for background methods
	// the numbers will be produced for all categories individually
	categoryFlag = "NoCategory";

	// run BG methods using DY-MC (to estimate DY yield)
	useEmbedding = false;

	// run Skim always using MC for WJets BG
	wJetsBGSource = "MC";

	// don't use QCD shape for background estimation
	qcdShapeFromData = false;
}

HToTaumuTauhBackgrounds::~HToTaumuTauhBackgrounds() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "HToTaumuTauhBackgrounds::~HToTaumuTauhBackgrounds Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  std::cout << "HToTaumuTauhBackgrounds::~HToTaumuTauhBackgrounds()" << std::endl;
}

void HToTaumuTauhBackgrounds::Setup(){
	if (verbose) std::cout << "HToTaumuTauhBackgrounds::Setup()" << std::endl;

	Cat0JetLowMt = HConfig.GetTH1D(Name + "_Cat0JetLowMt", "Cat0JetLowMt", 125, 0., 250., "0JL: m_{T}/GeV");
	Cat0JetLowMtSideband = HConfig.GetTH1D(Name + "_Cat0JetLowMtSideband", "Cat0JetLowMtSideband", 90, 70., 250., "0JL: m_{T}/GeV");
	Cat0JetLowMtExtrapolation = HConfig.GetTH1D(Name + "_Cat0JetLowMtExtrapolation", "Cat0JetLowMtExtrapolation", 2, 0.5, 2.5, "0JL: m_{T} signal and sideband");
	Cat0JetHighMt = HConfig.GetTH1D(Name + "_Cat0JetHighMt", "Cat0JetHighMt", 125, 0., 250., "0JH: m_{T}/GeV");
	Cat0JetHighMtSideband = HConfig.GetTH1D(Name + "_Cat0JetHighMtSideband", "Cat0JetHighMtSideband", 90, 70., 250., "0JH: m_{T}/GeV");
	Cat0JetHighMtExtrapolation = HConfig.GetTH1D(Name + "_Cat0JetHighMtExtrapolation", "Cat0JetHighMtExtrapolation", 2, 0.5, 2.5, "0JH: m_{T} signal and sideband");
	Cat1JetLowMt = HConfig.GetTH1D(Name + "_Cat1JetLowMt", "Cat1JetLowMt", 125, 0., 250., "1JL: m_{T}/GeV");
	Cat1JetLowMtSideband = HConfig.GetTH1D(Name + "_Cat1JetLowMtSideband", "Cat1JetLowMtSideband", 90, 70., 250., "1JL: m_{T}/GeV");
	Cat1JetLowMtExtrapolation = HConfig.GetTH1D(Name + "_Cat1JetLowMtExtrapolation", "Cat1JetLowMtExtrapolation", 2, 0.5, 2.5, "1JL: m_{T} signal and sideband");
	Cat1JetHighMt = HConfig.GetTH1D(Name + "_Cat1JetHighMt", "Cat1JetHighMt", 125, 0., 250., "1JH: m_{T}/GeV");
	Cat1JetHighMtSideband = HConfig.GetTH1D(Name + "_Cat1JetHighMtSideband", "Cat1JetHighMtSideband", 90, 70., 250., "1JH: m_{T}/GeV");
	Cat1JetHighMtExtrapolation = HConfig.GetTH1D(Name + "_Cat1JetHighMtExtrapolation", "Cat1JetHighMtExtrapolation", 2, 0.5, 2.5, "1JH: m_{T} signal and sideband");
	Cat1JetBoostMt = HConfig.GetTH1D(Name + "_Cat1JetBoostMt", "Cat1JetBoostMt", 125, 0., 250., "1JB: m_{T}/GeV");
	Cat1JetBoostMtSideband = HConfig.GetTH1D(Name + "_Cat1JetBoostMtSideband", "Cat1JetBoostMtSideband", 90, 70., 250., "1JB: m_{T}/GeV");
	Cat1JetBoostMtExtrapolation = HConfig.GetTH1D(Name + "_Cat1JetBoostMtExtrapolation", "Cat1JetBoostMtExtrapolation", 2, 0.5, 2.5, "1JB: m_{T} signal and sideband");
	CatVBFLooseMt = HConfig.GetTH1D(Name + "_CatVBFLooseMt", "CatVBFLooseMt", 125, 0., 250., "VBFL: m_{T}/GeV");
	CatVBFLooseMtSideband = HConfig.GetTH1D(Name + "_CatVBFLooseMtSideband", "CatVBFLooseMtSideband", 30, 60., 120., "VBFL: m_{T}/GeV");
	CatVBFLooseRelaxMt = HConfig.GetTH1D(Name + "_CatVBFLooseRelaxMt", "CatVBFLooseRelaxMt", 125, 0., 250., "VBFLRelax: m_{T}/GeV");
	CatVBFLooseRelaxMtExtrapolation = HConfig.GetTH1D(Name + "_CatVBFLooseRelaxMtExtrapolation", "CatVBFLooseRelaxMtExtrapolation", 2, 0.5, 2.5, "VBFLRelax: m_{T} signal and sideband");
	CatVBFTightMt = HConfig.GetTH1D(Name + "_CatVBFTightMt", "CatVBFTightMt", 125, 0., 250., "VBFT: m_{T}/GeV");
	CatVBFTightMtSideband = HConfig.GetTH1D(Name + "_CatVBFTightMtSideband", "CatVBFTightMtSideband", 30, 60., 120., "VBFT: m_{T}/GeV");
	CatVBFTightRelaxMt = HConfig.GetTH1D(Name + "_CatVBFTightRelaxMt", "CatVBFTightRelaxMt", 125, 0., 250., "VBFTRelax: m_{T}/GeV");
	CatVBFTightRelaxMtExtrapolation = HConfig.GetTH1D(Name + "_CatVBFTightRelaxMtExtrapolation", "CatVBFTightRelaxMtExtrapolation", 2, 0.5, 2.5, "VBFTRelax: m_{T} signal and sideband");
	CatInclusiveMt = HConfig.GetTH1D(Name + "_CatInclusiveMt", "CatInclusiveMt", 125, 0., 250., "Incl: m_{T}/GeV");
	CatInclusiveMtSideband = HConfig.GetTH1D(Name + "_CatInclusiveMtSideband", "CatInclusiveMtSideband", 90, 70., 250., "Incl: m_{T}/GeV");
	CatInclusiveMtExtrapolation = HConfig.GetTH1D(Name + "_CatInclusiveMtExtrapolation", "CatInclusiveMtExtrapolation", 2, 0.5, 2.5, "Incl: m_{T} signal and sideband");

	Cat0JetLowMtSS = HConfig.GetTH1D(Name + "_Cat0JetLowMtSS", "Cat0JetLowMtSS", 125, 0., 250., "0JL SS: m_{T}/GeV");
	Cat0JetLowMtSidebandSS = HConfig.GetTH1D(Name + "_Cat0JetLowMtSidebandSS", "Cat0JetLowMtSidebandSS", 90, 70., 250., "0JL SS: m_{T}/GeV");
	Cat0JetLowMtExtrapolationSS = HConfig.GetTH1D(Name + "_Cat0JetLowMtExtrapolationSS", "Cat0JetLowMtExtrapolationSS", 2, 0.5, 2.5, "0JL SS: m_{T} signal and sideband");
	Cat0JetHighMtSS = HConfig.GetTH1D(Name + "_Cat0JetHighMtSS", "Cat0JetHighMtSS", 125, 0., 250., "0JH SS: m_{T}/GeV");
	Cat0JetHighMtSidebandSS = HConfig.GetTH1D(Name + "_Cat0JetHighMtSidebandSS", "Cat0JetHighMtSidebandSS", 90, 70., 250., "0JH SS: m_{T}/GeV");
	Cat0JetHighMtExtrapolationSS = HConfig.GetTH1D(Name + "_Cat0JetHighMtExtrapolationSS", "Cat0JetHighMtExtrapolationSS", 2, 0.5, 2.5, "0JH SS: m_{T} signal and sideband");
	Cat1JetLowMtSS = HConfig.GetTH1D(Name + "_Cat1JetLowMtSS", "Cat1JetLowMtSS", 125, 0., 250., "1JL SS: m_{T}/GeV");
	Cat1JetLowMtSidebandSS = HConfig.GetTH1D(Name + "_Cat1JetLowMtSidebandSS", "Cat1JetLowMtSidebandSS", 90, 70., 250., "1JL SS: m_{T}/GeV");
	Cat1JetLowMtExtrapolationSS = HConfig.GetTH1D(Name + "_Cat1JetLowMtExtrapolationSS", "Cat1JetLowMtExtrapolationSS", 2, 0.5, 2.5, "1JL SS: m_{T} signal and sideband");
	Cat1JetHighMtSS = HConfig.GetTH1D(Name + "_Cat1JetHighMtSS", "Cat1JetHighMtSS", 125, 0., 250., "1JH SS: m_{T}/GeV");
	Cat1JetHighMtSidebandSS = HConfig.GetTH1D(Name + "_Cat1JetHighMtSidebandSS", "Cat1JetHighMtSidebandSS", 90, 70., 250., "1JH SS: m_{T}/GeV");
	Cat1JetHighMtExtrapolationSS = HConfig.GetTH1D(Name + "_Cat1JetHighMtExtrapolationSS", "Cat1JetHighMtExtrapolationSS", 2, 0.5, 2.5, "1JH SS: m_{T} signal and sideband");
	Cat1JetBoostMtSS = HConfig.GetTH1D(Name + "_Cat1JetBoostMtSS", "Cat1JetBoostMtSS", 125, 0., 250., "1JB SS: m_{T}/GeV");
	Cat1JetBoostMtSidebandSS = HConfig.GetTH1D(Name + "_Cat1JetBoostMtSidebandSS", "Cat1JetBoostMtSidebandSS", 90, 70., 250., "1JB SS: m_{T}/GeV");
	Cat1JetBoostMtExtrapolationSS = HConfig.GetTH1D(Name + "_Cat1JetBoostMtExtrapolationSS", "Cat1JetBoostMtExtrapolationSS", 2, 0.5, 2.5, "1JB SS: m_{T} signal and sideband");
	CatVBFLooseMtSS = HConfig.GetTH1D(Name + "_CatVBFLooseMtSS", "CatVBFLooseMtSS", 125, 0., 250., "VBFL SS: m_{T}/GeV");
	CatVBFLooseMtSidebandSS = HConfig.GetTH1D(Name + "_CatVBFLooseMtSidebandSS", "CatVBFLooseMtSidebandSS", 30, 60., 120., "VBFL SS: m_{T}/GeV");
	CatVBFTightMtSS = HConfig.GetTH1D(Name + "_CatVBFTightMtSS", "CatVBFTightMtSS", 125, 0., 250., "VBFT SS: m_{T}/GeV");
	CatVBFTightMtSidebandSS = HConfig.GetTH1D(Name + "_CatVBFTightMtSidebandSS", "CatVBFTightMtSidebandSS", 30, 60., 120., "VBFT SS: m_{T}/GeV");
	CatInclusiveMtSS = HConfig.GetTH1D(Name + "_CatInclusiveMtSS", "CatInclusiveMtSS", 125, 0., 250., "Incl SS: m_{T}/GeV");
	CatInclusiveMtSidebandSS = HConfig.GetTH1D(Name + "_CatInclusiveMtSidebandSS", "CatInclusiveMtSidebandSS", 90, 70., 250., "Incl SS: m_{T}/GeV");
	CatInclusiveMtExtrapolationSS = HConfig.GetTH1D(Name + "_CatInclusiveMtExtrapolationSS", "CatInclusiveMtExtrapolationSS", 2, 0.5, 2.5, "Incl SS: m_{T} signal and sideband");

	Cat0JetLowMtAntiIso = HConfig.GetTH1D(Name + "_Cat0JetLowMtAntiIso", "Cat0JetLowMtAntiIso", 100, 0., 150., "0JL noiso: m_{T}/GeV");
	Cat0JetHighMtAntiIso = HConfig.GetTH1D(Name + "_Cat0JetHighMtAntiIso", "Cat0JetHighMtAntiIso", 100, 0., 150., "0JH noiso: m_{T}/GeV");
	Cat1JetLowMtAntiIso = HConfig.GetTH1D(Name + "_Cat1JetLowMtAntiIso", "Cat1JetLowMtAntiIso", 100, 0., 150., "1JL noiso: m_{T}/GeV");
	Cat1JetHighMtAntiIso = HConfig.GetTH1D(Name + "_Cat1JetHighMtAntiIso", "Cat1JetHighMtAntiIso", 100, 0., 150., "1JH noiso: m_{T}/GeV");
	Cat1JetBoostMtAntiIso = HConfig.GetTH1D(Name + "_Cat1JetBoostMtAntiIso", "Cat1JetBoostMtAntiIso", 100, 0., 150., "1JB noiso: m_{T}/GeV");
	CatVBFLooseMtAntiIso = HConfig.GetTH1D(Name + "_CatVBFLooseMtAntiIso", "CatVBFLooseMtAntiIso", 100, 0., 150., "VBFL noiso: m_{T}/GeV");
	CatVBFTightMtAntiIso = HConfig.GetTH1D(Name + "_CatVBFTightMtAntiIso", "CatVBFTightMtAntiIso", 100, 0., 150., "VBFT noiso: m_{T}/GeV");
	CatInclusiveMtAntiIso = HConfig.GetTH1D(Name + "_CatInclusiveMtAntiIso", "CatInclusiveMtAntiIso", 100, 0., 150., "Incl noiso: m_{T}/GeV");
	Cat0JetLowMtAntiIsoSS = HConfig.GetTH1D(Name + "_Cat0JetLowMtAntiIsoSS", "Cat0JetLowMtAntiIsoSS", 100, 0., 150., "0JL noiso SS: m_{T}/GeV");
	Cat0JetHighMtAntiIsoSS = HConfig.GetTH1D(Name + "_Cat0JetHighMtAntiIsoSS", "Cat0JetHighMtAntiIsoSS", 100, 0., 150., "0JH noiso SS: m_{T}/GeV");
	Cat1JetLowMtAntiIsoSS = HConfig.GetTH1D(Name + "_Cat1JetLowMtAntiIsoSS", "Cat1JetLowMtAntiIsoSS", 100, 0., 150., "1JL noiso SS: m_{T}/GeV");
	Cat1JetHighMtAntiIsoSS = HConfig.GetTH1D(Name + "_Cat1JetHighMtAntiIsoSS", "Cat1JetHighMtAntiIsoSS", 100, 0., 150., "1JH noiso SS: m_{T}/GeV");
	Cat1JetBoostMtAntiIsoSS = HConfig.GetTH1D(Name + "_Cat1JetBoostMtAntiIsoSS", "Cat1JetBoostMtAntiIsoSS", 100, 0., 150., "1JB noiso SS: m_{T}/GeV");
	CatVBFLooseMtAntiIsoSS = HConfig.GetTH1D(Name + "_CatVBFLooseMtAntiIsoSS", "CatVBFLooseMtAntiIsoSS", 100, 0., 150., "VBFL noiso SS: m_{T}/GeV");
	CatVBFTightMtAntiIsoSS = HConfig.GetTH1D(Name + "_CatVBFTightMtAntiIsoSS", "CatVBFTightMtAntiIsoSS", 100, 0., 150., "VBFT noiso SS: m_{T}/GeV");
	CatInclusiveMtAntiIsoSS = HConfig.GetTH1D(Name + "_CatInclusiveMtAntiIsoSS", "CatInclusiveMtAntiIsoSS", 100, 0., 150., "Incl noiso SS: m_{T}/GeV");

	Cat0JetLowQcdAbcd = HConfig.GetTH1D(Name+"_Cat0JetLowQcdAbcd","Cat0JetLowQcdAbcd",5,-0.5,4.5,"0JL: ABCD");
	Cat0JetHighQcdAbcd = HConfig.GetTH1D(Name+"_Cat0JetHighQcdAbcd","Cat0JetHighQcdAbcd",5,-0.5,4.5,"0JH: ABCD");
	Cat1JetLowQcdAbcd = HConfig.GetTH1D(Name+"_Cat1JetLowQcdAbcd","Cat1JetLowQcdAbcd",5,-0.5,4.5,"1JL: ABCD");
	Cat1JetHighQcdAbcd = HConfig.GetTH1D(Name+"_Cat1JetHighQcdAbcd","Cat1JetHighQcdAbcd",5,-0.5,4.5,"1JH: ABCD");
	Cat1JetBoostQcdAbcd = HConfig.GetTH1D(Name+"_Cat1JetBoostQcdAbcd","Cat1JetBoostQcdAbcd",5,-0.5,4.5,"1JB: ABCD");
	CatVBFLooseQcdAbcd = HConfig.GetTH1D(Name+"_CatVBFLooseQcdAbcd","CatVBFLooseQcdAbcd",5,-0.5,4.5,"VBFL: ABCD");
	CatVBFTightQcdAbcd = HConfig.GetTH1D(Name+"_CatVBFTightQcdAbcd","CatVBFTightQcdAbcd",5,-0.5,4.5,"VBFT: ABCD");
	CatInclusiveQcdAbcd = HConfig.GetTH1D(Name+"_CatInclusiveQcdAbcd","CatInclusiveQcdAbcd",5,-0.5,4.5,"Incl: ABCD");

	Cat0JetLowQcdOSMuIso = HConfig.GetTH1D(Name + "_Cat0JetLowQcdOSMuIso", "Cat0JetLowQcdOSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat0JetLowQcdOSTauIso = HConfig.GetTH1D(Name + "_Cat0JetLowQcdOSTauIso", "Cat0JetLowQcdOSTauIso", 50, 0., 20., "iso(#tau_{h})");
	Cat0JetLowQcdSSMuIso = HConfig.GetTH1D(Name + "_Cat0JetLowQcdSSMuIso", "Cat0JetLowQcdSSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat0JetLowQcdSSTauIso = HConfig.GetTH1D(Name + "_Cat0JetLowQcdSSTauIso", "Cat0JetLowQcdSSTauIso", 50, 0., 20., "iso(#tau_{h})");
	Cat0JetHighQcdOSMuIso = HConfig.GetTH1D(Name + "_Cat0JetHighQcdOSMuIso", "Cat0JetHighQcdOSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat0JetHighQcdOSTauIso = HConfig.GetTH1D(Name + "_Cat0JetHighQcdOSTauIso", "Cat0JetHighQcdOSTauIso", 50, 0., 20., "iso(#tau_{h})");
	Cat0JetHighQcdSSMuIso = HConfig.GetTH1D(Name + "_Cat0JetHighQcdSSMuIso", "Cat0JetHighQcdSSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat0JetHighQcdSSTauIso = HConfig.GetTH1D(Name + "_Cat0JetHighQcdSSTauIso", "Cat0JetHighQcdSSTauIso", 50, 0., 20., "iso(#tau_{h})");
	Cat1JetLowQcdOSMuIso = HConfig.GetTH1D(Name + "_Cat1JetLowQcdOSMuIso", "Cat1JetLowQcdOSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat1JetLowQcdOSTauIso = HConfig.GetTH1D(Name + "_Cat1JetLowQcdOSTauIso", "Cat1JetLowQcdOSTauIso", 50, 0., 20., "iso(#tau_{h})");
	Cat1JetLowQcdSSMuIso = HConfig.GetTH1D(Name + "_Cat1JetLowQcdSSMuIso", "Cat1JetLowQcdSSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat1JetLowQcdSSTauIso = HConfig.GetTH1D(Name + "_Cat1JetLowQcdSSTauIso", "Cat1JetLowQcdSSTauIso", 50, 0., 20., "iso(#tau_{h})");
	Cat1JetHighQcdOSMuIso = HConfig.GetTH1D(Name + "_Cat1JetHighQcdOSMuIso", "Cat1JetHighQcdOSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat1JetHighQcdOSTauIso = HConfig.GetTH1D(Name + "_Cat1JetHighQcdOSTauIso", "Cat1JetHighQcdOSTauIso", 50, 0., 20., "iso(#tau_{h})");
	Cat1JetHighQcdSSMuIso = HConfig.GetTH1D(Name + "_Cat1JetHighQcdSSMuIso", "Cat1JetHighQcdSSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat1JetHighQcdSSTauIso = HConfig.GetTH1D(Name + "_Cat1JetHighQcdSSTauIso", "Cat1JetHighQcdSSTauIso", 50, 0., 20., "iso(#tau_{h})");
	Cat1JetBoostQcdOSMuIso = HConfig.GetTH1D(Name + "_Cat1JetBoostQcdOSMuIso", "Cat1JetBoostQcdOSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat1JetBoostQcdOSTauIso = HConfig.GetTH1D(Name + "_Cat1JetBoostQcdOSTauIso", "Cat1JetBoostQcdOSTauIso", 50, 0., 20., "iso(#tau_{h})");
	Cat1JetBoostQcdSSMuIso = HConfig.GetTH1D(Name + "_Cat1JetBoostQcdSSMuIso", "Cat1JetBoostQcdSSMuIso", 50, 0., 1., "relIso(#mu)");
	Cat1JetBoostQcdSSTauIso = HConfig.GetTH1D(Name + "_Cat1JetBoostQcdSSTauIso", "Cat1JetBoostQcdSSTauIso", 50, 0., 20., "iso(#tau_{h})");
	CatVBFLooseQcdOSMuIso = HConfig.GetTH1D(Name + "_CatVBFLooseQcdOSMuIso", "CatVBFLooseQcdOSMuIso", 50, 0., 1., "relIso(#mu)");
	CatVBFLooseQcdOSTauIso = HConfig.GetTH1D(Name + "_CatVBFLooseQcdOSTauIso", "CatVBFLooseQcdOSTauIso", 50, 0., 20., "iso(#tau_{h})");
	CatVBFLooseQcdSSMuIso = HConfig.GetTH1D(Name + "_CatVBFLooseQcdSSMuIso", "CatVBFLooseQcdSSMuIso", 50, 0., 1., "relIso(#mu)");
	CatVBFLooseQcdSSTauIso = HConfig.GetTH1D(Name + "_CatVBFLooseQcdSSTauIso", "CatVBFLooseQcdSSTauIso", 50, 0., 20., "iso(#tau_{h})");
	CatVBFTightQcdOSMuIso = HConfig.GetTH1D(Name + "_CatVBFTightQcdOSMuIso", "CatVBFTightQcdOSMuIso", 50, 0., 1., "relIso(#mu)");
	CatVBFTightQcdOSTauIso = HConfig.GetTH1D(Name + "_CatVBFTightQcdOSTauIso", "CatVBFTightQcdOSTauIso", 50, 0., 20., "iso(#tau_{h})");
	CatVBFTightQcdSSMuIso = HConfig.GetTH1D(Name + "_CatVBFTightQcdSSMuIso", "CatVBFTightQcdSSMuIso", 50, 0., 1., "relIso(#mu)");
	CatVBFTightQcdSSTauIso = HConfig.GetTH1D(Name + "_CatVBFTightQcdSSTauIso", "CatVBFTightQcdSSTauIso", 50, 0., 20., "iso(#tau_{h})");
	CatInclusiveQcdOSMuIso = HConfig.GetTH1D(Name + "_CatInclusiveQcdOSMuIso", "CatInclusiveQcdOSMuIso", 50, 0., 1., "relIso(#mu)");
	CatInclusiveQcdOSTauIso = HConfig.GetTH1D(Name + "_CatInclusiveQcdOSTauIso", "CatInclusiveQcdOSTauIso", 50, 0., 20., "iso(#tau_{h})");
	CatInclusiveQcdSSMuIso = HConfig.GetTH1D(Name + "_CatInclusiveQcdSSMuIso", "CatInclusiveQcdSSMuIso", 50, 0., 1., "relIso(#mu)");
	CatInclusiveQcdSSTauIso = HConfig.GetTH1D(Name + "_CatInclusiveQcdSSTauIso", "CatInclusiveQcdSSTauIso", 50, 0., 20., "iso(#tau_{h})");

	CatInclusiveMtOSChargeSum = HConfig.GetTH1D(Name + "_CatInclusiveMtOSChargeSum", "CatInclusiveMtOSChargeSum", 5,-2.5,2.5, "q(#mu)+q(#tau)");
	CatInclusiveMtSSChargeSum = HConfig.GetTH1D(Name + "_CatInclusiveMtSSChargeSum", "CatInclusiveMtSSChargeSum", 5,-2.5,2.5, "q(#mu)+q(#tau)");

	CatVBFLooseQcdEff  = HConfig.GetTH1D(Name + "_CatVBFLooseQcdEff",  "CatVBFLooseQcdEff",  2,-0.5,1.5, "passed VBFL selection");
	CatVBFTightQcdEff  = HConfig.GetTH1D(Name + "_CatVBFTightQcdEff",  "CatVBFTightQcdEff",  2,-0.5,1.5, "passed VBFT selection");
	Cat1JetBoostQcdEff = HConfig.GetTH1D(Name + "_Cat1JetBoostQcdEff", "Cat1JetBoostQcdEff", 2,-0.5,1.5, "passed 1JB selection");
}

void HToTaumuTauhBackgrounds::Configure(){
	if (verbose) std::cout << "HToTaumuTauhBackgrounds::Configure()" << std::endl;
	HToTaumuTauh::Setup();
	Setup();
	Selection::ConfigureHistograms();
	HConfig.GetHistoInfo(types, CrossSectionandAcceptance, legend, colour);
}

void HToTaumuTauhBackgrounds::Store_ExtraDist(){
	if (verbose) std::cout << "HToTaumuTauhBackgrounds::Store_ExtraDist()" << std::endl;
	HToTaumuTauh::Store_ExtraDist();

	Extradist1d.push_back(&Cat0JetLowMt);
	Extradist1d.push_back(&Cat0JetLowMtSideband);
	Extradist1d.push_back(&Cat0JetLowMtExtrapolation);
	Extradist1d.push_back(&Cat0JetHighMt);
	Extradist1d.push_back(&Cat0JetHighMtSideband);
	Extradist1d.push_back(&Cat0JetHighMtExtrapolation);
	Extradist1d.push_back(&Cat1JetLowMt);
	Extradist1d.push_back(&Cat1JetLowMtSideband);
	Extradist1d.push_back(&Cat1JetLowMtExtrapolation);
	Extradist1d.push_back(&Cat1JetHighMt);
	Extradist1d.push_back(&Cat1JetHighMtSideband);
	Extradist1d.push_back(&Cat1JetHighMtExtrapolation);
	Extradist1d.push_back(&Cat1JetBoostMt);
	Extradist1d.push_back(&Cat1JetBoostMtSideband);
	Extradist1d.push_back(&Cat1JetBoostMtExtrapolation);
	Extradist1d.push_back(&CatVBFLooseMt);
	Extradist1d.push_back(&CatVBFLooseMtSideband);
	Extradist1d.push_back(&CatVBFLooseRelaxMt);
	Extradist1d.push_back(&CatVBFLooseRelaxMtExtrapolation);
	Extradist1d.push_back(&CatVBFTightMt);
	Extradist1d.push_back(&CatVBFTightMtSideband);
	Extradist1d.push_back(&CatVBFTightRelaxMt);
	Extradist1d.push_back(&CatVBFTightRelaxMtExtrapolation);
	Extradist1d.push_back(&CatInclusiveMt);
	Extradist1d.push_back(&CatInclusiveMtSideband);
	Extradist1d.push_back(&CatInclusiveMtExtrapolation);

	Extradist1d.push_back(&Cat0JetLowMtSS);
	Extradist1d.push_back(&Cat0JetLowMtSidebandSS);
	Extradist1d.push_back(&Cat0JetLowMtExtrapolationSS);
	Extradist1d.push_back(&Cat0JetHighMtSS);
	Extradist1d.push_back(&Cat0JetHighMtSidebandSS);
	Extradist1d.push_back(&Cat0JetHighMtExtrapolationSS);
	Extradist1d.push_back(&Cat1JetLowMtSS);
	Extradist1d.push_back(&Cat1JetLowMtSidebandSS);
	Extradist1d.push_back(&Cat1JetLowMtExtrapolationSS);
	Extradist1d.push_back(&Cat1JetHighMtSS);
	Extradist1d.push_back(&Cat1JetHighMtSidebandSS);
	Extradist1d.push_back(&Cat1JetHighMtExtrapolationSS);
	Extradist1d.push_back(&Cat1JetBoostMtSS);
	Extradist1d.push_back(&Cat1JetBoostMtSidebandSS);
	Extradist1d.push_back(&Cat1JetBoostMtExtrapolationSS);
	Extradist1d.push_back(&CatVBFLooseMtSS);
	Extradist1d.push_back(&CatVBFLooseMtSidebandSS);
	Extradist1d.push_back(&CatVBFTightMtSS);
	Extradist1d.push_back(&CatVBFTightMtSidebandSS);
	Extradist1d.push_back(&CatInclusiveMtSS);
	Extradist1d.push_back(&CatInclusiveMtSidebandSS);
	Extradist1d.push_back(&CatInclusiveMtExtrapolationSS);

	Extradist1d.push_back(&Cat0JetLowMtAntiIso);
	Extradist1d.push_back(&Cat0JetHighMtAntiIso);
	Extradist1d.push_back(&Cat1JetLowMtAntiIso);
	Extradist1d.push_back(&Cat1JetHighMtAntiIso);
	Extradist1d.push_back(&Cat1JetBoostMtAntiIso);
	Extradist1d.push_back(&CatVBFLooseMtAntiIso);
	Extradist1d.push_back(&CatVBFTightMtAntiIso);
	Extradist1d.push_back(&CatInclusiveMtAntiIso);
	Extradist1d.push_back(&Cat0JetLowMtAntiIsoSS);
	Extradist1d.push_back(&Cat0JetHighMtAntiIsoSS);
	Extradist1d.push_back(&Cat1JetLowMtAntiIsoSS);
	Extradist1d.push_back(&Cat1JetHighMtAntiIsoSS);
	Extradist1d.push_back(&Cat1JetBoostMtAntiIsoSS);
	Extradist1d.push_back(&CatVBFLooseMtAntiIsoSS);
	Extradist1d.push_back(&CatVBFTightMtAntiIsoSS);
	Extradist1d.push_back(&CatInclusiveMtAntiIsoSS);

	Extradist1d.push_back(&Cat0JetLowQcdAbcd);
	Extradist1d.push_back(&Cat0JetHighQcdAbcd);
	Extradist1d.push_back(&Cat1JetLowQcdAbcd);
	Extradist1d.push_back(&Cat1JetHighQcdAbcd);
	Extradist1d.push_back(&Cat1JetBoostQcdAbcd);
	Extradist1d.push_back(&CatVBFLooseQcdAbcd);
	Extradist1d.push_back(&CatVBFTightQcdAbcd);
	Extradist1d.push_back(&CatInclusiveQcdAbcd);

	Extradist1d.push_back(&Cat0JetLowQcdOSMuIso);
	Extradist1d.push_back(&Cat0JetLowQcdOSTauIso);
	Extradist1d.push_back(&Cat0JetLowQcdSSMuIso);
	Extradist1d.push_back(&Cat0JetLowQcdSSTauIso);
	Extradist1d.push_back(&Cat0JetHighQcdOSMuIso);
	Extradist1d.push_back(&Cat0JetHighQcdOSTauIso);
	Extradist1d.push_back(&Cat0JetHighQcdSSMuIso);
	Extradist1d.push_back(&Cat0JetHighQcdSSTauIso);
	Extradist1d.push_back(&Cat1JetLowQcdOSMuIso);
	Extradist1d.push_back(&Cat1JetLowQcdOSTauIso);
	Extradist1d.push_back(&Cat1JetLowQcdSSMuIso);
	Extradist1d.push_back(&Cat1JetLowQcdSSTauIso);
	Extradist1d.push_back(&Cat1JetHighQcdOSMuIso);
	Extradist1d.push_back(&Cat1JetHighQcdOSTauIso);
	Extradist1d.push_back(&Cat1JetHighQcdSSMuIso);
	Extradist1d.push_back(&Cat1JetHighQcdSSTauIso);
	Extradist1d.push_back(&Cat1JetBoostQcdOSMuIso);
	Extradist1d.push_back(&Cat1JetBoostQcdOSTauIso);
	Extradist1d.push_back(&Cat1JetBoostQcdSSMuIso);
	Extradist1d.push_back(&Cat1JetBoostQcdSSTauIso);
	Extradist1d.push_back(&CatVBFLooseQcdOSMuIso);
	Extradist1d.push_back(&CatVBFLooseQcdOSTauIso);
	Extradist1d.push_back(&CatVBFLooseQcdSSMuIso);
	Extradist1d.push_back(&CatVBFLooseQcdSSTauIso);
	Extradist1d.push_back(&CatVBFTightQcdOSMuIso);
	Extradist1d.push_back(&CatVBFTightQcdOSTauIso);
	Extradist1d.push_back(&CatVBFTightQcdSSMuIso);
	Extradist1d.push_back(&CatVBFTightQcdSSTauIso);
	Extradist1d.push_back(&CatInclusiveQcdOSMuIso);
	Extradist1d.push_back(&CatInclusiveQcdOSTauIso);
	Extradist1d.push_back(&CatInclusiveQcdSSMuIso);
	Extradist1d.push_back(&CatInclusiveQcdSSTauIso);
	Extradist1d.push_back(&CatInclusiveMtOSChargeSum);
	Extradist1d.push_back(&CatInclusiveMtSSChargeSum);

	Extradist1d.push_back(&CatVBFLooseQcdEff);
	Extradist1d.push_back(&CatVBFTightQcdEff);
	Extradist1d.push_back(&Cat1JetBoostQcdEff);
}

void HToTaumuTauhBackgrounds::doEvent() {
	if (verbose)
		std::cout << "HToTaumuTauhBackgrounds::doEvent()" << std::endl;
	HToTaumuTauh::doEvent();

	////// W+Jets Background estimation
	if (passedFullInclusiveSelNoMtNoOS && passed_ZeroJetLow) {
		if (pass.at(OppCharge)) {
			Cat0JetLowMt.at(t).Fill(value.at(MT), w);
			Cat0JetLowMtSideband.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat0JetLowMtExtrapolation.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat0JetLowMtExtrapolation.at(t).Fill(2, w);
			}
		} else {
			Cat0JetLowMtSS.at(t).Fill(value.at(MT), w);
			Cat0JetLowMtSidebandSS.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat0JetLowMtExtrapolationSS.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat0JetLowMtExtrapolationSS.at(t).Fill(2, w);
			}
		}
	}
	if (passedFullInclusiveSelNoMtNoOS && passed_ZeroJetHigh) {
		if (pass.at(OppCharge)) {
			Cat0JetHighMt.at(t).Fill(value.at(MT), w);
			Cat0JetHighMtSideband.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat0JetHighMtExtrapolation.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat0JetHighMtExtrapolation.at(t).Fill(2, w);
			}
		} else {
			Cat0JetHighMtSS.at(t).Fill(value.at(MT), w);
			Cat0JetHighMtSidebandSS.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat0JetHighMtExtrapolationSS.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat0JetHighMtExtrapolationSS.at(t).Fill(2, w);
			}
		}
	}
	if (passedFullInclusiveSelNoMtNoOS && passed_OneJetLow) {
		if (pass.at(OppCharge)) {
			Cat1JetLowMt.at(t).Fill(value.at(MT), w);
			Cat1JetLowMtSideband.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat1JetLowMtExtrapolation.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat1JetLowMtExtrapolation.at(t).Fill(2, w);
			}
		} else {
			Cat1JetLowMtSS.at(t).Fill(value.at(MT), w);
			Cat1JetLowMtSidebandSS.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat1JetLowMtExtrapolationSS.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat1JetLowMtExtrapolationSS.at(t).Fill(2, w);
			}
		}
	}
	if (passedFullInclusiveSelNoMtNoOS && passed_OneJetHigh) {
		if (pass.at(OppCharge)) {
			Cat1JetHighMt.at(t).Fill(value.at(MT), w);
			Cat1JetHighMtSideband.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat1JetHighMtExtrapolation.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat1JetHighMtExtrapolation.at(t).Fill(2, w);
			}
		} else {
			Cat1JetHighMtSS.at(t).Fill(value.at(MT), w);
			Cat1JetHighMtSidebandSS.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat1JetHighMtExtrapolationSS.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat1JetHighMtExtrapolationSS.at(t).Fill(2, w);
			}
		}
	}
	if (passedFullInclusiveSelNoMtNoOS && passed_OneJetBoost) {
		if (pass.at(OppCharge)) {
			Cat1JetBoostMt.at(t).Fill(value.at(MT), w);
			Cat1JetBoostMtSideband.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat1JetBoostMtExtrapolation.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat1JetBoostMtExtrapolation.at(t).Fill(2, w);
			}
		} else {
			Cat1JetBoostMtSS.at(t).Fill(value.at(MT), w);
			Cat1JetBoostMtSidebandSS.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					Cat1JetBoostMtExtrapolationSS.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					Cat1JetBoostMtExtrapolationSS.at(t).Fill(2, w);
			}
		}
	}
	if (passedFullInclusiveSelNoMtNoOS) {
		if (passed_VBFLoose) {
			if (pass.at(OppCharge)) {
				CatVBFLooseMt.at(t).Fill(value.at(MT), w);
				CatVBFLooseMtSideband.at(t).Fill(value.at(MT), w);
			} else {
				CatVBFLooseMtSS.at(t).Fill(value.at(MT), w);
				CatVBFLooseMtSidebandSS.at(t).Fill(value.at(MT), w);
			}
		}
		if (passed_VBFLooseRelaxed) {
			// VBFLoose: Do not apply OS cut for mT extrapolation factor
			CatVBFLooseRelaxMt.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					CatVBFLooseRelaxMtExtrapolation.at(t).Fill(1, w);
				if (value.at(MT) > 60. && value.at(MT) < 120.)
					CatVBFLooseRelaxMtExtrapolation.at(t).Fill(2, w);
			}
		}
		if (passed_VBFTight) {
			if (pass.at(OppCharge)) {
				CatVBFTightMt.at(t).Fill(value.at(MT), w);
				CatVBFTightMtSideband.at(t).Fill(value.at(MT), w);
			} else {
				CatVBFTightMtSS.at(t).Fill(value.at(MT), w);
				CatVBFTightMtSidebandSS.at(t).Fill(value.at(MT), w);
			}
		}
		if (passed_VBFTightRelaxed) {
			// VBFTight: Do not apply OS cut for mT extrapolation factor
			CatVBFTightRelaxMt.at(t).Fill(value.at(MT), w);
			if (isWJetMC) {
				if (pass.at(MT))
					CatVBFTightRelaxMtExtrapolation.at(t).Fill(1, w);
				if (value.at(MT) > 60. && value.at(MT) < 120.)
					CatVBFTightRelaxMtExtrapolation.at(t).Fill(2, w);
			}
		}
	}
	if (passedFullInclusiveSelNoMtNoOS) {
		if (pass.at(OppCharge)) {
			CatInclusiveMt.at(t).Fill(value.at(MT), w);
			CatInclusiveMtSideband.at(t).Fill(value.at(MT), w);
			CatInclusiveMtOSChargeSum.at(t).Fill(Ntp->Muon_Charge(selMuon) + Ntp->PFTau_Charge(selTau));
			if (isWJetMC) {
				if (pass.at(MT))
					CatInclusiveMtExtrapolation.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					CatInclusiveMtExtrapolation.at(t).Fill(2, w);
			}
		} else {
			CatInclusiveMtSS.at(t).Fill(value.at(MT), w);
			CatInclusiveMtSidebandSS.at(t).Fill(value.at(MT), w);
			CatInclusiveMtSSChargeSum.at(t).Fill(Ntp->Muon_Charge(selMuon) + Ntp->PFTau_Charge(selTau));
			if (isWJetMC) {
				if (pass.at(MT))
					CatInclusiveMtExtrapolationSS.at(t).Fill(1, w);
				if (value.at(MT) > 70.)
					CatInclusiveMtExtrapolationSS.at(t).Fill(2, w);
			}
		}
	}

	////// QCD Background estimation
	//     OS/SS
	//       ^
	//    C  |  D
	//   ---------> relIso(mu)
	//    A  |  B
	if (verbose)
		std::cout << "	QCD Background plots (ABCD)" << std::endl;
	if (passedFullInclusiveNoTauNoMuNoCharge) {
		// veto events with signal muon AND antiIsoMuon, as in these cases mT etc. are calculated using the signal muon
		bool isA = pass.at(OppCharge) && passedObjects;
		bool isB = pass.at(OppCharge) && !passedMu && hasRelaxedIsoTau && hasAntiIsoMuon;
		bool isC = !pass.at(OppCharge) && passedObjects;
		bool isD = !pass.at(OppCharge) && !passedMu && hasRelaxedIsoTau && hasAntiIsoMuon;
		if (isA + isB + isC + isD > 1)
			printf("WARNING: Event %i enters more than 1 ABCD region! (A%i, B%i, C%i, D%i)\n", Ntp->EventNumber(), isA, isB, isC, isD);
		//	  if (isA+isB+isC+isD == 0) {
		//		  printf("ATTENTION: Event %9d enters no ABCD region! Sum(q) = %d, passedMu = %d, passedTau = %d\n", Ntp->EventNumber(), value.at(OppCharge), passedMu, passedTau);
		//		  printf("       		                                    hasRelTau = %d, hasAntiIsoMu = %d\n", hasRelaxedIsoTau, hasAntiIsoMuon);
		//	  }
		// save ABCD information in a 1D plot
		int abcd(0);
		if (isA)
			abcd = 1;
		if (isB)
			abcd = 2;
		if (isC)
			abcd = 3;
		if (isD)
			abcd = 4;

		CatInclusiveQcdAbcd.at(t).Fill(abcd, w);
		if (abcd != 0) {
			if (pass.at(OppCharge)) {
				CatInclusiveQcdOSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
				CatInclusiveQcdOSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
			}
			if (!pass.at(OppCharge)) {
				CatInclusiveQcdSSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
				CatInclusiveQcdSSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
			}
			if (passed_ZeroJetLow) {
				Cat0JetLowQcdAbcd.at(t).Fill(abcd, w);
				if (pass.at(OppCharge)) {
					Cat0JetLowQcdOSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat0JetLowQcdOSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
				if (!pass.at(OppCharge)) {
					Cat0JetLowQcdSSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat0JetLowQcdSSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
			}
			if (passed_ZeroJetHigh) {
				Cat0JetHighQcdAbcd.at(t).Fill(abcd, w);
				if (pass.at(OppCharge)) {
					Cat0JetHighQcdOSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat0JetHighQcdOSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
				if (!pass.at(OppCharge)) {
					Cat0JetHighQcdSSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat0JetHighQcdSSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
			}
			if (passed_OneJetLow) {
				Cat1JetLowQcdAbcd.at(t).Fill(abcd, w);
				if (pass.at(OppCharge)) {
					Cat1JetLowQcdOSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat1JetLowQcdOSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
				if (!pass.at(OppCharge)) {
					Cat1JetLowQcdSSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat1JetLowQcdSSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
			}
			if (passed_OneJetHigh) {
				Cat1JetHighQcdAbcd.at(t).Fill(abcd, w);
				if (pass.at(OppCharge)) {
					Cat1JetHighQcdOSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat1JetHighQcdOSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
				if (!pass.at(OppCharge)) {
					Cat1JetHighQcdSSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat1JetHighQcdSSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
			}
			if (passed_OneJetBoost) {
				Cat1JetBoostQcdAbcd.at(t).Fill(abcd, w);
				if (pass.at(OppCharge)) {
					Cat1JetBoostQcdOSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat1JetBoostQcdOSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
				if (!pass.at(OppCharge)) {
					Cat1JetBoostQcdSSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					Cat1JetBoostQcdSSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
			}
			if (passed_VBFLoose) {
				CatVBFLooseQcdAbcd.at(t).Fill(abcd, w);
				if (pass.at(OppCharge)) {
					CatVBFLooseQcdOSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					CatVBFLooseQcdOSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
				if (!pass.at(OppCharge)) {
					CatVBFLooseQcdSSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					CatVBFLooseQcdSSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
			}
			if (passed_VBFTight) {
				CatVBFTightQcdAbcd.at(t).Fill(abcd, w);
				if (pass.at(OppCharge)) {
					CatVBFTightQcdOSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					CatVBFTightQcdOSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
				if (!pass.at(OppCharge)) {
					CatVBFTightQcdSSMuIso.at(t).Fill(Ntp->Muon_RelIso(selMuon), w);
					CatVBFTightQcdSSTauIso.at(t).Fill(Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau), w);
				}
			}
		}

		// === QCD efficiency method for VBF loose, VBF tight and 1Jet Boost categories ===
		// VBF loose: efficiency from sideband with same sign and anti-iso muon
		if (pass.at(NTauKin) && !pass.at(OppCharge) && !passedMu && hasAntiIsoMuon) {
			CatVBFLooseQcdEff.at(t).Fill(0., w);
			if (passed_VBFLoose) CatVBFLooseQcdEff.at(t).Fill(1., w);
		}
		// VBF tight: efficiency from sideband with same sign and anti-iso muon and relaxed tau iso
		if (!pass.at(OppCharge) && !passedMu && hasRelaxedIsoTau && hasAntiIsoMuon) {
			CatVBFTightQcdEff.at(t).Fill(0., w);
			if (passed_VBFTight) CatVBFTightQcdEff.at(t).Fill(1., w);
		}
		// 1Jet boost: efficiency from sideband with anti-iso muon and relaxed tau iso
		if (pass.at(OppCharge) && !passedMu && hasRelaxedIsoTau && hasAntiIsoMuon) {
			Cat1JetBoostQcdEff.at(t).Fill(0., w);
			if (passed_OneJetBoost) Cat1JetBoostQcdEff.at(t).Fill(1., w);
		}
	}

	if (verbose)
			std::cout << "	Mt in AntiIso region" << std::endl;
	if(!passedMu && hasRelaxedIsoTau && hasAntiIsoMuon){
		if(pass.at(OppCharge)){
			CatInclusiveMtAntiIso.at(t).Fill(value.at(MT), w);
			if (passed_ZeroJetLow) Cat0JetLowMtAntiIso.at(t).Fill(value.at(MT), w);
			if (passed_ZeroJetHigh) Cat0JetHighMtAntiIso.at(t).Fill(value.at(MT), w);
			if (passed_OneJetLow) Cat1JetLowMtAntiIso.at(t).Fill(value.at(MT), w);
			if (passed_OneJetHigh) Cat1JetHighMtAntiIso.at(t).Fill(value.at(MT), w);
			if (passed_OneJetBoost) Cat1JetBoostMtAntiIso.at(t).Fill(value.at(MT), w);
			if (passed_VBFLoose) CatVBFLooseMtAntiIso.at(t).Fill(value.at(MT), w);
			if (passed_VBFTight) CatVBFTightMtAntiIso.at(t).Fill(value.at(MT), w);
		}
		else{
			CatInclusiveMtAntiIsoSS.at(t).Fill(value.at(MT), w);
			if (passed_ZeroJetLow) Cat0JetLowMtAntiIsoSS.at(t).Fill(value.at(MT), w);
			if (passed_ZeroJetHigh) Cat0JetHighMtAntiIsoSS.at(t).Fill(value.at(MT), w);
			if (passed_OneJetLow) Cat1JetLowMtAntiIsoSS.at(t).Fill(value.at(MT), w);
			if (passed_OneJetHigh) Cat1JetHighMtAntiIsoSS.at(t).Fill(value.at(MT), w);
			if (passed_OneJetBoost) Cat1JetBoostMtAntiIsoSS.at(t).Fill(value.at(MT), w);
			if (passed_VBFLoose) CatVBFLooseMtAntiIsoSS.at(t).Fill(value.at(MT), w);
			if (passed_VBFTight) CatVBFTightMtAntiIsoSS.at(t).Fill(value.at(MT), w);
		}
	}

}

void HToTaumuTauhBackgrounds::Finish() {
	if (verbose) std::cout << "HToTaumuTauhBackgrounds::Finish()" << std::endl;

	if(wJetsBGSource != "MC"){
		std::cout << "Please set wJetsBGSource = \"MC\" to obtain background yields. Abort...";
		return;
	}

	// suppress embedding completely, if available
	if(HConfig.hasID(DataMCType::DY_mutau_embedded)) ScaleAllHistOfType(HConfig.GetType(DataMCType::DY_mutau_embedded), 0.0);

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

//	// print MC scales before scaling
//	for (unsigned id = 2; id < 100; id++){
//			if (HConfig.GetHisto(false,id,histo)){
//				double scale = Lumi * HConfig.GetCrossSection(id) / Npassed.at(histo).GetBinContent(0);
//				printf("ID = %2i will be scaled by %4f \n", id, scale);
//			}
//	}

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
		if (id == DataMCType::QCD || id == DataMCType::DY_mutau_embedded ) continue;
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
	std::vector<TH1D> catQcdABCD(nCat,TH1D());
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
		catOsSsRatio.at(icat) = catQcdABCD.at(icat).GetBinContent(catQcdABCD.at(icat).FindFixBin(2)) / catQcdABCD.at(icat).GetBinContent(catQcdABCD.at(icat).FindFixBin(4));
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
			int bin = Cat0JetLowQcdAbcd.at(histo).FindFixBin(3);
			catQcdSSYieldMCBG.at(0) += Cat0JetLowQcdAbcd.at(histo).GetBinContent(bin);
			catQcdSSYieldMCBG.at(1) += Cat0JetHighQcdAbcd.at(histo).GetBinContent(bin);
			catQcdSSYieldMCBG.at(2) += Cat1JetLowQcdAbcd.at(histo).GetBinContent(bin);
			catQcdSSYieldMCBG.at(3) += Cat1JetHighQcdAbcd.at(histo).GetBinContent(bin);
			catQcdSSYieldMCBG.at(4) += Cat1JetBoostQcdAbcd.at(histo).GetBinContent(bin);
			catQcdSSYieldMCBG.at(5) += CatVBFLooseQcdAbcd.at(histo).GetBinContent(bin);
			catQcdSSYieldMCBG.at(6) += CatVBFTightQcdAbcd.at(histo).GetBinContent(bin);
			catQcdSSYieldMCBG.at(7) += CatInclusiveQcdAbcd.at(histo).GetBinContent(bin);
		}
	}
	for (unsigned int icat = 0; icat < nCat; icat++){
		catQcdSSYieldData.at(icat) = catQcdABCD.at(icat).GetBinContent(catQcdABCD.at(icat).FindFixBin(3));
		catQcdSSYieldWJets.at(icat) = catWJetsYieldSS.at(icat);

		catQcdSSYieldBGCleaned.at(icat) = catQcdSSYieldData.at(icat) - catQcdSSYieldWJets.at(icat) - catQcdSSYieldMCBG.at(icat);
		catQcdOSYield.at(icat) = catQcdSSYieldBGCleaned.at(icat) * catOsSsRatio.at(icat);
	}

	// QCD efficiency method
	std::vector<double> catQCDEffNum(nCat, 0.0);
	std::vector<double> catQCDEffDenom(nCat, 0.0);
	std::vector<double> catQCDEfficiency(nCat, 0.0);
	std::vector<double> catQCDEffYield(nCat, 0.0);
	if (HConfig.GetHisto(true,1,histo)){
		// 1 jet boost
		unsigned int i_cat = 4;
		catQCDEffNum.at(i_cat)   = Cat1JetBoostQcdEff.at(histo).GetBinContent(2);
		catQCDEffDenom.at(i_cat) = Cat1JetBoostQcdEff.at(histo).GetBinContent(1);
		if (catQCDEffDenom.at(i_cat) != 0) catQCDEfficiency.at(i_cat) = catQCDEffNum.at(i_cat) / catQCDEffDenom.at(i_cat);
		else catQCDEfficiency.at(4) = -999;
		catQCDEffYield.at(i_cat) = catQCDEfficiency.at(i_cat) * catQcdOSYield.at(7); // efficiency * ABCD yield inclusive category
		// VBF loose
		i_cat = 5;
		catQCDEffNum.at(i_cat)   = CatVBFLooseQcdEff.at(histo).GetBinContent(2);
		catQCDEffDenom.at(i_cat) = CatVBFLooseQcdEff.at(histo).GetBinContent(1);
		if (catQCDEffDenom.at(i_cat) != 0) catQCDEfficiency.at(i_cat) = catQCDEffNum.at(i_cat) / catQCDEffDenom.at(i_cat);
		else catQCDEfficiency.at(4) = -999;
		catQCDEffYield.at(i_cat) = catQCDEfficiency.at(i_cat) * catQcdOSYield.at(7); // efficiency * ABCD yield inclusive category
		// VBF tight
		i_cat = 6;
		catQCDEffNum.at(i_cat)   = CatVBFTightQcdEff.at(histo).GetBinContent(2);
		catQCDEffDenom.at(i_cat) = CatVBFTightQcdEff.at(histo).GetBinContent(1);
		if (catQCDEffDenom.at(i_cat) != 0) catQCDEfficiency.at(i_cat) = catQCDEffNum.at(i_cat) / catQCDEffDenom.at(i_cat);
		else catQCDEfficiency.at(4) = -999;
		catQCDEffYield.at(i_cat) = catQCDEfficiency.at(i_cat) * catQcdOSYield.at(7); // efficiency * ABCD yield inclusive category
	}

	std::cout << "  ############# QCD: OS/SS ratio #######################" << std::endl;
	printf("%12s  %12s / %12s = %12s\n", "Category", "N(OS)", "N(SS)", "OS/SS ratio");
	format = "%12s  %12.1f / %12.1f = %12f\n";
	for (unsigned int icat = 0; icat < nCat; icat++){
		double os = catQcdABCD.at(icat).GetBinContent(catQcdABCD.at(icat).FindFixBin(2));
		double ss = catQcdABCD.at(icat).GetBinContent(catQcdABCD.at(icat).FindFixBin(4));
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

	std::cout << "  ############# QCD: Efficiency Method #######################" << std::endl;
	printf("%12s  %12s / %12s = %12s => %12s\n", "Category", "N(category)", "N(inclusive)", "efficiency", "yield");
	format = "%12s  %12.1f / %12.1f = %12f => %12f\n";
	for (unsigned int icat = 0; icat < nCat; icat++){
		if (catQCDEfficiency.at(icat) != 0.0)
			printf(format, catNames.at(icat).Data(), catQCDEffNum.at(icat), catQCDEffDenom.at(icat), catQCDEfficiency.at(icat), catQCDEffYield.at(icat));
	}

	std::cout << "  ##########################################################\n" << std::endl;
	printf("Please copy the following numbers in order to use the data driven WJets yield:\n");
	for (unsigned int icat = 0; icat < nCat; icat++){
		printf("%12s : %14.8f\n", catNames.at(icat).Data(), catWJetsYield.at(icat));
	}
	printf("Please copy the following numbers in order to use the data driven QCD yield (ABCD method):\n");
	for (unsigned int icat = 0; icat < nCat; icat++){
		printf("%12s : %14.8f\n", catNames.at(icat).Data(), catQcdOSYield.at(icat));
	}
	printf("Please copy the following numbers in order to use the data driven QCD yield (efficiency method):\n");
	for (unsigned int icat = 0; icat < nCat; icat++){
		if (catQCDEfficiency.at(icat) != 0.0)
			printf("%12s : %14.8f\n", catNames.at(icat).Data(), catQCDEffYield.at(icat));
	}

}
