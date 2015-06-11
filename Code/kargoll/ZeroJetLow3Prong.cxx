/*
 * ZeroJetLow3Prong.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "ZeroJetLow3Prong.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "../Tools.h"

#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"

ZeroJetLow3Prong::ZeroJetLow3Prong(TString Name_, TString id_):
	HToTaumuTauh(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class ZeroJetLow3Prong" << std::endl;
	// run ZeroJetLow3Prong category
	Logger(Logger::Info) << "Use ZeroJetLow category selection" << std::endl;
	categoryFlag = "ZeroJetLow";

	// run Categories using embedding
	useEmbedding = false; //todo

	// run Categories using data-driven WJets BG
	wJetsBGSource = "MC"; //todo

	// run Categories using data-driven QCD BG
	qcdShapeFromData = false; //todo
}

ZeroJetLow3Prong::~ZeroJetLow3Prong() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

void  ZeroJetLow3Prong::Configure(){
  Logger(Logger::Verbose) << std::endl;
  HToTaumuTauh::Setup();
  ZeroJetLow3Prong::Setup();

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types, CrossSectionandAcceptance, legend, colour);
}

void ZeroJetLow3Prong::Setup(){
	Logger(Logger::Verbose) << std::endl;
	// Setup Cut Values
	// start from CatCut3, as 1&2 are already in use by ZeroJetLow catetgory
	for(unsigned i = CatCut3; i< NCuts; i++){
	  if(i==CatCut3)	cut.at(CatCut3)	= 10; // decay mode finding
	  if(i==CatCut4)	cut.at(CatCut4)	= 3.0; // SV significance
	}

	TString hlabel;
	TString htitle;
	TString c;

	title.at(CatCut3)="decayMode $==$";
	title.at(CatCut3)+=cut.at(CatCut3);
	htitle=title.at(CatCut3);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="#tau_{h} decay mode";
	c="_Cut_";c+=CatCut3;
	Nminus1.at(CatCut3) = HConfig.GetTH1D(Name+c+"_Nminus1_CatCut3_",htitle,13,-0.5,12.5,hlabel,"Events");
	Nminus0.at(CatCut3) = HConfig.GetTH1D(Name+c+"_Nminus0_CatCut3_",htitle,13,-0.5,12.5,hlabel,"Events");

	title.at(CatCut4)="$\\sigma(SV) >=$";
	title.at(CatCut4)+=cut.at(CatCut4);
	htitle=title.at(CatCut4);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="\\sigma(SV)";
	c="_Cut_";c+=CatCut4;
	Nminus1.at(CatCut4) = HConfig.GetTH1D(Name+c+"_Nminus1_CatCut4_",htitle,80,-10.,30.,hlabel,"Events");
	Nminus0.at(CatCut4) = HConfig.GetTH1D(Name+c+"_Nminus0_CatCut4_",htitle,80,-10.,30.,hlabel,"Events");

	// additional histograms
	Tau3p_Plus_Pt	= HConfig.GetTH1D(Name+"_Tau3p_Plus_Pt","Tau3p_Pt"	,50,0.,100.,"p_{T}(#tau_{3p}^{+ sol.})/GeV");
	Tau3p_Plus_Eta	= HConfig.GetTH1D(Name+"_Tau3p_Plus_Eta","Tau3p_Eta"	,50,-2.5,2.5,"#eta(#tau_{3p}^{+ sol.})");
	Tau3p_Plus_Phi	= HConfig.GetTH1D(Name+"_Tau3p_Plus_Phi","Tau3p_Phi"	,50,-3.14159,3.14159,"#phi(#tau_{3p}^{+ sol.})");
	Tau3p_Plus_E	= HConfig.GetTH1D(Name+"_Tau3p_Plus_E","Tau3p_E"		,50,0.,200.,"E(#tau_{3p}^{+ sol.})/GeV");
	Tau3p_Minus_Pt	= HConfig.GetTH1D(Name+"_Tau3p_Minus_Pt","Tau3p_Pt"	,50,0.,100.,"p_{T}(#tau_{3p}^{- sol.})/GeV");
	Tau3p_Minus_Eta	= HConfig.GetTH1D(Name+"_Tau3p_Minus_Eta","Tau3p_Eta"	,50,-2.5,2.5,"#eta(#tau_{3p}^{- sol.})");
	Tau3p_Minus_Phi	= HConfig.GetTH1D(Name+"_Tau3p_Minus_Phi","Tau3p_Phi"	,50,-3.14159,3.14159,"#phi(#tau_{3p}^{- sol.})");
	Tau3p_Minus_E	= HConfig.GetTH1D(Name+"_Tau3p_Minus_E","Tau3p_E"		,50,0.,200.,"E(#tau_{3p}^{- sol.})/GeV");
	Tau3p_Zero_Pt	= HConfig.GetTH1D(Name+"_Tau3p_Zero_Pt","Tau3p_Pt"	,50,0.,100.,"p_{T}(#tau_{3p}^{0 sol.})/GeV");
	Tau3p_Zero_Eta	= HConfig.GetTH1D(Name+"_Tau3p_Zero_Eta","Tau3p_Eta"	,50,-2.5,2.5,"#eta(#tau_{3p}^{0 sol.})");
	Tau3p_Zero_Phi	= HConfig.GetTH1D(Name+"_Tau3p_Zero_Phi","Tau3p_Phi"	,50,-3.14159,3.14159,"#phi(#tau_{3p}^{0 sol.})");
	Tau3p_Zero_E	= HConfig.GetTH1D(Name+"_Tau3p_Zero_E","Tau3p_E"		,50,0.,200.,"E(#tau_{3p}^{0 sol.})/GeV");

	Tau3p_fit_Pt	= HConfig.GetTH1D(Name+"_Tau3p_fit_Pt","Tau3p_fit_Pt"  ,50,0.,100.,"p_{T}(#tau_{3p}^{fit})/GeV");
	Tau3p_fit_Eta	= HConfig.GetTH1D(Name+"_Tau3p_fit_Eta","Tau3p_fit_Eta",50,-2.5,2.5,"#eta(#tau_{3p}^{fit})");
	Tau3p_fit_Phi	= HConfig.GetTH1D(Name+"_Tau3p_fit_Phi","Tau3p_fit_Phi",50,-3.14159,3.14159,"#phi(#tau_{3p}^{fit})");
	Tau3p_fit_E		= HConfig.GetTH1D(Name+"_Tau3p_fit_E","Tau3p_fit_E"    ,50,0.,200.,"E(#tau_{3p}^{fit})/GeV");

	TauMu_fit_Pt	= HConfig.GetTH1D(Name+"_TauMu_fit_Pt","TauMu_fit_Pt"  ,50,0.,100.,"p_{T}(#tau_{#mu}^{fit})/GeV");
	TauMu_fit_Eta	= HConfig.GetTH1D(Name+"_TauMu_fit_Eta","TauMu_fit_Eta",50,-2.5,2.5,"#eta(#tau_{#mu}^{fit})");
	TauMu_fit_Phi	= HConfig.GetTH1D(Name+"_TauMu_fit_Phi","TauMu_fit_Phi",50,-3.14159,3.14159,"#phi(#tau_{#mu}^{fit})");
	TauMu_fit_E		= HConfig.GetTH1D(Name+"_TauMu_fit_E","TauMu_fit_E"    ,50,0.,200.,"E(#tau_{#mu}^{fit})/GeV");

	EventFit_converged	= HConfig.GetTH1D(Name+"_EventFit_converged","EventFit_converged", 2, -0.5, 1.5, "Fit: converged");
	EventFit_Chi2		= HConfig.GetTH1D(Name+"_EventFit_Chi2","EventFit_Chi2",60,0.,20.,"Fit: #chi^{2}");
	EventFit_nIterations= HConfig.GetTH1D(Name+"_EventFit_nIterations","EventFit_nIterations",100,0.,100.,"Fit: N(iterations)");
	EventFit_deltaCSum	= HConfig.GetTH1D(Name+"_EventFit_deltaCSum","EventFit_deltaCSum", 100, 0., 200., "Fit: #Sigma #Delta c");

	EventFit_Res_Pt	= HConfig.GetTH1D(Name+"_EventFit_Res_Pt","EventFit_Res_Pt"   ,50,0.,100.,"p_{T}(X^{fit})/GeV");
	EventFit_Res_Eta= HConfig.GetTH1D(Name+"_EventFit_Res_Eta","EventFit_Res_Eta" ,50,-2.5,2.5,"#eta(X^{fit})");
	EventFit_Res_Phi= HConfig.GetTH1D(Name+"_EventFit_Res_Phi","EventFit_Res_Phi" ,50,-3.14159,3.14159,"#phi(X^{fit})");
	EventFit_Res_E	= HConfig.GetTH1D(Name+"_EventFit_Res_E","EventFit_Res_E"     ,50,0.,200.,"E(X^{fit})/GeV");

	svFit3pMass_PlusSol		= HConfig.GetTH1D(Name+"_svFit3pMass_PlusSol",	"svFit3pMass_PlusSol",	50, 0., 250., "m_{SVfit}(#tau_{3p}^{+},#mu)/GeV");
	svFit3pMass_MinusSol	= HConfig.GetTH1D(Name+"_svFit3pMass_MinusSol",	"svFit3pMass_MinusSol",	50, 0., 250., "m_{SVfit}(#tau_{3p}^{-},#mu)/GeV");
	svFit3pMass_ZeroSol		= HConfig.GetTH1D(Name+"_svFit3pMass_ZeroSol",	"svFit3pMass_ZeroSol",	50, 0., 250., "m_{SVfit}(#tau_{3p}^{0},#mu)/GeV");
	svFit3pMass_TrueSol		= HConfig.GetTH1D(Name+"_svFit3pMass_TrueSol",	"svFit3pMass_TrueSol",	50, 0., 250., "m_{SVfit}(#tau_{3p}^{true},#mu)/GeV");
	svFit3pMass_TruePlusSol		= HConfig.GetTH1D(Name+"_svFit3pMass_TruePlusSol",	"svFit3pMass_TruePlusSol",	50, 0., 250., "m_{SVfit}(#tau_{3p}^{+true},#mu)/GeV");
	svFit3pMass_TrueMinusSol	= HConfig.GetTH1D(Name+"_svFit3pMass_TrueMinusSol",	"svFit3pMass_TrueMinusSol",	50, 0., 250., "m_{SVfit}(#tau_{3p}^{-true},#mu)/GeV");
	svFit3pMassResol_PlusSol	= HConfig.GetTH1D(Name+"_svFit3pMassResol_PlusSol",	"svFit3pMassResol_PlusSol",	50, -1., 1., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{+},#mu)");
	svFit3pMassResol_MinusSol	= HConfig.GetTH1D(Name+"_svFit3pMassResol_MinusSol","svFit3pMassResol_MinusSol",50, -1., 1., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{-},#mu)");
	svFit3pMassResol_ZeroSol	= HConfig.GetTH1D(Name+"_svFit3pMassResol_ZeroSol",	"svFit3pMassResol_ZeroSol",	50, -1., 1., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{0},#mu)");
	svFit3pMassResol_TrueSol	= HConfig.GetTH1D(Name+"_svFit3pMassResol_TrueSol",	"svFit3pMassResol_TrueSol",	50, -1., 1., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{true},#mu)");

	Tau_Pt_3pTrueVsSVfit3pTrueSol	= HConfig.GetTH2D(Name+"_Tau_Pt_3pTrueVsSVfit3pTrueSol","Tau_Pt_3pTrueVsSVfit3pTrueSol", 50, 0., 200., 50, 0., 200., "p_{T}(#tau_{3p}^{true})/GeV", "p_{T}(#tau_{SVFit}^{3p true})/GeV");
	Tau_Pt_3pTrueSVfit3pTrueSolResol	= HConfig.GetTH1D(Name+"_Tau_Pt_3pTrueSVfit3pTrueSolResol","Tau_Pt_3pTrueSVfit3pTrueSolResol",50, -1., 1., "#frac{p_{T}(#tau_{3p}^{true}) - p_{T}(#tau_{SVFit}^{3p true})}{p_{T}(#tau_{3p}^{true}}");
	Tau_Pt_3pTrueVsSVfit3pZeroSol = HConfig.GetTH2D(Name+"_Tau_Pt_3pTrueVsSVfit3pZeroSol","Tau_Pt_3pTrueVsSVfit3pZeroSol", 50, 0., 200., 50, 0., 200., "p_{T}(#tau_{3p}^{true})/GeV", "p_{T}(#tau_{SVFit}^{3p zero})/GeV");
	Tau_Pt_3pTrueSVfit3pZeroSolResol	= HConfig.GetTH1D(Name+"_Tau_Pt_3pTrueSVfit3pZeroSolResol","Tau_Pt_3pTrueSVfit3pZeroSolResol",50, -1., 1., "#frac{p_{T}(#tau_{3p}^{true}) - p_{T}(#tau_{SVFit}^{3p zero})}{p_{T}(#tau_{3p}^{true}}");

	vis3pMass_PlusSol		= HConfig.GetTH1D(Name+"_vis3pMass_PlusSol",	"vis3pMass_PlusSol",	50, 0., 250., "m_{vis}(#tau_{3p}^{+},#mu)/GeV");
	vis3pMass_MinusSol		= HConfig.GetTH1D(Name+"_vis3pMass_MinusSol",	"vis3pMass_MinusSol",	50, 0., 250., "m_{vis}(#tau_{3p}^{-},#mu)/GeV");
	vis3pMass_ZeroSol		= HConfig.GetTH1D(Name+"_vis3pMass_ZeroSol",	"vis3pMass_ZeroSol",	50, 0., 250., "m_{vis}(#tau_{3p}^{0},#mu)/GeV");
	vis3pMass_TrueSol		= HConfig.GetTH1D(Name+"_vis3pMass_TrueSol",	"vis3pMass_TrueSol",	50, 0., 250., "m_{vis}(#tau_{3p}^{true},#mu)/GeV");
	vis3pMassResol_PlusSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_PlusSol",	"vis3pMassResol_PlusSol",	50, -1., 1., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{+},#mu)");
	vis3pMassResol_MinusSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_MinusSol","vis3pMassResol_MinusSol",50, -1., 1., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{-},#mu)");
	vis3pMassResol_ZeroSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_ZeroSol",	"vis3pMassResol_ZeroSol",	50, -1., 1., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{0},#mu)");
	vis3pMassResol_TrueSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_TrueSol",	"vis3pMassResol_TrueSol",	50, -1., 1., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{true},#mu)");

	Tau3p_True_Pt  = HConfig.GetTH1D(Name+"_Tau3p_True_Pt","Tau3p_True_Pt"	,50,0.,100.,"p_{T}(#tau_{3p}^{true sol.})/GeV");
	Tau3p_True_Eta = HConfig.GetTH1D(Name+"_Tau3p_True_Eta","Tau3p_True_Eta",50,-2.5,2.5,"#eta(#tau_{3p}^{true sol.})");
	Tau3p_True_Phi = HConfig.GetTH1D(Name+"_Tau3p_True_Phi","Tau3p_True_Phi",50,-3.14159,3.14159,"#phi(#tau_{3p}^{true sol.})");
	Tau3p_True_E   = HConfig.GetTH1D(Name+"_Tau3p_True_E","Tau3p_True_E"	,50,0.,200.,"E(#tau_{3p}^{true sol.})/GeV");
	Tau3p_True_Pt_Resol  = HConfig.GetTH1D(Name+"_Tau3p_True_Pt_Resol","Tau3p_True_Pt_Resol"	,50, -1., 1.,"p_{T} resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}");
	Tau3p_True_Eta_Resol = HConfig.GetTH1D(Name+"_Tau3p_True_Eta_Resol","Tau3p_True_Eta_Resol"	,50, -1., 1.,"#eta resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}");
	Tau3p_True_Phi_Resol = HConfig.GetTH1D(Name+"_Tau3p_True_Phi_Resol","Tau3p_True_Phi_Resol"	,50, -1., 1.,"#phi resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}");
	Tau3p_True_E_Resol   = HConfig.GetTH1D(Name+"_Tau3p_True_E_Resol","Tau3p_True_E_Resol"		,50, -1., 1.,"E resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}");

	TrueSignalMET           = HConfig.GetTH1D(Name+"_TrueSignalMET","TrueSignalMET", 50, 0., 100., "gen. E_{T}^{miss}/GeV");
	MetPtResol              = HConfig.GetTH1D(Name+"_MetPtResol","MetPtResol", 50, -2., 2., "#frac{E_{T,reco}^{miss} - E_{T,gen}^{miss}}{E_{T,gen}^{miss}}");
	MetPhiResol             = HConfig.GetTH1D(Name+"_MetPhiResol","MetPhiResol", 50, -1.5, 1.5, "#phi(E_{T,reco}^{miss}) - #phi(E_{T,gen}^{miss})");
	MetPxResol              = HConfig.GetTH1D(Name+"_MetXResol","MetXResol", 50, -2., 2., "#frac{E_{x,reco}^{miss} - E_{x,gen}^{miss}}{E_{x,gen}^{miss}}");
	MetPyResol              = HConfig.GetTH1D(Name+"_MetYResol","MetYResol", 50, -2., 2., "#frac{E_{y,reco}^{miss} - E_{y,gen}^{miss}}{E_{y,gen}^{miss}}");
	Tau3p_Neutrino_PtResol  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_PtResol","Tau3p_Neutrino_PtResol", 50, -1., 1., "p_{T} resol. #frac{#nu_{3p}^{reco} - #nu_{3p}^{gen}}{#nu_{3p}^{gen}}");
	Tau3p_Neutrino_PhiResol = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_PhiResol","Tau3p_Neutrino_PhiResol", 50, -1., 1., "#phi resol. #frac{#nu_{3p}^{reco} - #nu_{3p}^{gen}}{#nu_{3p}^{gen}}");
	Tau3p_Neutrino_XResol  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_XResol","Tau3p_Neutrino_XResol", 50, -1., 1., "x resol. #frac{#nu_{3p}^{reco} - #nu_{3p}^{gen}}{#nu_{3p}^{gen}}");
	Tau3p_Neutrino_YResol  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_YResol","Tau3p_Neutrino_YResol", 50, -1., 1., "y resol. #frac{#nu_{3p}^{reco} - #nu_{3p}^{gen}}{#nu_{3p}^{gen}}");
	MetMinus3pNeutrino_PtResol = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PtResol","MetMinus3pNeutrino_PtResol", 50, -2., 2., "p_{T} resol. #frac{#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen}}{#nu_{#tau,#mu}^{gen}}");
	MetMinus3pNeutrino_PhiResol = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PhiResol","MetMinus3pNeutrino_PhiResol", 50, -2., 2., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen})");
	MetMinus3pNeutrino_XResol = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_XResol","MetMinus3pNeutrino_XResol", 50, -2., 2., "x resol. #frac{#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen}}{#nu_{#tau,#mu}^{gen}}");
	MetMinus3pNeutrino_YResol = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_YResol","MetMinus3pNeutrino_YResol", 50, -2., 2., "y resol. #frac{#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen}}{#nu_{#tau,#mu}^{gen}}");

	svFitMass_Default = HConfig.GetTH1D(Name+"_svFitMass_Default",	"svFitMass_Default",	50, 0., 250., "m_{SVfit}(#tau_{h},#mu)/GeV");
	svFit3pMass = HConfig.GetTH1D(Name+"_svFit3pMass",	"svFit3pMass",	50, 0., 250., "m_{SVfit}(#tau_{3p},#mu)/GeV");;
	svFit3pMassResol = HConfig.GetTH1D(Name+"_svFit3pMassResol",	"svFit3pMassResol",	50, -1., 1., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p},#mu)");

	MetPxPull					= HConfig.GetTH1D(Name+"_MetPxPull",					"MetPxPull", 40, -3., 3., "E_{x}^{miss} pull");
	MetPyPull					= HConfig.GetTH1D(Name+"_MetPyPull",					"MetPyPull", 40, -3., 3., "E_{y}^{miss} pull");
	MetPxPyPull					= HConfig.GetTH2D(Name+"_MetPxPyPull",					"MetPxPyPull", 30, -3., 3., 30, -3., 3., "E_{x}^{miss} pull", "E_{y}^{miss} pull");
	MetMinus3pNeutrino_PxPull	= HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PxPull",	"MetMinus3pNeutrino_PxPull", 40, -3., 3., "E_{x}(#nu_{#tau,#mu}) pull");
	MetMinus3pNeutrino_PyPull	= HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PyPull",	"MetMinus3pNeutrino_PyPull", 40, -3., 3., "E_{y}(#nu_{#tau,#mu}) pull");
	MetMinus3pNeutrino_PxPyPull	= HConfig.GetTH2D(Name+"_MetMinus3pNeutrino_PxPyPull",	"MetMinus3pNeutrino_PxPyPull", 30, -3., 3.,30, -3., 3., "E_{x}(#nu_{#tau,#mu}) pull", "E_{y}(#nu_{#tau,#mu}) pull");
}

void ZeroJetLow3Prong::Store_ExtraDist(){
	Logger(Logger::Verbose) << std::endl;
	HToTaumuTauh::Store_ExtraDist();

	Extradist1d.push_back(&Tau3p_Plus_Pt);
	Extradist1d.push_back(&Tau3p_Plus_Eta);
	Extradist1d.push_back(&Tau3p_Plus_Phi);
	Extradist1d.push_back(&Tau3p_Plus_E);
	Extradist1d.push_back(&Tau3p_Minus_Pt);
	Extradist1d.push_back(&Tau3p_Minus_Eta);
	Extradist1d.push_back(&Tau3p_Minus_Phi);
	Extradist1d.push_back(&Tau3p_Minus_E);
	Extradist1d.push_back(&Tau3p_Zero_Pt);
	Extradist1d.push_back(&Tau3p_Zero_Eta);
	Extradist1d.push_back(&Tau3p_Zero_Phi);
	Extradist1d.push_back(&Tau3p_Zero_E);

	Extradist1d.push_back(&Tau3p_fit_Pt);
	Extradist1d.push_back(&Tau3p_fit_Eta);
	Extradist1d.push_back(&Tau3p_fit_Phi);
	Extradist1d.push_back(&Tau3p_fit_E);

	Extradist1d.push_back(&TauMu_fit_Pt);
	Extradist1d.push_back(&TauMu_fit_Eta);
	Extradist1d.push_back(&TauMu_fit_Phi);
	Extradist1d.push_back(&TauMu_fit_E);

	Extradist1d.push_back(&EventFit_converged);
	Extradist1d.push_back(&EventFit_Chi2);
	Extradist1d.push_back(&EventFit_nIterations);
	Extradist1d.push_back(&EventFit_deltaCSum);

	Extradist1d.push_back(&EventFit_Res_Pt);
	Extradist1d.push_back(&EventFit_Res_Eta);
	Extradist1d.push_back(&EventFit_Res_Phi);
	Extradist1d.push_back(&EventFit_Res_E);

	Extradist1d.push_back(&svFit3pMass_PlusSol);
	Extradist1d.push_back(&svFit3pMass_MinusSol);
	Extradist1d.push_back(&svFit3pMass_ZeroSol);
	Extradist1d.push_back(&svFit3pMass_TrueSol);
	Extradist1d.push_back(&svFit3pMass_TruePlusSol);
	Extradist1d.push_back(&svFit3pMass_TrueMinusSol);
	Extradist1d.push_back(&svFit3pMassResol_PlusSol);
	Extradist1d.push_back(&svFit3pMassResol_MinusSol);
	Extradist1d.push_back(&svFit3pMassResol_ZeroSol);
	Extradist1d.push_back(&svFit3pMassResol_TrueSol);

	Extradist2d.push_back(&Tau_Pt_3pTrueVsSVfit3pTrueSol);
	Extradist1d.push_back(&Tau_Pt_3pTrueSVfit3pTrueSolResol);
	Extradist2d.push_back(&Tau_Pt_3pTrueVsSVfit3pZeroSol);
	Extradist1d.push_back(&Tau_Pt_3pTrueSVfit3pZeroSolResol);

	Extradist1d.push_back(&vis3pMass_PlusSol);
	Extradist1d.push_back(&vis3pMass_MinusSol);
	Extradist1d.push_back(&vis3pMass_ZeroSol);
	Extradist1d.push_back(&vis3pMass_TrueSol);
	Extradist1d.push_back(&vis3pMassResol_PlusSol);
	Extradist1d.push_back(&vis3pMassResol_MinusSol);
	Extradist1d.push_back(&vis3pMassResol_ZeroSol);
	Extradist1d.push_back(&vis3pMassResol_TrueSol);

	Extradist1d.push_back(&Tau3p_True_Pt);
	Extradist1d.push_back(&Tau3p_True_Eta);
	Extradist1d.push_back(&Tau3p_True_Phi);
	Extradist1d.push_back(&Tau3p_True_E);
	Extradist1d.push_back(&Tau3p_True_Pt_Resol);
	Extradist1d.push_back(&Tau3p_True_Eta_Resol);
	Extradist1d.push_back(&Tau3p_True_Phi_Resol);
	Extradist1d.push_back(&Tau3p_True_E_Resol);

	Extradist1d.push_back(&TrueSignalMET);
	Extradist1d.push_back(&MetPtResol);
	Extradist1d.push_back(&MetPhiResol);
	Extradist1d.push_back(&MetPxResol);
	Extradist1d.push_back(&MetPyResol);
	Extradist1d.push_back(&Tau3p_Neutrino_PtResol);
	Extradist1d.push_back(&Tau3p_Neutrino_PhiResol);
	Extradist1d.push_back(&Tau3p_Neutrino_XResol);
	Extradist1d.push_back(&Tau3p_Neutrino_YResol);
	Extradist1d.push_back(&MetMinus3pNeutrino_PtResol);
	Extradist1d.push_back(&MetMinus3pNeutrino_PhiResol);
	Extradist1d.push_back(&MetMinus3pNeutrino_XResol);
	Extradist1d.push_back(&MetMinus3pNeutrino_YResol);

	Extradist1d.push_back(&svFitMass_Default);
	Extradist1d.push_back(&svFit3pMass);
	Extradist1d.push_back(&svFit3pMassResol);

	Extradist1d.push_back(&MetPxPull);
	Extradist1d.push_back(&MetPyPull);
	Extradist2d.push_back(&MetPxPyPull);
	Extradist1d.push_back(&MetMinus3pNeutrino_PxPull);
	Extradist1d.push_back(&MetMinus3pNeutrino_PyPull);
	Extradist2d.push_back(&MetMinus3pNeutrino_PxPyPull);
}

void ZeroJetLow3Prong::doEvent(){
	Logger(Logger::Verbose) << std::endl;
	HToTaumuTauh::doEvent(false); // run HToTaumuTauh::doEvent without filling framework plots

	// todo: split analysis part from plotting part in HToTaumuTauh

	Logger(Logger::Debug) << "Cut: Tau decay mode" << std::endl;
	if (selTau < 0){
		// set cuts to true for nice N-0 and N-1 plots
		value.at(CatCut3) = -10;
		pass.at(CatCut3) = true;
	}
	else{
		bool foundSignalTauh(false);
		bool foundSignalTaumu(false);
		for(int i_tau = 0; i_tau < Ntp->NMCTaus(); i_tau++){
			if(Ntp->MCTau_JAK(i_tau) == 5) foundSignalTauh = true;
			if(Ntp->MCTau_JAK(i_tau) == 2) foundSignalTaumu = true;
		}
		if( Ntp->PFTau_hpsDecayMode(selTau) == 10 && (not foundSignalTauh || not foundSignalTaumu)){
			value.at(CatCut3) = 7; // remove all events which are not tau->mu tau->3pi
		}
		else value.at(CatCut3)=Ntp->PFTau_hpsDecayMode(selTau);
		pass.at(CatCut3)=(value.at(CatCut3)==cut.at(CatCut3));
	}

	Logger(Logger::Debug) << "Cut: SV significance" << std::endl;
	if (selTau < 0 || not pass.at(CatCut3)) {
		value.at(CatCut4) = -999;
		pass.at(CatCut4) = true;
	}
	else if ( not Ntp->PFTau_TIP_hassecondaryVertex(selTau) ){
		value.at(CatCut4) = -9;
		pass.at(CatCut4) = false;
	}
	else {
		// set sign of flight length significance by projection of tau momentum direction
		// on fitted PV-SV direction
		int sign = ( Ntp->PFTau_FlightLength3d(selTau).Dot( Ntp->PFTau_3PS_A1_LV(selTau).Vect()) > 0 ) ? +1 : -1;
		value.at(CatCut4) = sign * Ntp->PFTau_FlightLength_significance(selTau);
		pass.at(CatCut4) = ( value.at(CatCut4) >= cut.at(CatCut4) );
	}

	double wobs = 1;
	bool status = AnalysisCuts(t,w,wobs); // true only if full selection passed

	TPTRObject TPResults;
	if (status){
		LorentzVectorParticle A1 = Ntp->PFTau_a1_lvp(selTau);
		TrackParticle MuonTP = Ntp->Muon_TrackParticle(selMuon);
		TVector3 PV = Ntp->PFTau_TIP_primaryVertex_pos(selTau);
		TMatrixTSym<double> PVCov = Ntp->PFTau_TIP_primaryVertex_cov(selTau);

		TLorentzVector Recoil;
		for(unsigned i; i<Ntp->Vtx_nTrk(selVertex); i++){
			Recoil += Ntp->Vtx_TracksP4(selVertex, i);
		}
		Recoil -= Ntp->Muon_p4(selMuon);
		Recoil -= Ntp->PFTau_p4(selTau);
		double Phi_Res = (Recoil.Phi() > 0) ? Recoil.Phi() - TMath::Pi() : Recoil.Phi() + TMath::Pi();

		//Logger::Instance()->setLevelForClass("GlobalEventFit", Logger::Debug);
		GlobalEventFit GEF(MuonTP, A1, Phi_Res, PV, PVCov);
		TPResults = GEF.getTPTRObject();
		LorentzVectorParticle test = TPResults.getA1();
		Logger(Logger::Debug) << test.LV().E() << std::endl;
		GEFObject Results = GEF.Fit();
		Logger(Logger::Debug) << Results.getChi2() << std::endl;

		// fill plots
		if (TPResults.isAmbiguous()) {
			Tau3p_Plus_Pt.at(t).Fill(TPResults.getTauPlus().LV().Pt(), w);
			Tau3p_Plus_Eta.at(t).Fill(TPResults.getTauPlus().LV().Eta(), w);
			Tau3p_Plus_Phi.at(t).Fill(TPResults.getTauPlus().LV().Phi(), w);
			Tau3p_Plus_E.at(t).Fill(TPResults.getTauPlus().LV().E(), w);
			Tau3p_Minus_Pt.at(t).Fill(TPResults.getTauMinus().LV().Pt(), w);
			Tau3p_Minus_Eta.at(t).Fill(TPResults.getTauMinus().LV().Eta(), w);
			Tau3p_Minus_Phi.at(t).Fill(TPResults.getTauMinus().LV().Phi(), w);
			Tau3p_Minus_E.at(t).Fill(TPResults.getTauMinus().LV().E(), w);
		}
		else {
			Tau3p_Zero_Pt.at(t).Fill(TPResults.getTauZero().LV().Pt(), w);
			Tau3p_Zero_Eta.at(t).Fill(TPResults.getTauZero().LV().Eta(), w);
			Tau3p_Zero_Phi.at(t).Fill(TPResults.getTauZero().LV().Phi(), w);
			Tau3p_Zero_E.at(t).Fill(TPResults.getTauZero().LV().E(), w);
		}

		Tau3p_fit_Pt.at(t).Fill(Results.getTauH().LV().Pt(), w);
		Tau3p_fit_Eta.at(t).Fill(Results.getTauH().LV().Eta(), w);
		Tau3p_fit_Phi.at(t).Fill(Results.getTauH().LV().Phi(), w);
		Tau3p_fit_E.at(t).Fill(Results.getTauH().LV().E(), w);

		TauMu_fit_Pt.at(t).Fill(Results.getTauMu().LV().Pt(), w);
		TauMu_fit_Eta.at(t).Fill(Results.getTauMu().LV().Eta(), w);
		TauMu_fit_Phi.at(t).Fill(Results.getTauMu().LV().Phi(), w);
		TauMu_fit_E.at(t).Fill(Results.getTauMu().LV().E(), w);

		EventFit_converged.at(t).Fill(Results.Fitconverged(), w);
		EventFit_Chi2.at(t).Fill(Results.getChi2(), w);
		EventFit_nIterations.at(t).Fill(Results.getNiterations(), w);
		EventFit_deltaCSum.at(t).Fill(Results.getCsum(), w);

		EventFit_Res_Pt.at(t).Fill(Results.getResonance().LV().Pt(), w);
		EventFit_Res_Eta.at(t).Fill(Results.getResonance().LV().Eta(), w);
		EventFit_Res_Phi.at(t).Fill(Results.getResonance().LV().Phi(), w);
		EventFit_Res_E.at(t).Fill(Results.getResonance().LV().E(), w);
	}

	// calculate SVFit stuff for comparison
	if (status) {
		double trueHiggsMass = Ntp->getResonanceMassFromGenInfo(false, true, false);
		int i_matchedMCTau = Ntp->matchTauTruth(selTau, true);
		int i_otherTau(-1);
		if (i_matchedMCTau == 0) i_otherTau = 1;
		else if (i_matchedMCTau == 1) i_otherTau = 0;
		else Logger(Logger::Warning) << "Cannot determine index of tau_mu." << std::endl;
		int true3ProngAmbig = (i_matchedMCTau >= 0) ? Ntp->MCTau_true3prongAmbiguity( i_matchedMCTau ) : -9;
		bool isSignalTauDecay = (i_matchedMCTau >= 0 && true3ProngAmbig >= 0);
		double ptSVFittrue(-9), ptSVFitZero(-9);
		TLorentzVector trueTauP4 = isSignalTauDecay ? Ntp->MCTau_p4(i_matchedMCTau) : TLorentzVector(0,0,0,0);

		TLorentzVector trueMet(0,0,0,0);
		for(int i_tau = 0; i_tau < Ntp->NMCTaus(); i_tau++){
			if(i_tau == 2) Logger(Logger::Error) << "More than 2 signal taus in event!" << std::endl;
			trueMet += Ntp->MCTau_invisiblePart(i_tau);
		}

		TrueSignalMET.at(t).Fill(trueMet.Pt(), w);
		MetPtResol.at(t).Fill( (Ntp->MET_CorrMVAMuTau_et() - trueMet.Pt())/trueMet.Pt() , w);
		MetPhiResol.at(t).Fill( Tools::DeltaPhi(Ntp->MET_CorrMVAMuTau_phi(),trueMet.Phi()) , w);
		MetPxResol.at(t).Fill( (Ntp->MET_CorrMVAMuTau_ex() - trueMet.Px())/trueMet.Px() , w);
		MetPyResol.at(t).Fill( (Ntp->MET_CorrMVAMuTau_ey() - trueMet.Py())/trueMet.Py() , w);
		if(Ntp->MET_CorrMVAMuTau_significance_xx() != 0 && Ntp->MET_CorrMVAMuTau_significance_yy() != 0){
			double pxPull = (Ntp->MET_CorrMVAMuTau_ex() - trueMet.Px()) / sqrt(Ntp->MET_CorrMVAMuTau_significance_xx());
			double pyPull = (Ntp->MET_CorrMVAMuTau_ey() - trueMet.Py()) / sqrt(Ntp->MET_CorrMVAMuTau_significance_yy());
			MetPxPull.at(t).Fill( pxPull , w);
			MetPyPull.at(t).Fill( pyPull , w);
			MetPxPyPull.at(t).Fill( pxPull, pyPull, w);
		}
		else
			Logger(Logger::Warning) << "MET covariance is zero!" << std::endl;

		SVFitObject *svfObjDefault = Ntp->getSVFitResult_MuTauh(svfitstorage, "CorrMVAMuTau", selMuon, selTau, 50000);
		svFitMass_Default.at(t).Fill(svfObjDefault->get_mass(), w);

		if (TPResults.isAmbiguous()) {
			SVFitObject *svfObjPlus  = Ntp->getSVFitResult_MuTau3p(svfitstor_Tau3pPlus, "CorrMVAMuTau", selMuon, TPResults.getTauPlus().LV(), TPResults.getNeutrinoPlus(), "3ProngRecoPlus");
			SVFitObject *svfObjMinus = Ntp->getSVFitResult_MuTau3p(svfitstor_Tau3pMinus, "CorrMVAMuTau", selMuon, TPResults.getTauMinus().LV(), TPResults.getNeutrinoMinus(), "3ProngRecoMinus");

			svFit3pMass_PlusSol.at(t).Fill(svfObjPlus->get_mass(), w);
			svFit3pMass_MinusSol.at(t).Fill(svfObjMinus->get_mass(), w);
			svFit3pMassResol_PlusSol.at(t).Fill( (svfObjPlus->get_mass() - trueHiggsMass) / trueHiggsMass, w);
			svFit3pMassResol_MinusSol.at(t).Fill((svfObjMinus->get_mass() - trueHiggsMass) / trueHiggsMass, w);
			double visMassPlus = (Ntp->Muon_p4(selMuon) + TPResults.getTauPlus().LV()).M();
			double visMassMinus = (Ntp->Muon_p4(selMuon) + TPResults.getTauMinus().LV()).M();
			vis3pMass_PlusSol.at(t).Fill( visMassPlus, w);
			vis3pMass_MinusSol.at(t).Fill( visMassMinus, w);
			vis3pMassResol_PlusSol.at(t).Fill( (visMassPlus - trueHiggsMass)/trueHiggsMass, w);
			vis3pMassResol_MinusSol.at(t).Fill( (visMassMinus - trueHiggsMass)/trueHiggsMass, w);
			objects::MET metPlus(Ntp, "CorrMVAMuTau");
			metPlus.subtractNeutrino(TPResults.getNeutrinoPlus());
			objects::MET metMinus(Ntp, "CorrMVAMuTau");
			metMinus.subtractNeutrino(TPResults.getNeutrinoMinus());
			switch (true3ProngAmbig) {
				case MultiProngTauSolver::plus:
					svFit3pMass_TrueSol.at(t).Fill(svfObjPlus->get_mass(), w);
					svFit3pMass_TruePlusSol.at(t).Fill(svfObjPlus->get_mass(), w);
					svFit3pMassResol_TrueSol.at(t).Fill((svfObjPlus->get_mass() - trueHiggsMass) / trueHiggsMass, w);
					vis3pMass_TrueSol.at(t).Fill((Ntp->Muon_p4(selMuon) + TPResults.getTauPlus().LV()).M(), w);
					vis3pMassResol_TrueSol.at(t).Fill((visMassPlus - trueHiggsMass)/trueHiggsMass, w);
					svFit3pMass.at(t).Fill(svfObjPlus->get_mass(), w);
					svFit3pMassResol.at(t).Fill((svfObjPlus->get_mass() - trueHiggsMass) / trueHiggsMass, w);
					if(svfObjPlus->get_fittedTauLeptons().size() > 1) ptSVFittrue = svfObjPlus->get_fittedTauLeptons().at(1).Pt();

					Tau3p_True_Pt.at(t).Fill(TPResults.getTauPlus().LV().Pt(), w);
					Tau3p_True_Eta.at(t).Fill(TPResults.getTauPlus().LV().Eta(), w);
					Tau3p_True_Phi.at(t).Fill(TPResults.getTauPlus().LV().Phi(), w);
					Tau3p_True_E.at(t).Fill(TPResults.getTauPlus().LV().E(), w);
					if (isSignalTauDecay){
						Tau3p_True_Pt_Resol.at(t).Fill( (TPResults.getTauPlus().LV().Pt() - trueTauP4.Pt())/trueTauP4.Pt() , w);
						Tau3p_True_Eta_Resol.at(t).Fill((TPResults.getTauPlus().LV().Eta() - trueTauP4.Eta())/trueTauP4.Eta(), w);
						Tau3p_True_Phi_Resol.at(t).Fill((TPResults.getTauPlus().LV().Phi() - trueTauP4.Phi())/trueTauP4.Phi(), w);
						Tau3p_True_E_Resol.at(t).Fill((TPResults.getTauPlus().LV().E() - trueTauP4.E())/trueTauP4.E(), w);
					}
					if (Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt() > 0.001){
						Tau3p_Neutrino_PtResol.at(t).Fill( (TPResults.getNeutrinoPlus().LV().Pt() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt(), w);
						Tau3p_Neutrino_XResol.at(t).Fill( (TPResults.getNeutrinoPlus().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Px(), w);
						Tau3p_Neutrino_YResol.at(t).Fill( (TPResults.getNeutrinoPlus().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Py(), w);
						Tau3p_Neutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoPlus().LV().Phi(), Ntp->MCTau_invisiblePart(i_matchedMCTau).Phi()), w);
					}
					if (Ntp->MCTau_invisiblePart(i_otherTau).Pt() > 0.001){
						MetMinus3pNeutrino_PtResol.at(t).Fill( (metPlus.et() - Ntp->MCTau_invisiblePart(i_otherTau).Pt())/Ntp->MCTau_invisiblePart(i_otherTau).Pt(), w);
						MetMinus3pNeutrino_XResol.at(t).Fill( (metPlus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px())/Ntp->MCTau_invisiblePart(i_otherTau).Px(), w);
						MetMinus3pNeutrino_YResol.at(t).Fill( (metPlus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py())/Ntp->MCTau_invisiblePart(i_otherTau).Py(), w);
						MetMinus3pNeutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
					}
					if(metPlus.significanceXX() != 0 && metPlus.significanceYY() != 0){
						double pxPull = (metPlus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px()) / sqrt(metPlus.significanceXX());
						double pyPull = (metPlus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py()) / sqrt(metPlus.significanceYY());
						MetMinus3pNeutrino_PxPull.at(t).Fill( pxPull , w);
						MetMinus3pNeutrino_PyPull.at(t).Fill( pyPull , w);
						MetMinus3pNeutrino_PxPyPull.at(t).Fill( pxPull, pyPull, w);
					}
					else
						Logger(Logger::Warning) << "MET-plus covariance is zero!" << std::endl;
					break;
				case MultiProngTauSolver::minus:
					svFit3pMass_TrueSol.at(t).Fill(svfObjMinus->get_mass(), w);
					svFit3pMass_TrueMinusSol.at(t).Fill(svfObjMinus->get_mass(), w);
					svFit3pMassResol_TrueSol.at(t).Fill((svfObjMinus->get_mass() - trueHiggsMass) / trueHiggsMass, w);
					vis3pMass_TrueSol.at(t).Fill((Ntp->Muon_p4(selMuon) + TPResults.getTauMinus().LV()).M(), w);
					vis3pMassResol_TrueSol.at(t).Fill((visMassMinus - trueHiggsMass)/trueHiggsMass, w);
					svFit3pMass.at(t).Fill(svfObjMinus->get_mass(), w);
					svFit3pMassResol.at(t).Fill((svfObjMinus->get_mass() - trueHiggsMass) / trueHiggsMass, w);
					if(svfObjPlus->get_fittedTauLeptons().size() > 1) ptSVFittrue = svfObjMinus->get_fittedTauLeptons().at(1).Pt();

					Tau3p_True_Pt.at(t).Fill(TPResults.getTauMinus().LV().Pt(), w);
					Tau3p_True_Eta.at(t).Fill(TPResults.getTauMinus().LV().Eta(), w);
					Tau3p_True_Phi.at(t).Fill(TPResults.getTauMinus().LV().Phi(), w);
					Tau3p_True_E.at(t).Fill(TPResults.getTauMinus().LV().E(), w);
					if(isSignalTauDecay){
						Tau3p_True_Pt_Resol.at(t).Fill( (TPResults.getTauMinus().LV().Pt() - trueTauP4.Pt())/trueTauP4.Pt() , w);
						Tau3p_True_Eta_Resol.at(t).Fill((TPResults.getTauMinus().LV().Eta() - trueTauP4.Eta())/trueTauP4.Eta(), w);
						Tau3p_True_Phi_Resol.at(t).Fill((TPResults.getTauMinus().LV().Phi() - trueTauP4.Phi())/trueTauP4.Phi(), w);
						Tau3p_True_E_Resol.at(t).Fill((TPResults.getTauMinus().LV().E() - trueTauP4.E())/trueTauP4.E(), w);
					}
					if (Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt() > 0.001){
						Tau3p_Neutrino_PtResol.at(t).Fill( (TPResults.getNeutrinoMinus().LV().Pt() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt(), w);
						Tau3p_Neutrino_XResol.at(t).Fill( (TPResults.getNeutrinoMinus().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Px(), w);
						Tau3p_Neutrino_YResol.at(t).Fill( (TPResults.getNeutrinoMinus().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Py(), w);
						Tau3p_Neutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoMinus().LV().Phi(), Ntp->MCTau_invisiblePart(i_matchedMCTau).Phi()), w);
					}
					if (Ntp->MCTau_invisiblePart(i_otherTau).Pt() > 0.001){
						MetMinus3pNeutrino_PtResol.at(t).Fill( (metMinus.et() - Ntp->MCTau_invisiblePart(i_otherTau).Pt())/Ntp->MCTau_invisiblePart(i_otherTau).Pt(), w);
						MetMinus3pNeutrino_XResol.at(t).Fill( (metMinus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px())/Ntp->MCTau_invisiblePart(i_otherTau).Px(), w);
						MetMinus3pNeutrino_YResol.at(t).Fill( (metMinus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py())/Ntp->MCTau_invisiblePart(i_otherTau).Py(), w);
						MetMinus3pNeutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
					}
					if(metMinus.significanceXX() != 0 && metMinus.significanceYY() != 0){
						double pxPull = (metMinus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px()) / sqrt(metMinus.significanceXX());
						double pyPull = (metMinus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py()) / sqrt(metMinus.significanceYY());
						MetMinus3pNeutrino_PxPull.at(t).Fill( pxPull , w);
						MetMinus3pNeutrino_PyPull.at(t).Fill( pyPull , w);
						MetMinus3pNeutrino_PxPyPull.at(t).Fill( pxPull, pyPull, w);
					}
					else
						Logger(Logger::Warning) << "MET-minus covariance is zero!" << std::endl;
					break;
				case -9:
					Logger(Logger::Verbose) << "HPS tau not matched to signal tau decay." << std::endl;
					break;
				default:
					Logger(Logger::Warning) << "True ambiguity value is " << true3ProngAmbig << ", which does not make sense." << std::endl;
					break;
			}

		}
		else {
			SVFitObject *svfObjZero  = Ntp->getSVFitResult_MuTau3p(svfitstor_Tau3pZero, "CorrMVAMuTau", selMuon, TPResults.getTauZero().LV(), TPResults.getNeutrinoZero(), "3ProngRecoZero");

			svFit3pMass_ZeroSol.at(t).Fill(svfObjZero->get_mass(), w);
			svFit3pMassResol_ZeroSol.at(t).Fill((svfObjZero->get_mass() - trueHiggsMass) / trueHiggsMass, w);
			double visMassZero = (Ntp->Muon_p4(selMuon) + TPResults.getTauZero().LV()).M();
			vis3pMass_ZeroSol.at(t).Fill(visMassZero, w);
			vis3pMassResol_ZeroSol.at(t).Fill((visMassZero - trueHiggsMass)/trueHiggsMass, w);
			svFit3pMass.at(t).Fill(svfObjZero->get_mass(), w);
			svFit3pMassResol.at(t).Fill((svfObjZero->get_mass() - trueHiggsMass) / trueHiggsMass, w);

			if(svfObjZero->get_fittedTauLeptons().size() > 1) ptSVFitZero = svfObjZero->get_fittedTauLeptons().at(1).Pt();\

			if (Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt() > 0.001)
				Tau3p_Neutrino_PtResol.at(t).Fill( (TPResults.getNeutrinoZero().LV().Pt() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt(), w);
			Tau3p_Neutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoZero().LV().Phi(), Ntp->MCTau_invisiblePart(i_matchedMCTau).Phi()), w);
			objects::MET metZero(Ntp, "CorrMVAMuTau");
			metZero.subtractNeutrino(TPResults.getNeutrinoZero());
			if (Ntp->MCTau_invisiblePart(i_otherTau).Pt() > 0.001){
				MetMinus3pNeutrino_PtResol.at(t).Fill( (metZero.et() - Ntp->MCTau_invisiblePart(i_otherTau).Pt())/Ntp->MCTau_invisiblePart(i_otherTau).Pt(), w);
				MetMinus3pNeutrino_XResol.at(t).Fill( (metZero.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px())/Ntp->MCTau_invisiblePart(i_otherTau).Px(), w);
				MetMinus3pNeutrino_YResol.at(t).Fill( (metZero.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py())/Ntp->MCTau_invisiblePart(i_otherTau).Py(), w);
				MetMinus3pNeutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
			}
			if(metZero.significanceXX() != 0 && metZero.significanceYY() != 0){
				double pxPull = (metZero.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px()) / sqrt(metZero.significanceXX());
				double pyPull = (metZero.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py()) / sqrt(metZero.significanceYY());
				MetMinus3pNeutrino_PxPull.at(t).Fill( pxPull , w);
				MetMinus3pNeutrino_PyPull.at(t).Fill( pyPull , w);
				MetMinus3pNeutrino_PxPyPull.at(t).Fill( pxPull, pyPull, w);
			}
			else
				Logger(Logger::Warning) << "MET-zero covariance is zero!" << std::endl;
		}

		if(isSignalTauDecay){
			if(ptSVFittrue > -1){
				Tau_Pt_3pTrueVsSVfit3pTrueSol.at(t).Fill(trueTauP4.Pt(), ptSVFittrue, w);
				Tau_Pt_3pTrueSVfit3pTrueSolResol.at(t).Fill( (trueTauP4.Pt() - ptSVFittrue)/ptSVFittrue, w);
			}
			if(ptSVFitZero > -1){
				Tau_Pt_3pTrueVsSVfit3pZeroSol.at(t).Fill(trueTauP4.Pt(), ptSVFitZero, w);
				Tau_Pt_3pTrueSVfit3pZeroSolResol.at(t).Fill( (trueTauP4.Pt() - ptSVFitZero)/ptSVFitZero, w);
			}
		}
	}

}
