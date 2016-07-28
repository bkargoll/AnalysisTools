/*
 * ZeroJet3Prong.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "ZeroJet3Prong.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "../Tools.h"

#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/GlobalEventFit.h"

ZeroJet3Prong::ZeroJet3Prong(TString Name_, TString id_):
	Category(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class ZeroJet3Prong" << std::endl;
	// run ZeroJet3Prong category
	categoryFlag = "ZeroJet3Prong";

	runSVFit_ = false;
}

ZeroJet3Prong::~ZeroJet3Prong() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

void  ZeroJet3Prong::categoryConfiguration(){
	// Setup cut Values for this category
	cut.at(NJet)		= 0;
	cut.at(DecayMode)	= 10;	// decay mode finding
	cut.at(SigmaSV)		= 3.0;	// SV significance

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(NJet)="Number Jets $<=$";
	title.at(NJet)+=cut.at(NJet);
	htitle=title.at(NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=NJet;
	Nminus1.at(NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(DecayMode)="decayMode $==$";
	title.at(DecayMode)+=cut.at(DecayMode);
	htitle=title.at(DecayMode);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="#tau_{h} decay mode";
	c="_Cut_";c+=DecayMode;
	Nminus1.at(DecayMode) = HConfig.GetTH1D(Name+c+"_Nminus1_DecayMode_",htitle,13,-0.5,12.5,hlabel,"Events");
	Nminus0.at(DecayMode) = HConfig.GetTH1D(Name+c+"_Nminus0_DecayMode_",htitle,13,-0.5,12.5,hlabel,"Events");

	title.at(SigmaSV)="$\\sigma(SV) >=$";
	title.at(SigmaSV)+=cut.at(SigmaSV);
	htitle=title.at(SigmaSV);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="\\sigma(SV)";
	c="_Cut_";c+=SigmaSV;
	Nminus1.at(SigmaSV) = HConfig.GetTH1D(Name+c+"_Nminus1_SigmaSV_",htitle,80,-10.,30.,hlabel,"Events");
	Nminus0.at(SigmaSV) = HConfig.GetTH1D(Name+c+"_Nminus0_SigmaSV_",htitle,80,-10.,30.,hlabel,"Events");

	// setup additional histograms
	svFit3pMass_PlusSol		= HConfig.GetTH1D(Name+"_svFit3pMass_PlusSol",	"svFit3pMass_PlusSol",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{+},#mu)/GeV");
	svFit3pMass_MinusSol	= HConfig.GetTH1D(Name+"_svFit3pMass_MinusSol",	"svFit3pMass_MinusSol",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{-},#mu)/GeV");
	svFit3pMass_ZeroSol		= HConfig.GetTH1D(Name+"_svFit3pMass_ZeroSol",	"svFit3pMass_ZeroSol",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{0},#mu)/GeV");
	svFit3pMass_TrueSol		= HConfig.GetTH1D(Name+"_svFit3pMass_TrueSol",	"svFit3pMass_TrueSol",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{true},#mu)/GeV");
	svFit3pMass_TruePlusSol		= HConfig.GetTH1D(Name+"_svFit3pMass_TruePlusSol",	"svFit3pMass_TruePlusSol",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{+true},#mu)/GeV");
	svFit3pMass_TrueMinusSol	= HConfig.GetTH1D(Name+"_svFit3pMass_TrueMinusSol",	"svFit3pMass_TrueMinusSol",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{-true},#mu)/GeV");
	svFit3pMassResol_PlusSol	= HConfig.GetTH1D(Name+"_svFit3pMassResol_PlusSol",	"svFit3pMassResol_PlusSol",	50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{+},#mu)");
	svFit3pMassResol_MinusSol	= HConfig.GetTH1D(Name+"_svFit3pMassResol_MinusSol","svFit3pMassResol_MinusSol",50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{-},#mu)");
	svFit3pMassResol_ZeroSol	= HConfig.GetTH1D(Name+"_svFit3pMassResol_ZeroSol",	"svFit3pMassResol_ZeroSol",	50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{0},#mu)");
	svFit3pMassResol_TrueSol	= HConfig.GetTH1D(Name+"_svFit3pMassResol_TrueSol",	"svFit3pMassResol_TrueSol",	50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{true},#mu)");

	svFit3pMass_PlusSol_notCleaned		= HConfig.GetTH1D(Name+"_svFit3pMass_PlusSol_notCleaned",	"svFit3pMass_PlusSol_notCleaned",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{+},#mu)/GeV");
	svFit3pMass_MinusSol_notCleaned		= HConfig.GetTH1D(Name+"_svFit3pMass_MinusSol_notCleaned",	"svFit3pMass_MinusSol_notCleaned",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{-},#mu)/GeV");
	svFit3pMass_ZeroSol_notCleaned		= HConfig.GetTH1D(Name+"_svFit3pMass_ZeroSol_notCleaned",	"svFit3pMass_ZeroSol_notCleaned",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{0},#mu)/GeV");
	svFit3pMass_TrueSol_notCleaned		= HConfig.GetTH1D(Name+"_svFit3pMass_TrueSol_notCleaned",	"svFit3pMass_TrueSol_notCleaned",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{true},#mu)/GeV");
	svFit3pMass_TruePlusSol_notCleaned	= HConfig.GetTH1D(Name+"_svFit3pMass_TruePlusSol_notCleaned",	"svFit3pMass_TruePlusSol_notCleaned",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{+true},#mu)/GeV");
	svFit3pMass_TrueMinusSol_notCleaned	= HConfig.GetTH1D(Name+"_svFit3pMass_TrueMinusSol_notCleaned",	"svFit3pMass_TrueMinusSol_notCleaned",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{-true},#mu)/GeV");
	svFit3pMassResol_PlusSol_notCleaned	= HConfig.GetTH1D(Name+"_svFit3pMassResol_PlusSol_notCleaned",	"svFit3pMassResol_PlusSol_notCleaned",	50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{+},#mu)");
	svFit3pMassResol_MinusSol_notCleaned= HConfig.GetTH1D(Name+"_svFit3pMassResol_MinusSol_notCleaned","svFit3pMassResol_MinusSol_notCleaned",50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{-},#mu)");
	svFit3pMassResol_ZeroSol_notCleaned	= HConfig.GetTH1D(Name+"_svFit3pMassResol_ZeroSol_notCleaned",	"svFit3pMassResol_ZeroSol_notCleaned",	50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{0},#mu)");
	svFit3pMassResol_TrueSol_notCleaned	= HConfig.GetTH1D(Name+"_svFit3pMassResol_TrueSol_notCleaned",	"svFit3pMassResol_TrueSol_notCleaned",	50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{true},#mu)");

	svFit3pMass_PlusSol_invalid			= HConfig.GetTH1D(Name+"_svFit3pMass_PlusSol_invalid",	"svFit3pMass_PlusSol_invalid",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{+},#mu)/GeV");
	svFit3pMass_MinusSol_invalid		= HConfig.GetTH1D(Name+"_svFit3pMass_MinusSol_invalid",	"svFit3pMass_MinusSol_invalid",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{-},#mu)/GeV");
	svFit3pMass_ZeroSol_invalid			= HConfig.GetTH1D(Name+"_svFit3pMass_ZeroSol_invalid",	"svFit3pMass_ZeroSol_invalid",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{0},#mu)/GeV");
	svFit3pMass_TrueSol_invalid			= HConfig.GetTH1D(Name+"_svFit3pMass_TrueSol_invalid",	"svFit3pMass_TrueSol_invalid",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{true},#mu)/GeV");
	svFit3pMass_TruePlusSol_invalid		= HConfig.GetTH1D(Name+"_svFit3pMass_TruePlusSol_invalid",	"svFit3pMass_TruePlusSol_invalid",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{+true},#mu)/GeV");
	svFit3pMass_TrueMinusSol_invalid	= HConfig.GetTH1D(Name+"_svFit3pMass_TrueMinusSol_invalid",	"svFit3pMass_TrueMinusSol_invalid",	80, 0., 400., "m_{SVfit}(#tau_{3p}^{-true},#mu)/GeV");
	svFit3pMassResol_PlusSol_invalid	= HConfig.GetTH1D(Name+"_svFit3pMassResol_PlusSol_invalid",	"svFit3pMassResol_PlusSol_invalid",	50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{+},#mu)");
	svFit3pMassResol_MinusSol_invalid	= HConfig.GetTH1D(Name+"_svFit3pMassResol_MinusSol_invalid","svFit3pMassResol_MinusSol_invalid",50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{-},#mu)");
	svFit3pMassResol_ZeroSol_invalid	= HConfig.GetTH1D(Name+"_svFit3pMassResol_ZeroSol_invalid",	"svFit3pMassResol_ZeroSol_invalid",	50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{0},#mu)");
	svFit3pMassResol_TrueSol_invalid	= HConfig.GetTH1D(Name+"_svFit3pMassResol_TrueSol_invalid",	"svFit3pMassResol_TrueSol_invalid",	50, -1., 3., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p}^{true},#mu)");

	svFit3p_MassLMax_TrueSol= HConfig.GetTH1D(Name+"_svFit3p_MassLMax_TrueSol",	"svFit3p_MassLMax_TrueSol",	100, 0., 1000000., "L^{m}_{max} (true sol.)");
	svFit3p_MassLMax_WrongSol= HConfig.GetTH1D(Name+"_svFit3p_MassLMax_WrongSol",	"svFit3p_MassLMax_WrongSol", 100, 0., 1000000., "L^{m}_{max} (false sol.)");
	svFit3p_MassLMax_ZeroSol= HConfig.GetTH1D(Name+"_svFit3p_MassLMax_ZeroSol",	"svFit3p_MassLMax_ZeroSol",	100, 0., 1000000., "L^{m}_{max} (zero sol.)");
	svFit3p_MassLMax_PlusSol= HConfig.GetTH1D(Name+"_svFit3p_MassLMax_PlusSol",	"svFit3p_MassLMax_PlusSol",	100, 0., 1000000., "L^{m}_{max} (plus sol.)");
	svFit3p_MassLMax_MinusSol= HConfig.GetTH1D(Name+"_svFit3p_MassLMax_MinusSol",	"svFit3p_MassLMax_MinusSol", 100, 0., 1000000., "L^{m}_{max} (minus sol.)");
	svFit3p_TrueSolByHigherLMax= HConfig.GetTH1D(Name+"_svFit3p_TrueSolByHigherLMax",	"svFit3p_TrueSolByHigherLMax", 2, -0.5, 1.5, "correct by L^{m}_{max}");
	svFit3p_SolutionMatrixByLMax = HConfig.GetTH2D(Name+"_svFit3p_SolutionMatrixByLMax", "svFit3p_SolutionMatrixByLMax", 3, -0.5, 2.5, 3, -0.5, 2.5, "True amb.", "higher L^{max}(SVfit)");

	invalidSvFit3p_Tau3pNeutrinoXPull        = HConfig.GetTH1D(Name+"_invalidSvFit3p_Tau3pNeutrinoXPull", "invalidSvFit3p_Tau3pNeutrinoXPull", 50, -5., 5., "x(#nu_{#tau}) pull");
	invalidSvFit3p_Tau3pNeutrinoYPull        = HConfig.GetTH1D(Name+"_invalidSvFit3p_Tau3pNeutrinoYPull", "invalidSvFit3p_Tau3pNeutrinoYPull", 50, -5., 5., "y(#nu_{#tau}) pull");
	invalidSvFit3p_Tau3pNeutrinoPhiPull      = HConfig.GetTH1D(Name+"_invalidSvFit3p_Tau3pNeutrinoPhiPull", "invalidSvFit3p_Tau3pNeutrinoPhiPull", 50, -5., 5., "#phi(#nu_{#tau}) pull");
	invalidSvFit3p_MetMinus3pNeutrinoXPull   = HConfig.GetTH1D(Name+"_invalidSvFit3p_MetMinus3pNeutrinoXPull", "invalidSvFit3p_MetMinus3pNeutrinoXPull", 50, -5., 5., "E_{x}(#nu_{#tau,#mu}) pull");
	invalidSvFit3p_MetMinus3pNeutrinoYPull   = HConfig.GetTH1D(Name+"_invalidSvFit3p_MetMinus3pNeutrinoYPull", "invalidSvFit3p_MetMinus3pNeutrinoYPull", 50, -5., 5., "E_{y}(#nu_{#tau,#mu}) pull");
	invalidSvFit3p_MetMinus3pNeutrinoPhiPull = HConfig.GetTH1D(Name+"_invalidSvFit3p_MetMinus3pNeutrinoPhiPull", "invalidSvFit3p_MetMinus3pNeutrinoPhiPull", 50, -5., 5., "#phi(E_{T}(#nu_{#tau,#mu})) pull");

	vis3pMass_PlusSol		= HConfig.GetTH1D(Name+"_vis3pMass_PlusSol",	"vis3pMass_PlusSol",	80, 0., 400., "m_{vis}(#tau_{3p}^{+},#mu)/GeV");
	vis3pMass_MinusSol		= HConfig.GetTH1D(Name+"_vis3pMass_MinusSol",	"vis3pMass_MinusSol",	80, 0., 400., "m_{vis}(#tau_{3p}^{-},#mu)/GeV");
	vis3pMass_ZeroSol		= HConfig.GetTH1D(Name+"_vis3pMass_ZeroSol",	"vis3pMass_ZeroSol",	80, 0., 400., "m_{vis}(#tau_{3p}^{0},#mu)/GeV");
	vis3pMass_TrueSol		= HConfig.GetTH1D(Name+"_vis3pMass_TrueSol",	"vis3pMass_TrueSol",	80, 0., 400., "m_{vis}(#tau_{3p}^{true},#mu)/GeV");
	vis3pMassResol_PlusSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_PlusSol",	"vis3pMassResol_PlusSol",	50, -1., 3., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{+},#mu)");
	vis3pMassResol_MinusSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_MinusSol","vis3pMassResol_MinusSol",	50, -1., 3., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{-},#mu)");
	vis3pMassResol_ZeroSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_ZeroSol",	"vis3pMassResol_ZeroSol",	50, -1., 3., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{0},#mu)");
	vis3pMassResol_TrueSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_TrueSol",	"vis3pMassResol_TrueSol",	50, -1., 3., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{true},#mu)");

	Tau3p_Neutrino_PtResol  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_PtResol","Tau3p_Neutrino_PtResol", 50, -2., 2., "p_{T} resol. #frac{#nu_{3p}^{reco} - #nu_{3p}^{gen}}{#nu_{3p}^{gen}}");
	Tau3p_Neutrino_PhiResol = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_PhiResol","Tau3p_Neutrino_PhiResol", 50, -1., 1., "#phi resol. #nu_{3p}^{reco} - #nu_{3p}^{gen}");
	Tau3p_Neutrino_XResol  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_XResol","Tau3p_Neutrino_XResol", 50, -50., 50., "x resol. #nu_{3p}^{reco} - #nu_{3p}^{gen}");
	Tau3p_Neutrino_YResol  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_YResol","Tau3p_Neutrino_YResol", 50, -50., 50., "y resol. #nu_{3p}^{reco} - #nu_{3p}^{gen}");
	Tau3p_Neutrino_EPull  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_PtPull","Tau3p_Neutrino_PtPull", 50, -5., 5., "p_{T}(#nu_{#tau}) pull");
	Tau3p_Neutrino_XPull  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_XPull","Tau3p_Neutrino_XPull", 50, -5., 5., "x(#nu_{#tau}) pull");
	Tau3p_Neutrino_YPull  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_YPull","Tau3p_Neutrino_YPull", 50, -5., 5., "y(#nu_{#tau}) pull");
	Tau3p_Neutrino_PhiPull  = HConfig.GetTH1D(Name+"_Tau3p_Neutrino_PhiPull","Tau3p_Neutrino_PhiPull", 50, -5., 5., "#phi(#nu_{#tau}) pull");
	MetMinus3pNeutrino_PtResol = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PtResol","MetMinus3pNeutrino_PtResol", 50, -2., 2., "p_{T} resol. #frac{#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen}}{#nu_{#tau,#mu}^{gen}}");
	MetMinus3pNeutrino_PhiResol = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PhiResol","MetMinus3pNeutrino_PhiResol", 50, -2., 2., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen})");
	MetMinus3pNeutrino_XResol = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_XResol","MetMinus3pNeutrino_XResol", 50, -50., 50., "x resol. #nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen}");
	MetMinus3pNeutrino_YResol = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_YResol","MetMinus3pNeutrino_YResol", 50, -50., 50., "y resol. #nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen}");
	MetMinus3pNeutrino_PhiPtResol  = HConfig.GetTH2D(Name+"_MetMinus3pNeutrino_PhiPtResol","MetMinus3pNeutrino_PhiPtResol", 50, -2., 2.,  50, -2., 2., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen})", "p_{T} resol. #frac{#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen}}{#nu_{#tau,#mu}^{gen}}");
	MetMinus3pNeutrino_PhiResol_FailedSVFit = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PhiResol_FailedSVFit","MetMinus3pNeutrino_PhiResol_FailedSVFit", 50, -2., 2., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen}) (SVFit failed)");
	MetMinus3pNeutrino_DeltaPhiMuon = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_DeltaPhiMuon","MetMinus3pNeutrino_DeltaPhiMuon", 50, -3.14159, 3.14159, "#Delta#phi(#nu_{#tau,#mu}^{reco},#mu_{sel})");
	MetMinus3pNeutrino_DeltaPhiMuon_FailedSVFit = HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_DeltaPhiMuon_FailedSVFit","MetMinus3pNeutrino_DeltaPhiMuon_FailedSVFit", 50, -3.14159, 3.14159, "#Delta#phi(#nu_{#tau,#mu}^{reco},#mu_{sel}) (SVFit failed)");

	svFitMass_Default = HConfig.GetTH1D(Name+"_svFitMass_Default",	"svFitMass_Default",	50, 0., 250., "m_{SVfit}(#tau_{h},#mu)/GeV");
	svFit3pMass = HConfig.GetTH1D(Name+"_svFit3pMass",	"svFit3pMass",	50, 0., 250., "m_{SVfit}(#tau_{3p},#mu)/GeV");;
	svFit3pMassResol = HConfig.GetTH1D(Name+"_svFit3pMassResol",	"svFit3pMassResol",	50, -1., 1., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p},#mu)");

	svFit3pMass_invalid 	= HConfig.GetTH1D(Name+"_svFit3pMass_invalid",	"svFit3pMass_invalid",	50, 0., 250., "m_{SVfit}(#tau_{3p},#mu)/GeV");;
	svFit3pMassResol_invalid= HConfig.GetTH1D(Name+"_svFit3pMassResol_invalid",	"svFit3pMassResol_invalid",	50, -1., 1., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p},#mu)");
	svFit3pMass_notCleaned	= HConfig.GetTH1D(Name+"_svFit3pMass_notCleaned",	"svFit3pMass_notCleaned",	50, 0., 250., "m_{SVfit}(#tau_{3p},#mu)/GeV");;
	svFit3pMassResol_notCleaned	= HConfig.GetTH1D(Name+"_svFit3pMassResol_notCleaned",	"svFit3pMassResol_notCleaned",	50, -1., 1., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p},#mu)");

	MetPxPull					= HConfig.GetTH1D(Name+"_MetPxPull",					"MetPxPull", 50, -5., 5., "E_{x}^{miss} pull");
	MetPyPull					= HConfig.GetTH1D(Name+"_MetPyPull",					"MetPyPull", 50, -5., 5., "E_{y}^{miss} pull");
	MetPxPyPull					= HConfig.GetTH2D(Name+"_MetPxPyPull",					"MetPxPyPull", 50, -5., 5., 50, -5., 5., "E_{x}^{miss} pull", "E_{y}^{miss} pull");
	MetPhiPull					= HConfig.GetTH1D(Name+"_MetPhiPull",					"MetPhiPull", 50, -5., 5., "#phi(E_{T}^{miss}) pull");
	MetMinus3pNeutrino_PxPull	= HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PxPull",	"MetMinus3pNeutrino_PxPull", 50, -5., 5., "E_{x}(#nu_{#tau,#mu}) pull");
	MetMinus3pNeutrino_PyPull	= HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PyPull",	"MetMinus3pNeutrino_PyPull", 50, -5., 5., "E_{y}(#nu_{#tau,#mu}) pull");
	MetMinus3pNeutrino_PxPyPull	= HConfig.GetTH2D(Name+"_MetMinus3pNeutrino_PxPyPull",	"MetMinus3pNeutrino_PxPyPull", 50, -5., 5.,50, -5., 5., "E_{x}(#nu_{#tau,#mu}) pull", "E_{y}(#nu_{#tau,#mu}) pull");
	MetMinus3pNeutrino_PhiPull	= HConfig.GetTH1D(Name+"_MetMinus3pNeutrino_PhiPull",	"MetMinus3pNeutrino_PhiPull", 50, -5., 5., "#phi(E_{T}(#nu_{#tau,#mu})) pull");

	mtt_calculateZnu_PlusSol = 			HConfig.GetTH1D(Name+"_mtt_calculateZnu_PlusSol", "mtt_calculateZnu_PlusSol", 50, 0., 250., "calc. z_{#nu} m_{#tau#tau}/GeV");
	mtt_calculateZnu_MinusSol = 		HConfig.GetTH1D(Name+"_mtt_calculateZnu_MinusSol", "mtt_calculateZnu_MinusSol", 50, 0., 250., "calc. z_{#nu} m_{#tau#tau}/GeV");
	mtt_calculateZnu_ZeroSol = 			HConfig.GetTH1D(Name+"_mtt_calculateZnu_ZeroSol", "mtt_calculateZnu_ZeroSol", 50, 0., 250., "calc. z_{#nu} m_{#tau#tau}/GeV");
	mtt_calculateZnu_TrueSol =			HConfig.GetTH1D(Name+"_mtt_calculateZnu_TrueSol", "mtt_calculateZnu_TrueSol", 50, 0., 250., "calc. z_{#nu} m_{#tau#tau}/GeV");
	mtt_calculateZnuNoAmb_PlusSol = 	HConfig.GetTH1D(Name+"_mtt_calculateZnuNoAmb_PlusSol", "mtt_calculateZnuNoAmb_PlusSol", 50, 0., 250., "calc. z_{#nu} m_{#tau#tau}/GeV");
	mtt_calculateZnuNoAmb_MinusSol = 	HConfig.GetTH1D(Name+"_mtt_calculateZnuNoAmb_MinusSol", "mtt_calculateZnuNoAmb_MinusSol", 50, 0., 250., "calc. z_{#nu} m_{#tau#tau}/GeV");
	mtt_calculateZnuNoAmb_ZeroSol = 	HConfig.GetTH1D(Name+"_mtt_calculateZnuNoAmb_ZeroSol", "mtt_calculateZnuNoAmb_ZeroSol", 50, 0., 250., "calc. z_{#nu} m_{#tau#tau}/GeV");
	mtt_calculateZnuNoAmb_TrueSol =		HConfig.GetTH1D(Name+"_mtt_calculateZnuNoAmb_TrueSol", "mtt_calculateZnuNoAmb_TrueSol", 50, 0., 250., "calc. z_{#nu} m_{#tau#tau}/GeV");
	mtt_3prongMuonRawMET_PlusSol = 		HConfig.GetTH1D(Name+"_mtt_3prongMuonRawMET_PlusSol", "mtt_3prongMuonRawMET_PlusSol", 50, 0., 250., "raw MET m_{#tau#tau}/GeV");
	mtt_3prongMuonRawMET_MinusSol = 	HConfig.GetTH1D(Name+"_mtt_3prongMuonRawMET_MinusSol", "mtt_3prongMuonRawMET_MinusSol", 50, 0., 250., "raw MET m_{#tau#tau}/GeV");
	mtt_3prongMuonRawMET_ZeroSol = 		HConfig.GetTH1D(Name+"_mtt_3prongMuonRawMET_ZeroSol", "mtt_3prongMuonRawMET_ZeroSol", 50, 0., 250., "raw MET m_{#tau#tau}/GeV");
	mtt_3prongMuonRawMET_TrueSol = 		HConfig.GetTH1D(Name+"_mtt_3prongMuonRawMET_TrueSol", "mtt_3prongMuonRawMET_TrueSol", 50, 0., 250., "raw MET m_{#tau#tau}/GeV");
	mtt_projectMetOnMu_PlusSol = 		HConfig.GetTH1D(Name+"_mtt_projectMetOnMu_PlusSol", "mtt_projectMetOnMu_PlusSol", 50, 0., 250., "Projection m_{#tau#tau}/GeV");
	mtt_projectMetOnMu_MinusSol = 		HConfig.GetTH1D(Name+"_mtt_projectMetOnMu_MinusSol", "mtt_projectMetOnMu_MinusSol", 50, 0., 250., "Projection m_{#tau#tau}/GeV");
	mtt_projectMetOnMu_ZeroSol = 		HConfig.GetTH1D(Name+"_mtt_projectMetOnMu_ZeroSol", "mtt_projectMetOnMu_ZeroSol", 50, 0., 250., "Projection m_{#tau#tau}/GeV");
	mtt_projectMetOnMu_TrueSol = 		HConfig.GetTH1D(Name+"_mtt_projectMetOnMu_TrueSol", "mtt_projectMetOnMu_TrueSol", 50, 0., 250., "Projection m_{#tau#tau}/GeV");
	mtt_projectMetOnMuRotate_PlusSol = 	HConfig.GetTH1D(Name+"_mtt_projectMetOnMuRotate_PlusSol", "mtt_projectMetOnMuRotate_PlusSol", 50, 0., 250., "Project+Rotate m_{#tau#tau}/GeV");
	mtt_projectMetOnMuRotate_MinusSol = HConfig.GetTH1D(Name+"_mtt_projectMetOnMuRotate_MinusSol", "mtt_projectMetOnMuRotate_MinusSol", 50, 0., 250., "Project+Rotate m_{#tau#tau}/GeV");
	mtt_projectMetOnMuRotate_ZeroSol = 	HConfig.GetTH1D(Name+"_mtt_projectMetOnMuRotate_ZeroSol", "mtt_projectMetOnMuRotate_ZeroSol", 50, 0., 250., "Project+Rotate m_{#tau#tau}/GeV");
	mtt_projectMetOnMuRotate_TrueSol = 	HConfig.GetTH1D(Name+"_mtt_projectMetOnMuRotate_TrueSol", "mtt_projectMetOnMuRotate_TrueSol", 50, 0., 250., "Project+Rotate m_{#tau#tau}/GeV");

	deltaMtt_calculateZnu = HConfig.GetTH1D(Name+"_deltaMtt_calculateZnu", "deltaMtt_calculateZnu", 50, -25., 25., "#Deltam_{#tau#tau}^{amb.}/GeV");

	svFit3pMassResol_vs_MetPhiResol = HConfig.GetTH2D(Name+"_svFit3pMassResol_vs_MetPhiResol", "svFit3pMassResol_vs_MetPhiResol", 30, 0., 3., 50, -100., 100., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen})", "m_{SVfit} - m_{true}(#tau_{3p},#mu) /GeV");
	mttCalculateZnuResol_vs_MetPhiResol = HConfig.GetTH2D(Name+"_mttCalculateZnuResol_vs_MetPhiResol", "mttCalculateZnuResol_vs_MetPhiResol", 30, 0., 3., 50, -100., 100., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen})", "(calc. z_{#nu} m_{#tau#tau} - m_{true})/GeV");
	mttCalculateZnuNoAmbResol_vs_MetPhiResol = HConfig.GetTH2D(Name+"_mttCalculateZnuNoAmbResol_vs_MetPhiResol", "mttCalculateZnuNoAmbResol_vs_MetPhiResol", 30, 0., 3., 50, -100., 100., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen})", "(calc. z_{#nu} m_{#tau#tau} - m_{true})/GeV");
	mtt3prongMuonRawMETResol_vs_MetPhiResol = HConfig.GetTH2D(Name+"_mtt3prongMuonRawMETResol_vs_MetPhiResol", "mtt3prongMuonRawMETResol_vs_MetPhiResol", 30, 0., 3., 50, -100., 100., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen})", "(raw MET m_{#tau#tau} - m_{true})/GeV");
	mttprojectMetOnMuRotateResol_vs_MetPhiResol = HConfig.GetTH2D(Name+"_mttprojectMetOnMuRotateResol_vs_MetPhiResol", "mttprojectMetOnMuRotateResol_vs_MetPhiResol", 30, 0., 3., 50, -100., 100., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen})", "(Project+Rotate m_{#tau#tau} - m_{true})/GeV");
	mttprojectMetOnMuResol_vs_MetPhiResol = HConfig.GetTH2D(Name+"_mttprojectMetOnMuResol_vs_MetPhiResol", "mttprojectMetOnMuResol_vs_MetPhiResol", 30, 0., 3., 50, -100., 100., "#phi resol. (#nu_{#tau,#mu}^{reco} - #nu_{#tau,#mu}^{gen})", "(Projection m_{#tau#tau} - m_{true})/GeV");

	svFit3pMassResol_vs_Tau3pEResol = HConfig.GetTH2D(Name+"_svFit3pMassResol_vs_Tau3pEResol", "svFit3pMassResol_vs_Tau3pEResol", 50, -1., 1., 50, -100., 100., "E resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}", "m_{SVfit} - m_{true}(#tau_{3p},#mu) /GeV");
	mttCalculateZnuResol_vs_Tau3pEResol = HConfig.GetTH2D(Name+"_mttCalculateZnuResol_vs_Tau3pEResol", "mttCalculateZnuResol_vs_Tau3pEResol", 50, -1., 1., 50, -100., 100., "E resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}", "m_{SVfit} - m_{true}(#tau_{3p},#mu) /GeV");
	mttCalculateZnuNoAmbResol_vs_Tau3pEResol = HConfig.GetTH2D(Name+"_mttCalculateZnuNoAmbResol_vs_Tau3pEResol", "mttCalculateZnuNoAmbResol_vs_Tau3pEResol", 50, -1., 1., 50, -100., 100., "E resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}", "m_{SVfit} - m_{true}(#tau_{3p},#mu) /GeV");
	mtt3prongMuonRawMETResol_vs_Tau3pEResol = HConfig.GetTH2D(Name+"_mtt3prongMuonRawMETResol_vs_Tau3pEResol", "mtt3prongMuonRawMETResol_vs_Tau3pEResol", 50, -1., 1., 50, -100., 100., "E resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}", "m_{SVfit} - m_{true}(#tau_{3p},#mu) /GeV");
	mttprojectMetOnMuRotateResol_vs_Tau3pEResol = HConfig.GetTH2D(Name+"_mttprojectMetOnMuRotateResol_vs_Tau3pEResol", "mttprojectMetOnMuRotateResol_vs_Tau3pEResol", 50, -1., 1., 50, -100., 100., "E resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}", "m_{SVfit} - m_{true}(#tau_{3p},#mu) /GeV");
	mttprojectMetOnMuResol_vs_Tau3pEResol = HConfig.GetTH2D(Name+"_mttprojectMetOnMuResol_vs_Tau3pEResol", "mttprojectMetOnMuResol_vs_Tau3pEResol", 50, -1., 1., 50, -100., 100., "E resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}", "m_{SVfit} - m_{true}(#tau_{3p},#mu) /GeV");
}

void ZeroJet3Prong::categoryExtradist(){
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

	Extradist1d.push_back(&svFit3pMass_PlusSol_notCleaned);
	Extradist1d.push_back(&svFit3pMass_MinusSol_notCleaned);
	Extradist1d.push_back(&svFit3pMass_ZeroSol_notCleaned);
	Extradist1d.push_back(&svFit3pMass_TrueSol_notCleaned);
	Extradist1d.push_back(&svFit3pMass_TruePlusSol_notCleaned);
	Extradist1d.push_back(&svFit3pMass_TrueMinusSol_notCleaned);
	Extradist1d.push_back(&svFit3pMassResol_PlusSol_notCleaned);
	Extradist1d.push_back(&svFit3pMassResol_MinusSol_notCleaned);
	Extradist1d.push_back(&svFit3pMassResol_ZeroSol_notCleaned);
	Extradist1d.push_back(&svFit3pMassResol_TrueSol_notCleaned);

	Extradist1d.push_back(&svFit3pMass_PlusSol_invalid);
	Extradist1d.push_back(&svFit3pMass_MinusSol_invalid);
	Extradist1d.push_back(&svFit3pMass_ZeroSol_invalid);
	Extradist1d.push_back(&svFit3pMass_TrueSol_invalid);
	Extradist1d.push_back(&svFit3pMass_TruePlusSol_invalid);
	Extradist1d.push_back(&svFit3pMass_TrueMinusSol_invalid);
	Extradist1d.push_back(&svFit3pMassResol_PlusSol_invalid);
	Extradist1d.push_back(&svFit3pMassResol_MinusSol_invalid);
	Extradist1d.push_back(&svFit3pMassResol_ZeroSol_invalid);
	Extradist1d.push_back(&svFit3pMassResol_TrueSol_invalid);

	Extradist1d.push_back(&svFit3p_MassLMax_TrueSol);
	Extradist1d.push_back(&svFit3p_MassLMax_WrongSol);
	Extradist1d.push_back(&svFit3p_MassLMax_ZeroSol);
	Extradist1d.push_back(&svFit3p_MassLMax_PlusSol);
	Extradist1d.push_back(&svFit3p_MassLMax_MinusSol);
	Extradist1d.push_back(&svFit3p_TrueSolByHigherLMax);
	Extradist2d.push_back(&svFit3p_SolutionMatrixByLMax);

	Extradist1d.push_back(&invalidSvFit3p_Tau3pNeutrinoXPull);
	Extradist1d.push_back(&invalidSvFit3p_Tau3pNeutrinoYPull);
	Extradist1d.push_back(&invalidSvFit3p_Tau3pNeutrinoPhiPull);
	Extradist1d.push_back(&invalidSvFit3p_MetMinus3pNeutrinoXPull);
	Extradist1d.push_back(&invalidSvFit3p_MetMinus3pNeutrinoYPull);
	Extradist1d.push_back(&invalidSvFit3p_MetMinus3pNeutrinoPhiPull);

	Extradist1d.push_back(&vis3pMass_PlusSol);
	Extradist1d.push_back(&vis3pMass_MinusSol);
	Extradist1d.push_back(&vis3pMass_ZeroSol);
	Extradist1d.push_back(&vis3pMass_TrueSol);
	Extradist1d.push_back(&vis3pMassResol_PlusSol);
	Extradist1d.push_back(&vis3pMassResol_MinusSol);
	Extradist1d.push_back(&vis3pMassResol_ZeroSol);
	Extradist1d.push_back(&vis3pMassResol_TrueSol);

	Extradist1d.push_back(&Tau3p_Neutrino_PtResol);
	Extradist1d.push_back(&Tau3p_Neutrino_PhiResol);
	Extradist1d.push_back(&Tau3p_Neutrino_XResol);
	Extradist1d.push_back(&Tau3p_Neutrino_YResol);
	Extradist1d.push_back(&Tau3p_Neutrino_EPull);
	Extradist1d.push_back(&Tau3p_Neutrino_XPull);
	Extradist1d.push_back(&Tau3p_Neutrino_YPull);
	Extradist1d.push_back(&Tau3p_Neutrino_PhiPull);
	Extradist1d.push_back(&MetMinus3pNeutrino_PtResol);
	Extradist1d.push_back(&MetMinus3pNeutrino_PhiResol);
	Extradist1d.push_back(&MetMinus3pNeutrino_XResol);
	Extradist1d.push_back(&MetMinus3pNeutrino_YResol);
	Extradist2d.push_back(&MetMinus3pNeutrino_PhiPtResol);

	Extradist1d.push_back(&MetMinus3pNeutrino_PhiResol_FailedSVFit);
	Extradist1d.push_back(&MetMinus3pNeutrino_DeltaPhiMuon);
	Extradist1d.push_back(&MetMinus3pNeutrino_DeltaPhiMuon_FailedSVFit);

	Extradist1d.push_back(&svFitMass_Default);
	Extradist1d.push_back(&svFit3pMass);
	Extradist1d.push_back(&svFit3pMassResol);

	Extradist1d.push_back(&svFit3pMass_notCleaned);
	Extradist1d.push_back(&svFit3pMassResol_notCleaned);
	Extradist1d.push_back(&svFit3pMass_invalid);
	Extradist1d.push_back(&svFit3pMassResol_invalid);

	Extradist1d.push_back(&MetPxPull);
	Extradist1d.push_back(&MetPyPull);
	Extradist2d.push_back(&MetPxPyPull);
	Extradist1d.push_back(&MetPhiPull);
	Extradist1d.push_back(&MetMinus3pNeutrino_PxPull);
	Extradist1d.push_back(&MetMinus3pNeutrino_PyPull);
	Extradist2d.push_back(&MetMinus3pNeutrino_PxPyPull);
	Extradist1d.push_back(&MetMinus3pNeutrino_PhiPull);

	Extradist1d.push_back(&mtt_calculateZnu_PlusSol);
	Extradist1d.push_back(&mtt_calculateZnu_MinusSol);
	Extradist1d.push_back(&mtt_calculateZnu_ZeroSol);
	Extradist1d.push_back(&mtt_calculateZnu_TrueSol);
	Extradist1d.push_back(&mtt_calculateZnuNoAmb_PlusSol);
	Extradist1d.push_back(&mtt_calculateZnuNoAmb_MinusSol);
	Extradist1d.push_back(&mtt_calculateZnuNoAmb_ZeroSol);
	Extradist1d.push_back(&mtt_calculateZnuNoAmb_TrueSol);
	Extradist1d.push_back(&mtt_3prongMuonRawMET_PlusSol);
	Extradist1d.push_back(&mtt_3prongMuonRawMET_MinusSol);
	Extradist1d.push_back(&mtt_3prongMuonRawMET_ZeroSol);
	Extradist1d.push_back(&mtt_3prongMuonRawMET_TrueSol);
	Extradist1d.push_back(&mtt_projectMetOnMu_PlusSol);
	Extradist1d.push_back(&mtt_projectMetOnMu_MinusSol);
	Extradist1d.push_back(&mtt_projectMetOnMu_ZeroSol);
	Extradist1d.push_back(&mtt_projectMetOnMu_TrueSol);
	Extradist1d.push_back(&mtt_projectMetOnMuRotate_PlusSol);
	Extradist1d.push_back(&mtt_projectMetOnMuRotate_MinusSol);
	Extradist1d.push_back(&mtt_projectMetOnMuRotate_ZeroSol);
	Extradist1d.push_back(&mtt_projectMetOnMuRotate_TrueSol);

	Extradist1d.push_back(&deltaMtt_calculateZnu);

	Extradist2d.push_back(&svFit3pMassResol_vs_MetPhiResol);
	Extradist2d.push_back(&mttCalculateZnuResol_vs_MetPhiResol);
	Extradist2d.push_back(&mttCalculateZnuNoAmbResol_vs_MetPhiResol);
	Extradist2d.push_back(&mtt3prongMuonRawMETResol_vs_MetPhiResol);
	Extradist2d.push_back(&mttprojectMetOnMuRotateResol_vs_MetPhiResol);
	Extradist2d.push_back(&mttprojectMetOnMuResol_vs_MetPhiResol);

	Extradist2d.push_back(&svFit3pMassResol_vs_Tau3pEResol);
	Extradist2d.push_back(&mttCalculateZnuResol_vs_Tau3pEResol);
	Extradist2d.push_back(&mttCalculateZnuNoAmbResol_vs_Tau3pEResol);
	Extradist2d.push_back(&mtt3prongMuonRawMETResol_vs_Tau3pEResol);
	Extradist2d.push_back(&mttprojectMetOnMuRotateResol_vs_Tau3pEResol);
	Extradist2d.push_back(&mttprojectMetOnMuResol_vs_Tau3pEResol);
}

bool ZeroJet3Prong::categorySelection(){
	bool categoryPass = true;
	std::vector<float> value_ZeroJet3Prong(NCuts,-10);
	std::vector<float> pass_ZeroJet3Prong(NCuts,false);

	Logger(Logger::Debug) << "Cut: Number of Jets" << std::endl;
	value_ZeroJet3Prong.at(NJet) = nJets_;
	pass_ZeroJet3Prong.at(NJet) = ( value_ZeroJet3Prong.at(NJet) <= cut.at(NJet) );

	if (selTau == -1){
		// decayMode cut is set to true for nice N-0 and N-1 plots
		value_ZeroJet3Prong.at(DecayMode) = -10;
		pass_ZeroJet3Prong.at(DecayMode) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		/* enable and adapt this to only select mu+3prong tau events
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
		*/
		Logger(Logger::Debug) << "Cut: Tau decay mode" << std::endl;
		value_ZeroJet3Prong.at(DecayMode)	= Ntp->PFTau_hpsDecayMode(selTau);
		pass_ZeroJet3Prong.at(DecayMode)		= (value_ZeroJet3Prong.at(DecayMode)==cut.at(DecayMode));
	}

	if (selTau == -1 || not pass_ZeroJet3Prong.at(DecayMode)) {
		value_ZeroJet3Prong.at(SigmaSV) = -999;
		pass_ZeroJet3Prong.at(SigmaSV) = true;
	}
	else if ( not Ntp->PFTau_TIP_hassecondaryVertex(selTau) ){
		value_ZeroJet3Prong.at(SigmaSV) = -9;
		pass_ZeroJet3Prong.at(SigmaSV) = false;
	}
	else {
		Logger(Logger::Debug) << "Cut: SV significance" << std::endl;
		value_ZeroJet3Prong.at(SigmaSV) = flightLengthSig_;
		pass_ZeroJet3Prong.at(SigmaSV) = ( value_ZeroJet3Prong.at(SigmaSV) >= cut.at(SigmaSV) );
	}

	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("ZeroJet3Prong",value_ZeroJet3Prong, pass_ZeroJet3Prong,NCuts) && categoryPass;
	return categoryPass;
}

void ZeroJet3Prong::categoryPlotting(){
	TPTRObject TPResults;

	if (status){
		bool isMC = (not Ntp->isData()) || (idStripped == DataMCType::DY_mutau_embedded);

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
		GEF.SetCorrectPt(false);
		TPResults = GEF.getTPTRObject();

		// calculate SVFit stuff for comparison
		SVFitObject *svfObjDefault = NULL;
		if (runSVFit_){
			svfObjDefault = Ntp->getSVFitResult_MuTauh(svfitstorage, "CorrMVAMuTau", selMuon, selTau, 50000);
			svFitMass_Default.at(t).Fill(svfObjDefault->get_mass(), w);
		}

		double trueResonanceMass = (isMC) ? Ntp->getResonanceMassFromGenInfo(true, true, false) : -999;
		int i_matchedMCTau = (isMC) ? Ntp->matchTauTruth(selTau, true) : -999;
		int i_otherTau(-1);
		if (i_matchedMCTau == 0) i_otherTau = 1;
		else if (i_matchedMCTau == 1) i_otherTau = 0;
		else if (!isMC) i_otherTau = -999;
		else Logger(Logger::Warning) << "Cannot determine index of tau_mu." << std::endl;
		bool isSignalTauDecay = i_matchedMCTau >= 0 && Ntp->MCTau_JAK(i_matchedMCTau) == 5;
		int true3ProngAmbig = isSignalTauDecay ? Ntp->MCTau_true3prongAmbiguity( i_matchedMCTau ) : -999;
		TLorentzVector trueTauP4 = isSignalTauDecay ? Ntp->MCTau_p4(i_matchedMCTau) : TLorentzVector(0,0,0,0);

		TLorentzVector trueMet(0,0,0,0);
		for(int i_tau = 0; i_tau < Ntp->NMCTaus(); i_tau++){
			if(i_tau == 2) Logger(Logger::Error) << "More than 2 signal taus in event!" << std::endl;
			trueMet += Ntp->MCTau_invisiblePart(i_tau);
		}

		if (isMC) {
			if(Ntp->MET_CorrMVAMuTau_significance_xx() != 0 && Ntp->MET_CorrMVAMuTau_significance_yy() != 0){
				double pxPull = (Ntp->MET_CorrMVAMuTau_ex() - trueMet.Px()) / sqrt(Ntp->MET_CorrMVAMuTau_significance_xx());
				double pyPull = (Ntp->MET_CorrMVAMuTau_ey() - trueMet.Py()) / sqrt(Ntp->MET_CorrMVAMuTau_significance_yy());
				double phiPull = Tools::DeltaPhi(Ntp->MET_CorrMVAMuTau_phi(), trueMet.Phi()) /
						Tools::phiUncertainty(Ntp->MET_CorrMVAMuTau_ex(), Ntp->MET_CorrMVAMuTau_ey(), Ntp->MET_CorrMVAMuTau_significance_xx(), Ntp->MET_CorrMVAMuTau_significance_yy(), Ntp->MET_CorrMVAMuTau_significance_xy());
				MetPxPull.at(t).Fill( pxPull , w);
				MetPyPull.at(t).Fill( pyPull , w);
				MetPxPyPull.at(t).Fill( pxPull, pyPull, w);
				MetPhiPull.at(t).Fill( phiPull , w);
			}
			else
				Logger(Logger::Warning) << "MET covariance is zero!" << std::endl;
		}

		if (TPResults.isAmbiguous()) {
			SVFitObject *svfObjPlus  = NULL;
			SVFitObject *svfObjMinus = NULL;
			bool plusLmaxIsHigher;
			if (runSVFit_){
				svfObjPlus = Ntp->getSVFitResult_MuTau3p(svfitstor_Tau3pPlus, "CorrMVAMuTau", selMuon, TPResults.getTauPlus().LV(), TPResults.getNeutrinoPlus(), "3ProngRecoPlus");
				svfObjMinus = Ntp->getSVFitResult_MuTau3p(svfitstor_Tau3pMinus, "CorrMVAMuTau", selMuon, TPResults.getTauMinus().LV(), TPResults.getNeutrinoMinus(), "3ProngRecoMinus");
				svFit3pMass_PlusSol_notCleaned.at(t).Fill(svfObjPlus->get_mass(), w);
				if (svfObjPlus->isValid())	svFit3pMass_PlusSol.at(t).Fill(svfObjPlus->get_mass(), w);
				else						svFit3pMass_PlusSol_invalid.at(t).Fill(svfObjPlus->get_mass(), w);
				svFit3pMass_MinusSol_notCleaned.at(t).Fill(svfObjMinus->get_mass(), w);
				if (svfObjMinus->isValid())	svFit3pMass_MinusSol.at(t).Fill(svfObjMinus->get_mass(), w);
				else						svFit3pMass_MinusSol_invalid.at(t).Fill(svfObjMinus->get_mass(), w);
				plusLmaxIsHigher = (svfObjPlus->get_massLmax() > svfObjMinus->get_massLmax());

				svFit3p_MassLMax_PlusSol.at(t).Fill(svfObjPlus->get_massLmax(), w);
				svFit3p_MassLMax_MinusSol.at(t).Fill(svfObjMinus->get_massLmax(), w);
			}

			double visMassPlus = (Ntp->Muon_p4(selMuon) + TPResults.getTauPlus().LV()).M();
			double visMassMinus = (Ntp->Muon_p4(selMuon) + TPResults.getTauMinus().LV()).M();
			vis3pMass_PlusSol.at(t).Fill( visMassPlus, w);
			vis3pMass_MinusSol.at(t).Fill( visMassMinus, w);

			objects::MET metPlus(Ntp, "CorrMVAMuTau");
			metPlus.subtractNeutrino(TPResults.getNeutrinoPlus());
			objects::MET metMinus(Ntp, "CorrMVAMuTau");
			metMinus.subtractNeutrino(TPResults.getNeutrinoMinus());
			//TLorentzVector genMet = Ntp->MCTau_invisiblePart(i_otherTau);
			//objects::MET metPlus("gen", genMet.Px(), genMet.Py());
			//objects::MET metMinus(metPlus);

			mtt_calculateZnu_PlusSol.at(t).Fill(mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);
			mtt_calculateZnu_MinusSol.at(t).Fill(mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);
			mtt_calculateZnuNoAmb_PlusSol.at(t).Fill(mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);
			mtt_calculateZnuNoAmb_MinusSol.at(t).Fill(mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);
			mtt_3prongMuonRawMET_PlusSol.at(t).Fill(mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);
			mtt_3prongMuonRawMET_MinusSol.at(t).Fill(mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);
			mtt_projectMetOnMu_PlusSol.at(t).Fill(mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);
			mtt_projectMetOnMu_MinusSol.at(t).Fill(mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);
			mtt_projectMetOnMuRotate_PlusSol.at(t).Fill(mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);
			mtt_projectMetOnMuRotate_MinusSol.at(t).Fill(mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);

			if (isMC) {
				if (runSVFit_) {
					svFit3pMassResol_PlusSol_notCleaned.at(t).Fill( (svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
					if (svfObjPlus->isValid())	svFit3pMassResol_PlusSol.at(t).Fill( (svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
					else						svFit3pMassResol_PlusSol_invalid.at(t).Fill( (svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
					svFit3pMassResol_MinusSol_notCleaned.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
					if (svfObjMinus->isValid())	svFit3pMassResol_MinusSol.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
					else 						svFit3pMassResol_MinusSol_invalid.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);

				}
				vis3pMassResol_PlusSol.at(t).Fill( (visMassPlus - trueResonanceMass)/trueResonanceMass, w);
				vis3pMassResol_MinusSol.at(t).Fill( (visMassMinus - trueResonanceMass)/trueResonanceMass, w);

				double phiUnc(-1);
				double tau3pEresol(-1);
				switch (true3ProngAmbig) {
					case MultiProngTauSolver::plus:
						if (runSVFit_) {
							svFit3pMass_TrueSol_notCleaned.at(t).Fill(svfObjPlus->get_mass(), w);
							svFit3pMass_TruePlusSol_notCleaned.at(t).Fill(svfObjPlus->get_mass(), w);
							svFit3pMassResol_TrueSol_notCleaned.at(t).Fill((svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
							svFit3pMass_notCleaned.at(t).Fill(svfObjPlus->get_mass(), w);
							svFit3pMassResol_notCleaned.at(t).Fill((svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);

							if(svfObjPlus->isValid()){
								svFit3pMass_TrueSol.at(t).Fill(svfObjPlus->get_mass(), w);
								svFit3pMass_TruePlusSol.at(t).Fill(svfObjPlus->get_mass(), w);
								svFit3pMassResol_TrueSol.at(t).Fill((svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
								svFit3pMass.at(t).Fill(svfObjPlus->get_mass(), w);
								svFit3pMassResol.at(t).Fill((svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);

								svFit3pMassResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), (svfObjPlus->get_mass() - trueResonanceMass), w);
								svFit3pMassResol_vs_Tau3pEResol.at(t).Fill((TPResults.getTauPlus().LV().E() - trueTauP4.E())/trueTauP4.E(), (svfObjPlus->get_mass() - trueResonanceMass), w);
							}
							else{
								svFit3pMass_TrueSol_invalid.at(t).Fill(svfObjPlus->get_mass(), w);
								svFit3pMass_TruePlusSol_invalid.at(t).Fill(svfObjPlus->get_mass(), w);
								svFit3pMassResol_TrueSol_invalid.at(t).Fill((svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
								svFit3pMass_invalid.at(t).Fill(svfObjPlus->get_mass(), w);
								svFit3pMassResol_invalid.at(t).Fill((svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);

								if (TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px) > 0.001)
									invalidSvFit3p_Tau3pNeutrinoXPull.at(t).Fill( (TPResults.getNeutrinoPlus().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px()) / TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px), w);
								if (TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py) > 0.001)
									invalidSvFit3p_Tau3pNeutrinoYPull.at(t).Fill( (TPResults.getNeutrinoPlus().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py()) / TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py), w);
								double phiUnc = Tools::phiUncertainty(TPResults.getNeutrinoPlus().LV().Px(), TPResults.getNeutrinoPlus().LV().Py(), TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px),
										TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py));
								if (phiUnc > 0.001)
									invalidSvFit3p_Tau3pNeutrinoPhiPull.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoPlus().LV(), Ntp->MCTau_invisiblePart(i_matchedMCTau)) / phiUnc, w);

								if(metPlus.significanceXX() != 0 && metPlus.significanceYY() != 0){
									double pxPull = (metPlus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px()) / sqrt(metPlus.significanceXX());
									double pyPull = (metPlus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py()) / sqrt(metPlus.significanceYY());
									double phiPull = Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()) / metPlus.phiUncertainty();
									invalidSvFit3p_MetMinus3pNeutrinoXPull.at(t).Fill( pxPull , w);
									invalidSvFit3p_MetMinus3pNeutrinoYPull.at(t).Fill( pyPull , w);
									invalidSvFit3p_MetMinus3pNeutrinoPhiPull.at(t).Fill( phiPull , w);
								}
							}

							MetMinus3pNeutrino_DeltaPhiMuon.at(t).Fill( Tools::DeltaPhi(metPlus.phi(), Ntp->Muon_p4(selMuon).Phi()) );
							if (not svfObjPlus->isValid()) {
								MetMinus3pNeutrino_PhiResol_FailedSVFit.at(t).Fill( Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
								MetMinus3pNeutrino_DeltaPhiMuon_FailedSVFit.at(t).Fill( Tools::DeltaPhi(metPlus.phi(), Ntp->Muon_p4(selMuon).Phi()) );
							}

							svFit3p_MassLMax_TrueSol.at(t).Fill(svfObjPlus->get_massLmax(), w);
							svFit3p_MassLMax_WrongSol.at(t).Fill(svfObjMinus->get_massLmax(), w);
							svFit3p_TrueSolByHigherLMax.at(t).Fill(plusLmaxIsHigher,w);
							svFit3p_SolutionMatrixByLMax.at(t).Fill(true3ProngAmbig, plusLmaxIsHigher+1, w);
						}
						vis3pMass_TrueSol.at(t).Fill((Ntp->Muon_p4(selMuon) + TPResults.getTauPlus().LV()).M(), w);
						vis3pMassResol_TrueSol.at(t).Fill((visMassPlus - trueResonanceMass)/trueResonanceMass, w);

						mtt_calculateZnu_TrueSol.at(t).Fill(mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);
						mtt_calculateZnuNoAmb_TrueSol.at(t).Fill(mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);
						mtt_3prongMuonRawMET_TrueSol.at(t).Fill(mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);
						mtt_projectMetOnMu_TrueSol.at(t).Fill(mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);
						mtt_projectMetOnMuRotate_TrueSol.at(t).Fill(mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::plus, selMuon, metPlus), w);

						mttCalculateZnuResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), 			mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);
						mttCalculateZnuNoAmbResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), 		mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);
						mtt3prongMuonRawMETResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), 		mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);
						mttprojectMetOnMuRotateResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), 	mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);
						mttprojectMetOnMuResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), 			mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);

						tau3pEresol = (TPResults.getTauPlus().LV().E() - trueTauP4.E())/trueTauP4.E();
						mttCalculateZnuResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol, 		mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);
						mttCalculateZnuNoAmbResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,	mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);
						mtt3prongMuonRawMETResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,		mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);
						mttprojectMetOnMuRotateResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,	mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);
						mttprojectMetOnMuResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,		mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::plus, selMuon, metPlus) - trueResonanceMass, w);

						if (Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt() > 0.001){
							Tau3p_Neutrino_PtResol.at(t).Fill( (TPResults.getNeutrinoPlus().LV().Pt() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt(), w);
							Tau3p_Neutrino_XResol.at(t).Fill( TPResults.getNeutrinoPlus().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px(), w);
							Tau3p_Neutrino_YResol.at(t).Fill( TPResults.getNeutrinoPlus().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py(), w);
							Tau3p_Neutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoPlus().LV().Phi(), Ntp->MCTau_invisiblePart(i_matchedMCTau).Phi()), w);
						}
						if (TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::E, LorentzVectorParticle::E) > 0.001)
							Tau3p_Neutrino_EPull.at(t).Fill( (TPResults.getNeutrinoPlus().LV().E() - Ntp->MCTau_invisiblePart(i_matchedMCTau).E())/TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::E, LorentzVectorParticle::E), w);
						if (TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px) > 0.001)
							Tau3p_Neutrino_XPull.at(t).Fill( (TPResults.getNeutrinoPlus().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px()) / TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px), w);
						if (TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py) > 0.001)
							Tau3p_Neutrino_YPull.at(t).Fill( (TPResults.getNeutrinoPlus().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py()) / TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py), w);
						phiUnc = Tools::phiUncertainty(TPResults.getNeutrinoPlus().LV().Px(), TPResults.getNeutrinoPlus().LV().Py(), TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px),
								TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), TPResults.getNeutrinoPlus().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py));
						if (phiUnc > 0.001)
							Tau3p_Neutrino_PhiPull.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoPlus().LV(), Ntp->MCTau_invisiblePart(i_matchedMCTau)) / phiUnc, w);

						if (Ntp->MCTau_invisiblePart(i_otherTau).Pt() > 0.001){
							MetMinus3pNeutrino_PtResol.at(t).Fill( (metPlus.et() - Ntp->MCTau_invisiblePart(i_otherTau).Pt())/Ntp->MCTau_invisiblePart(i_otherTau).Pt(), w);
							MetMinus3pNeutrino_XResol.at(t).Fill( metPlus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px(), w);
							MetMinus3pNeutrino_YResol.at(t).Fill( metPlus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py(), w);
							MetMinus3pNeutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
							MetMinus3pNeutrino_PhiPtResol.at(t).Fill( Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()),(metPlus.et() - Ntp->MCTau_invisiblePart(i_otherTau).Pt())/Ntp->MCTau_invisiblePart(i_otherTau).Pt(), w);
						}
						if(metPlus.significanceXX() != 0 && metPlus.significanceYY() != 0){
							double pxPull = (metPlus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px()) / sqrt(metPlus.significanceXX());
							double pyPull = (metPlus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py()) / sqrt(metPlus.significanceYY());
							double phiPull = Tools::DeltaPhi(metPlus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()) / metPlus.phiUncertainty();
							MetMinus3pNeutrino_PxPull.at(t).Fill( pxPull , w);
							MetMinus3pNeutrino_PyPull.at(t).Fill( pyPull , w);
							MetMinus3pNeutrino_PxPyPull.at(t).Fill( pxPull, pyPull, w);
							MetMinus3pNeutrino_PhiPull.at(t).Fill( phiPull , w);
						}
						else
							Logger(Logger::Warning) << "MET-plus covariance is zero!" << std::endl;
						break;
					case MultiProngTauSolver::minus:
						if (runSVFit_) {
							svFit3pMass_TrueSol_notCleaned.at(t).Fill(svfObjMinus->get_mass(), w);
							svFit3pMass_TrueMinusSol_notCleaned.at(t).Fill(svfObjMinus->get_mass(), w);
							svFit3pMassResol_TrueSol_notCleaned.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
							svFit3pMass_notCleaned.at(t).Fill(svfObjMinus->get_mass(), w);
							svFit3pMassResol_notCleaned.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);

							if (svfObjMinus->isValid()) {
								svFit3pMass_TrueSol.at(t).Fill(svfObjMinus->get_mass(), w);
								svFit3pMass_TrueMinusSol.at(t).Fill(svfObjMinus->get_mass(), w);
								svFit3pMassResol_TrueSol.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
								svFit3pMass.at(t).Fill(svfObjMinus->get_mass(), w);
								svFit3pMassResol.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);

								svFit3pMassResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), (svfObjMinus->get_mass() - trueResonanceMass), w);
								svFit3pMassResol_vs_Tau3pEResol.at(t).Fill((TPResults.getTauMinus().LV().E() - trueTauP4.E())/trueTauP4.E(), (svfObjMinus->get_mass() - trueResonanceMass), w);
							}
							else {
								svFit3pMass_TrueSol_invalid.at(t).Fill(svfObjMinus->get_mass(), w);
								svFit3pMass_TrueMinusSol_invalid.at(t).Fill(svfObjMinus->get_mass(), w);
								svFit3pMassResol_TrueSol_invalid.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
								svFit3pMass_invalid.at(t).Fill(svfObjMinus->get_mass(), w);
								svFit3pMassResol_invalid.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);

								if (TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px) > 0.001)
									invalidSvFit3p_Tau3pNeutrinoXPull.at(t).Fill( (TPResults.getNeutrinoMinus().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px()) / TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px), w);
								if (TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py) > 0.001)
									invalidSvFit3p_Tau3pNeutrinoXPull.at(t).Fill( (TPResults.getNeutrinoMinus().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py()) / TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py), w);
								phiUnc = Tools::phiUncertainty(TPResults.getNeutrinoMinus().LV().Px(), TPResults.getNeutrinoMinus().LV().Py(), TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px),
										TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py));
								if (phiUnc > 0.001)
									invalidSvFit3p_Tau3pNeutrinoPhiPull.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoMinus().LV(), Ntp->MCTau_invisiblePart(i_matchedMCTau)) / phiUnc, w);

								if(metMinus.significanceXX() != 0 && metMinus.significanceYY() != 0){
									double pxPull = (metMinus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px()) / sqrt(metMinus.significanceXX());
									double pyPull = (metMinus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py()) / sqrt(metMinus.significanceYY());
									double phiPull = Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()) / metMinus.phiUncertainty();
									invalidSvFit3p_MetMinus3pNeutrinoXPull.at(t).Fill( pxPull , w);
									invalidSvFit3p_MetMinus3pNeutrinoYPull.at(t).Fill( pyPull , w);
									invalidSvFit3p_MetMinus3pNeutrinoPhiPull.at(t).Fill( phiPull , w);
								}
							}

							MetMinus3pNeutrino_DeltaPhiMuon.at(t).Fill( Tools::DeltaPhi(metMinus.phi(), Ntp->Muon_p4(selMuon).Phi()) );
							if (not svfObjMinus->isValid()) {
								MetMinus3pNeutrino_PhiResol_FailedSVFit.at(t).Fill( Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
								MetMinus3pNeutrino_DeltaPhiMuon_FailedSVFit.at(t).Fill( Tools::DeltaPhi(metMinus.phi(), Ntp->Muon_p4(selMuon).Phi()) );
							}

							svFit3p_MassLMax_TrueSol.at(t).Fill(svfObjMinus->get_massLmax(), w);
							svFit3p_MassLMax_WrongSol.at(t).Fill(svfObjPlus->get_massLmax(), w);
							svFit3p_TrueSolByHigherLMax.at(t).Fill(!plusLmaxIsHigher,w);
							svFit3p_SolutionMatrixByLMax.at(t).Fill(true3ProngAmbig, plusLmaxIsHigher+1, w);
						}
						vis3pMass_TrueSol.at(t).Fill((Ntp->Muon_p4(selMuon) + TPResults.getTauMinus().LV()).M(), w);
						vis3pMassResol_TrueSol.at(t).Fill((visMassMinus - trueResonanceMass)/trueResonanceMass, w);

						mtt_calculateZnu_TrueSol.at(t).Fill(mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);
						mtt_calculateZnuNoAmb_TrueSol.at(t).Fill(mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);
						mtt_3prongMuonRawMET_TrueSol.at(t).Fill(mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);
						mtt_projectMetOnMu_TrueSol.at(t).Fill(mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);
						mtt_projectMetOnMuRotate_TrueSol.at(t).Fill(mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::minus, selMuon, metMinus), w);

						mttCalculateZnuResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);
						mttCalculateZnuNoAmbResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);
						mtt3prongMuonRawMETResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);
						mttprojectMetOnMuRotateResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);
						mttprojectMetOnMuResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);

						tau3pEresol = (TPResults.getTauMinus().LV().E() - trueTauP4.E())/trueTauP4.E();
						mttCalculateZnuResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol, 		mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);
						mttCalculateZnuNoAmbResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,	mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);
						mtt3prongMuonRawMETResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,		mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);
						mttprojectMetOnMuRotateResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,	mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);
						mttprojectMetOnMuResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,		mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::minus, selMuon, metMinus) - trueResonanceMass, w);

						if (Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt() > 0.001){
							Tau3p_Neutrino_PtResol.at(t).Fill( (TPResults.getNeutrinoMinus().LV().Pt() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt(), w);
							Tau3p_Neutrino_XResol.at(t).Fill( TPResults.getNeutrinoMinus().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px(), w);
							Tau3p_Neutrino_YResol.at(t).Fill( TPResults.getNeutrinoMinus().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py(), w);
							Tau3p_Neutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoMinus().LV().Phi(), Ntp->MCTau_invisiblePart(i_matchedMCTau).Phi()), w);
						}
						if (TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::E, LorentzVectorParticle::E) > 0.001)
							Tau3p_Neutrino_EPull.at(t).Fill( (TPResults.getNeutrinoMinus().LV().E() - Ntp->MCTau_invisiblePart(i_matchedMCTau).E())/TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::E, LorentzVectorParticle::E), w);
						if (TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px) > 0.001)
							Tau3p_Neutrino_XPull.at(t).Fill( (TPResults.getNeutrinoMinus().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px()) / TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px), w);
						if (TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py) > 0.001)
							Tau3p_Neutrino_YPull.at(t).Fill( (TPResults.getNeutrinoMinus().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py()) / TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py), w);
						phiUnc = Tools::phiUncertainty(TPResults.getNeutrinoMinus().LV().Px(), TPResults.getNeutrinoMinus().LV().Py(), TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px),
								TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), TPResults.getNeutrinoMinus().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py));
						if (phiUnc > 0.001)
							Tau3p_Neutrino_PhiPull.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoMinus().LV(), Ntp->MCTau_invisiblePart(i_matchedMCTau)) / phiUnc, w);

						if (Ntp->MCTau_invisiblePart(i_otherTau).Pt() > 0.001){
							MetMinus3pNeutrino_PtResol.at(t).Fill( (metMinus.et() - Ntp->MCTau_invisiblePart(i_otherTau).Pt())/Ntp->MCTau_invisiblePart(i_otherTau).Pt(), w);
							MetMinus3pNeutrino_XResol.at(t).Fill( metMinus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px(), w);
							MetMinus3pNeutrino_YResol.at(t).Fill( metMinus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py(), w);
							MetMinus3pNeutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
							MetMinus3pNeutrino_PhiPtResol.at(t).Fill( Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()),(metMinus.et() - Ntp->MCTau_invisiblePart(i_otherTau).Pt())/Ntp->MCTau_invisiblePart(i_otherTau).Pt(), w);
						}
						if(metMinus.significanceXX() != 0 && metMinus.significanceYY() != 0){
							double pxPull = (metMinus.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px()) / sqrt(metMinus.significanceXX());
							double pyPull = (metMinus.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py()) / sqrt(metMinus.significanceYY());
							double phiPull = Tools::DeltaPhi(metMinus.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()) / metMinus.phiUncertainty();
							MetMinus3pNeutrino_PxPull.at(t).Fill( pxPull , w);
							MetMinus3pNeutrino_PyPull.at(t).Fill( pyPull , w);
							MetMinus3pNeutrino_PxPyPull.at(t).Fill( pxPull, pyPull, w);
							MetMinus3pNeutrino_PhiPull.at(t).Fill( phiPull , w);
						}
						else
							Logger(Logger::Warning) << "MET-minus covariance is zero!" << std::endl;
						break;
					case -9:
					case -999:
						Logger(Logger::Verbose) << "HPS tau not matched to signal tau decay." << std::endl;
						break;
					default:
						Logger(Logger::Warning) << "True ambiguity value is " << true3ProngAmbig << ", which does not make sense." << std::endl;
						break;
				}
			}
		}
		else {
			SVFitObject *svfObjZero  = NULL;
			if (runSVFit_) {
				svfObjZero = Ntp->getSVFitResult_MuTau3p(svfitstor_Tau3pZero, "CorrMVAMuTau", selMuon, TPResults.getTauZero().LV(), TPResults.getNeutrinoZero(), "3ProngRecoZero");
				svFit3pMass_ZeroSol_notCleaned.at(t).Fill(svfObjZero->get_mass(), w);
				if (svfObjZero->isValid())	svFit3pMass_ZeroSol.at(t).Fill(svfObjZero->get_mass(), w);
				else						svFit3pMass_ZeroSol_invalid.at(t).Fill(svfObjZero->get_mass(), w);
			}

			double visMassZero = (Ntp->Muon_p4(selMuon) + TPResults.getTauZero().LV()).M();
			vis3pMass_ZeroSol.at(t).Fill(visMassZero, w);

			if (isMC) {
				objects::MET metZero(Ntp, "CorrMVAMuTau");
				metZero.subtractNeutrino(TPResults.getNeutrinoZero());
				//TLorentzVector genMet = Ntp->MCTau_invisiblePart(i_otherTau);
				//objects::MET metZero("gen", genMet.Px(), genMet.Py());
				if (runSVFit_) {
					svFit3pMassResol_ZeroSol_notCleaned.at(t).Fill((svfObjZero->get_mass() - trueResonanceMass) / trueResonanceMass, w);
					svFit3pMass_notCleaned.at(t).Fill(svfObjZero->get_mass(), w);
					svFit3pMassResol_notCleaned.at(t).Fill((svfObjZero->get_mass() - trueResonanceMass) / trueResonanceMass, w);
					if (svfObjZero->isValid()){
						svFit3pMassResol_ZeroSol.at(t).Fill((svfObjZero->get_mass() - trueResonanceMass) / trueResonanceMass, w);
						svFit3pMass.at(t).Fill(svfObjZero->get_mass(), w);
						svFit3pMassResol.at(t).Fill((svfObjZero->get_mass() - trueResonanceMass) / trueResonanceMass, w);

						svFit3pMassResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), (svfObjZero->get_mass() - trueResonanceMass), w);
						svFit3pMassResol_vs_Tau3pEResol.at(t).Fill((TPResults.getTauZero().LV().E() - trueTauP4.E())/trueTauP4.E(), (svfObjZero->get_mass() - trueResonanceMass), w);
					}
					else{
						svFit3pMassResol_ZeroSol_invalid.at(t).Fill((svfObjZero->get_mass() - trueResonanceMass) / trueResonanceMass, w);
						svFit3pMass_invalid.at(t).Fill(svfObjZero->get_mass(), w);
						svFit3pMassResol_invalid.at(t).Fill((svfObjZero->get_mass() - trueResonanceMass) / trueResonanceMass, w);

						if (TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px) > 0.001)
							invalidSvFit3p_Tau3pNeutrinoXPull.at(t).Fill( (TPResults.getNeutrinoZero().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px()) / TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px), w);
						if (TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py) > 0.001)
							invalidSvFit3p_Tau3pNeutrinoYPull.at(t).Fill( (TPResults.getNeutrinoZero().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py()) / TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py), w);
						double phiUnc = Tools::phiUncertainty(TPResults.getNeutrinoZero().LV().Px(), TPResults.getNeutrinoZero().LV().Py(), TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px),
								TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py));
						if (phiUnc > 0.001)
							invalidSvFit3p_Tau3pNeutrinoPhiPull.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoZero().LV(), Ntp->MCTau_invisiblePart(i_matchedMCTau)) / phiUnc, w);

						if(metZero.significanceXX() != 0 && metZero.significanceYY() != 0){
							double pxPull = (metZero.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px()) / sqrt(metZero.significanceXX());
							double pyPull = (metZero.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py()) / sqrt(metZero.significanceYY());
							double phiPull = Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()) / metZero.phiUncertainty();
							invalidSvFit3p_MetMinus3pNeutrinoXPull.at(t).Fill( pxPull , w);
							invalidSvFit3p_MetMinus3pNeutrinoXPull.at(t).Fill( pyPull , w);
							invalidSvFit3p_MetMinus3pNeutrinoPhiPull.at(t).Fill( phiPull , w);
						}
					}

					MetMinus3pNeutrino_DeltaPhiMuon.at(t).Fill( Tools::DeltaPhi(metZero.phi(), Ntp->Muon_p4(selMuon).Phi()) );
					if (not svfObjZero->isValid()) {
						MetMinus3pNeutrino_PhiResol_FailedSVFit.at(t).Fill( Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
						MetMinus3pNeutrino_DeltaPhiMuon_FailedSVFit.at(t).Fill( Tools::DeltaPhi(metZero.phi(), Ntp->Muon_p4(selMuon).Phi()) );
					}

					svFit3p_MassLMax_ZeroSol.at(t).Fill(svfObjZero->get_massLmax(), w);
					svFit3p_SolutionMatrixByLMax.at(t).Fill(true3ProngAmbig, 0., w);
				}
				vis3pMassResol_ZeroSol.at(t).Fill((visMassZero - trueResonanceMass)/trueResonanceMass, w);

				mtt_calculateZnu_ZeroSol.at(t).Fill(mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::zero, selMuon, metZero), w);
				mtt_calculateZnuNoAmb_ZeroSol.at(t).Fill(mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::zero, selMuon, metZero), w);
				mtt_3prongMuonRawMET_ZeroSol.at(t).Fill(mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::zero, selMuon, metZero), w);
				mtt_projectMetOnMu_ZeroSol.at(t).Fill(mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::zero, selMuon, metZero), w);
				mtt_projectMetOnMuRotate_ZeroSol.at(t).Fill(mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::zero, selMuon, metZero), w);

				mttCalculateZnuResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);
				mttCalculateZnuNoAmbResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);
				mtt3prongMuonRawMETResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);
				mttprojectMetOnMuRotateResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);
				mttprojectMetOnMuResol_vs_MetPhiResol.at(t).Fill(fabs(Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi())), mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);

				double tau3pEresol = (TPResults.getTauZero().LV().E() - trueTauP4.E())/trueTauP4.E();
				mttCalculateZnuResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol, 		mtt_mu3prong_calculateZnu(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);
				mttCalculateZnuNoAmbResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,	mtt_mu3prong_calculateZnuNoAmb(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);
				mtt3prongMuonRawMETResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,		mtt_mu3prong_UseRaw3prongMuonMET(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);
				mttprojectMetOnMuRotateResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,	mtt_mu3prong_projectMetOnMuAndRotate(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);
				mttprojectMetOnMuResol_vs_Tau3pEResol.at(t).Fill(tau3pEresol,		mtt_mu3prong_projectMetOnMu(TPResults, MultiProngTauSolver::zero, selMuon, metZero) - trueResonanceMass, w);

				if (Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt() > 0.001)
					Tau3p_Neutrino_PtResol.at(t).Fill( (TPResults.getNeutrinoZero().LV().Pt() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt())/Ntp->MCTau_invisiblePart(i_matchedMCTau).Pt(), w);
				Tau3p_Neutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoZero().LV().Phi(), Ntp->MCTau_invisiblePart(i_matchedMCTau).Phi()), w);

				if (TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::E, LorentzVectorParticle::E) > 0.001)
					Tau3p_Neutrino_EPull.at(t).Fill( (TPResults.getNeutrinoZero().LV().E() - Ntp->MCTau_invisiblePart(i_matchedMCTau).E())/TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::E, LorentzVectorParticle::E), w);
				if (TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px) > 0.001)
					Tau3p_Neutrino_XPull.at(t).Fill( (TPResults.getNeutrinoZero().LV().Px() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Px()) / TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::px, LorentzVectorParticle::px), w);
				if (TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py) > 0.001)
					Tau3p_Neutrino_YPull.at(t).Fill( (TPResults.getNeutrinoZero().LV().Py() - Ntp->MCTau_invisiblePart(i_matchedMCTau).Py()) / TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::py, LorentzVectorParticle::py), w);
				double phiUnc = Tools::phiUncertainty(TPResults.getNeutrinoZero().LV().Px(), TPResults.getNeutrinoZero().LV().Py(), TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::px),
						TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::py,LorentzVectorParticle::py), TPResults.getNeutrinoZero().Covariance(LorentzVectorParticle::px,LorentzVectorParticle::py));
				if (phiUnc > 0.001)
					Tau3p_Neutrino_PhiPull.at(t).Fill( Tools::DeltaPhi(TPResults.getNeutrinoZero().LV(), Ntp->MCTau_invisiblePart(i_matchedMCTau)) / phiUnc, w);

				if (Ntp->MCTau_invisiblePart(i_otherTau).Pt() > 0.001){
					MetMinus3pNeutrino_PtResol.at(t).Fill( (metZero.et() - Ntp->MCTau_invisiblePart(i_otherTau).Pt())/Ntp->MCTau_invisiblePart(i_otherTau).Pt(), w);
					MetMinus3pNeutrino_XResol.at(t).Fill( metZero.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px(), w);
					MetMinus3pNeutrino_YResol.at(t).Fill( metZero.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py(), w);
					MetMinus3pNeutrino_PhiResol.at(t).Fill( Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
					MetMinus3pNeutrino_PhiPtResol.at(t).Fill( Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()),(metZero.et() - Ntp->MCTau_invisiblePart(i_otherTau).Pt())/Ntp->MCTau_invisiblePart(i_otherTau).Pt(), w);

					if (runSVFit_ && not svfObjZero->isValid())
						MetMinus3pNeutrino_PhiResol_FailedSVFit.at(t).Fill( Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()), w);
					MetMinus3pNeutrino_DeltaPhiMuon.at(t).Fill( Tools::DeltaPhi(metZero.phi(), Ntp->Muon_p4(selMuon).Phi()) );
				}
				if(metZero.significanceXX() != 0 && metZero.significanceYY() != 0){
					double pxPull = (metZero.ex() - Ntp->MCTau_invisiblePart(i_otherTau).Px()) / sqrt(metZero.significanceXX());
					double pyPull = (metZero.ey() - Ntp->MCTau_invisiblePart(i_otherTau).Py()) / sqrt(metZero.significanceYY());
					double phiPull = Tools::DeltaPhi(metZero.phi(), Ntp->MCTau_invisiblePart(i_otherTau).Phi()) / metZero.phiUncertainty();
					MetMinus3pNeutrino_PxPull.at(t).Fill( pxPull , w);
					MetMinus3pNeutrino_PyPull.at(t).Fill( pyPull , w);
					MetMinus3pNeutrino_PxPyPull.at(t).Fill( pxPull, pyPull, w);
					MetMinus3pNeutrino_PhiPull.at(t).Fill( phiPull , w);
				}
				else
					Logger(Logger::Warning) << "MET-zero covariance is zero!" << std::endl;
			}
		}
	}

}

TLorentzVector ZeroJet3Prong::getTauLV(const TPTRObject& tau, int tau_ambiguity){
	TLorentzVector lv_tau;
	switch (tau_ambiguity) {
		case MultiProngTauSolver::plus:
			lv_tau = tau.getTauPlus().LV();
			break;
		case MultiProngTauSolver::minus:
			lv_tau = tau.getTauMinus().LV();
			break;
		case MultiProngTauSolver::zero:
			lv_tau = tau.getTauZero().LV();
			break;
		default:
			Logger(Logger::Warning) << tau_ambiguity << " is an invalid tau ambiguity."<< std::endl;
			break;
	}
	return lv_tau;
}

double ZeroJet3Prong::mtt_mu3prong_projectMetOnMu(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met){
	const TLorentzVector lv_tau = getTauLV(tau, tau_ambiguity);
	const TLorentzVector lv_mu	= Ntp->Muon_p4(i_mu);
	const TLorentzVector lv_met(met.ex(), met.ey(), 0., met.et());

	// project MET on muon in transverse plane to obtain transverse part of neutrino vector
	double s = (lv_met.Vect().Dot(lv_mu.Vect())) / lv_mu.Pt();

	TLorentzVector p_nu;
	p_nu.SetVect( (s / lv_mu.Pt()) * lv_mu.Vect());
	p_nu.SetE( p_nu.P() );

	return (lv_tau + lv_mu + p_nu).M();
}

double ZeroJet3Prong::mtt_mu3prong_projectMetOnMuAndRotate(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met){
	const TLorentzVector lv_tau = getTauLV(tau, tau_ambiguity);
	const TLorentzVector lv_mu	= Ntp->Muon_p4(i_mu);
	const TLorentzVector lv_met(met.ex(), met.ey(), 0., met.et());

	// project MET on muon in transverse plane to obtain transverse part of neutrino vector
	double s = (lv_met.Vect().Dot(lv_mu.Vect())) / lv_mu.Pt();

	TLorentzVector p_nu;
	p_nu.SetVect( (s / lv_mu.P()) * lv_mu.Vect()); // difference to other version is usage of P() instead of Pt() here
	p_nu.SetE( p_nu.P() );

	return (lv_tau + lv_mu + p_nu).M();
}

double ZeroJet3Prong:: mtt_mu3prong_UseRaw3prongMuonMET(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met){
	const TLorentzVector lv_tau = getTauLV(tau, tau_ambiguity);
	const TLorentzVector lv_mu	= Ntp->Muon_p4(i_mu);
	const TLorentzVector lv_met(met.ex(), met.ey(), 0., met.et());

	return (lv_tau + lv_mu + lv_met).M();
}

double ZeroJet3Prong::mtt_mu3prong_calculateZnu(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met){
	const TLorentzVector lv_tau = getTauLV(tau, tau_ambiguity);
	const TLorentzVector lv_mu	= Ntp->Muon_p4(i_mu);
	TVector2 p_met(met.ex(), met.ey());

/*
	Logger(Logger::Info) << "Muon: " << lv_mu.X() << " " <<  lv_mu.Y() << " " <<  lv_mu.Z() << std::endl;
	Logger(Logger::Info) << "MET:  " << met.ex() << " " <<  met.ey() <<  std::endl;
	Logger(Logger::Info) << "deltaPhi = " << Tools::DeltaPhi(lv_mu.Phi(), met.phi()) <<  std::endl;
*/

	double m_mu 	= 0.10565837;
	double m_tau	= PDGInfo::tau_mass();
	double dm2		= pow(m_tau, 2) - pow(m_mu, 2);

	double mtt = 0;

	// solve quadratic equation for z component of neutrino: C*z^2 - A*z + B = 0
	double A = lv_mu.Z() * (dm2 + 2*lv_mu.X()*p_met.X() + 2*lv_mu.Y()*p_met.Y());
	double B = pow(lv_mu.E(),2)*( pow(p_met.X(),2) + pow(p_met.Y(),2) ) - pow(dm2,2) - pow(lv_mu.X()*p_met.X(),2) - pow(lv_mu.Y()*p_met.Y(),2)
			   - dm2*(lv_mu.X()*p_met.X() + lv_mu.Y()*p_met.Y()) - 2*lv_mu.X()*lv_mu.Y()*p_met.X()*p_met.Y();
	double C = pow(lv_mu.E(),2) - pow(lv_mu.Z(),2);

	double p1 = 0.5*A/C;
	double p2 = B/C;

//	Logger(Logger::Info) << p1 << " +/- sqrt " << pow(p1,2) - p2 << std::endl;

	// two solutions
	if( pow(p1,2) >= p2){
//		Logger(Logger::Info) << "Two solutions" << std::endl;
		std::pair<double, double> z;
		z.first		= p1 + sqrt(pow(p1,2) - p2);
		z.second	= p1 - sqrt(pow(p1,2) - p2);

		std::pair<TVector3, TVector3> p3_nu;
		std::pair<TLorentzVector, TLorentzVector> p_nu;

		p3_nu.first = TVector3(p_met.X(), p_met.Y(), z.first);
		p3_nu.second = TVector3(p_met.X(), p_met.Y(), z.second);
		p_nu.first = TLorentzVector(p3_nu.first, p3_nu.first.Mag());
		p_nu.second = TLorentzVector(p3_nu.second, p3_nu.second.Mag());

		// select solution where angle to muon is smallest
		std::pair<double,double> dR = std::make_pair(lv_mu.DeltaR(p_nu.first), lv_mu.DeltaR(p_nu.second));
		std::pair<double,double> invM = std::make_pair( (lv_mu + p_nu.first + lv_tau).M(), (lv_mu + p_nu.second + lv_tau).M());
		deltaMtt_calculateZnu.at(t).Fill(invM.first - invM.second, w);
		if ( dR.first < dR.second ){
			//Logger(Logger::Info) << "Selected + solution (dR = " << dR.first << ", mtt = " << invM.first << ") over - solution (dR = "
			//		<< dR.second << ", mtt = " << invM.second  << ")" << std::endl;
			mtt = invM.first;
		}
		else if ( dR.first > dR.second ) {
			//Logger(Logger::Info) << "Selected - solution (dR = " << dR.second << ", mtt = " << invM.second << ") over + solution (dR = "
			//				<< dR.first << ", mtt = " << invM.first  << ")" << std::endl;
			mtt = invM.second;
		}
		else { // happens when p2 == 0 (should never occur using double precision)
			Logger(Logger::Warning) << "Both solutions are identical = " << dR.first << std::endl;
		}

		if ( fabs((lv_mu + p_nu.first).M() - m_tau) > 0.001 || fabs((lv_mu + p_nu.second).M() - m_tau) > 0.001 )
			Logger(Logger::Warning) << "Unphysical tau masses calculated: " << (lv_mu + p_nu.first).M() << " " << (lv_mu + p_nu.second).M() << std::endl;
	}
	else {
//		Logger(Logger::Info) << "No solution" << std::endl;
	// no solution: assume sqrt to be 0
/*
		double z = p1;
		TVector3 p3_nu(p_met.X(), p_met.Y(), z);
		TLorentzVector p_nu(p3_nu, p3_nu.Mag());

		mtt = (lv_mu + p_nu + lv_tau).M();

		if ( fabs((lv_mu + p_nu).M() - m_tau) > 0.001)
				Logger(Logger::Warning) << "Unphysical tau mass calculated: " << (lv_mu + p_nu).M() << std::endl;
*/

		//no solution: discard event
		mtt = 0;
	}

	return mtt;
}

double ZeroJet3Prong::mtt_mu3prong_calculateZnuNoAmb(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met){
	const TLorentzVector lv_tau = getTauLV(tau, tau_ambiguity);
	const TLorentzVector lv_mu	= Ntp->Muon_p4(i_mu);
	TVector2 p_met(met.ex(), met.ey());

	double m_mu 	= 0.10565837;
	double m_tau	= PDGInfo::tau_mass();
	double dm2		= pow(m_tau, 2) - pow(m_mu, 2);

	double mtt = 0;

	// solve quadratic equation for z component of neutrino: C*z^2 - A*z + B = 0
	double A = lv_mu.Z() * (dm2 + 2*lv_mu.X()*p_met.X() + 2*lv_mu.Y()*p_met.Y());
//	double B = pow(lv_mu.E(),2)*( pow(p_met.X(),2) + pow(p_met.Y(),2) ) - pow(dm2,2) - pow(lv_mu.X()*p_met.X(),2) - pow(lv_mu.Y()*p_met.Y(),2)
//			   - dm2*(lv_mu.X()*p_met.X() + lv_mu.Y()*p_met.Y()) - 2*lv_mu.X()*lv_mu.Y()*p_met.X()*p_met.Y();
	double C = pow(lv_mu.E(),2) - pow(lv_mu.Z(),2);

	double p1 = 0.5*A/C;
	// double p2 = B/C; // completely ignore the squareroot

	double z = p1;

	TVector3 p3_nu(p_met.X(), p_met.Y(), z);
	TLorentzVector p_nu(p3_nu, p3_nu.Mag());

	mtt = (lv_mu + p_nu+ lv_tau).M();

	return mtt;
}
