/*
 * ZeroJetLow3Prong.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "ZeroJetLow3Prong.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

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
}

void ZeroJetLow3Prong::doEvent(){
	Logger(Logger::Verbose) << std::endl;
	HToTaumuTauh::doEvent();
	// todo: split analysis part from plotting part in HToTaumuTauh

	Logger(Logger::Debug) << "Cut: Tau decay mode" << std::endl;
	if (selTau < 0){
		// set cuts to true for nice N-0 and N-1 plots
		value.at(CatCut3) = -10;
		pass.at(CatCut3) = true;
	}
	else{
		value.at(CatCut3)=Ntp->PFTau_hpsDecayMode(selTau);
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

		GlobalEventFit GEF(MuonTP, A1, Phi_Res, PV, PVCov);
		TPTRObject TPResults = GEF.getTPTRObject();
		LorentzVectorParticle test = TPResults.getA1();
		Logger(Logger::Debug) << test.LV().E() << std::endl;
		GEFObject Results = GEF.Fit();
		Logger(Logger::Debug) << Results.getChi2() << std::endl;

		// fill plots
		if (not TPResults.isValid() ) Logger(Logger::Error) << "3-prong tau object is not valid!" << std::endl;

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


}
