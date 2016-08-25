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

	runSVFit_ = true;
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

	svFit3p_TrueSolByHigherLMax= HConfig.GetTH1D(Name+"_svFit3p_TrueSolByHigherLMax",	"svFit3p_TrueSolByHigherLMax", 2, -0.5, 1.5, "correct by L^{m}_{max}");
	svFit3p_SolutionMatrixByLMax = HConfig.GetTH2D(Name+"_svFit3p_SolutionMatrixByLMax", "svFit3p_SolutionMatrixByLMax", 3, -0.5, 2.5, 3, -0.5, 2.5, "True amb.", "higher L^{max}(SVfit)");

	vis3pMass_PlusSol		= HConfig.GetTH1D(Name+"_vis3pMass_PlusSol",	"vis3pMass_PlusSol",	80, 0., 400., "m_{vis}(#tau_{3p}^{+},#mu)/GeV");
	vis3pMass_MinusSol		= HConfig.GetTH1D(Name+"_vis3pMass_MinusSol",	"vis3pMass_MinusSol",	80, 0., 400., "m_{vis}(#tau_{3p}^{-},#mu)/GeV");
	vis3pMass_ZeroSol		= HConfig.GetTH1D(Name+"_vis3pMass_ZeroSol",	"vis3pMass_ZeroSol",	80, 0., 400., "m_{vis}(#tau_{3p}^{0},#mu)/GeV");
	vis3pMass_TrueSol		= HConfig.GetTH1D(Name+"_vis3pMass_TrueSol",	"vis3pMass_TrueSol",	80, 0., 400., "m_{vis}(#tau_{3p}^{true},#mu)/GeV");
	vis3pMassResol_PlusSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_PlusSol",	"vis3pMassResol_PlusSol",	50, -1., 3., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{+},#mu)");
	vis3pMassResol_MinusSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_MinusSol","vis3pMassResol_MinusSol",	50, -1., 3., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{-},#mu)");
	vis3pMassResol_ZeroSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_ZeroSol",	"vis3pMassResol_ZeroSol",	50, -1., 3., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{0},#mu)");
	vis3pMassResol_TrueSol	= HConfig.GetTH1D(Name+"_vis3pMassResol_TrueSol",	"vis3pMassResol_TrueSol",	50, -1., 3., "#frac{m_{vis} - m_{true}}{m_{true}}(#tau_{3p}^{true},#mu)");

	svFit3pMass = HConfig.GetTH1D(Name+"_svFit3pMass",	"svFit3pMass",	50, 0., 250., "m_{SVfit}(#tau_{3p},#mu)/GeV");;
	svFit3pMassResol = HConfig.GetTH1D(Name+"_svFit3pMassResol",	"svFit3pMassResol",	50, -1., 1., "#frac{m_{SVfit} - m_{true}}{m_{true}}(#tau_{3p},#mu)");

	svFit3pMassResol_vs_Tau3pEResol = HConfig.GetTH2D(Name+"_svFit3pMassResol_vs_Tau3pEResol", "svFit3pMassResol_vs_Tau3pEResol", 50, -1., 1., 50, -100., 100., "E resol. #frac{#tau_{3p}^{true sol.} - #tau_{gen}}{#tau_{gen}}", "m_{SVfit} - m_{true}(#tau_{3p},#mu) /GeV");
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

	Extradist1d.push_back(&svFit3p_TrueSolByHigherLMax);
	Extradist2d.push_back(&svFit3p_SolutionMatrixByLMax);

	Extradist1d.push_back(&vis3pMass_PlusSol);
	Extradist1d.push_back(&vis3pMass_MinusSol);
	Extradist1d.push_back(&vis3pMass_ZeroSol);
	Extradist1d.push_back(&vis3pMass_TrueSol);
	Extradist1d.push_back(&vis3pMassResol_PlusSol);
	Extradist1d.push_back(&vis3pMassResol_MinusSol);
	Extradist1d.push_back(&vis3pMassResol_ZeroSol);
	Extradist1d.push_back(&vis3pMassResol_TrueSol);

	Extradist1d.push_back(&svFit3pMass);
	Extradist1d.push_back(&svFit3pMassResol);

	Extradist2d.push_back(&svFit3pMassResol_vs_Tau3pEResol);
}

bool ZeroJet3Prong::categorySelection(){
	bool categoryPass = true;
	std::vector<float> value_ZeroJet3Prong(NCuts,-10);
	std::vector<float> pass_ZeroJet3Prong(NCuts,true); // disable all selections in category (was false)

	Logger(Logger::Debug) << "Cut: Number of Jets" << std::endl;
	value_ZeroJet3Prong.at(NJet) = nJets_;
	//pass_ZeroJet3Prong.at(NJet) = ( value_ZeroJet3Prong.at(NJet) <= cut.at(NJet) );

	if (selTau == -1){
		// decayMode cut is set to true for nice N-0 and N-1 plots
		value_ZeroJet3Prong.at(DecayMode) = -10;
		//pass_ZeroJet3Prong.at(DecayMode) = true;
		// whole category is failing selection, to avoid NCat > 1
		//categoryPass = false;
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
		//pass_ZeroJet3Prong.at(DecayMode)		= (value_ZeroJet3Prong.at(DecayMode)==cut.at(DecayMode));
	}

	if (selTau == -1 || not pass_ZeroJet3Prong.at(DecayMode)) {
		value_ZeroJet3Prong.at(SigmaSV) = -999;
		//pass_ZeroJet3Prong.at(SigmaSV) = true;
	}
	else if ( not Ntp->PFTau_TIP_hassecondaryVertex(selTau) ){
		value_ZeroJet3Prong.at(SigmaSV) = -9;
		//pass_ZeroJet3Prong.at(SigmaSV) = false;
	}
	else {
		Logger(Logger::Debug) << "Cut: SV significance" << std::endl;
		value_ZeroJet3Prong.at(SigmaSV) = flightLengthSig_;
		//pass_ZeroJet3Prong.at(SigmaSV) = ( value_ZeroJet3Prong.at(SigmaSV) >= cut.at(SigmaSV) );
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

		if (TPResults.isAmbiguous()) {
			SVFitObject *svfObjPlus  = NULL;
			SVFitObject *svfObjMinus = NULL;
			bool plusLmaxIsHigher;
			if (runSVFit_){
				svfObjPlus = Ntp->getSVFitResult_MuTau3p(svfitstor_Tau3pPlus, "CorrMVAMuTau", selMuon, TPResults.getTauPlus().LV(), TPResults.getNeutrinoPlus(), "3ProngRecoPlus");
				svfObjMinus = Ntp->getSVFitResult_MuTau3p(svfitstor_Tau3pMinus, "CorrMVAMuTau", selMuon, TPResults.getTauMinus().LV(), TPResults.getNeutrinoMinus(), "3ProngRecoMinus");
				if (svfObjPlus->isValid())	svFit3pMass_PlusSol.at(t).Fill(svfObjPlus->get_mass(), w);
				if (svfObjMinus->isValid())	svFit3pMass_MinusSol.at(t).Fill(svfObjMinus->get_mass(), w);
				plusLmaxIsHigher = (svfObjPlus->get_massLmax() > svfObjMinus->get_massLmax());
			}

			double visMassPlus = (Ntp->Muon_p4(selMuon) + TPResults.getTauPlus().LV()).M();
			double visMassMinus = (Ntp->Muon_p4(selMuon) + TPResults.getTauMinus().LV()).M();
			vis3pMass_PlusSol.at(t).Fill( visMassPlus, w);
			vis3pMass_MinusSol.at(t).Fill( visMassMinus, w);

			objects::MET metPlus(Ntp, "CorrMVAMuTau");
			metPlus.subtractNeutrino(TPResults.getNeutrinoPlus());
			objects::MET metMinus(Ntp, "CorrMVAMuTau");
			metMinus.subtractNeutrino(TPResults.getNeutrinoMinus());


			if (isMC) {
				if (runSVFit_) {
					if (svfObjPlus->isValid())	svFit3pMassResol_PlusSol.at(t).Fill( (svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
					if (svfObjMinus->isValid())	svFit3pMassResol_MinusSol.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
				}
				vis3pMassResol_PlusSol.at(t).Fill( (visMassPlus - trueResonanceMass)/trueResonanceMass, w);
				vis3pMassResol_MinusSol.at(t).Fill( (visMassMinus - trueResonanceMass)/trueResonanceMass, w);

				double phiUnc(-1);
				switch (true3ProngAmbig) {
					case MultiProngTauSolver::plus:
						if (runSVFit_) {
							if(svfObjPlus->isValid()){
								svFit3pMass_TrueSol.at(t).Fill(svfObjPlus->get_mass(), w);
								svFit3pMass_TruePlusSol.at(t).Fill(svfObjPlus->get_mass(), w);
								svFit3pMassResol_TrueSol.at(t).Fill((svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
								svFit3pMass.at(t).Fill(svfObjPlus->get_mass(), w);
								svFit3pMassResol.at(t).Fill((svfObjPlus->get_mass() - trueResonanceMass) / trueResonanceMass, w);

								svFit3pMassResol_vs_Tau3pEResol.at(t).Fill((TPResults.getTauPlus().LV().E() - trueTauP4.E())/trueTauP4.E(), (svfObjPlus->get_mass() - trueResonanceMass), w);
							}
							svFit3p_TrueSolByHigherLMax.at(t).Fill(plusLmaxIsHigher,w);
							svFit3p_SolutionMatrixByLMax.at(t).Fill(true3ProngAmbig, plusLmaxIsHigher+1, w);
						}
						vis3pMass_TrueSol.at(t).Fill((Ntp->Muon_p4(selMuon) + TPResults.getTauPlus().LV()).M(), w);
						vis3pMassResol_TrueSol.at(t).Fill((visMassPlus - trueResonanceMass)/trueResonanceMass, w);

						break;

					case MultiProngTauSolver::minus:
						if (runSVFit_) {
							if (svfObjMinus->isValid()) {
								svFit3pMass_TrueSol.at(t).Fill(svfObjMinus->get_mass(), w);
								svFit3pMass_TrueMinusSol.at(t).Fill(svfObjMinus->get_mass(), w);
								svFit3pMassResol_TrueSol.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);
								svFit3pMass.at(t).Fill(svfObjMinus->get_mass(), w);
								svFit3pMassResol.at(t).Fill((svfObjMinus->get_mass() - trueResonanceMass) / trueResonanceMass, w);

								svFit3pMassResol_vs_Tau3pEResol.at(t).Fill((TPResults.getTauMinus().LV().E() - trueTauP4.E())/trueTauP4.E(), (svfObjMinus->get_mass() - trueResonanceMass), w);
							}
							svFit3p_TrueSolByHigherLMax.at(t).Fill(!plusLmaxIsHigher,w);
							svFit3p_SolutionMatrixByLMax.at(t).Fill(true3ProngAmbig, plusLmaxIsHigher+1, w);
						}
						vis3pMass_TrueSol.at(t).Fill((Ntp->Muon_p4(selMuon) + TPResults.getTauMinus().LV()).M(), w);
						vis3pMassResol_TrueSol.at(t).Fill((visMassMinus - trueResonanceMass)/trueResonanceMass, w);

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
				if (svfObjZero->isValid())	svFit3pMass_ZeroSol.at(t).Fill(svfObjZero->get_mass(), w);
			}

			double visMassZero = (Ntp->Muon_p4(selMuon) + TPResults.getTauZero().LV()).M();
			vis3pMass_ZeroSol.at(t).Fill(visMassZero, w);

			if (isMC) {
				objects::MET metZero(Ntp, "CorrMVAMuTau");
				metZero.subtractNeutrino(TPResults.getNeutrinoZero());
				//TLorentzVector genMet = Ntp->MCTau_invisiblePart(i_otherTau);
				//objects::MET metZero("gen", genMet.Px(), genMet.Py());
				if (runSVFit_) {
					if (svfObjZero->isValid()){
						svFit3pMassResol_ZeroSol.at(t).Fill((svfObjZero->get_mass() - trueResonanceMass) / trueResonanceMass, w);
						svFit3pMass.at(t).Fill(svfObjZero->get_mass(), w);
						svFit3pMassResol.at(t).Fill((svfObjZero->get_mass() - trueResonanceMass) / trueResonanceMass, w);

						svFit3pMassResol_vs_Tau3pEResol.at(t).Fill((TPResults.getTauZero().LV().E() - trueTauP4.E())/trueTauP4.E(), (svfObjZero->get_mass() - trueResonanceMass), w);
					}

					svFit3p_SolutionMatrixByLMax.at(t).Fill(true3ProngAmbig, 0., w);
				}
				vis3pMassResol_ZeroSol.at(t).Fill((visMassZero - trueResonanceMass)/trueResonanceMass, w);

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
