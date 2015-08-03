/*
 * ZeroJetLow3Prong.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ZeroJetLow3Prong_H_
#define ZeroJetLow3Prong_H_

#include "Category.h"

class ZeroJetLow3Prong: public Category {
public:
	ZeroJetLow3Prong(TString Name_, TString id_);
	virtual ~ZeroJetLow3Prong();

	// enumerator for the cuts in this category
	enum cuts_ZeroJetLow3Prong {
		NJet = HToTaumuTauh::CatCut1,
		TauPt,
		DecayMode,
		SigmaSV,
		NCuts
	};

protected:
	void categoryConfiguration();
	bool categorySelection();
	void categoryExtradist();
	void categoryPlotting();

	// svfit
	SVFitStorage svfitstor_Tau3pPlus;
	SVFitStorage svfitstor_Tau3pMinus;
	SVFitStorage svfitstor_Tau3pZero;

	// histograms
	std::vector<TH1D> Tau3p_Plus_Pt;
	std::vector<TH1D> Tau3p_Plus_Eta;
	std::vector<TH1D> Tau3p_Plus_Phi;
	std::vector<TH1D> Tau3p_Plus_E;
	std::vector<TH1D> Tau3p_Minus_Pt;
	std::vector<TH1D> Tau3p_Minus_Eta;
	std::vector<TH1D> Tau3p_Minus_Phi;
	std::vector<TH1D> Tau3p_Minus_E;
	std::vector<TH1D> Tau3p_Zero_Pt;
	std::vector<TH1D> Tau3p_Zero_Eta;
	std::vector<TH1D> Tau3p_Zero_Phi;
	std::vector<TH1D> Tau3p_Zero_E;

	std::vector<TH1D> Tau3p_fit_Pt;
	std::vector<TH1D> Tau3p_fit_Eta;
	std::vector<TH1D> Tau3p_fit_Phi;
	std::vector<TH1D> Tau3p_fit_E;

	std::vector<TH1D> TauMu_fit_Pt;
	std::vector<TH1D> TauMu_fit_Eta;
	std::vector<TH1D> TauMu_fit_Phi;
	std::vector<TH1D> TauMu_fit_E;

	std::vector<TH1D> EventFit_converged;
	std::vector<TH1D> EventFit_Chi2;
	std::vector<TH1D> EventFit_nIterations;
	std::vector<TH1D> EventFit_deltaCSum;

	std::vector<TH1D> EventFit_Res_Pt;
	std::vector<TH1D> EventFit_Res_Eta;
	std::vector<TH1D> EventFit_Res_Phi;
	std::vector<TH1D> EventFit_Res_E;

	std::vector<TH1D> svFit3pMass_PlusSol;
	std::vector<TH1D> svFit3pMass_MinusSol;
	std::vector<TH1D> svFit3pMass_ZeroSol;
	std::vector<TH1D> svFit3pMass_TrueSol;
	std::vector<TH1D> svFit3pMass_TruePlusSol;
	std::vector<TH1D> svFit3pMass_TrueMinusSol;
	std::vector<TH1D> svFit3pMassResol_PlusSol;
	std::vector<TH1D> svFit3pMassResol_MinusSol;
	std::vector<TH1D> svFit3pMassResol_ZeroSol;
	std::vector<TH1D> svFit3pMassResol_TrueSol;

	std::vector<TH2D> Tau_Pt_3pTrueVsSVfit3pTrueSol;
	std::vector<TH1D> Tau_Pt_3pTrueSVfit3pTrueSolResol;
	std::vector<TH2D> Tau_Pt_3pTrueVsSVfit3pZeroSol;
	std::vector<TH1D> Tau_Pt_3pTrueSVfit3pZeroSolResol;

	std::vector<TH1D> vis3pMass_PlusSol;
	std::vector<TH1D> vis3pMass_MinusSol;
	std::vector<TH1D> vis3pMass_ZeroSol;
	std::vector<TH1D> vis3pMass_TrueSol;
	std::vector<TH1D> vis3pMassResol_PlusSol;
	std::vector<TH1D> vis3pMassResol_MinusSol;
	std::vector<TH1D> vis3pMassResol_ZeroSol;
	std::vector<TH1D> vis3pMassResol_TrueSol;

	std::vector<TH1D> Tau3p_True_Pt;
	std::vector<TH1D> Tau3p_True_Eta;
	std::vector<TH1D> Tau3p_True_Phi;
	std::vector<TH1D> Tau3p_True_E;
	std::vector<TH1D> Tau3p_True_Pt_Resol;
	std::vector<TH1D> Tau3p_True_Eta_Resol;
	std::vector<TH1D> Tau3p_True_Phi_Resol;
	std::vector<TH1D> Tau3p_True_E_Resol;

	std::vector<TH1D> TrueSignalMET;
	std::vector<TH1D> MetPtResol;
	std::vector<TH1D> MetPhiResol;
	std::vector<TH1D> MetPxResol;
	std::vector<TH1D> MetPyResol;
	std::vector<TH1D> Tau3p_Neutrino_PtResol;
	std::vector<TH1D> Tau3p_Neutrino_PhiResol;
	std::vector<TH1D> Tau3p_Neutrino_XResol;
	std::vector<TH1D> Tau3p_Neutrino_YResol;
	std::vector<TH1D> MetMinus3pNeutrino_PtResol;
	std::vector<TH1D> MetMinus3pNeutrino_PhiResol;
	std::vector<TH1D> MetMinus3pNeutrino_XResol;
	std::vector<TH1D> MetMinus3pNeutrino_YResol;

	std::vector<TH1D> svFitMass_Default;
	std::vector<TH1D> svFit3pMass;
	std::vector<TH1D> svFit3pMassResol;

	std::vector<TH1D> MetPxPull;
	std::vector<TH1D> MetPyPull;
	std::vector<TH2D> MetPxPyPull;
	std::vector<TH1D> MetMinus3pNeutrino_PxPull;
	std::vector<TH1D> MetMinus3pNeutrino_PyPull;
	std::vector<TH2D> MetMinus3pNeutrino_PxPyPull;

};

#endif /* ZeroJetLow3Prong_H_ */
