/*
 * ZeroJet3Prong.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ZeroJet3Prong_H_
#define ZeroJet3Prong_H_

#include "Category.h"
class TPTRObject;

class ZeroJet3Prong: public Category {
public:
	ZeroJet3Prong(TString Name_, TString id_);
	virtual ~ZeroJet3Prong();

	// enumerator for the cuts in this category
	enum cuts_ZeroJet3Prong {
		NJet = HToTaumuTauh::CatCut1,
		DecayMode,
		SigmaSV,
		NCuts
	};

	TLorentzVector getTauLV(const TPTRObject& tau, int tau_ambiguity);
	double mtt_mu3prong_projectMetOnMu(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met);
	double mtt_mu3prong_UseRaw3prongMuonMET(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met);
	double mtt_mu3prong_projectMetOnMuAndRotate(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met);
	double mtt_mu3prong_calculateZnu(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met);
	double mtt_mu3prong_calculateZnuNoAmb(const TPTRObject& tau, int tau_ambiguity, unsigned i_mu, const objects::MET& met);

protected:
	void categoryConfiguration();
	bool categorySelection();
	void categoryExtradist();
	void categoryPlotting();

	// flag to enable/disable SVFit calculation
	bool runSVFit_;

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
	std::vector<TH1D> Tau3p_Zero_RotSignificance;
	std::vector<TH1D> Tau3p_FlightLengthSignificanceDefPV;

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

	std::vector<TH1D> svFit3pMass_PlusSol_notCleaned;
	std::vector<TH1D> svFit3pMass_MinusSol_notCleaned;
	std::vector<TH1D> svFit3pMass_ZeroSol_notCleaned;
	std::vector<TH1D> svFit3pMass_TrueSol_notCleaned;
	std::vector<TH1D> svFit3pMass_TruePlusSol_notCleaned;
	std::vector<TH1D> svFit3pMass_TrueMinusSol_notCleaned;
	std::vector<TH1D> svFit3pMassResol_PlusSol_notCleaned;
	std::vector<TH1D> svFit3pMassResol_MinusSol_notCleaned;
	std::vector<TH1D> svFit3pMassResol_ZeroSol_notCleaned;
	std::vector<TH1D> svFit3pMassResol_TrueSol_notCleaned;

	std::vector<TH1D> svFit3pMass_PlusSol_invalid;
	std::vector<TH1D> svFit3pMass_MinusSol_invalid;
	std::vector<TH1D> svFit3pMass_ZeroSol_invalid;
	std::vector<TH1D> svFit3pMass_TrueSol_invalid;
	std::vector<TH1D> svFit3pMass_TruePlusSol_invalid;
	std::vector<TH1D> svFit3pMass_TrueMinusSol_invalid;
	std::vector<TH1D> svFit3pMassResol_PlusSol_invalid;
	std::vector<TH1D> svFit3pMassResol_MinusSol_invalid;
	std::vector<TH1D> svFit3pMassResol_ZeroSol_invalid;
	std::vector<TH1D> svFit3pMassResol_TrueSol_invalid;

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
	std::vector<TH1D> Tau3p_Neutrino_EPull;
	std::vector<TH1D> Tau3p_Neutrino_XPull;
	std::vector<TH1D> Tau3p_Neutrino_YPull;
	std::vector<TH1D> Tau3p_Neutrino_PhiPull;

	std::vector<TH1D> MetMinus3pNeutrino_PtResol;
	std::vector<TH1D> MetMinus3pNeutrino_PhiResol;
	std::vector<TH1D> MetMinus3pNeutrino_XResol;
	std::vector<TH1D> MetMinus3pNeutrino_YResol;
	std::vector<TH2D> MetMinus3pNeutrino_PhiPtResol;

	std::vector<TH1D> MetMinus3pNeutrino_PhiResol_FailedSVFit;
	std::vector<TH1D> MetMinus3pNeutrino_DeltaPhiMuon;
	std::vector<TH1D> MetMinus3pNeutrino_DeltaPhiMuon_FailedSVFit;

	std::vector<TH1D> svFitMass_Default;
	std::vector<TH1D> svFit3pMass;
	std::vector<TH1D> svFit3pMassResol;

	std::vector<TH1D> svFit3pMass_notCleaned;
	std::vector<TH1D> svFit3pMassResol_notCleaned;
	std::vector<TH1D> svFit3pMass_invalid;
	std::vector<TH1D> svFit3pMassResol_invalid;

	std::vector<TH1D> svFit3p_MassLMax_TrueSol;
	std::vector<TH1D> svFit3p_MassLMax_WrongSol;
	std::vector<TH1D> svFit3p_MassLMax_ZeroSol;
	std::vector<TH1D> svFit3p_MassLMax_PlusSol;
	std::vector<TH1D> svFit3p_MassLMax_MinusSol;
	std::vector<TH1D> svFit3p_TrueSolByHigherLMax;
	std::vector<TH2D> svFit3p_SolutionMatrixByLMax;

	std::vector<TH1D> invalidSvFit3p_Tau3pNeutrinoXPull;
	std::vector<TH1D> invalidSvFit3p_Tau3pNeutrinoYPull;
	std::vector<TH1D> invalidSvFit3p_Tau3pNeutrinoPhiPull;
	std::vector<TH1D> invalidSvFit3p_MetMinus3pNeutrinoXPull;
	std::vector<TH1D> invalidSvFit3p_MetMinus3pNeutrinoYPull;
	std::vector<TH1D> invalidSvFit3p_MetMinus3pNeutrinoPhiPull;

	std::vector<TH1D> MetPxPull;
	std::vector<TH1D> MetPyPull;
	std::vector<TH2D> MetPxPyPull;
	std::vector<TH1D> MetPhiPull;
	std::vector<TH1D> MetMinus3pNeutrino_PxPull;
	std::vector<TH1D> MetMinus3pNeutrino_PyPull;
	std::vector<TH2D> MetMinus3pNeutrino_PxPyPull;
	std::vector<TH1D> MetMinus3pNeutrino_PhiPull;

	std::vector<TH1D> mtt_calculateZnu_PlusSol;
    std::vector<TH1D> mtt_calculateZnu_MinusSol;
	std::vector<TH1D> mtt_calculateZnu_ZeroSol;
	std::vector<TH1D> mtt_calculateZnu_TrueSol;
	std::vector<TH1D> mtt_calculateZnuNoAmb_PlusSol;
    std::vector<TH1D> mtt_calculateZnuNoAmb_MinusSol;
	std::vector<TH1D> mtt_calculateZnuNoAmb_ZeroSol;
	std::vector<TH1D> mtt_calculateZnuNoAmb_TrueSol;
	std::vector<TH1D> mtt_3prongMuonRawMET_PlusSol;
	std::vector<TH1D> mtt_3prongMuonRawMET_MinusSol;
	std::vector<TH1D> mtt_3prongMuonRawMET_ZeroSol;
	std::vector<TH1D> mtt_3prongMuonRawMET_TrueSol;
	std::vector<TH1D> mtt_projectMetOnMu_PlusSol;
	std::vector<TH1D> mtt_projectMetOnMu_MinusSol;
	std::vector<TH1D> mtt_projectMetOnMu_ZeroSol;
	std::vector<TH1D> mtt_projectMetOnMu_TrueSol;
	std::vector<TH1D> mtt_projectMetOnMuRotate_PlusSol;
	std::vector<TH1D> mtt_projectMetOnMuRotate_MinusSol;
	std::vector<TH1D> mtt_projectMetOnMuRotate_ZeroSol;
	std::vector<TH1D> mtt_projectMetOnMuRotate_TrueSol;

	std::vector<TH1D> deltaMtt_calculateZnu;

	std::vector<TH2D> svFit3pMassResol_vs_MetPhiResol;
	std::vector<TH2D> mttCalculateZnuResol_vs_MetPhiResol;
	std::vector<TH2D> mttCalculateZnuNoAmbResol_vs_MetPhiResol;
	std::vector<TH2D> mtt3prongMuonRawMETResol_vs_MetPhiResol;
	std::vector<TH2D> mttprojectMetOnMuRotateResol_vs_MetPhiResol;
	std::vector<TH2D> mttprojectMetOnMuResol_vs_MetPhiResol;

	std::vector<TH2D> svFit3pMassResol_vs_Tau3pEResol;
	std::vector<TH2D> mttCalculateZnuResol_vs_Tau3pEResol;
	std::vector<TH2D> mttCalculateZnuNoAmbResol_vs_Tau3pEResol;
	std::vector<TH2D> mtt3prongMuonRawMETResol_vs_Tau3pEResol;
	std::vector<TH2D> mttprojectMetOnMuRotateResol_vs_Tau3pEResol;
	std::vector<TH2D> mttprojectMetOnMuResol_vs_Tau3pEResol;
};

#endif /* ZeroJet3Prong_H_ */
