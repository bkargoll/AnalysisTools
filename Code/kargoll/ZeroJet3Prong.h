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

	std::vector<TH1D> vis3pMass_PlusSol;
	std::vector<TH1D> vis3pMass_MinusSol;
	std::vector<TH1D> vis3pMass_ZeroSol;
	std::vector<TH1D> vis3pMass_TrueSol;
	std::vector<TH1D> vis3pMassResol_PlusSol;
	std::vector<TH1D> vis3pMassResol_MinusSol;
	std::vector<TH1D> vis3pMassResol_ZeroSol;
	std::vector<TH1D> vis3pMassResol_TrueSol;

	std::vector<TH1D> svFit3pMass;
	std::vector<TH1D> svFit3pMassResol;

	std::vector<TH1D> svFit3p_TrueSolByHigherLMax;
	std::vector<TH2D> svFit3p_SolutionMatrixByLMax;

	std::vector<TH2D> svFit3pMassResol_vs_Tau3pEResol;
};

#endif /* ZeroJet3Prong_H_ */
