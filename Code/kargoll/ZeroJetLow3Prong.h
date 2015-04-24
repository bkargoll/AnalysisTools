/*
 * ZeroJetLow3Prong.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ZeroJetLow3Prong_H_
#define ZeroJetLow3Prong_H_

#include "HToTaumuTauh.h"

class ZeroJetLow3Prong: public HToTaumuTauh {
public:
	ZeroJetLow3Prong(TString Name_, TString id_);
	virtual ~ZeroJetLow3Prong();

	virtual void Configure();

protected:
	virtual void Setup();
	virtual void doEvent();
	virtual void Store_ExtraDist();

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
};

#endif /* ZeroJetLow3Prong_H_ */
