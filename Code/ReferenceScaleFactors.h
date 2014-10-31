/*
 * ReferenceScaleFactors.h
 *
 *  Created on: Aug 13, 2014
 *      Author: nehrkorn
 */

#ifndef REFERENCESCALEFACTORS_H_
#define REFERENCESCALEFACTORS_H_

#include "TFile.h"
#include "TH2F.h"
#include "TLorentzVector.h"

class ReferenceScaleFactors {

public:

	ReferenceScaleFactors(int runType);
	virtual ~ReferenceScaleFactors();

	// Muon Id scale factors
	double MuonIdTight2012(TLorentzVector vect);
	double MuonIdUncTight2012(TLorentzVector vect);
	double MuonIsoTight2012(TLorentzVector vect);
	double MuonIsoUncTight2012(TLorentzVector vect);
	double TrackingEfficiency2012(TLorentzVector vect);
	double TrackingEfficiencyUnc2012(TLorentzVector vect);
	double HiggsEMuId_Mu(TLorentzVector vect);
	double HiggsEMuIdUnc_Mu(TLorentzVector vect);

	// Electron Id scale factors
	double ElectronReconstruction2012(double Et, double Eta);
	double ElectronReconstructionUnc2012(double Et, double Eta);
	double ElectronIdTrig2012(double Et, double Eta);
	double ElectronIdTrigUnc2012(double Et, double Eta);
	double ElectronIdNonTrig2012(double Et, double Eta);
	double ElectronIdNonTrigUnc2012(double Et, double Eta);
	double ElectronEmbedding2012(double Et, double Eta);
	double HiggsTauTau_EMu_Id_E(double Et, double Eta);
	double HiggsTauTau_EMu_IdUnc_E(double Et, double Eta);

	// Tau scale factors

	//
	// Trigger efficiency scale factors
	//

	// ++++++++++ Final state: e+mu ++++++++++
	// use this for full trigger efficiency
	double HiggsWW_EMu_Trigger(TLorentzVector mu_vect, double e_et, double e_eta, TString path);
	// efficiencies of individual trigger legs
	double HiggsWW_EMu_SingleEle(double Et, double Eta);
	double HiggsWW_EMu_DoubleEleLeading(double Et, double Eta);
	double HiggsWW_EMu_DoubleEleTrailing(double Et, double Eta);
	double HiggsWW_EMu_SingleMu(TLorentzVector vect);
	double HiggsWW_EMu_DoubleMuLeading(TLorentzVector vect);
	double HiggsWW_EMu_DoubleMuTrailing(TLorentzVector vect);
	double HiggsTauTau_EMu_Trigger_Mu(TLorentzVector vect);
	double HiggsTauTau_EMu_TriggerUnc_Mu(TLorentzVector vect);
	double HiggsTauTau_EMu_Trigger_E(double Et, double Eta);
	double HiggsTauTau_EMu_TriggerUnc_E(double Et, double Eta);

	// Final state: mu+tau

private:
	//
	// Root files for scale factors
	//

	// Electron Id's
	TFile* ETrigIdEffFile;
	TFile* ENonTrigIdEffFile;
	TFile* ERecoEffFile;
	// Trigger efficiencies
	TFile* HWW_TriggerEfficiencies;

	//
	// Histograms for scale factors
	//

	// Electron Id's
	TH2D* ElectronTrigEff;
	TH2D* ElectronNonTrigEff;
	TH2D* ElectronRecoEff;
	// Trigger efficiencies
	TH1D* HiggsWW_EMu_SingleEle15;
	TH1D* HiggsWW_EMu_SingleEle25;
	TH1D* HiggsWW_EMu_DoubleEleLead15;
	TH1D* HiggsWW_EMu_DoubleEleLead25;
	TH1D* HiggsWW_EMu_DoubleEleTrail15;
	TH1D* HiggsWW_EMu_DoubleEleTrail25;
	TH1D* HiggsWW_EMu_SingleMu08;
	TH1D* HiggsWW_EMu_SingleMu12;
	TH1D* HiggsWW_EMu_SingleMu21;
	TH1D* HiggsWW_EMu_SingleMu25;
	TH1D* HiggsWW_EMu_DoubleMuLead12;
	TH1D* HiggsWW_EMu_DoubleMuLead21;
	TH1D* HiggsWW_EMu_DoubleMuLead25;
	TH1D* HiggsWW_EMu_DoubleMuTrail12;
	TH1D* HiggsWW_EMu_DoubleMuTrail21;
	TH1D* HiggsWW_EMu_DoubleMuTrail25;

};



#endif /* REFERENCESCALEFACTORS_H_ */
