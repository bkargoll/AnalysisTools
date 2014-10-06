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

	ReferenceScaleFactors(TString basedir);
	virtual ~ReferenceScaleFactors();

//protected:

	// Muon scale factors
	double MuonIdTight2012(TLorentzVector vect);
	double MuonIdUncTight2012(TLorentzVector vect);
	double MuonIsoTight2012(TLorentzVector vect);
	double MuonIsoUncTight2012(TLorentzVector vect);
	double TrackingEfficiency2012(TLorentzVector vect);

	// Electron scale factors
	double ElectronReconstruction2012(double Et, double Eta);
	double ElectronReconstructionUnc2012(double Et, double Eta);
	double ElectronIdTrig2012(double Et, double Eta);
	double ElectronIdTrigUnc2012(double Et, double Eta);
	double ElectronIdNonTrig2012(double Et, double Eta);
	double ElectronIdNonTrigUnc2012(double Et, double Eta);

	// Tau scale factors

private:
	// Root files for scale factors
	TFile* ETrigIdEffFile;
	TFile* ENonTrigIdEffFile;
	TFile* ERecoEffFile;

	// Histograms for scale factors
	TH2D* ElectronTrigEff;
	TH2D* ElectronNonTrigEff;
	TH2D* ElectronRecoEff;

};



#endif /* REFERENCESCALEFACTORS_H_ */
