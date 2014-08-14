/*
 * ReferenceScaleFactors.cxx
 *
 *  Created on: Aug 13, 2014
 *      Author: nehrkorn
 */

#include "ReferenceScaleFactors.h"

///////////////////////////
//
// Constructor
//

ReferenceScaleFactors::ReferenceScaleFactors(){
	// Read root files
	TString base = "CommonFiles/";// todo: make sure that base is always known
	ETrigIdEffFile = new TFile(base+"ElectronEfficiencies_Run2012ReReco_53X_Trig.root");
	ENonTrigIdEffFile = new TFile(base+"ElectronEfficiencies_Run2012ReReco_53X_NonTrig.root");
	ERecoEffFile = new TFile(base+"Electrons_ScaleFactors_Reco_8TeV.root");
	// Get histograms
	ElectronTrigEff = (TH2D*)(ETrigIdEffFile->Get("electronsDATAMCratio_FO_ID_ISO"));
	ElectronNonTrigEff = (TH2D*)(ENonTrigIdEffFile->Get("h_electronScaleFactor_IdIsoSip"));
	ElectronRecoEff = (TH2D*)(ERecoEffFile->Get("h_electronScaleFactor_RECO"));
}

ReferenceScaleFactors::~ReferenceScaleFactors(){}

///////////////////////////
//
// Muon scale factors
//

double ReferenceScaleFactors::MuonIdTight2012(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=10. && pt<20.){
		if(eta<0.9) eff = 0.970;
		if(eta>=0.9 && eta<1.2) eff = 1.002;
		if(eta>=1.2 && eta<2.1) eff = 1.018;
		if(eta>=2.1 && eta<2.4) eff = 1.005;
	}
	if(pt>=20. && pt<25.){
		if(eta<0.9) eff = 0.989;
		if(eta>=0.9 && eta<1.2) eff = 0.994;
		if(eta>=1.2 && eta<2.1) eff = 1.000;
		if(eta>=2.1 && eta<2.4) eff = 0.998;
	}
	if(pt>=25. && pt<30.){
		if(eta<0.9) eff = 0.992;
		if(eta>=0.9 && eta<1.2) eff = 0.995;
		if(eta>=1.2 && eta<2.1) eff = 0.998;
		if(eta>=2.1 && eta<2.4) eff = 0.996;
	}
	if(pt>=30. && pt<35.){
		if(eta<0.9) eff = 0.993;
		if(eta>=0.9 && eta<1.2) eff = 0.993;
		if(eta>=1.2 && eta<2.1) eff = 0.997;
		if(eta>=2.1 && eta<2.4) eff = 1.001;
	}
	if(pt>=35. && pt<40.){
		if(eta<0.9) eff = 0.994;
		if(eta>=0.9 && eta<1.2) eff = 0.992;
		if(eta>=1.2 && eta<2.1) eff = 0.996;
		if(eta>=2.1 && eta<2.4) eff = 0.993;
	}
	if(pt>=40. && pt<50.){
		if(eta<0.9) eff = 0.992;
		if(eta>=0.9 && eta<1.2) eff = 0.992;
		if(eta>=1.2 && eta<2.1) eff = 0.996;
		if(eta>=2.1 && eta<2.4) eff = 0.995;
	}
	if(pt>=50. && pt<60.){
		if(eta<0.9) eff = 0.991;
		if(eta>=0.9 && eta<1.2) eff = 0.995;
		if(eta>=1.2 && eta<2.1) eff = 0.995;
		if(eta>=2.1 && eta<2.4) eff = 0.994;
	}
	if(pt>=60. && pt<90.){
		if(eta<0.9) eff = 0.989;
		if(eta>=0.9 && eta<1.2) eff = 0.990;
		if(eta>=1.2 && eta<2.1) eff = 0.992;
		if(eta>=2.1 && eta<2.4) eff = 0.989;
	}
	if(pt>=90. && pt<140.){
		if(eta<0.9) eff = 1.004;
		if(eta>=0.9 && eta<1.2) eff = 1.009;
		if(eta>=1.2 && eta<2.1) eff = 1.023;
		if(eta>=2.1 && eta<2.4) eff = 1.060;
	}
	if(pt>=140. && pt<300.){
		if(eta<0.9) eff = 1.019;
		if(eta>=0.9 && eta<1.2) eff = 1.011;
		if(eta>=1.2 && eta<2.1) eff = 0.975;
		if(eta>=2.1 && eta<2.4) eff = 0.891;
	}
	return eff;
}

double ReferenceScaleFactors::MuonIdUncTight2012(TLorentzVector vect){ // todo: find correct uncertainties
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double err = 0.;
	if(pt>=10. && pt<20.){
		if(eta<0.9) err = 0.005;
		if(eta>=0.9 && eta<1.2) err = 0.007;
		if(eta>=1.2 && eta<2.1) err = 0.004;
		if(eta>=2.1 && eta<2.4) err = 0.007;
	}
	if(pt>=20. && pt<25.){
		if(eta<0.9) err = 0.0016;
		if(eta>=0.9 && eta<1.2) err = 0.0025;
		if(eta>=1.2 && eta<2.1) err = 0.0014;
		if(eta>=2.1 && eta<2.4) err = 0.003;
	}
	if(pt>=25. && pt<30.){
		if(eta<0.9) err = 0.0008;
		if(eta>=0.9 && eta<1.2) err = 0.0014;
		if(eta>=1.2 && eta<2.1) err = 0.0009;
		if(eta>=2.1 && eta<2.4) err = 0.0018;
	}
	if(pt>=30. && pt<35.){
		if(eta<0.9) err = 0.0005;
		if(eta>=0.9 && eta<1.2) err = 0.001;
		if(eta>=1.2 && eta<2.1) err = 0.0007;
		if(eta>=2.1 && eta<2.4) err = 0.0015;
	}
	if(pt>=35. && pt<40.){
		if(eta<0.9) err = 0.0004;
		if(eta>=0.9 && eta<1.2) err = 0.0008;
		if(eta>=1.2 && eta<2.1) err = 0.0006;
		if(eta>=2.1 && eta<2.4) err = 0.0013;
	}
	if(pt>=40. && pt<50.){
		if(eta<0.9) err = 0.0003;
		if(eta>=0.9 && eta<1.2) err = 0.0005;
		if(eta>=1.2 && eta<2.1) err = 0.00024;
		if(eta>=2.1 && eta<2.4) err = 0.001;
	}
	if(pt>=50. && pt<60.){
		if(eta<0.9) err = 0.0007;
		if(eta>=0.9 && eta<1.2) err = 0.0012;
		if(eta>=1.2 && eta<2.1) err = 0.0009;
		if(eta>=2.1 && eta<2.4) err = 0.0026;
	}
	if(pt>=60. && pt<90.){
		if(eta<0.9) err = 0.001;
		if(eta>=0.9 && eta<1.2) err = 0.0019;
		if(eta>=1.2 && eta<2.1) err = 0.0015;
		if(eta>=2.1 && eta<2.4) err = 0.005;
	}
	if(pt>=90. && pt<140.){
		if(eta<0.9) err = 0.003;
		if(eta>=0.9 && eta<1.2) err = 0.0062;
		if(eta>=1.2 && eta<2.1) err = 0.0054;
		if(eta>=2.1 && eta<2.4) err = 0.013;
	}
	if(pt>=140. && pt<300.){
		if(eta<0.9) err = 0.017;
		if(eta>=0.9 && eta<1.2) err = 0.03;
		if(eta>=1.2 && eta<2.1) err = 0.03;
		if(eta>=2.1 && eta<2.4) err = 0.14;
	}
	return err;
}

double ReferenceScaleFactors::MuonIsoTight2012(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=10. && pt<20.){
		if(eta<0.9) eff = 0.940;
		if(eta>=0.9 && eta<1.2) eff = 0.948;
		if(eta>=1.2 && eta<2.1) eff = 0.972;
		if(eta>=2.1 && eta<2.4) eff = 1.117;
	}
	if(pt>=20. && pt<25.){
		if(eta<0.9) eff = 0.977;
		if(eta>=0.9 && eta<1.2) eff = 0.986;
		if(eta>=1.2 && eta<2.1) eff = 0.990;
		if(eta>=2.1 && eta<2.4) eff = 1.116;
	}
	if(pt>=25. && pt<30.){
		if(eta<0.9) eff = 0.996;
		if(eta>=0.9 && eta<1.2) eff = 1.000;
		if(eta>=1.2 && eta<2.1) eff = 1.003;
		if(eta>=2.1 && eta<2.4) eff = 1.097;
	}
	if(pt>=30. && pt<35.){
		if(eta<0.9) eff = 0.993;
		if(eta>=0.9 && eta<1.2) eff = 1.000;
		if(eta>=1.2 && eta<2.1) eff = 1.004;
		if(eta>=2.1 && eta<2.4) eff = 1.075;
	}
	if(pt>=35. && pt<40.){
		if(eta<0.9) eff = 0.994;
		if(eta>=0.9 && eta<1.2) eff = 0.999;
		if(eta>=1.2 && eta<2.1) eff = 1.002;
		if(eta>=2.1 && eta<2.4) eff = 1.061;
	}
	if(pt>=40. && pt<50.){
		if(eta<0.9) eff = 0.994;
		if(eta>=0.9 && eta<1.2) eff = 0.999;
		if(eta>=1.2 && eta<2.1) eff = 1.001;
		if(eta>=2.1 && eta<2.4) eff = 1.034;
	}
	if(pt>=50. && pt<60.){
		if(eta<0.9) eff = 0.996;
		if(eta>=0.9 && eta<1.2) eff = 0.998;
		if(eta>=1.2 && eta<2.1) eff = 1.000;
		if(eta>=2.1 && eta<2.4) eff = 1.025;
	}
	if(pt>=60. && pt<90.){
		if(eta<0.9) eff = 0.999;
		if(eta>=0.9 && eta<1.2) eff = 0.999;
		if(eta>=1.2 && eta<2.1) eff = 1.001;
		if(eta>=2.1 && eta<2.4) eff = 1.015;
	}
	if(pt>=90. && pt<140.){
		if(eta<0.9) eff = 1.000;
		if(eta>=0.9 && eta<1.2) eff = 1.001;
		if(eta>=1.2 && eta<2.1) eff = 0.999;
		if(eta>=2.1 && eta<2.4) eff = 1.008;
	}
	if(pt>=140. && pt<300.){
		if(eta<0.9) eff = 0.999;
		if(eta>=0.9 && eta<1.2) eff = 1.002;
		if(eta>=1.2 && eta<2.1) eff = 0.996;
		if(eta>=2.1 && eta<2.4) eff = 1.011;
	}
	return eff;
}

double ReferenceScaleFactors::MuonIsoUncTight2012(TLorentzVector vect){ // todo: find correct uncertainties
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double err = 0.;
	if(pt>=10. && pt<20.){
		if(eta<0.9) err = 0.0043;
		if(eta>=0.9 && eta<1.2) err = 0.0048;
		if(eta>=1.2 && eta<2.1) err = 0.0023;
		if(eta>=2.1 && eta<2.4) err = 0.0055;
	}
	if(pt>=20. && pt<25.){
		if(eta<0.9) err = 0.0022;
		if(eta>=0.9 && eta<1.2) err = 0.0031;
		if(eta>=1.2 && eta<2.1) err = 0.0015;
		if(eta>=2.1 && eta<2.4) err = 0.0041;
	}
	if(pt>=25. && pt<30.){
		if(eta<0.9) err = 0.0012;
		if(eta>=0.9 && eta<1.2) err = 0.002;
		if(eta>=1.2 && eta<2.1) err = 0.001;
		if(eta>=2.1 && eta<2.4) err = 0.0028;
	}
	if(pt>=30. && pt<35.){
		if(eta<0.9) err = 0.0008;
		if(eta>=0.9 && eta<1.2) err = 0.0014;
		if(eta>=1.2 && eta<2.1) err = 0.0008;
		if(eta>=2.1 && eta<2.4) err = 0.0022;
	}
	if(pt>=35. && pt<40.){
		if(eta<0.9) err = 0.00054;
		if(eta>=0.9 && eta<1.2) err = 0.00092;
		if(eta>=1.2 && eta<2.1) err = 0.00056;
		if(eta>=2.1 && eta<2.4) err = 0.00017;
	}
	if(pt>=40. && pt<50.){
		if(eta<0.9) err = 0.00027;
		if(eta>=0.9 && eta<1.2) err = 0.00044;
		if(eta>=1.2 && eta<2.1) err = 0.00026;
		if(eta>=2.1 && eta<2.4) err = 0.00092;
	}
	if(pt>=50. && pt<60.){
		if(eta<0.9) err = 0.0005;
		if(eta>=0.9 && eta<1.2) err = 0.0008;
		if(eta>=1.2 && eta<2.1) err = 0.00048;
		if(eta>=2.1 && eta<2.4) err = 0.00017;
	}
	if(pt>=60. && pt<90.){
		if(eta<0.9) err = 0.00058;
		if(eta>=0.9 && eta<1.2) err = 0.001;
		if(eta>=1.2 && eta<2.1) err = 0.0006;
		if(eta>=2.1 && eta<2.4) err = 0.0021;
	}
	if(pt>=90. && pt<140.){
		if(eta<0.9) err = 0.0011;
		if(eta>=0.9 && eta<1.2) err = 0.0018;
		if(eta>=1.2 && eta<2.1) err = 0.0012;
		if(eta>=2.1 && eta<2.4) err = 0.0045;
	}
	if(pt>=140. && pt<300.){
		if(eta<0.9) err = 0.0026;
		if(eta>=0.9 && eta<1.2) err = 0.005;
		if(eta>=1.2 && eta<2.1) err = 0.0028;
		if(eta>=2.1 && eta<2.4) err = 0.018;
	}
	return err;
}

double ReferenceScaleFactors::TrackingEfficiency2012(TLorentzVector vect){
	double eta = vect.Eta();
	double eff = 1.;
	if(eta>=-2.4 && eta<-2.2) eff = 0.9861;
	if(eta>=-2.2 && eta<-2.0) eff = 0.9901;
	if(eta>=-2.0 && eta<-1.8) eff = 0.9946;
	if(eta>=-1.8 && eta<-1.6) eff = 0.9966;
	if(eta>=-1.6 && eta<-1.4) eff = 0.9964;
	if(eta>=-1.4 && eta<-1.2) eff = 0.9970;
	if(eta>=-1.2 && eta<-1.0) eff = 0.9973;
	if(eta>=-1.0 && eta<-0.8) eff = 0.9977;
	if(eta>=-0.8 && eta<-0.6) eff = 0.9980;
	if(eta>=-0.6 && eta<-0.4) eff = 0.9981;
	if(eta>=-0.4 && eta<-0.2) eff = 0.9974;
	if(eta>=-0.2 && eta<0.0) eff = 0.9967;
	if(eta>=0.0 && eta<0.2) eff = 0.9960;
	if(eta>=0.2 && eta<0.4) eff = 0.9978;
	if(eta>=0.4 && eta<0.6) eff = 0.9977;
	if(eta>=0.6 && eta<0.8) eff = 0.9975;
	if(eta>=0.8 && eta<1.0) eff = 0.9976;
	if(eta>=1.0 && eta<1.2) eff = 0.9964;
	if(eta>=1.2 && eta<1.4) eff = 0.9970;
	if(eta>=1.4 && eta<1.6) eff = 0.9948;
	if(eta>=1.6 && eta<1.8) eff = 0.9977;
	if(eta>=1.8 && eta<2.0) eff = 0.9974;
	if(eta>=2.0 && eta<2.2) eff = 0.9917;
	if(eta>=2.2 && eta<2.4) eff = 0.9806;
	return eff;
}

///////////////////////////
//
// Electron scale factors
//

// Electrons need eta from supercluster -> not using four vector as argument
double ReferenceScaleFactors::ElectronReconstruction2012(double Et, double Eta){
	double et = Et;
	double eff = 1.;
	if(fabs(Eta)<2.5){
		if(et>=200) et = 199.;
		eff = ElectronRecoEff->GetBinContent(ElectronRecoEff->FindFixBin(et,Eta));
	}
	return eff;
}

double ReferenceScaleFactors::ElectronReconstructionUnc2012(double Et, double Eta){
	double et = Et;
	double err = 0.;
	if(fabs(Eta)<2.5){
		if(et>=200) et = 199.;
		err = ElectronRecoEff->GetBinError(ElectronRecoEff->FindFixBin(et,Eta));
	}
	return err;
}

double ReferenceScaleFactors::ElectronIdTrig2012(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(eta<2.5){
		if(et>=200) et = 199;
		eff = ElectronTrigEff->GetBinContent(ElectronTrigEff->FindFixBin(eta,et));
	}
	return eff;
}

double ReferenceScaleFactors::ElectronIdTrigUnc2012(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double err = 0.;
	if(eta<2.5){
		if(et>=200) et = 199;
		err = ElectronTrigEff->GetBinError(ElectronTrigEff->FindFixBin(eta,et));
	}
	return err;
}

double ReferenceScaleFactors::ElectronIdNonTrig2012(double Et, double Eta){
	double et = Et;
	double eff = 1.;
	if(fabs(Eta)<2.5){
		if(et>=200) et = 199;
		eff = ElectronNonTrigEff->GetBinContent(ElectronNonTrigEff->FindFixBin(et,Eta));
	}
	return eff;
}

double ReferenceScaleFactors::ElectronIdNonTrigUnc2012(double Et, double Eta){
	double et = Et;
	double err = 0.;
	if(fabs(Eta)<2.5){
		if(et>=200) et = 199;
		err = ElectronNonTrigEff->GetBinError(ElectronNonTrigEff->FindFixBin(et,Eta));
	}
	return err;
}

///////////////////////////
//
// Tau scale factors
//
