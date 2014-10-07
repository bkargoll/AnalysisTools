/*
 * ReferenceScaleFactors.cxx
 *
 *  Created on: Aug 13, 2014
 *      Author: nehrkorn
 */

#include "ReferenceScaleFactors.h"
#include "Selection_Base.h"

///////////////////////////
//
// Constructor
//

ReferenceScaleFactors::ReferenceScaleFactors(int runType){
	TString basedir = "";
	if(runType==Selection_Base::GRID){
		basedir = (TString)std::getenv("PWD")+"/Code/CommonFiles/";
	}
	else if(runType==Selection_Base::Local){
		basedir = (TString)std::getenv("workdir")+"/Code/CommonFiles/";
	}
	//
	// Read root files
	//

	// Electron Id's
	ETrigIdEffFile = new TFile(basedir+"ElectronEfficiencies_Run2012ReReco_53X_Trig.root");
	ENonTrigIdEffFile = new TFile(basedir+"ElectronEfficiencies_Run2012ReReco_53X_NonTrig.root");
	ERecoEffFile = new TFile(basedir+"Electrons_ScaleFactors_Reco_8TeV.root");
	// Trigger efficiencies
	HWW_TriggerEfficiencies = new TFile(basedir+"TriggerEfficienciesWW_TH1D.root");

	//
	// Get histograms
	//

	// Electron Id's
	ElectronTrigEff = (TH2D*)(ETrigIdEffFile->Get("electronsDATAMCratio_FO_ID_ISO"));
	ElectronNonTrigEff = (TH2D*)(ENonTrigIdEffFile->Get("h_electronScaleFactor_IdIsoSip"));
	ElectronRecoEff = (TH2D*)(ERecoEffFile->Get("h_electronScaleFactor_RECO"));
	// Trigger efficiencies
	HiggsWW_EMu_SingleEle15 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleEle15"));
	HiggsWW_EMu_SingleEle25 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleEle25"));
	HiggsWW_EMu_DoubleEleLead15 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleEleLead15"));
	HiggsWW_EMu_DoubleEleLead25 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleEleLead25"));
	HiggsWW_EMu_DoubleEleTrail15 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleEleTrail15"));
	HiggsWW_EMu_DoubleEleTrail25 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleEleTrail25"));
	HiggsWW_EMu_SingleMu08 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleMu08"));
	HiggsWW_EMu_SingleMu12 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleMu12"));
	HiggsWW_EMu_SingleMu21 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleMu21"));
	HiggsWW_EMu_SingleMu25 = (TH1D*)(HWW_TriggerEfficiencies->Get("SingleMu25"));
	HiggsWW_EMu_DoubleMuLead12 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuLead12"));
	HiggsWW_EMu_DoubleMuLead21 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuLead21"));
	HiggsWW_EMu_DoubleMuLead25 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuLead25"));
	HiggsWW_EMu_DoubleMuTrail12 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuTrail12"));
	HiggsWW_EMu_DoubleMuTrail21 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuTrail21"));
	HiggsWW_EMu_DoubleMuTrail25 = (TH1D*)(HWW_TriggerEfficiencies->Get("DoubleMuTrail25"));

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

double ReferenceScaleFactors::MuonIdUncTight2012(TLorentzVector vect){
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

double ReferenceScaleFactors::MuonIsoUncTight2012(TLorentzVector vect){
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

double ReferenceScaleFactors::HiggsEMuId_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			eff = 0.9771;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9746;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9644;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9891;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			eff = 0.9548;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9701;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9766;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9892;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			eff = 0.9648;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9836;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9820;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9909;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			eff = 0.9676;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9817;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9886;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9883;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			eff = 0.9883;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9833;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9910;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9900;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			eff = 0.9826;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9841;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9900;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9886;
		}
	}
	return eff;
}

double ReferenceScaleFactors::HiggsEMuIdUnc_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double err = 0.;
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			err = 0.0107;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0091;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0078;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0080;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			err = 0.0046;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0049;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0049;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0049;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			err = 0.0023;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0030;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0029;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0030;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			err = 0.0012;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0018;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0018;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0019;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			err = 0.0008;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0012;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0013;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0016;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			err = 0.0005;
		}else if(eta>=0.8 && eta<1.2){
			err = 0.0004;
		}else if(eta>=1.2 && eta<1.6){
			err = 0.0003;
		}else if(eta>=1.6 && eta<2.1){
			err = 0.0004;
		}
	}
	return err;
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

double ReferenceScaleFactors::ElectronEmbedding2012(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	double xAxis[10] = {10,15,20,25,30,40,55,70,100,200};
	double yAxis[4] = {0,0.8,1.479,2.5};
	TH2D hPtEtaSFL("hPtEtaSFL","",9,xAxis,3,yAxis);
	hPtEtaSFL.SetBinContent(12,0.81);
	hPtEtaSFL.SetBinContent(13,0.91);
	hPtEtaSFL.SetBinContent(14,0.95);
	hPtEtaSFL.SetBinContent(15,0.96);
	hPtEtaSFL.SetBinContent(16,0.97);
	hPtEtaSFL.SetBinContent(17,0.98);
	hPtEtaSFL.SetBinContent(18,0.99);
	hPtEtaSFL.SetBinContent(19,0.98);
	hPtEtaSFL.SetBinContent(20,0.99);
	hPtEtaSFL.SetBinContent(21,0.98);
	hPtEtaSFL.SetBinContent(23,0.78);
	hPtEtaSFL.SetBinContent(24,0.89);
	hPtEtaSFL.SetBinContent(25,0.92);
	hPtEtaSFL.SetBinContent(26,0.94);
	hPtEtaSFL.SetBinContent(27,0.94);
	hPtEtaSFL.SetBinContent(28,0.97);
	hPtEtaSFL.SetBinContent(29,0.97);
	hPtEtaSFL.SetBinContent(30,0.99);
	hPtEtaSFL.SetBinContent(31,1.00);
	hPtEtaSFL.SetBinContent(32,1.00);
	hPtEtaSFL.SetBinContent(34,0.46);
	hPtEtaSFL.SetBinContent(35,0.66);
	hPtEtaSFL.SetBinContent(36,0.73);
	hPtEtaSFL.SetBinContent(37,0.80);
	hPtEtaSFL.SetBinContent(38,0.83);
	hPtEtaSFL.SetBinContent(39,0.86);
	hPtEtaSFL.SetBinContent(40,0.88);
	hPtEtaSFL.SetBinContent(41,0.91);
	hPtEtaSFL.SetBinContent(42,0.93);
	hPtEtaSFL.SetBinContent(43,1.00);

	if(et>199.99) et = 199.9;
	if(eta>2.49) eta = 2.49;
	if(et<10.) eff = 0.;
	eff = hPtEtaSFL.GetBinContent(hPtEtaSFL.FindFixBin(et,eta));

	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_Id_E(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(et>10 && et<=15){
		if(eta>=0 && eta<0.8){
			eff = 0.7654;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.7693;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.5719;
		}
	}else if(et>15 && et<=20){
		if(eta>=0 && eta<0.8){
			eff = 0.8394;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.8457;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.7024;
		}
	}else if(et>20 && et<=25){
		if(eta>=0 && eta<0.8){
			eff = 0.8772;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.8530;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.7631;
		}
	}else if(et>25 && et<=30){
		if(eta>=0 && eta<0.8){
			eff = 0.9006;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.8874;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.8092;
		}
	}else if(et>30 && et<=35){
		if(eta>=0 && eta<0.8){
			eff = 0.9261;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.9199;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.8469;
		}
	}else if(et>35){
		if(eta>=0 && eta<0.8){
			eff = 0.9514;
		}else if(eta>=0.8 && eta<1.5){
			eff = 0.9445;
		}else if(eta>=1.5 && eta<2.3){
			eff = 0.9078;
		}
	}
	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_IdUnc_E(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double err = 0.;
	if(et>10 && et<=15){
		if(eta>=0 && eta<0.8){
			err = 0.0149;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0164;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0131;
		}
	}else if(et>15 && et<=20){
		if(eta>=0 && eta<0.8){
			err = 0.0045;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0061;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0075;
		}
	}else if(et>20 && et<=25){
		if(eta>=0 && eta<0.8){
			err = 0.0023;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0039;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0061;
		}
	}else if(et>25 && et<=30){
		if(eta>=0 && eta<0.8){
			err = 0.0018;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0017;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0024;
		}
	}else if(et>30 && et<=35){
		if(eta>=0 && eta<0.8){
			err = 0.0007;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0010;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0024;
		}
	}else if(et>35){
		if(eta>=0 && eta<0.8){
			err = 0.0002;
		}else if(eta>=0.8 && eta<1.5){
			err = 0.0003;
		}else if(eta>=1.5 && eta<2.3){
			err = 0.0007;
		}
	}
	return err;
}

///////////////////////////
//
// Tau scale factors
//

///////////////////////////
//
// Trigger scale factors
//

// Final state e+mu

double ReferenceScaleFactors::HiggsWW_EMu_Trigger(TLorentzVector mu_vect, double e_et, double e_eta, TString path){
	double eff = 1.;
	if(path.Contains("Mu17_Ele8")){
		eff = HiggsWW_EMu_SingleMu(mu_vect) + (1-HiggsWW_EMu_SingleMu(mu_vect))*HiggsWW_EMu_SingleEle(e_et, e_eta)
				+ (1-HiggsWW_EMu_SingleMu(mu_vect))*(1-HiggsWW_EMu_SingleEle(e_et, e_eta))*
				(HiggsWW_EMu_DoubleMuLeading(mu_vect)*HiggsWW_EMu_DoubleEleTrailing(e_et, e_eta) + (1-HiggsWW_EMu_DoubleMuLeading(mu_vect)*HiggsWW_EMu_DoubleEleTrailing(e_et, e_eta))*HiggsWW_EMu_DoubleEleLeading(e_et, e_eta)*HiggsWW_EMu_DoubleMuTrailing(mu_vect));
	}
	if(path.Contains("Mu8_Ele17")){
		eff = HiggsWW_EMu_SingleEle(e_et, e_eta) + (1-HiggsWW_EMu_SingleEle(e_et, e_eta))*HiggsWW_EMu_SingleMu(mu_vect)
				+ (1-HiggsWW_EMu_SingleEle(e_et, e_eta))*(1-HiggsWW_EMu_SingleMu(mu_vect))*
				(HiggsWW_EMu_DoubleEleLeading(e_et, e_eta)*HiggsWW_EMu_DoubleMuTrailing(mu_vect) + (1-HiggsWW_EMu_DoubleEleLeading(e_et, e_eta)*HiggsWW_EMu_DoubleMuTrailing(mu_vect))*HiggsWW_EMu_DoubleMuLeading(mu_vect)*HiggsWW_EMu_DoubleEleTrailing(e_et, e_eta));
	}
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_SingleEle(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(et>=60.)et=55.;
	if(eta<1.5) eff = HiggsWW_EMu_SingleEle15->GetBinContent(HiggsWW_EMu_SingleEle15->FindFixBin(et));
	if(eta>=1.5 && eta<2.5) eff = HiggsWW_EMu_SingleEle25->GetBinContent(HiggsWW_EMu_SingleEle25->FindFixBin(et));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_DoubleEleLeading(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(et>=60.)et=55.;
	if(eta<1.5) eff = HiggsWW_EMu_DoubleEleLead15->GetBinContent(HiggsWW_EMu_DoubleEleLead15->FindFixBin(et));
	if(eta>=1.5 && eta<2.5) eff = HiggsWW_EMu_DoubleEleLead25->GetBinContent(HiggsWW_EMu_DoubleEleLead25->FindFixBin(et));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_DoubleEleTrailing(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	double eff = 1.;
	if(et>=60.)et=55.;
	if(eta<1.5) eff = HiggsWW_EMu_DoubleEleTrail15->GetBinContent(HiggsWW_EMu_DoubleEleTrail15->FindFixBin(et));
	if(eta>=1.5 && eta<2.5) eff = HiggsWW_EMu_DoubleEleTrail25->GetBinContent(HiggsWW_EMu_DoubleEleTrail25->FindFixBin(et));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_SingleMu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=60.)pt=55.;
	if(eta<0.8) eff = HiggsWW_EMu_SingleMu08->GetBinContent(HiggsWW_EMu_SingleMu08->FindFixBin(pt));
	if(eta>=0.8 && eta<1.2) eff = HiggsWW_EMu_SingleMu12->GetBinContent(HiggsWW_EMu_SingleMu12->FindFixBin(pt));
	if(eta>=1.2 && eta<2.1) eff = HiggsWW_EMu_SingleMu21->GetBinContent(HiggsWW_EMu_SingleMu21->FindFixBin(pt));
	if(eta>=2.1 && eta<2.5) eff = HiggsWW_EMu_SingleMu25->GetBinContent(HiggsWW_EMu_SingleMu25->FindFixBin(pt));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_DoubleMuLeading(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=60.)pt=55.;
	if(eta<1.2) eff = HiggsWW_EMu_DoubleMuLead12->GetBinContent(HiggsWW_EMu_DoubleMuLead12->FindFixBin(pt));
	if(eta>=1.2 && eta<2.1) eff = HiggsWW_EMu_DoubleMuLead21->GetBinContent(HiggsWW_EMu_DoubleMuLead21->FindFixBin(pt));
	if(eta>=2.1 && eta<2.5) eff = HiggsWW_EMu_DoubleMuLead25->GetBinContent(HiggsWW_EMu_DoubleMuLead25->FindFixBin(pt));
	return eff;
}

double ReferenceScaleFactors::HiggsWW_EMu_DoubleMuTrailing(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	double eff = 1.;
	if(pt>=60.)pt=55.;
	if(eta<1.2) eff = HiggsWW_EMu_DoubleMuTrail12->GetBinContent(HiggsWW_EMu_DoubleMuTrail12->FindFixBin(pt));
	if(eta>=1.2 && eta<2.1) eff = HiggsWW_EMu_DoubleMuTrail21->GetBinContent(HiggsWW_EMu_DoubleMuTrail21->FindFixBin(pt));
	if(eta>=2.1 && eta<2.5) eff = HiggsWW_EMu_DoubleMuTrail25->GetBinContent(HiggsWW_EMu_DoubleMuTrail25->FindFixBin(pt));
	return eff;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_Trigger_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.9829;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9745;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9943;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9158;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.9850;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9852;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9743;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9333;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.9951;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9610;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9716;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9459;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.9869;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9779;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9665;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9501;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.9959;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9881;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9932;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9391;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.9986;
		}else if(eta>=0.8 && eta<1.2){
			return 0.9540;
		}else if(eta>=1.2 && eta<1.6){
			return 0.9549;
		}else if(eta>=1.6 && eta<2.1){
			return 0.9386;
		}
	}
	return 1.;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_TriggerUnc_Mu(TLorentzVector vect){
	double pt = vect.Pt();
	double eta = fabs(vect.Eta());
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			return 0.0058;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0124;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0164;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0176;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			return 0.0056;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0171;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0179;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0162;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			return 0.0060;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0116;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0141;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0159;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			return 0.0074;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0187;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0184;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0251;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			return 0.0085;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0227;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0271;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0307;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			return 0.0087;
		}else if(eta>=0.8 && eta<1.2){
			return 0.0165;
		}else if(eta>=1.2 && eta<1.6){
			return 0.0211;
		}else if(eta>=1.6 && eta<2.1){
			return 0.0209;
		}
	}
	return 0.;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_Trigger_E(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	if(et>10 && et<=15){
		if(eta>=0 && eta<0.8){
			return 0.9548;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9015;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9017;
		}
	}else if(et>15 && et<=20){
		if(eta>=0 && eta<0.8){
			return 0.9830;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9672;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9463;
		}
	}else if(et>20 && et<=25){
		if(eta>=0 && eta<0.8){
			return 0.9707;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9731;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9691;
		}
	}else if(et>25 && et<=30){
		if(eta>=0 && eta<0.8){
			return 0.9768;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9870;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9727;
		}
	}else if(et>30 && et<=35){
		if(eta>=0 && eta<0.8){
			return 1.0047;
		}else if(eta>=0.8 && eta<1.5){
			return 0.9891;
		}else if(eta>=1.5 && eta<2.3){
			return 0.9858;
		}
	}else if(et>35){
		if(eta>=0 && eta<0.8){
			return 1.0063;
		}else if(eta>=0.8 && eta<1.5){
			return 1.0047;
		}else if(eta>=1.5 && eta<2.3){
			return 1.0015;
		}
	}
	return 1.;
}

double ReferenceScaleFactors::HiggsTauTau_EMu_TriggerUnc_E(double Et, double Eta){
	double et = Et;
	double eta = fabs(Eta);
	if(et>10 && et<=15){
		if(eta>=0 && eta<0.8){
			return 0.0197;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0205;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0470;
		}
	}else if(et>15 && et<=20){
		if(eta>=0 && eta<0.8){
			return 0.0115;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0113;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0212;
		}
	}else if(et>20 && et<=25){
		if(eta>=0 && eta<0.8){
			return 0.0087;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0083;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0149;
		}
	}else if(et>25 && et<=30){
		if(eta>=0 && eta<0.8){
			return 0.0084;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0083;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0162;
		}
	}else if(et>30 && et<=35){
		if(eta>=0 && eta<0.8){
			return 0.0100;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0111;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0112;
		}
	}else if(et>35){
		if(eta>=0 && eta<0.8){
			return 0.0078;
		}else if(eta>=0.8 && eta<1.5){
			return 0.0073;
		}else if(eta>=1.5 && eta<2.3){
			return 0.0135;
		}
	}
	return 0.;
}
