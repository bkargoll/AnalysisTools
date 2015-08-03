/*
 * VBFLoose.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "VBFLoose.h"
#include "VBFLooseStandalone.h"
#include "VBFTightStandalone.h"
#include "RelaxedVBFLooseStandalone.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

VBFLoose::VBFLoose(TString Name_, TString id_):
	Category(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class VBFLoose" << std::endl;
	// run VBFLoose category
	categoryFlag = "VBFLoose";
}

VBFLoose::~VBFLoose() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

void VBFLoose::categoryConfiguration(){
	// set cut values to be the cut values of this category
	VBFLooseStandalone vbfl;
	std::vector<float> vbflCuts = vbfl.get_cutValues();
	cut.at(NJet) 		= vbflCuts.at(VBFLooseStandalone::NJet);
	cut.at(DeltaEta)	= vbflCuts.at(VBFLooseStandalone::DeltaEta);
	cut.at(NJetRapGap)	= vbflCuts.at(VBFLooseStandalone::NJetRapGap);
	cut.at(JetInvM)		= vbflCuts.at(VBFLooseStandalone::JetInvM);
	cut.at(NotVbfTight)	= vbflCuts.at(VBFLooseStandalone::NotVbfTight);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(NJet)="Number VBF Jets $>=$";
	title.at(NJet)+=cut.at(NJet);
	htitle=title.at(NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jet_{VBF}";
	c="_Cut_";c+=NJet;
	Nminus1.at(NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_NJet_",htitle,6,-0.5,5.5,hlabel,"Events");
	Nminus0.at(NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_NJet_",htitle,6,-0.5,5.5,hlabel,"Events");

	title.at(DeltaEta)="$\\Delta\\eta(jj) >$";
	title.at(DeltaEta)+=cut.at(DeltaEta);
	htitle=title.at(DeltaEta);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="#Delta#eta(Jet_{VBF}^{1},Jet_{VBF}^{2})";
	c="_Cut_";c+=DeltaEta;
	Nminus1.at(DeltaEta) = HConfig.GetTH1D(Name+c+"_Nminus1_DeltaEta_",htitle,32,-8.,8.,hlabel,"Events");
	Nminus0.at(DeltaEta) = HConfig.GetTH1D(Name+c+"_Nminus0_DeltaEta_",htitle,32,-8.,8.,hlabel,"Events");

	title.at(NJetRapGap)="Number Jets in $\\eta$ gap $<=$";
	title.at(NJetRapGap)+=cut.at(NJetRapGap);
	htitle=title.at(NJetRapGap);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jet_{VBF} in rapidity gap";
	c="_Cut_";c+=NJetRapGap;
	Nminus1.at(NJetRapGap) = HConfig.GetTH1D(Name+c+"_Nminus1_NJetRapGap_",htitle,6,-0.5,5.5,hlabel,"Events");
	Nminus0.at(NJetRapGap) = HConfig.GetTH1D(Name+c+"_Nminus0_NJetRapGap_",htitle,6,-0.5,5.5,hlabel,"Events");

	title.at(JetInvM)="$m_{jj}(VBF) >$";
	title.at(JetInvM)+=cut.at(JetInvM);
	title.at(JetInvM)+=" GeV";
	htitle=title.at(JetInvM);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="m_{inv}(jj) of VBF-jets";
	c="_Cut_";c+=JetInvM;
	Nminus1.at(JetInvM) = HConfig.GetTH1D(Name+c+"_Nminus1_JetInvM_",htitle,20,0.,2000.,hlabel,"Events");
	Nminus0.at(JetInvM) = HConfig.GetTH1D(Name+c+"_Nminus0_JetInvM_",htitle,20,0.,2000.,hlabel,"Events");

	title.at(NotVbfTight)="Not VBFTight $==$";
	title.at(NotVbfTight)+=cut.at(NotVbfTight);
	htitle=title.at(NotVbfTight);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Did not pass VBFTight cat.";
	c="_Cut_";c+=NotVbfTight;
	Nminus1.at(NotVbfTight) = HConfig.GetTH1D(Name+c+"_Nminus1_NotVbfTight_",htitle,2,-0.5,1.5,hlabel,"Events");
	Nminus0.at(NotVbfTight) = HConfig.GetTH1D(Name+c+"_Nminus0_NotVbfTight_",htitle,2,-0.5,1.5,hlabel,"Events");
}

bool VBFLoose::categorySelection(){
	// changes to default selection in this category
	// loose relaxed Jet selection for QCD shape
	if(isQCDShapeEvent) calculateJetVariables(selectedLooseJets);

	// categorisation
	std::vector<float> value_VBFLoose(NCuts,-10);
	std::vector<float> pass_VBFLoose(NCuts,false);

	std::vector<float> vbfl_values;
	std::vector<bool>  vbfl_pass;

	// run VBFTight category to veto against that passed it
	VBFTightStandalone vbft(nJets_, jetdEta_, nJetsInGap_, mjj_, higgsPt_);
	vbft.run();

	VBFLooseStandalone vbfl(nJets_, jetdEta_, nJetsInGap_, mjj_, !vbft.passed());
	vbfl.run();
	RelaxedVBFLooseStandalone rvbfl(nJets_, jetdEta_, nJetsInGap_, mjj_);
	rvbfl.run();
	passed_VBFLooseRelaxed = rvbfl.passed();

	if (isWJetShapeEvent || isQCDShapeEvent){
		// use relaxed category for QCD and WJets shape in VBFLoose
		vbfl_values	= rvbfl.get_eventValues();
		vbfl_pass 	= rvbfl.get_passCut();
	}
	else{
		vbfl_values	= vbfl.get_eventValues();
		vbfl_pass 	= vbfl.get_passCut();
	}

	for (unsigned i_cut = 0; i_cut < vbfl_values.size(); i_cut++){
		value_VBFLoose.at(HToTaumuTauh::CatCut1 + i_cut) = vbfl_values.at(i_cut);
		pass_VBFLoose.at(HToTaumuTauh::CatCut1 + i_cut) = vbfl_pass.at(i_cut);
	}

	// migrate into main analysis if this is chosen category
	return migrateCategoryIntoMain("VBFLoose",value_VBFLoose, pass_VBFLoose,NCuts);
}

void VBFLoose::categoryPlotting(){
	// === QCD efficiency method ===
	// VBF loose: efficiency from sideband with same sign and anti-iso muon
	if (getStatusBoolean(FullInclusiveNoTauNoMuNoCharge, originalPass) && originalPass.at(NTauKin) && !originalPass.at(OppCharge) && !getStatusBoolean(VtxMu, originalPass) && hasAntiIsoMuon) {
		// take care of events in QCD shape region: set t back to Data temporarily
		if (isQCDShapeEvent){
			if (!HConfig.GetHisto(true, DataMCType::Data, t)){
				Logger(Logger::Error) << "failed to find id " << DataMCType::Data << std::endl;
				return;
			}
		}

		// QCD efficiency has to be calculated using default cut values
		// -> switch jet collection temporarily to default
		calculateJetVariables(selectedJets);

		VBFTightStandalone vbft(nJets_, jetdEta_, nJetsInGap_, mjj_, higgsPt_);
		vbft.run();
		VBFLooseStandalone vbfl(nJets_, jetdEta_, nJetsInGap_, mjj_, !vbft.passed());
		vbfl.run();

		h_BGM_QcdEff.at(t).Fill(0., w);
		if (vbfl.passed()) h_BGM_QcdEff.at(t).Fill(1., w);

		// switch jet collection back
		if(isQCDShapeEvent) calculateJetVariables(selectedLooseJets);

		// take care of events in QCD shape region: set t back to QCD
		if (isQCDShapeEvent){
			if (!HConfig.GetHisto(false, DataMCType::QCD, t)){
				Logger(Logger::Error) << "failed to find id " << DataMCType::QCD << std::endl;
				return;
			}
		}
	}
}
