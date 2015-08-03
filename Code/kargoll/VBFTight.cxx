/*
 * VBFTight.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "VBFTight.h"
#include "VBFTightStandalone.h"
#include "RelaxedVBFTightStandalone.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

VBFTight::VBFTight(TString Name_, TString id_):
	Category(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class VBFTight" << std::endl;
	// run VBFTight category
	categoryFlag = "VBFTight";
}

VBFTight::~VBFTight() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

void VBFTight::categoryConfiguration(){
	// set cut values to be the cut values of this category
	VBFTightStandalone vbft;
	std::vector<float> vbftCuts = vbft.get_cutValues();
	cut.at(NJet) 		= vbftCuts.at(VBFTightStandalone::NJet);
	cut.at(DeltaEta)	= vbftCuts.at(VBFTightStandalone::DeltaEta);
	cut.at(NJetRapGap)	= vbftCuts.at(VBFTightStandalone::NJetRapGap);
	cut.at(JetInvM)		= vbftCuts.at(VBFTightStandalone::JetInvM);
	cut.at(HiggsPt)		= vbftCuts.at(VBFTightStandalone::HiggsPt);

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

	title.at(HiggsPt)="$p_{T}(H) >$";
	title.at(HiggsPt)+=cut.at(HiggsPt);
	title.at(HiggsPt)+=" GeV";
	htitle=title.at(HiggsPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T} of Higgs candidate";
	c="_Cut_";c+=HiggsPt;
	Nminus1.at(HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus1_HiggsPt_",htitle,25,0.,250.,hlabel,"Events");
	Nminus0.at(HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus0_HiggsPt_",htitle,25,0.,250.,hlabel,"Events");
}

bool VBFTight::categorySelection(){
	// changes to default selection in this category
	// loose relaxed Jet selection for QCD shape
	if(isQCDShapeEvent) calculateJetVariables(selectedLooseJets);

	// categorisation
	std::vector<float> value_VBFTight(NCuts,-10);
	std::vector<float> pass_VBFTight(NCuts,false);

	std::vector<float> vbft_values;
	std::vector<bool>  vbft_pass;

	VBFTightStandalone vbft(nJets_, jetdEta_, nJetsInGap_, mjj_, higgsPt_);
	vbft.run();
	RelaxedVBFTightStandalone rvbft(nJets_, jetdEta_, nJetsInGap_, mjj_, higgsPt_);
	rvbft.run();
	passed_VBFTightRelaxed = rvbft.passed();

	if (isQCDShapeEvent){
		// use relaxed category for QCD shape in VBFTight
		vbft_values	= rvbft.get_eventValues();
		vbft_pass 	= rvbft.get_passCut();
	}
	else {
		vbft_values	= vbft.get_eventValues();
		vbft_pass 	= vbft.get_passCut();
	}

	for (unsigned i_cut = 0; i_cut < vbft_values.size(); i_cut++){
		value_VBFTight.at(HToTaumuTauh::CatCut1 + i_cut) = vbft_values.at(i_cut);
		pass_VBFTight.at(HToTaumuTauh::CatCut1 + i_cut) = vbft_pass.at(i_cut);
	}

	// migrate into main analysis if this is chosen category
	return migrateCategoryIntoMain("VBFTight",value_VBFTight, pass_VBFTight,NCuts);
}

void VBFTight::categoryPlotting(){
	// === QCD efficiency method ===
	// VBF tight: efficiency from sideband with same sign and anti-iso muon and relaxed tau iso
	if (getStatusBoolean(FullInclusiveNoTauNoMuNoCharge, originalPass) && !originalPass.at(OppCharge) && !getStatusBoolean(VtxMu, originalPass) && hasRelaxedIsoTau && hasAntiIsoMuon) {
		if (isQCDShapeEvent && hasRelaxedIsoTau){
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

		h_BGM_QcdEff.at(t).Fill(0., w);
		if (vbft.passed()) h_BGM_QcdEff.at(t).Fill(1., w);

		// switch jet collection back
		if(isQCDShapeEvent) calculateJetVariables(selectedLooseJets);

		// take care of events in QCD shape region: set t back to QCD
		if (isQCDShapeEvent && hasRelaxedIsoTau){
			if (!HConfig.GetHisto(false, DataMCType::QCD, t)){
				Logger(Logger::Error) << "failed to find id " << DataMCType::QCD << std::endl;
				return;
			}
		}
	}
}
