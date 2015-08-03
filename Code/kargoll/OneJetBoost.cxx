/*
 * OneJetBoost.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "OneJetBoost.h"
#include "VBFTightStandalone.h"
#include "VBFLooseStandalone.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

OneJetBoost::OneJetBoost(TString Name_, TString id_):
	Category(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class OneJetBoost" << std::endl;
	// run OneJetBoost category
	categoryFlag = "OneJetBoost";
}

OneJetBoost::~OneJetBoost() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

void OneJetBoost::categoryConfiguration(){
	// set cut values to be the cut values of this category
	cut.at(NJet)	= 1;
	cut.at(NotVbf)	= true;
	cut.at(TauPt)	= cCat_splitTauPt;
	cut.at(HiggsPt)	= 100.0;

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(NJet)="Number Jets $>=$";
	title.at(NJet)+=cut.at(NJet);
	htitle=title.at(NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=NJet;
	Nminus1.at(NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(NotVbf)="No VBF $==$";
	title.at(NotVbf)+=cut.at(NotVbf);
	htitle=title.at(NotVbf);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Did not pass VBF cat.";
	c="_Cut_";c+=NotVbf;
	Nminus1.at(NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus1_NotVbf_",htitle,2,-0.5,1.5,hlabel,"Events");
	Nminus0.at(NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus0_NotVbf_",htitle,2,-0.5,1.5,hlabel,"Events");

	title.at(TauPt)="$p_{T}(\\tau_{h}) >=$";
	title.at(TauPt)+=cut.at(TauPt);
	title.at(TauPt)+=" GeV";
	htitle=title.at(TauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=TauPt;
	Nminus1.at(TauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_TauPt_",htitle,40,0.,200.,hlabel,"Events");
	Nminus0.at(TauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_TauPt_",htitle,40,0.,200.,hlabel,"Events");

	title.at(HiggsPt)="$p_{T}(H) >=$";
	title.at(HiggsPt)+=cut.at(HiggsPt);
	title.at(HiggsPt)+=" GeV";
	htitle=title.at(HiggsPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T} of Higgs candidate";
	c="_Cut_";c+=HiggsPt;
	Nminus1.at(HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus1_HiggsPt_",htitle,40,0.,200.,hlabel,"Events");
	Nminus0.at(HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus0_HiggsPt_",htitle,40,0.,200.,hlabel,"Events");
}

bool OneJetBoost::categorySelection(){
	bool categoryPass = true;
	std::vector<float> value_OneJetBoost(NCuts,-10);
	std::vector<float> pass_OneJetBoost(NCuts,false);

	value_OneJetBoost.at(NJet) = nJets_;
	pass_OneJetBoost.at(NJet) = ( value_OneJetBoost.at(NJet) >= cut.at(NJet) );

	VBFTightStandalone vbft(nJets_, jetdEta_, nJetsInGap_, mjj_, higgsPt_);
	vbft.run();
	VBFLooseStandalone vbfl(nJets_, jetdEta_, nJetsInGap_, mjj_, !vbft.passed());
	vbfl.run();

	value_OneJetBoost.at(NotVbf) = (not vbft.passed()) && (not vbfl.passed());
	pass_OneJetBoost.at(NotVbf) = ( value_OneJetBoost.at(NotVbf) == cut.at(NotVbf) );

	if (selTau == -1){
		// TauPt cut is set to true for nice N-0 and N-1 plots
		value_OneJetBoost.at(TauPt) = -10.;
		pass_OneJetBoost.at(TauPt) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		value_OneJetBoost.at(TauPt) = tauPt_;
		pass_OneJetBoost.at(TauPt) = ( value_OneJetBoost.at(TauPt) >= cut.at(TauPt) );
	}

	value_OneJetBoost.at(HiggsPt) = higgsPt_;
	pass_OneJetBoost.at(HiggsPt) = (value_OneJetBoost.at(HiggsPt) >= cut.at(HiggsPt));

	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("OneJetBoost",value_OneJetBoost, pass_OneJetBoost,NCuts) && categoryPass;
	return categoryPass;
}


void OneJetBoost::categoryPlotting(){
	// === QCD efficiency method ===
	// 1Jet boost: efficiency from sideband with anti-iso muon and relaxed tau iso
	if (getStatusBoolean(FullInclusiveNoTauNoMuNoCharge, originalPass) && originalPass.at(OppCharge) && !getStatusBoolean(VtxMu, originalPass) && hasRelaxedIsoTau && hasAntiIsoMuon){
		if (isQCDShapeEvent && hasRelaxedIsoTau){
			if (!HConfig.GetHisto(true, DataMCType::Data, t)){
				Logger(Logger::Error) << "failed to find id " << DataMCType::Data << std::endl;
				return;
			}
		}

		h_BGM_QcdEff.at(t).Fill(0., w);
		if (catPassed) h_BGM_QcdEff.at(t).Fill(1., w);

		// take care of events in QCD shape region: set t back to QCD
		if (isQCDShapeEvent && hasRelaxedIsoTau){
			if (!HConfig.GetHisto(false, DataMCType::QCD, t)){
				Logger(Logger::Error) << "failed to find id " << DataMCType::QCD << std::endl;
				return;
			}
		}
	}
}

