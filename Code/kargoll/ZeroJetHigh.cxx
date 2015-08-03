/*
 * ZeroJetHigh.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "ZeroJetHigh.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

ZeroJetHigh::ZeroJetHigh(TString Name_, TString id_):
	Category(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class ZeroJetHigh" << std::endl;
	// run ZeroJetHigh category
	categoryFlag = "ZeroJetHigh";
}

ZeroJetHigh::~ZeroJetHigh() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

void ZeroJetHigh::categoryConfiguration(){
	// set cut values to be the cut values of this category
	cut.at(NJet)	= 0;
	cut.at(TauPt)	= cCat_splitTauPt;

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(NJet)="Number Jets $<=$";
	title.at(NJet)+=cut.at(NJet);
	htitle=title.at(NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=NJet;
	Nminus1.at(NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(TauPt)="$p_{T}(\\tau_{h}) >=$";
	title.at(TauPt)+=cut.at(TauPt);
	title.at(TauPt)+=" GeV";
	htitle=title.at(TauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=TauPt;
	Nminus1.at(TauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_TauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(TauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_TauPt_",htitle,50,0.,200.,hlabel,"Events");
}

bool ZeroJetHigh::categorySelection(){
	bool categoryPass = true;
	std::vector<float> value_ZeroJetHigh(NCuts,-10);
	std::vector<float> pass_ZeroJetHigh(NCuts,false);

	value_ZeroJetHigh.at(NJet) = nJets_;
	pass_ZeroJetHigh.at(NJet) = ( value_ZeroJetHigh.at(NJet) <= cut.at(NJet) );

	if (selTau == -1){
		// TauPt cut is set to true for nice N-0 and N-1 plots
		value_ZeroJetHigh.at(TauPt) = -10.;
		pass_ZeroJetHigh.at(TauPt) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		value_ZeroJetHigh.at(TauPt) = tauPt_;
		pass_ZeroJetHigh.at(TauPt) = ( value_ZeroJetHigh.at(TauPt) >= cut.at(TauPt) );
	}

	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("ZeroJetHigh",value_ZeroJetHigh, pass_ZeroJetHigh,NCuts) && categoryPass;
	return categoryPass;
}
