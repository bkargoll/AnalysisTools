/*
 * ZeroJetLow.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "ZeroJetLow.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

ZeroJetLow::ZeroJetLow(TString Name_, TString id_):
	Category(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class ZeroJetLow" << std::endl;
	// run ZeroJetLow category
	categoryFlag = "ZeroJetLow";
}

ZeroJetLow::~ZeroJetLow() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

void ZeroJetLow::categoryConfiguration(){
	// set cut values for this category
	cut.at(NJet)		= 0;
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

	title.at(TauPt)="$p_{T}(\\tau_{h}) <$";
	title.at(TauPt)+=cut.at(TauPt);
	title.at(TauPt)+=" GeV";
	htitle=title.at(TauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=TauPt;
	Nminus1.at(TauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_TauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(TauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_TauPt_",htitle,40,0.,200.,hlabel,"Events");
}

bool ZeroJetLow::categorySelection(){
	bool categoryPass = true;
	std::vector<float> value_ZeroJetLow(NCuts,-10);
	std::vector<float> pass_ZeroJetLow(NCuts,false);

	value_ZeroJetLow.at(NJet) = nJets_;
	pass_ZeroJetLow.at(NJet) = ( value_ZeroJetLow.at(NJet) <= cut.at(NJet) );

	if (selTau == -1){
		// TauPt cut is set to true for nice N-0 and N-1 plots
		value_ZeroJetLow.at(TauPt) = -10.;
		pass_ZeroJetLow.at(TauPt) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		value_ZeroJetLow.at(TauPt) = tauPt_;
		pass_ZeroJetLow.at(TauPt) = ( value_ZeroJetLow.at(TauPt) < cut.at(TauPt) );
	}

	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("ZeroJetLow",value_ZeroJetLow, pass_ZeroJetLow,NCuts) && categoryPass;
	return categoryPass;

}
