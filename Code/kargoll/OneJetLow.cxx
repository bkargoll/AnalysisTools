/*
 * OneJetLow.cxx
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#include "OneJetLow.h"
#include "VBFTightStandalone.h"
#include "VBFLooseStandalone.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

OneJetLow::OneJetLow(TString Name_, TString id_):
	Category(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class OneJetLow" << std::endl;
	// run OneJetLow category
	categoryFlag = "OneJetLow";
}

OneJetLow::~OneJetLow() {
	  for(unsigned int j=0; j<Npassed.size(); j++){
	    Logger(Logger::Info) << "Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	  }
	  Logger(Logger::Info) << "Done." << std::endl;
}

void OneJetLow::categoryConfiguration(){
	// set cut values to be the cut values of this category
	cut.at(NJet) = 1;
	cut.at(NotVbf) = true;
	cut.at(TauPt) = cCat_splitTauPt;

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

	title.at(NotVbf)="Not VBFTight or VBFLoose $==$";
	title.at(NotVbf)+=cut.at(NotVbf);
	htitle=title.at(NotVbf);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Did not pass VBF cat.";
	c="_Cut_";c+=NotVbf;
	Nminus1.at(NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus1_NotVBF_",htitle,2,-0.5,1.5,hlabel,"Events");
	Nminus0.at(NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus0_NotVBF_",htitle,2,-0.5,1.5,hlabel,"Events");

	title.at(TauPt)="$p_{T}(\\tau_{h}) <$";
	title.at(TauPt)+=cut.at(TauPt);
	title.at(TauPt)+=" GeV";
	htitle=title.at(TauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=TauPt;
	Nminus1.at(TauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_TauPt_",htitle,40,0.,200.,hlabel,"Events");
	Nminus0.at(TauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_TauPt_",htitle,40,0.,200.,hlabel,"Events");
}

bool OneJetLow::categorySelection(){
	bool categoryPass = true;
	std::vector<float> value_OneJetLow(NCuts,-10);
	std::vector<float> pass_OneJetLow(NCuts,false);

	value_OneJetLow.at(NJet) = nJets_;
	pass_OneJetLow.at(NJet) = ( value_OneJetLow.at(NJet) >= cut.at(NJet) );

	VBFTightStandalone vbft(nJets_, jetdEta_, nJetsInGap_, mjj_, higgsPt_);
	vbft.run();
	VBFLooseStandalone vbfl(nJets_, jetdEta_, nJetsInGap_, mjj_, !vbft.passed());
	vbfl.run();

	value_OneJetLow.at(NotVbf) = (not vbft.passed()) && (not vbfl.passed());
	pass_OneJetLow.at(NotVbf) = ( value_OneJetLow.at(NotVbf) == cut.at(NotVbf) );

	if (selTau == -1){
		// TauPt cut is set to true for nice N-0 and N-1 plots
		value_OneJetLow.at(TauPt) = -10.;
		pass_OneJetLow.at(TauPt) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		value_OneJetLow.at(TauPt) = tauPt_;
		pass_OneJetLow.at(TauPt) = ( value_OneJetLow.at(TauPt) < cut.at(TauPt) );
	}

	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("OneJetLow",value_OneJetLow, pass_OneJetLow,NCuts) && categoryPass;
	return categoryPass;
}
