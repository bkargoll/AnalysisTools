#include "HToTaumuTauh.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

HToTaumuTauh::HToTaumuTauh(TString Name_, TString id_):
  Selection(Name_,id_),
  verbose(0),
  cMu_dxy(0.045),
  cMu_dz(0.2),
  cMu_relIso(0.1),
  cMu_pt(20.0),
  cMu_eta(2.1),
  cMu_dRHltMatch(0.5),
  cTau_pt(20.0),
  cTau_eta(2.3),
  cTau_rawIso(1.5),
  cMuTau_dR(0.5),
  cTau_dRHltMatch(0.5),
  cMuTriLep_pt(10.0),
  cMuTriLep_eta(2.4),
  cEleTriLep_pt(10.0),
  cEleTriLep_eta(2.5),
  cCat_jetPt(30.0),
  cCat_jetEta(4.7),
  cCat_bjetPt(20.0),
  cCat_bjetEta(2.4),
  cCat_btagDisc(0.679), // medium WP, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP#B_tagging_Operating_Points_for_5
  cCat_splitTauPt(45.0),
  cJetClean_dR(0.5)
{
	if (verbose) std::cout << "HToTaumuTauh::HToTaumuTauh(TString Name_, TString id_)" << std::endl;
	TString trigNames[] = {"HLT_IsoMu18_eta2p1_LooseIsoPFTau20","HLT_IsoMu17_eta2p1_LooseIsoPFTau20"};
	std::vector<TString> temp (trigNames, trigNames + sizeof(trigNames) / sizeof(TString) );
	cTriggerNames = temp;

	// implemented categories:
	// VBFTight, VBFLoose
	// OneJetHigh, OneJetLow, OneJetBoost
	// ZeroJetHigh, ZeroJetLow
	// Inclusive
	// NoCategory (for skimming)

	// Set it to "NoCategory" here.
	// For each category, there should be a special class inheriting from HToTaumuTauh
	categoryFlag = "NoCategory";

	// select which WJets Background source to use
	// chose between:
	// * "MC": use MC as given in Histo.txt
	// * "Data": use data driven method (make sure wJetsYieldScaleMap is filled correctly)
	wJetsBGSource = "MC";

	// this one is used to set the event yield for W+Jet
	wJetsYieldMap.insert(std::pair<TString,double>("ZeroJetLow",   6606.6) );
	wJetsYieldMap.insert(std::pair<TString,double>("ZeroJetHigh",  1144.5) );
	wJetsYieldMap.insert(std::pair<TString,double>("OneJetLow",    4842.3) );
	wJetsYieldMap.insert(std::pair<TString,double>("OneJetHigh",    668.2) );
	wJetsYieldMap.insert(std::pair<TString,double>("OneJetBoost",   153.9) );
	wJetsYieldMap.insert(std::pair<TString,double>("VBFLoose",       62.8) );
	wJetsYieldMap.insert(std::pair<TString,double>("VBFTight",        4.6) );
	wJetsYieldMap.insert(std::pair<TString,double>("Inclusive",   13307.8) );

	// this one is used for cross-check only
	// wJetsYieldScaleMap.insert(std::pair<TString,double>("ZeroJetLow",   0.80901627) );
	// wJetsYieldScaleMap.insert(std::pair<TString,double>("ZeroJetHigh",  0.68228370) );
	// wJetsYieldScaleMap.insert(std::pair<TString,double>("OneJetLow",    0.98325814) );
	// wJetsYieldScaleMap.insert(std::pair<TString,double>("OneJetHigh",   0.79677233) );
	// wJetsYieldScaleMap.insert(std::pair<TString,double>("OneJetBoost",  0.67725753) );
	// wJetsYieldScaleMap.insert(std::pair<TString,double>("VBFLoose",     0.14590955) );
	// wJetsYieldScaleMap.insert(std::pair<TString,double>("VBFTight",     0.04354578) );
	// wJetsYieldScaleMap.insert(std::pair<TString,double>("Inclusive",    0.83902189) );

}

HToTaumuTauh::~HToTaumuTauh(){
	if (verbose) std::cout << "HToTaumuTauh::~HToTaumuTauh()" << std::endl;
	for(int j=0; j<Npassed.size(); j++){
	std::cout << "HToTaumuTauh::~HToTaumuTauh Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
	}
	std::cout << "HToTaumuTauh::~HToTaumuTauh() done" << std::endl;
}

void  HToTaumuTauh::Configure(){
  if (verbose) std::cout << "HToTaumuTauh::Configure()" << std::endl;
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)    	cut.at(TriggerOk)=0;
    if(i==PrimeVtx)     	cut.at(PrimeVtx)=1;
    if(i==NMuId)			cut.at(NMuId)=1;
    if(i==NMuKin)			cut.at(NMuKin)=1;
    if(i==DiMuonVeto)		cut.at(DiMuonVeto)=0.15;
    if(i==NTauId)			cut.at(NTauId)=1;
    if(i==NTauIso)			cut.at(NTauIso)=1;
    if(i==NTauKin)			cut.at(NTauKin)=1;
    if(i==TriLeptonVeto)	cut.at(TriLeptonVeto)=0;
    if(i==OppCharge)		cut.at(OppCharge)=0;
    if(i==MT)				cut.at(MT)=30.0; // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#Event_Categories_SM
    if(i==BJetVeto)			cut.at(BJetVeto)=0;
    //category-specific values are set in the corresponding configure function
    // set them to dummy value -10.0 here
    cut_VBFTight.push_back(-10.);
    cut_VBFLoose.push_back(-10.);
    cut_OneJetHigh.push_back(-10.);
    cut_OneJetLow.push_back(-10.);
    cut_OneJetBoost.push_back(-10.);
    cut_ZeroJetHigh.push_back(-10.);
    cut_ZeroJetLow.push_back(-10.);
    cut_VBFTightRelaxed.push_back(-10.);
    cut_VBFLooseRelaxed.push_back(-10.);
    if(i>=CatCut1){
    	cut.at(i)=-10.0;
    }
  }

  // Setup Category Cut Values
  for(unsigned i = CatCut1; i< NCuts; i++){
	  if(i==VbfTight_NJet)		cut_VBFTight.at(VbfTight_NJet)		= 2;
	  if(i==VbfTight_DeltaEta)	cut_VBFTight.at(VbfTight_DeltaEta)	= 4.0;
	  if(i==VbfTight_NJetRapGap)cut_VBFTight.at(VbfTight_NJetRapGap)= 0;
	  if(i==VbfTight_JetInvM)	cut_VBFTight.at(VbfTight_JetInvM)	= 700.0;
	  if(i==VbfTight_HiggsPt)	cut_VBFTight.at(VbfTight_HiggsPt)	= 100.0;

	  if(i==VbfLoose_NJet)		cut_VBFLoose.at(VbfLoose_NJet)		= 2;
	  if(i==VbfLoose_DeltaEta)	cut_VBFLoose.at(VbfLoose_DeltaEta)	= 3.5;
	  if(i==VbfLoose_NJetRapGap)cut_VBFLoose.at(VbfLoose_NJetRapGap)= 0;
	  if(i==VbfLoose_JetInvM)	cut_VBFLoose.at(VbfLoose_JetInvM)	= 500.0;
	  if(i==VbfLoose_NotVbfTight)cut_VBFLoose.at(VbfLoose_NotVbfTight)	= true;

	  if(i==OneJetLow_NJet)     cut_OneJetLow.at(OneJetLow_NJet)    = 1;
	  if(i==OneJetLow_NotVbf)   cut_OneJetLow.at(OneJetLow_NotVbf)  = true;
	  if(i==OneJetLow_TauPt)    cut_OneJetLow.at(OneJetLow_TauPt)   = cCat_splitTauPt;

	  if(i==OneJetHigh_NJet)    cut_OneJetHigh.at(OneJetHigh_NJet)    = 1;
	  if(i==OneJetHigh_NotVbf)  cut_OneJetHigh.at(OneJetHigh_NotVbf)  = true;
	  if(i==OneJetHigh_TauPt)   cut_OneJetHigh.at(OneJetHigh_TauPt)   = cCat_splitTauPt;
	  if(i==OneJetHigh_HiggsPt) cut_OneJetHigh.at(OneJetHigh_HiggsPt) = 100.0;

	  if(i==OneJetBoost_NJet)   cut_OneJetBoost.at(OneJetBoost_NJet)    = 1;
	  if(i==OneJetBoost_NotVbf) cut_OneJetBoost.at(OneJetBoost_NotVbf)  = true;
	  if(i==OneJetBoost_TauPt)  cut_OneJetBoost.at(OneJetBoost_TauPt)   = cCat_splitTauPt;
	  if(i==OneJetBoost_HiggsPt)cut_OneJetBoost.at(OneJetBoost_HiggsPt) = 100.0;

	  if(i==ZeroJetHigh_NJet) 	cut_ZeroJetHigh.at(ZeroJetHigh_NJet) 	= 0;
	  if(i==ZeroJetHigh_TauPt)	cut_ZeroJetHigh.at(ZeroJetHigh_TauPt) 	= cCat_splitTauPt;

	  if(i==ZeroJetLow_NJet) 	cut_ZeroJetLow.at(ZeroJetLow_NJet)	 	= 0;
	  if(i==ZeroJetLow_TauPt)	cut_ZeroJetLow.at(ZeroJetLow_TauPt) 	= cCat_splitTauPt;

	  // relaxed categories
	  if(i==VbfTight_NJet)		cut_VBFTightRelaxed.at(VbfTight_NJet)		= 2;
	  if(i==VbfTight_DeltaEta)	cut_VBFTightRelaxed.at(VbfTight_DeltaEta)	= 2.0;
	  if(i==VbfTight_NJetRapGap)cut_VBFTightRelaxed.at(VbfTight_NJetRapGap)= 0;
	  if(i==VbfTight_JetInvM)	cut_VBFTightRelaxed.at(VbfTight_JetInvM)	= 200.0;
	  if(i==VbfTight_HiggsPt)	cut_VBFTightRelaxed.at(VbfTight_HiggsPt)	= 100.0;

	  if(i==VbfLoose_NJet)		cut_VBFLooseRelaxed.at(VbfLoose_NJet)		= 2;
	  if(i==VbfLoose_DeltaEta)	cut_VBFLooseRelaxed.at(VbfLoose_DeltaEta)	= 2.0;
	  if(i==VbfLoose_NJetRapGap)cut_VBFLooseRelaxed.at(VbfLoose_NJetRapGap)= 0;
	  if(i==VbfLoose_JetInvM)	cut_VBFLooseRelaxed.at(VbfLoose_JetInvM)	= 200.0;
	  if(i==VbfLoose_NotVbfTight)cut_VBFLooseRelaxed.at(VbfLoose_NotVbfTight)	= true; // disabled, set to true here
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i_cut=0; i_cut<NCuts; i_cut++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    TString c="_Cut_";c+=i_cut;
  
    if(i_cut==PrimeVtx){
      title.at(i_cut)="Number of Prime Vertices $(N>$";
      title.at(i_cut)+=cut.at(PrimeVtx);
      title.at(i_cut)+=")";
      htitle=title.at(i_cut);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,51,-0.5,50.5,hlabel,"Events"));
    }
    else if(i_cut==TriggerOk){
      title.at(i_cut)="Trigger ";
      hlabel="Trigger ";

      std::vector<TH1D> Nm1Temp = HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,cTriggerNames.size()+2,-1.5,cTriggerNames.size()+0.5,hlabel,"Events");
      std::vector<TH1D> Nm0Temp = HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,cTriggerNames.size()+2,-1.5,cTriggerNames.size()+0.5,hlabel,"Events");
      for (unsigned i_hist = 0; i_hist < Nm1Temp.size(); i_hist++){
    	  Nm1Temp.at(i_hist).GetXaxis()->SetBinLabel(1,"not fired");
    	  Nm0Temp.at(i_hist).GetXaxis()->SetBinLabel(1,"not fired");
    	  for (unsigned i_bin = 2; i_bin < cTriggerNames.size()+2; i_bin++){
    		  Nm1Temp.at(i_hist).GetXaxis()->SetBinLabel(i_bin,cTriggerNames.at(i_bin-2));
    		  Nm0Temp.at(i_hist).GetXaxis()->SetBinLabel(i_bin,cTriggerNames.at(i_bin-2));
    	  }
    	  Nm1Temp.at(i_hist).GetXaxis()->SetBinLabel(cTriggerNames.size()+2,"multiple fired");
    	  Nm0Temp.at(i_hist).GetXaxis()->SetBinLabel(cTriggerNames.size()+2,"multiple fired");
      }
      Nminus1.push_back(Nm1Temp);
      Nminus0.push_back(Nm0Temp);
    }
    else if(i_cut==NMuId){
    	title.at(i_cut)="Number $\\mu_{ID} >=$";
    	title.at(i_cut)+=cut.at(NMuId);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{ID}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuId_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==NMuKin){
    	title.at(i_cut)="Number $\\mu_{sel} >=$";
    	title.at(i_cut)+=cut.at(NMuKin);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #mu_{sel}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuKin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i_cut==DiMuonVeto){
    	title.at(i_cut)="$\\Delta R(\\mu_{veto}^{+},\\mu_{veto}^{-}) <$";
        char buffer[50];
        sprintf(buffer,"%5.1f",cut.at(DiMuonVeto));
    	title.at(i_cut)+=buffer;
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="#DeltaR(#mu_{veto}^{+},#mu_{veto}^{-})";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_DiMuonVeto_",htitle,100,0.,5.0,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_DiMuonVeto_",htitle,100,0.,5.0,hlabel,"Events"));
    }
    else if(i_cut==NTauId){
    	title.at(i_cut)="Number $\\tau_{ID} >=$";
    	title.at(i_cut)+=cut.at(NTauId);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{ID}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauId_",htitle,26,-0.5,25.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauId_",htitle,26,-0.5,25.5,hlabel,"Events"));
    }
    else if(i_cut==NTauIso){
    	title.at(i_cut)="Number $\\tau_{Iso} >=$";
    	title.at(i_cut)+=cut.at(NTauIso);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{Iso}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauIso_",htitle,16,-0.5,15.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauIso_",htitle,16,-0.5,15.5,hlabel,"Events"));
    }
    else if(i_cut==NTauKin){
    	title.at(i_cut)="Number $\\tau_{sel} >=$";
    	title.at(i_cut)+=cut.at(NTauKin);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of #tau_{sel}";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauKin_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut==TriLeptonVeto){
    	title.at(i_cut)="3 lepton veto: $N(\\mu)+N(e) =$";
    	title.at(i_cut)+=cut.at(TriLeptonVeto);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of tri-lepton veto leptons";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriLeptonVeto_",htitle,5,-0.5,4.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriLeptonVeto_",htitle,5,-0.5,4.5,hlabel,"Events"));
    }
    else if(i_cut==OppCharge){
    	title.at(i_cut)="$q(\\mu)+q(\\tau) =$";
    	title.at(i_cut)+=cut.at(OppCharge);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="q(#mu)+q(#tau)";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_OppCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_OppCharge_",htitle,5,-2.5,2.5,hlabel,"Events"));
    }
    else if(i_cut==MT){
    	title.at(i_cut)="$m_{T}(\\mu,E_{T}^{miss}) <$";
    	title.at(i_cut)+=cut.at(MT);
    	title.at(i_cut)+=" GeV";
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="m_{T}(#mu,E_{T}^{miss})/GeV";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MT_",htitle,50,0.,100.,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MT_",htitle,50,0.,100.,hlabel,"Events"));
    }
    else if(i_cut==BJetVeto){
    	title.at(i_cut)="Number b-Jets $<=$";
    	title.at(i_cut)+=cut.at(BJetVeto);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel="Number of b-Jets";
    	Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_BJetVeto_",htitle,11,-0.5,10.5,hlabel,"Events"));
    	Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_BJetVeto_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i_cut>=CatCut1){
    	// set histograms to dummy values
    	// will be overwritten in configure_{Category} method
    	title.at(i_cut) = "Category Dummy ";
    	title.at(i_cut)	+=(i_cut-CatCut1);
    	htitle=title.at(i_cut);
    	htitle.ReplaceAll("$","");
    	htitle.ReplaceAll("\\","#");
    	hlabel = title.at(i_cut);
    	TString n = Name+c+"_Nminus1_CatDummy";
    	n += (i_cut-CatCut1);
    	n += "_";
    	Nminus1.push_back(HConfig.GetTH1D(n,htitle,50,-50.,50.,hlabel,"Events"));
    	n.ReplaceAll("Nminus1","Nminus0");
    	Nminus0.push_back(HConfig.GetTH1D(n,htitle,50,-50.,50.,hlabel,"Events"));
    }
  } 
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  // Setup Extra Histograms
  NCatFired=HConfig.GetTH1D(Name+"_NCatFired","NCatFired",6,-0.5,5.5,"Num. of passed categories");
  CatFired=HConfig.GetTH1D(Name+"_CatFired","CatFired",8,-0.5,7.5,"Fired Categories");

  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"N(Vtx) before selection");
  NVtxFullSelection=HConfig.GetTH1D(Name+"_NVtxFullSelection","NVtxFullSelection",26,-0.5,25.5,"N(Vertex) after selection");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"N(good Vertex)");
  VtxZ=HConfig.GetTH1D(Name+"_VtxZ","VtxZ",50,-50.0,50.0,"z(Vtx)/cm");
  VtxRho=HConfig.GetTH1D(Name+"_VtxRho","VtxRho",50,0.0,2.0,"#rho(Vtx)/cm");
  VtxPhi=HConfig.GetTH1D(Name+"_VtxPhi","VtxPhi",50,0.0,3.2,"#phi(Vtx)");
  VtxNdof=HConfig.GetTH1D(Name+"_VtxNdof","VtxNdof",50,-0.5,49.5,"NDoF(Vtx Fit)");
  VtxIsfake=HConfig.GetTH1D(Name+"_VtxIsfake","VtxIsfake",2,-0.5,1.5,"IsFake(Vtx)");

  MuDxy=HConfig.GetTH1D(Name+"_MuDxy","MuDxy",60,-0.3,0.3,"d_{xy}(#mu,Vtx)/cm");
  MuDz=HConfig.GetTH1D(Name+"_MuDz","MuDz",60,-.6,.6,"d_{z}(#mu,Vtx)/cm");
  MuRelIso=HConfig.GetTH1D(Name+"_MuRelIso","MuRelIso",50,0.,1.,"relIso(#mu)");
  MuPt=HConfig.GetTH1D(Name+"_MuPt","MuPt",50,0.,200.,"p_{T}(#mu)/GeV");
  MuEta=HConfig.GetTH1D(Name+"_MuEta","MuEta",50,-2.5,2.5,"#eta(#mu)");
  MuPhi=HConfig.GetTH1D(Name+"_MuPhi","MuPhi",50,-3.14159,3.14159,"#phi(#mu)");

  MuSelPt=HConfig.GetTH1D(Name+"_MuSelPt","MuSelPt",50,0.,200.,"p_{T}(#mu_{sel})/GeV");
  MuSelEta=HConfig.GetTH1D(Name+"_MuSelEta","MuSelEta",50,-2.5,2.5,"#eta(#mu_{sel})");
  MuSelPhi=HConfig.GetTH1D(Name+"_MuSelPhi","MuSelPhi",50,-3.14159,3.14159,"#phi(#mu_{sel})");
  MuSelFakesTauID=HConfig.GetTH1D(Name+"_MuSelFakesTauID","MuSelFakesTauID",2,-0.5,1.5,"#mu_{sel} fakes #tau_{h}");
  MuSelDrHlt=HConfig.GetTH1D(Name+"_MuSelDrHlt","MuSelDrHLT",50,0.,1.,"#DeltaR(#mu_{sel},#mu_{HLT})");

  TauPt=HConfig.GetTH1D(Name+"_TauPt","TauPt",50,0.,200.,"p_{T}(#tau)/GeV");
  TauEta=HConfig.GetTH1D(Name+"_TauEta","TauEta",50,-2.5,2.5,"#eta(#tau)");
  TauPhi=HConfig.GetTH1D(Name+"_TauPhi","TauPhi",50,-3.14159,3.14159,"#phi(#tau)");

  TauSelPt=HConfig.GetTH1D(Name+"_TauSelPt","TauSelPt",50,0.,200.,"p_{T}(#tau_{sel})/GeV");
  TauSelEta=HConfig.GetTH1D(Name+"_TauSelEta","TauSelEta",50,-2.5,2.5,"#eta(#tau_{sel})");
  TauSelPhi=HConfig.GetTH1D(Name+"_TauSelPhi","TauSelPhi",50,-3.14159,3.14159,"#phi(#tau_{sel})");
  TauSelDrHlt=HConfig.GetTH1D(Name+"_TauSelDrHlt","TauSelDrHLT",50,0.,1.,"#DeltaR(#tau_{sel},#tau_{HLT})");
  TauSelDecayMode=HConfig.GetTH1D(Name+"_TauSelDecayMode","TauSelDecayMode",16,-0.5,15.5,"#tau_{sel} decay mode");

  MuVetoDPtSelMuon=HConfig.GetTH1D(Name+"_MuVetoDPtSelMuon","MuVetoDPtSelMuon",100,-100.,100.,"#Deltap_{T}(#mu_{veto},#mu)/GeV");
  MuVetoInvM=HConfig.GetTH1D(Name+"_MuVetoInvM","MuVetoInvM",100,0.,200,"m_{inv}(#mu_{veto}^{1},#mu_{veto}^{2})/GeV");
  MuVetoPtPositive=HConfig.GetTH1D(Name+"_MuVetoPtPositive","MuVetoPtPositive",50,0.,200.,"p_{T}(#mu_{veto}^{+})/GeV");
  MuVetoPtNegative=HConfig.GetTH1D(Name+"_MuVetoPtNegative","MuVetoPtNegative",50,0.,200.,"p_{T}(#mu_{veto}^{-})/GeV");
  MuVetoDRTau=HConfig.GetTH1D(Name+"_MuVetoDRTau","MuVetoDRTau",50,0.,5.,"#DeltaR(#mu_{veto},#tau_{h})");
  MuVetoDeltaR=HConfig.GetTH1D(Name+"_MuVetoDeltaR","MuVetoDeltaR",50,0.,5.,"#DeltaR(#mu^{+}_{veto},#mu^{-}_{veto})");

  NMuonTriLepVeto=HConfig.GetTH1D(Name+"_NMuonTriLepVeto","NMuonTriLepVeto",5,-0.5,4.5,"N(#mu_{3l veto})");
  NElecTriLepVeto=HConfig.GetTH1D(Name+"_NElecTriLepVeto","NElecTriLepVeto",5,-0.5,4.5,"N(e_{3l veto})");

  MuCharge=HConfig.GetTH1D(Name+"_MuCharge","MuCharge",3,-1.5,1.5,"q(#mu)/e");
  TauCharge=HConfig.GetTH1D(Name+"_TauCharge","TauCharge",3,-1.5,1.5,"q(#tau)/e");

  MuTauDR=HConfig.GetTH1D(Name+"_MuTauDR","MuTauDR",50,0.,5.,"#DeltaR(#mu,#tau_{h})");
  MuTauDPhi=HConfig.GetTH1D(Name+"_MuTauDPhi","MuTauDPhi",50,0.,3.2,"#Delta#phi(#mu,#tau_{h})");
  MuTauDEta=HConfig.GetTH1D(Name+"_MuTauDEta","MuTauDEta",100,-6.,6.,"#Delta#eta(#mu,#tau_{h})");
  MuTauDPt=HConfig.GetTH1D(Name+"_MuTauDPt","MuTauDPt",100,-100.,100.,"#Deltap_{T}(#mu,#tau_{h})/GeV");
  MuTauRelDPt=HConfig.GetTH1D(Name+"_MuTauRelDPt","MuTauRelDPt",100,-2.,2.,"#Deltap_{T}(#mu,#tau_{h})/p_{T}(#mu)");

  MetPt  = HConfig.GetTH1D(Name+"_MetPt","MetPt",50,0.,200.,"E_{T}^{miss}/GeV");
  MetPhi = HConfig.GetTH1D(Name+"_MetPhi","MetPhi",50,-3.14159,3.14159,"#phi(E_{T}^{miss})");

  NJetsKin = HConfig.GetTH1D(Name+"_NJetsKin","NJetsKin",11,-0.5,10.5,"N(j_{kin})");
  JetKin1Pt = HConfig.GetTH1D(Name+"_JetKin1Pt","JetKin1Pt",50,0.,200.,"p_{T}(j_{kin}^{1})/GeV");
  JetKin1Eta = HConfig.GetTH1D(Name+"_JetKin1Eta","JetKin1Eta",100,-5.0,5.0,"#eta(j_{kin}^{1})");
  JetKin1Phi = HConfig.GetTH1D(Name+"_JetKin1Phi","JetKin1Phi",50,-3.14159,3.14159,"#phi(j_{kin}^{1})");
  JetKin1IsLooseId = HConfig.GetTH1D(Name+"_JetKin1IsLooseId","JetKin1IsLooseId",2,-0.5,1.5,"isLoosePUJetID(j_{kin}^{1}");
  JetKin2IsLooseId = HConfig.GetTH1D(Name+"_JetKin2IsLooseId","JetKin2IsLooseId",2,-0.5,1.5,"isLoosePUJetID(j_{kin}^{2}");
  JetKin2Pt = HConfig.GetTH1D(Name+"_JetKin2Pt","JetKin2Pt",50,0.,200.,"p_{T}(j_{kin}^{2})/GeV");
  JetKin2Eta = HConfig.GetTH1D(Name+"_JetKin2Eta","JetKin2Eta",100,-5.0,5.0,"#eta(j_{kin}^{2})");
  JetKin2Phi = HConfig.GetTH1D(Name+"_JetKin2Phi","JetKin2Phi",50,-3.14159,3.14159,"#phi(j_{kin}^{2})");
  NJetsId = HConfig.GetTH1D(Name+"_NJetsId","NJetsId",11,-0.5,10.5,"N(jets)");
  Jet1Pt = HConfig.GetTH1D(Name+"_Jet1Pt","Jet1Pt",50,0.,200.,"p_{T}(j^{1})/GeV");
  Jet1Eta = HConfig.GetTH1D(Name+"_Jet1Eta","Jet1Eta",100,-5.0,5.0,"#eta(j^{1})");
  Jet1Phi = HConfig.GetTH1D(Name+"_Jet1Phi","Jet1Phi",50,-3.14159,3.14159,"#phi(j^{1})");
  Jet1IsB = HConfig.GetTH1D(Name+"_Jet1IsB","Jet1IsB",2,-0.5,1.5,"isBJet(j^{1})");
  Jet2Pt = HConfig.GetTH1D(Name+"_Jet2Pt","Jet2Pt",50,0.,200.,"p_{T}(j^{2})/GeV");
  Jet2Eta = HConfig.GetTH1D(Name+"_Jet2Eta","Jet2Eta",100,-5.0,5.0,"#eta(j^{2})");
  Jet2Phi = HConfig.GetTH1D(Name+"_Jet2Phi","Jet2Phi",50,-3.14159,3.14159,"#phi(j^{2})");
  Jet2IsB = HConfig.GetTH1D(Name+"_Jet2IsB","Jet2IsB",2,-0.5,1.5,"isBJet(j^{2})");

  NBJets = HConfig.GetTH1D(Name+"_NBJets","NBJets",6,-0.5,5.5,"N(bjets)");
  BJet1Pt = HConfig.GetTH1D(Name+"_BJet1Pt","BJet1Pt",50,0.,200.,"p_{T}(b^{1})/GeV");
  BJet1Eta = HConfig.GetTH1D(Name+"_BJet1Eta","BJet1Eta",100,-5.0,5.0,"#eta(b^{1})");
  BJet1Phi = HConfig.GetTH1D(Name+"_BJet1Phi","BJet1Phi",50,-3.14159,3.14159,"#phi(b^{1})");

  HiggsPt = HConfig.GetTH1D(Name+"_HiggsPt","HiggsPt",50,0.,200.,"p_{T}(H)/GeV");
  HiggsPhi = HConfig.GetTH1D(Name+"_HiggsPhi","HiggsPhi",50,-3.14159,3.14159,"#phi(H)");
  JetsDEta = HConfig.GetTH1D(Name+"_JetsDEta","JetsDEta",100,-10.,10.,"#Delta#eta(j^{1},j^{2})");
  JetsInEtaGap = HConfig.GetTH1D(Name+"_JetsInEtaGap","JetsInEtaGap",6,-0.5,5.5,"N(j in #eta gap)");
  JetsInvM = HConfig.GetTH1D(Name+"_JetsInvM","JetsInvM",100,0.,2000.,"m_{inv}(j^{1},j^{2})");

  Cat0JetLowMt = HConfig.GetTH1D(Name+"_Cat0JetLowMt","Cat0JetLowMt",125,0.,250.,"0JL: m_{T}/GeV");
  Cat0JetLowMtSideband = HConfig.GetTH1D(Name+"_Cat0JetLowMtSideband","Cat0JetLowMtSideband",90,70.,250.,"0JL: m_{T}/GeV");
  Cat0JetLowMtExtrapolation = HConfig.GetTH1D(Name+"_Cat0JetLowMtExtrapolation","Cat0JetLowMtExtrapolation",2,0.5,2.5,"0JL: m_{T} signal and sideband");
  Cat0JetHighMt = HConfig.GetTH1D(Name+"_Cat0JetHighMt","Cat0JetHighMt",125,0.,250.,"0JH: m_{T}/GeV");
  Cat0JetHighMtSideband = HConfig.GetTH1D(Name+"_Cat0JetHighMtSideband","Cat0JetHighMtSideband",90,70.,250.,"0JH: m_{T}/GeV");
  Cat0JetHighMtExtrapolation = HConfig.GetTH1D(Name+"_Cat0JetHighMtExtrapolation","Cat0JetHighMtExtrapolation",2,0.5,2.5,"0JH: m_{T} signal and sideband");
  Cat1JetLowMt = HConfig.GetTH1D(Name+"_Cat1JetLowMt","Cat1JetLowMt",125,0.,250.,"1JL: m_{T}/GeV");
  Cat1JetLowMtSideband = HConfig.GetTH1D(Name+"_Cat1JetLowMtSideband","Cat1JetLowMtSideband",90,70.,250.,"1JL: m_{T}/GeV");
  Cat1JetLowMtExtrapolation = HConfig.GetTH1D(Name+"_Cat1JetLowMtExtrapolation","Cat1JetLowMtExtrapolation",2,0.5,2.5,"1JL: m_{T} signal and sideband");
  Cat1JetHighMt = HConfig.GetTH1D(Name+"_Cat1JetHighMt","Cat1JetHighMt",125,0.,250.,"1JH: m_{T}/GeV");
  Cat1JetHighMtSideband = HConfig.GetTH1D(Name+"_Cat1JetHighMtSideband","Cat1JetHighMtSideband",90,70.,250.,"1JH: m_{T}/GeV");
  Cat1JetHighMtExtrapolation = HConfig.GetTH1D(Name+"_Cat1JetHighMtExtrapolation","Cat1JetHighMtExtrapolation",2,0.5,2.5,"1JH: m_{T} signal and sideband");
  Cat1JetBoostMt = HConfig.GetTH1D(Name+"_Cat1JetBoostMt","Cat1JetBoostMt",125,0.,250.,"1JB: m_{T}/GeV");
  Cat1JetBoostMtSideband = HConfig.GetTH1D(Name+"_Cat1JetBoostMtSideband","Cat1JetBoostMtSideband",90,70.,250.,"1JB: m_{T}/GeV");
  Cat1JetBoostMtExtrapolation = HConfig.GetTH1D(Name+"_Cat1JetBoostMtExtrapolation","Cat1JetBoostMtExtrapolation",2,0.5,2.5,"1JB: m_{T} signal and sideband");
  CatVBFLooseMt = HConfig.GetTH1D(Name+"_CatVBFLooseMt","CatVBFLooseMt",125,0.,250.,"VBFL: m_{T}/GeV");
  CatVBFLooseMtSideband = HConfig.GetTH1D(Name+"_CatVBFLooseMtSideband","CatVBFLooseMtSideband",30,60.,120.,"VBFL: m_{T}/GeV");
  CatVBFLooseRelaxMt = HConfig.GetTH1D(Name+"_CatVBFLooseRelaxMt","CatVBFLooseRelaxMt",125,0.,250.,"VBFLRelax: m_{T}/GeV");
  CatVBFLooseRelaxMtExtrapolation = HConfig.GetTH1D(Name+"_CatVBFLooseRelaxMtExtrapolation","CatVBFLooseRelaxMtExtrapolation",2,0.5,2.5,"VBFLRelax: m_{T} signal and sideband");
  CatVBFTightMt = HConfig.GetTH1D(Name+"_CatVBFTightMt","CatVBFTightMt",125,0.,250.,"VBFT: m_{T}/GeV");
  CatVBFTightMtSideband = HConfig.GetTH1D(Name+"_CatVBFTightMtSideband","CatVBFTightMtSideband",30,60.,120.,"VBFT: m_{T}/GeV");
  CatVBFTightRelaxMt = HConfig.GetTH1D(Name+"_CatVBFTightRelaxMt","CatVBFTightRelaxMt",125,0.,250.,"VBFTRelax: m_{T}/GeV");
  CatVBFTightRelaxMtExtrapolation = HConfig.GetTH1D(Name+"_CatVBFTightRelaxMtExtrapolation","CatVBFTightRelaxMtExtrapolation",2,0.5,2.5,"VBFTRelax: m_{T} signal and sideband");
  CatInclusiveMt = HConfig.GetTH1D(Name+"_CatInclusiveMt","CatInclusiveMt",125,0.,250.,"Incl: m_{T}/GeV");
  CatInclusiveMtSideband = HConfig.GetTH1D(Name+"_CatInclusiveMtSideband","CatInclusiveMtSideband",90,70.,250.,"Incl: m_{T}/GeV");
  CatInclusiveMtExtrapolation = HConfig.GetTH1D(Name+"_CatInclusiveMtExtrapolation","CatInclusiveMtExtrapolation",2,0.5,2.5,"Incl: m_{T} signal and sideband");

  Cat0JetLowQcdAbcd = HConfig.GetTH2D(Name+"_Cat0JetLowQcdAbcd","Cat0JetLowQcdAbcd",2,-1.,1.,2,-1.,1.,"0JL: relIso(#mu)","OS/SS");
  Cat0JetHighQcdAbcd = HConfig.GetTH2D(Name+"_Cat0JetHighQcdAbcd","Cat0JetHighQcdAbcd",2,-1.,1.,2,-1.,1.,"0JH:relIso(#mu)","OS/SS");
  Cat1JetLowQcdAbcd = HConfig.GetTH2D(Name+"_Cat1JetLowQcdAbcd","Cat1JetLowQcdAbcd",2,-1.,1.,2,-1.,1.,"1JL:relIso(#mu)","OS/SS");
  Cat1JetHighQcdAbcd = HConfig.GetTH2D(Name+"_Cat1JetHighQcdAbcd","Cat1JetHighQcdAbcd",2,-1.,1.,2,-1.,1.,"1JH:relIso(#mu)","OS/SS");
  Cat1JetBoostQcdAbcd = HConfig.GetTH2D(Name+"_Cat1JetBoostQcdAbcd","Cat1JetBoostQcdAbcd",2,-1.,1.,2,-1.,1.,"1JB:relIso(#mu)","OS/SS");
  CatVBFLooseQcdAbcd = HConfig.GetTH2D(Name+"_CatVBFLooseQcdAbcd","CatVBFLooseQcdAbcd",2,-1.,1.,2,-1.,1.,"VBFL:relIso(#mu)","OS/SS");
  CatVBFTightQcdAbcd = HConfig.GetTH2D(Name+"_CatVBFTightQcdAbcd","CatVBFTightQcdAbcd",2,-1.,1.,2,-1.,1.,"VBFT:relIso(#mu)","OS/SS");
  CatInclusiveQcdAbcd = HConfig.GetTH2D(Name+"_CatInclusiveQcdAbcd","CatInclusiveQcdAbcd",2,-1.,1.,2,-1.,1.,"Incl:relIso(#mu)","OS/SS");

  // configure category
  if (categoryFlag == "VBFTight")	configure_VBFTight();
  else if (categoryFlag == "VBFLoose")	configure_VBFLoose();
  else if (categoryFlag == "OneJetHigh")	configure_OneJetHigh();
  else if (categoryFlag == "OneJetLow")	configure_OneJetLow();
  else if (categoryFlag == "OneJetBoost")configure_OneJetBoost();
  else if (categoryFlag == "ZeroJetHigh")configure_ZeroJetHigh();
  else if (categoryFlag == "ZeroJetLow") configure_ZeroJetLow();
  else if (categoryFlag == "Inclusive") configure_NoCategory();
  else if (categoryFlag == "NoCategory")	configure_NoCategory();
  else{
	  std::cout << "WARNING: category " << categoryFlag << " does not exist. Using NoCategory instead." << std::endl;
	  configure_NoCategory();
  }

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types, CrossSectionandAcceptance, legend, colour);
}


void  HToTaumuTauh::Store_ExtraDist(){
 if (verbose) std::cout << "HToTaumuTauh::Store_ExtraDist()" << std::endl;
 Extradist1d.push_back(&NCatFired);
 Extradist1d.push_back(&CatFired);

 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NVtxFullSelection);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&VtxZ);
 Extradist1d.push_back(&VtxRho);
 Extradist1d.push_back(&VtxNdof);
 Extradist1d.push_back(&VtxIsfake);

 Extradist1d.push_back(&MuDxy);
 Extradist1d.push_back(&MuDz );
 Extradist1d.push_back(&MuRelIso);
 Extradist1d.push_back(&MuPt  );
 Extradist1d.push_back(&MuEta  );
 Extradist1d.push_back(&MuPhi  );

 Extradist1d.push_back(&MuSelPt  );
 Extradist1d.push_back(&MuSelEta  );
 Extradist1d.push_back(&MuSelPhi  );
 Extradist1d.push_back(&MuSelFakesTauID  );

 Extradist1d.push_back(&TauPt  );
 Extradist1d.push_back(&TauEta  );
 Extradist1d.push_back(&TauPhi  );

 Extradist1d.push_back(&TauSelPt  );
 Extradist1d.push_back(&TauSelEta  );
 Extradist1d.push_back(&TauSelPhi  );
 Extradist1d.push_back(&TauSelDecayMode  );

 Extradist1d.push_back(&MuVetoDPtSelMuon);
 Extradist1d.push_back(&MuVetoInvM);
 Extradist1d.push_back(&MuVetoPtPositive);
 Extradist1d.push_back(&MuVetoPtNegative);
 Extradist1d.push_back(&MuVetoDRTau);
 Extradist1d.push_back(&MuVetoDeltaR);

 Extradist1d.push_back(&NMuonTriLepVeto);
 Extradist1d.push_back(&NElecTriLepVeto);

 Extradist1d.push_back(&MuCharge  );
 Extradist1d.push_back(&TauCharge  );

 Extradist1d.push_back(&MuTauDR);
 Extradist1d.push_back(&MuTauDPhi);
 Extradist1d.push_back(&MuTauDEta);
 Extradist1d.push_back(&MuTauDPt);
 Extradist1d.push_back(&MuTauRelDPt);

 Extradist1d.push_back(&MetPt);
 Extradist1d.push_back(&MetPhi);

 Extradist1d.push_back(&NJetsKin);
 Extradist1d.push_back(&JetKin1Pt);
 Extradist1d.push_back(&JetKin1Eta);
 Extradist1d.push_back(&JetKin1Phi);
 Extradist1d.push_back(&JetKin1IsLooseId);
 Extradist1d.push_back(&JetKin2IsLooseId);
 Extradist1d.push_back(&JetKin2Pt);
 Extradist1d.push_back(&JetKin2Eta);
 Extradist1d.push_back(&JetKin2Phi);
 Extradist1d.push_back(&NJetsId);
 Extradist1d.push_back(&Jet1Pt);
 Extradist1d.push_back(&Jet1Eta);
 Extradist1d.push_back(&Jet1Phi);
 Extradist1d.push_back(&Jet1IsB);
 Extradist1d.push_back(&Jet2Pt);
 Extradist1d.push_back(&Jet2Eta);
 Extradist1d.push_back(&Jet2Phi);
 Extradist1d.push_back(&Jet2IsB);

 Extradist1d.push_back(&NBJets);
 Extradist1d.push_back(&BJet1Pt);
 Extradist1d.push_back(&BJet1Eta);
 Extradist1d.push_back(&BJet1Phi);

 Extradist1d.push_back(&HiggsPt);
 Extradist1d.push_back(&HiggsPhi);
 Extradist1d.push_back(&JetsDEta);
 Extradist1d.push_back(&JetsInEtaGap);
 Extradist1d.push_back(&JetsInvM);

 Extradist1d.push_back(&Cat0JetLowMt);
 Extradist1d.push_back(&Cat0JetLowMtSideband);
 Extradist1d.push_back(&Cat0JetLowMtExtrapolation);
 Extradist1d.push_back(&Cat0JetHighMt);
 Extradist1d.push_back(&Cat0JetHighMtSideband);
 Extradist1d.push_back(&Cat0JetHighMtExtrapolation);
 Extradist1d.push_back(&Cat1JetLowMt);
 Extradist1d.push_back(&Cat1JetLowMtSideband);
 Extradist1d.push_back(&Cat1JetLowMtExtrapolation);
 Extradist1d.push_back(&Cat1JetHighMt);
 Extradist1d.push_back(&Cat1JetHighMtSideband);
 Extradist1d.push_back(&Cat1JetHighMtExtrapolation);
 Extradist1d.push_back(&Cat1JetBoostMt);
 Extradist1d.push_back(&Cat1JetBoostMtSideband);
 Extradist1d.push_back(&Cat1JetBoostMtExtrapolation);
 Extradist1d.push_back(&CatVBFLooseMt);
 Extradist1d.push_back(&CatVBFLooseMtSideband);
 Extradist1d.push_back(&CatVBFLooseRelaxMt);
 Extradist1d.push_back(&CatVBFLooseRelaxMtExtrapolation);
 Extradist1d.push_back(&CatVBFTightMt);
 Extradist1d.push_back(&CatVBFTightMtSideband);
 Extradist1d.push_back(&CatVBFTightRelaxMt);
 Extradist1d.push_back(&CatVBFTightRelaxMtExtrapolation);
 Extradist1d.push_back(&CatInclusiveMt);
 Extradist1d.push_back(&CatInclusiveMtSideband);
 Extradist1d.push_back(&CatInclusiveMtExtrapolation);

 Extradist2d.push_back(&Cat0JetLowQcdAbcd);
 Extradist2d.push_back(&Cat0JetHighQcdAbcd);
 Extradist2d.push_back(&Cat1JetLowQcdAbcd);
 Extradist2d.push_back(&Cat1JetHighQcdAbcd);
 Extradist2d.push_back(&Cat1JetBoostQcdAbcd);
 Extradist2d.push_back(&CatVBFLooseQcdAbcd);
 Extradist2d.push_back(&CatVBFTightQcdAbcd);
 Extradist2d.push_back(&CatInclusiveQcdAbcd);

}

void  HToTaumuTauh::doEvent(){
  if (verbose) std::cout << "HToTaumuTauh::doEvent()" << std::endl;
  // set variables to hold selected objects to default values
  selVertex = -1;
  selMuon = -1;
  selTau = -1;
  selJets.clear();
  selBJets.clear();
  selMjj = -1;
  selJetdeta = -100;
  selNjetingap = -1;

  unsigned int t;
  int id(Ntp->GetMCID());
  //std::cout << "ID before = " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}
  
  double wobs=1;
  double w;
  if(!Ntp->isData()){w = Ntp->PUWeightFineBins();}
  else{w=1;}

  // Apply Selection

  // Vertex
  unsigned int nGoodVtx=0;
  for(unsigned int i_vtx=0;i_vtx<Ntp->NVtx();i_vtx++){
    if(Ntp->isGoodVtx(i_vtx)){
    	if(selVertex == -1) selVertex = i_vtx; // selected vertex = first vertex (highest sum[pT^2]) to fulfill vertex requirements
    	nGoodVtx++;
    }
  }

  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  
  // Trigger
  value.at(TriggerOk) = -1;
  for (std::vector<TString>::iterator it_trig = cTriggerNames.begin(); it_trig != cTriggerNames.end(); ++it_trig){
	  if(Ntp->TriggerAccept(*it_trig)){
		  if ( value.at(TriggerOk) == -1 )
			  value.at(TriggerOk) = it_trig - cTriggerNames.begin();
		  else // more than 1 trigger fired, save this separately
			  value.at(TriggerOk) = cTriggerNames.size();
	  }
  }
  pass.at(TriggerOk) = (value.at(TriggerOk) >= cut.at(TriggerOk));
  
  // Muon cuts
  std::vector<int> selectedMuonsId;
  selectedMuonsId.clear();
  for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
	  if( selectMuon_Id(i_mu,selVertex) ) {
		  selectedMuonsId.push_back(i_mu);
	  }
  }
  value.at(NMuId)=selectedMuonsId.size();
  pass.at(NMuId)=(value.at(NMuId)>=cut.at(NMuId));

  std::vector<int> selectedMuons;	// full selection: ID and Kinematics
  selectedMuons.clear();
  for(std::vector<int>::iterator it_mu = selectedMuonsId.begin(); it_mu != selectedMuonsId.end(); ++it_mu){
	  if( selectMuon_Kinematics(*it_mu)) {
		  selectedMuons.push_back(*it_mu);
	  }
  }
  value.at(NMuKin)=selectedMuons.size();
  pass.at(NMuKin)=(value.at(NMuKin)>=cut.at(NMuKin));

  // muons for QCD background method
  std::vector<int> antiIsoMuons;
  antiIsoMuons.clear();
  for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
	  if( selectMuon_antiIso(i_mu,selVertex) ) {
		  antiIsoMuons.push_back(i_mu);
	  }
  }
  bool hasAntiIsoMuon = (antiIsoMuons.size() > 0);

  if (selectedMuons.size() > 0)
	  selMuon = selectedMuons.at(0); // use signal muon
  if (selectedMuons.size() <= 0 && hasAntiIsoMuon)
	  selMuon = antiIsoMuons.at(0); // for background methods: use anti-iso muon

  std::vector<int> diMuonVetoMuonsPositive;	// muons selected for the dimuon veto
  diMuonVetoMuonsPositive.clear();
  std::vector<int> diMuonVetoMuonsNegative;	// muons selected for the dimuon veto
  diMuonVetoMuonsNegative.clear();
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if( selectMuon_diMuonVeto(i, selVertex) ) {
		  if (Ntp->Muon_Charge(i) == 1) {
			  diMuonVetoMuonsPositive.push_back(i);
		  }
		  else if (Ntp->Muon_Charge(i) == -1) {
			  diMuonVetoMuonsNegative.push_back(i);
		  }
	  }
  }
  if (diMuonVetoMuonsPositive.size() == 0 || diMuonVetoMuonsNegative.size() == 0){
	  value.at(DiMuonVeto) = 0.0;
  }
  else{
	  value.at(DiMuonVeto) = Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)).DeltaR( Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0)) );
  }
  pass.at(DiMuonVeto) = (value.at(DiMuonVeto) < cut.at(DiMuonVeto));

  // Tau cuts
  std::vector<int> selectedTausId;
  selectedTausId.clear();
  for(unsigned i_tau=0; i_tau < Ntp->NPFTaus(); i_tau++){
	  if ( selectPFTau_Id(i_tau,selectedMuonsId) ){
		  selectedTausId.push_back(i_tau);
	  }
  }
  value.at(NTauId)=selectedTausId.size();
  pass.at(NTauId)=(value.at(NTauId)>=cut.at(NTauId));

  std::vector<int> selectedTausIso;
  selectedTausIso.clear();
  for(std::vector<int>::iterator it_tau = selectedTausId.begin(); it_tau != selectedTausId.end(); ++it_tau){
	  if ( selectPFTau_Iso(*it_tau) ){
		  selectedTausIso.push_back(*it_tau);
	  }
  }
  value.at(NTauIso)=selectedTausIso.size();
  pass.at(NTauIso)=(value.at(NTauIso)>=cut.at(NTauIso));

  std::vector<int> selectedTaus;
  selectedTaus.clear();
  for(std::vector<int>::iterator it_tau = selectedTausIso.begin(); it_tau != selectedTausIso.end(); ++it_tau){
	  if ( selectPFTau_Kinematics(*it_tau) ){
		  selectedTaus.push_back(*it_tau);
	  }
  }
  value.at(NTauKin)=selectedTaus.size();
  pass.at(NTauKin)=(value.at(NTauKin)>=cut.at(NTauKin));

  // taus for QCD background method
  std::vector<int> relaxedIsoTaus;
  relaxedIsoTaus.clear();
  for(unsigned i_tau=0; i_tau < Ntp->NPFTaus(); i_tau++){
	  if ( selectPFTau_relaxedIso(i_tau,selectedMuonsId) ){
		  relaxedIsoTaus.push_back(i_tau);
	  }
  }
  bool hasRelaxedIsoTau = (relaxedIsoTaus.size() > 0);

  if(selectedTaus.size() > 0)
	  selTau = selectedTaus.at(0); // use signal tau
  if(selectedTaus.size() <= 0 && hasRelaxedIsoTau)
	  selTau = relaxedIsoTaus.at(0); // relaxed isolation tau

  // Tri-lepton veto
  std::vector<int> triLepVetoMuons;
  triLepVetoMuons.clear();
  for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
	  if( selectMuon_triLeptonVeto(i_mu,selMuon,selVertex) ) {
		  triLepVetoMuons.push_back(i_mu);
	  }
  }
  std::vector<int> triLepVetoElecs;
  triLepVetoElecs.clear();
  for(unsigned i_el=0;i_el<Ntp->NElectrons();i_el++){
	  if( selectElectron_triLeptonVeto(i_el,selVertex,selectedMuonsId) ) {
		  triLepVetoElecs.push_back(i_el);
	  }
  }
  value.at(TriLeptonVeto) = triLepVetoMuons.size() + triLepVetoElecs.size();
  pass.at(TriLeptonVeto) = (value.at(TriLeptonVeto) <= cut.at(TriLeptonVeto));

  // Opposite charge
  if (selMuon != -1 && selTau != -1){
	  value.at(OppCharge) = Ntp->Muon_Charge(selMuon) + Ntp->PFTau_Charge(selTau);
  }
  else {
	  value.at(OppCharge) = -9;
	  pass.at(OppCharge) = true;
  }
  if (cut.at(OppCharge) == 999) // set to 999 to disable oppcharge cut
	  pass.at(OppCharge) = true;
  else
	  pass.at(OppCharge) = (value.at(OppCharge) == cut.at(OppCharge));

  // Transverse mass
  if(selMuon == -1){ // no good muon in event: set MT to small dummy value -10 -> pass cut
	  value.at(MT) = -10.0;
  }
  else{
	  double pT 	= Ntp->Muon_p4(selMuon).Pt();
	  double phi	= Ntp->Muon_p4(selMuon).Phi();
	  double eTmiss = Ntp->MET_CorrMVAMuTau_et();
	  double eTmPhi = Ntp->MET_CorrMVAMuTau_phi();
	  value.at(MT)	= transverseMass(pT,phi,eTmiss,eTmPhi);
  }
  if (cut.at(MT) == 999) // set to 999 to disable mt cut
	  pass.at(MT) = true;
  else
	  pass.at(MT) = (value.at(MT) < cut.at(MT));

  // sort jets by corrected pt
  std::vector<int> sortedPFJets = sortPFjets();
  // select jets for categories
  // PFJet and bjet collections can have mutual elements!
  std::vector<int> selectedJetsClean;
  selectedJetsClean.clear();
  std::vector<int> selectedJetsKin;
  selectedJetsKin.clear();
  std::vector<int> selectedJets;
  selectedJets.clear();
  std::vector<int> selectedBJets;
  selectedBJets.clear();
  for (unsigned i_jet = 0; i_jet < Ntp->NPFJets(); i_jet++){
	  if ( selectPFJet_Cleaning(sortedPFJets.at(i_jet), selMuon, selTau)){
		  selectedJetsClean.push_back(sortedPFJets.at(i_jet));
		  if ( selectPFJet_Kinematics(sortedPFJets.at(i_jet)) ) {
			  selectedJetsKin.push_back(sortedPFJets.at(i_jet));
			  if ( selectPFJet_Id(sortedPFJets.at(i_jet)) ){
				  selectedJets.push_back(sortedPFJets.at(i_jet));
			  }
		  }
	  }
	  if ( selectBJet(sortedPFJets.at(i_jet), selMuon, selTau) ) {
		  selectedBJets.push_back(sortedPFJets.at(i_jet));
	  }
  }
  selJets = selectedJetsClean;
  selBJets = selectedBJets;

  // b-Jet veto
  value.at(BJetVeto) = selectedBJets.size();
  pass.at(BJetVeto) = (value.at(BJetVeto) <= cut.at(BJetVeto));

  // store pt of selected tau for categories
  double tauPt = -12;
  if (selTau != -1){
	  tauPt = Ntp->PFTau_p4(selTau).Pt();
  }

  // calculate pt of higgs candidatef
  double higgsPt = -10;
  double higgsPhi = -10;
  if (selMuon != -1 && selTau != -1){
	  TVector3 muon3Vec = Ntp->Muon_p4(selMuon).Vect();
	  TVector3 tau3Vec = Ntp->PFTau_p4(selTau).Vect();
	  TVector3 met3Vec = TVector3(Ntp->MET_CorrMVAMuTau_ex(),Ntp->MET_CorrMVAMuTau_ey(),0);

	  higgsPt = (muon3Vec + tau3Vec + met3Vec).Pt();
	  higgsPhi = (muon3Vec + tau3Vec + met3Vec).Phi();
  }

  // calculate jet-related variables used by categories
  unsigned nJets = selectedJets.size();

  if (nJets >= 2){
	  double vbfJetEta1 = Ntp->PFJet_p4(selectedJets.at(0)).Eta();
	  double vbfJetEta2 = Ntp->PFJet_p4(selectedJets.at(1)).Eta();
	  selJetdeta = vbfJetEta1 - vbfJetEta2;

	  int jetsInRapidityGap = 0;
	  for(std::vector<int>::iterator it_jet = selectedJets.begin()+2; it_jet != selectedJets.end(); ++it_jet){
		  double etaPos = ( selJetdeta >= 0) ? vbfJetEta1 : vbfJetEta2;
		  double etaNeg = ( selJetdeta >= 0) ? vbfJetEta2 : vbfJetEta1;
		  if (	Ntp->PFJet_p4(*it_jet).Eta() > etaNeg &&
				Ntp->PFJet_p4(*it_jet).Eta() < etaPos){
			  jetsInRapidityGap++;
		  }
	  }
	  selNjetingap = jetsInRapidityGap;

	  double invM = (Ntp->PFJet_p4(selectedJets.at(0)) + Ntp->PFJet_p4(selectedJets.at(1))).M();
	  selMjj = invM;
  }
  else{
	  selJetdeta = -100;
	  selNjetingap = -1;
	  selMjj = -1;
  }


  // define booleans for different stages of selection
  bool passedVertex = pass.at(TriggerOk) && pass.at(PrimeVtx);
  bool passedMuId = passedVertex && pass.at(NMuId);
  bool passedMu = passedMuId && pass.at(NMuKin);
  bool passedTauIdIso = passedVertex && pass.at(NTauId) && pass.at(NTauIso);
  bool passedTau = passedTauIdIso && pass.at(NTauKin);
  bool passedObjects = passedMu && passedTau;
  bool passedDiMuonVeto = passedObjects && pass.at(DiMuonVeto);
  bool passedFullInclusiveSelNoBVeto = passedDiMuonVeto && pass.at(TriLeptonVeto) && pass.at(OppCharge) && pass.at(MT);
  bool passedFullInclusiveSel = passedFullInclusiveSelNoBVeto && pass.at(BJetVeto);

  // define booleans for analysis stages needed for background methods
  bool passedFullInclusiveSelNoMt = passedObjects && pass.at(DiMuonVeto) && pass.at(TriLeptonVeto) && pass.at(OppCharge) && pass.at(BJetVeto);
  bool passedFullInclusiveSelNoMtNoOS = passedObjects && pass.at(DiMuonVeto) && pass.at(TriLeptonVeto) && pass.at(BJetVeto);
  bool passedFullInclusiveNoIsoNoCharge = passedVertex && pass.at(DiMuonVeto) && pass.at(TriLeptonVeto) && pass.at(MT) && pass.at(BJetVeto);

  // define booleans for analysis stages for additional plots
  bool passedObjectsFailDiMuonVeto = passedObjects && !pass.at(DiMuonVeto);

  // run categories
  bool passed_VBFTight		= category_VBFTight(nJets, selJetdeta, selNjetingap, selMjj, higgsPt);
  bool passed_VBFLoose		= category_VBFLoose(nJets, selJetdeta, selNjetingap, selMjj, passed_VBFTight);
  bool passed_VBF = passed_VBFTight || passed_VBFLoose;
  bool passed_OneJetHigh	= category_OneJetHigh(nJets, tauPt, higgsPt, passed_VBF);
  bool passed_OneJetLow		= category_OneJetLow(nJets, tauPt, passed_VBF);
  bool passed_OneJetBoost	= category_OneJetBoost(nJets, tauPt, higgsPt, passed_VBF);
  bool passed_ZeroJetHigh	= category_ZeroJetHigh(nJets, tauPt);
  bool passed_ZeroJetLow	= category_ZeroJetLow(nJets, tauPt);
  bool passed_NoCategory	= category_NoCategory();

  // run relaxed categories for background methods
  bool isWJetMC = (Ntp->GetMCID() >= 20) && (Ntp->GetMCID() <= 23);
  bool useRelaxedForPlots =  (wJetsBGSource == "Data") && isWJetMC; // overwrite pass-vector with relaxed VBF categories (for WJets shape) only if wanted
  bool passed_VBFTightRelaxed = helperCategory_VBFTightRelaxed(useRelaxedForPlots, nJets, selJetdeta, selNjetingap, selMjj, higgsPt);
  bool passed_VBFLooseRelaxed = helperCategory_VBFLooseRelaxed(useRelaxedForPlots, nJets, selJetdeta, selNjetingap, selMjj);

  // fill plot checking if multiple categories have passed
  // this should never happen (except for WJets MC events when running with the "Data" flag
  unsigned nCat = 0;
  if (passedFullInclusiveSel) {CatFired.at(t).Fill(0., w);}
  if (passedFullInclusiveSel && passed_VBFTight	) {nCat++; CatFired.at(t).Fill(1, w);}
  if (passedFullInclusiveSel && passed_VBFLoose	) {nCat++; CatFired.at(t).Fill(2, w);}
  if (passedFullInclusiveSel && passed_OneJetHigh	) {nCat++; CatFired.at(t).Fill(3, w);}
  if (passedFullInclusiveSel && passed_OneJetLow	) {nCat++; CatFired.at(t).Fill(4, w);}
  if (passedFullInclusiveSel && passed_OneJetBoost) {nCat++; CatFired.at(t).Fill(5, w);}
  if (passedFullInclusiveSel && passed_ZeroJetHigh) {nCat++; CatFired.at(t).Fill(6, w);}
  if (passedFullInclusiveSel && passed_ZeroJetLow	) {nCat++; CatFired.at(t).Fill(7, w);}

  NCatFired.at(t).Fill(nCat, w);

  if (passedFullInclusiveSel && nCat == 0){
	  std::cout << "                       Here comes a bad event" << std::endl;
	  const char* format = "%12s : %5s %5s %5s %5s %5s %5s %5s %5s \n";
	  printf(format,"Event", "NJets", "dEta", "CJV", "mjj", "pT(H)", "isVBT", "isVBF", "pT(t)");
	  format = "%12i : %5i %5.2f %5i %5.2f %5.2f %5i %5i %5,2f \n";
	  printf(format,Ntp->EventNumber(), nJets, selJetdeta, selNjetingap, selMjj, higgsPt, passed_VBFTight, passed_VBF, tauPt);
  }

  //if (!(passed_VBFTight || passed_VBFLoose || passed_OneJetHigh|| passed_OneJetLow || passed_OneJetBoost || passed_ZeroJetHigh || passed_ZeroJetLow))
	//	  std::cout << "************* NO CATEGORY PASSED! ****************" << std::endl;

  bool status=AnalysisCuts(t,w,wobs); // true only if full selection passed

  ///////////////////////////////////////////////////////////
  // Add plots
  ///////////////////////////////////////////////////////////

  //////// plots filled before any cuts
  // Vertex plots
  NVtx.at(t).Fill(Ntp->NVtx(),w);
  for(unsigned int i_vtx=0;i_vtx<Ntp->NVtx();i_vtx++){
	VtxZ.at(t).Fill(Ntp->Vtx(i_vtx).z(),w);
	VtxRho.at(t).Fill(sqrt(Ntp->Vtx(i_vtx).x()*Ntp->Vtx(i_vtx).x() + Ntp->Vtx(i_vtx).y()*Ntp->Vtx(i_vtx).y()), w);
	VtxNdof.at(t).Fill(Ntp->Vtx_ndof(i_vtx), w);
	VtxIsfake.at(t).Fill(Ntp->Vtx_isFake(i_vtx), w);
  }
  NGoodVtx.at(t).Fill(nGoodVtx,w);

  //////// plots filled after Vertex selection: Object selection

  if(passedVertex){
	  for(unsigned i_mu=0;i_mu<Ntp->NMuons();i_mu++){
		  if(	Ntp->isTightMuon(i_mu,selVertex) ){
			  MuDxy.at(t).Fill(dxy(Ntp->Muon_p4(i_mu),Ntp->Muon_Poca(i_mu),Ntp->Vtx(selVertex)), w);
			  MuDz.at(t).Fill(dz(Ntp->Muon_p4(i_mu),Ntp->Muon_Poca(i_mu),Ntp->Vtx(selVertex)), w);
			  MuRelIso.at(t).Fill(Ntp->Muon_RelIso(i_mu), w);
		  }
	  }
  }

  //////// plots filled after muon ID selection: Muon Kinematics
  if(passedMuId){
	  for(std::vector<int>::iterator it_mu = selectedMuonsId.begin();it_mu != selectedMuonsId.end(); ++it_mu){
		  MuPt.at(t).Fill(Ntp->Muon_p4(*it_mu).Pt(), w);
		  MuEta.at(t).Fill(Ntp->Muon_p4(*it_mu).Eta(), w);
		  MuPhi.at(t).Fill(Ntp->Muon_p4(*it_mu).Phi(), w);
	  }

	  //////// plots filled only with selected muon
	  if(passedMu){
		  MuSelPt.at(t).Fill(Ntp->Muon_p4(selMuon).Pt(), w);
		  MuSelEta.at(t).Fill(Ntp->Muon_p4(selMuon).Eta(), w);
		  MuSelPhi.at(t).Fill(Ntp->Muon_p4(selMuon).Phi(), w);

		  // Does the muon fake the tau_ID+Iso?
		  bool fakes = false;
		  for( unsigned  i_tau = 0; i_tau < Ntp->NPFTaus(); i_tau++){
			  if (	  selectPFTau_Id(i_tau) &&
					  selectPFTau_Iso(i_tau) &&
					  Ntp->Muon_p4(selMuon).DeltaR(Ntp->PFTau_p4(i_tau)) < cMuTau_dR){
				  fakes = true;
				  break;
			  }
		  }
		  MuSelFakesTauID.at(t).Fill(fakes, w);
	  }
  }

  //////// plots filled after tau ID + Iso selection: Tau Kinematics
  if(passedTauIdIso){
	  for(std::vector<int>::iterator it_tau = selectedTausIso.begin(); it_tau != selectedTausIso.end(); ++it_tau){
		  TauPt.at(t).Fill(Ntp->PFTau_p4(*it_tau).Pt(), w);
		  TauEta.at(t).Fill(Ntp->PFTau_p4(*it_tau).Eta(), w);
		  TauPhi.at(t).Fill(Ntp->PFTau_p4(*it_tau).Phi(), w);
	  }

	  //////// plots filled only with selected tau
	  if(passedTau){
		  TauSelPt.at(t).Fill(Ntp->PFTau_p4(selTau).Pt(), w);
		  TauSelEta.at(t).Fill(Ntp->PFTau_p4(selTau).Eta(), w);
		  TauSelPhi.at(t).Fill(Ntp->PFTau_p4(selTau).Phi(), w);
		  TauSelDecayMode.at(t).Fill(Ntp->PFTau_hpsDecayMode(selTau), w);
	  }
  }

  //////// plots filled after full muon and tau selection
  if(passedObjectsFailDiMuonVeto){
	  // Investigate events discarded by the DiMuon Veto
	  if (Ntp->Muon_Charge(selMuon) == 1){
		  MuVetoDPtSelMuon.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0)).Pt() - Ntp->Muon_p4(selMuon).Pt(), w );
		  MuVetoDRTau.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0)).DeltaR(Ntp->PFTau_p4(selTau)), w);
	  }
	  else if (Ntp->Muon_Charge(selMuon) == -1){
		  MuVetoDPtSelMuon.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)).Pt() - Ntp->Muon_p4(selMuon).Pt(), w );
		  MuVetoDRTau.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)).DeltaR(Ntp->PFTau_p4(selTau)), w);
	  }
	  MuVetoInvM.at(t).Fill( (Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)) + Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0))).M() , w);
	  MuVetoPtPositive.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)).Pt(), w);
	  MuVetoPtNegative.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0)).Pt(), w);
	  MuVetoDeltaR.at(t).Fill( Ntp->Muon_p4(diMuonVetoMuonsPositive.at(0)).DeltaR(Ntp->Muon_p4(diMuonVetoMuonsNegative.at(0))), w );
  }

  if(passedDiMuonVeto){
	  // Tri-lepton vetoes
	  NMuonTriLepVeto.at(t).Fill(triLepVetoMuons.size(), w);
	  NElecTriLepVeto.at(t).Fill(triLepVetoElecs.size(), w);
  }
  //////// plots filled after full selection without BJetVeto
  if (passedFullInclusiveSelNoBVeto){
	  NBJets.at(t).Fill( selectedBJets.size(), w);
	  if (selectedBJets.size() > 0){
		  BJet1Pt.at(t).Fill( Ntp->PFJet_p4(selectedBJets.at(0)).Pt(), w);
		  BJet1Eta.at(t).Fill( Ntp->PFJet_p4(selectedBJets.at(0)).Eta(), w);
		  BJet1Phi.at(t).Fill( Ntp->PFJet_p4(selectedBJets.at(0)).Phi(), w);
	  }
  }

  //////// plots filled after full selection (without categories)
  if (passedFullInclusiveSel){
	  // Mu-Tau correlations
	  MuTauDR    .at(t).Fill( Ntp->Muon_p4(selMuon).DeltaR(Ntp->PFTau_p4(selTau)), w );
	  MuTauDPhi  .at(t).Fill( Ntp->Muon_p4(selMuon).DeltaPhi(Ntp->PFTau_p4(selTau)), w );
	  MuTauDEta  .at(t).Fill( Ntp->Muon_p4(selMuon).Eta() - Ntp->PFTau_p4(selTau).Eta(), w );
	  MuTauDPt   .at(t).Fill( Ntp->Muon_p4(selMuon).Pt() - Ntp->PFTau_p4(selTau).Pt(), w );
	  MuTauRelDPt.at(t).Fill( (Ntp->Muon_p4(selMuon).Pt() - Ntp->PFTau_p4(selTau).Pt()) / Ntp->Muon_p4(selMuon).Pt() , w);

	  // lepton charge
	  MuCharge.at(t).Fill( Ntp->Muon_Charge(selMuon), w);
	  TauCharge.at(t).Fill( Ntp->PFTau_Charge(selTau), w);

	  // MET
	  MetPt.at(t).Fill( Ntp->MET_CorrMVAMuTau_et(), w);
	  MetPhi.at(t).Fill( Ntp->MET_CorrMVAMuTau_phi(), w);

	  // Jets
	  NJetsKin.at(t).Fill( selectedJetsKin.size(), w);
	  if (selectedJetsKin.size() > 0){
		  JetKin1Pt.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(0)).Pt(), w);
		  JetKin1Eta.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(0)).Eta(), w);
		  JetKin1Phi.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(0)).Phi(), w);
		  JetKin1IsLooseId.at(t).Fill( Ntp->PFJet_PUJetID_looseWP(selectedJetsKin.at(0)), w);
	  }
	  if (selectedJetsKin.size() > 1){
		  JetKin2IsLooseId.at(t).Fill( Ntp->PFJet_PUJetID_looseWP(selectedJetsKin.at(1)), w);
		  JetKin2Pt.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(1)).Pt(), w);
		  JetKin2Eta.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(1)).Eta(), w);
		  JetKin2Phi.at(t).Fill( Ntp->PFJet_p4(selectedJetsKin.at(1)).Phi(), w);
	  }
	  NJetsId.at(t).Fill( selectedJets.size(), w);
	  if (selectedJets.size() > 0){
		  Jet1Pt.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(0)).Pt(), w);
		  Jet1Eta.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(0)).Eta(), w);
		  Jet1Phi.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(0)).Phi(), w);
	  	  Jet1IsB.at(t).Fill( Ntp->PFJet_bDiscriminator(selectedJets.at(0)) > cCat_btagDisc , w);
	  }
	  if (selectedJets.size() > 1){
		  Jet2Pt.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(1)).Pt(), w);
		  Jet2Eta.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(1)).Eta(), w);
		  Jet2Phi.at(t).Fill( Ntp->PFJet_p4(selectedJets.at(1)).Phi(), w);
		  Jet2IsB.at(t).Fill( Ntp->PFJet_bDiscriminator(selectedJets.at(1)) > cCat_btagDisc, w);
	  }

	  // variables for categorization
	  HiggsPt.at(t).Fill(higgsPt , w);
	  HiggsPhi.at(t).Fill(higgsPhi , w);
	  JetsDEta.at(t).Fill(selJetdeta , w);
	  JetsInEtaGap.at(t).Fill(selNjetingap , w);
	  JetsInvM.at(t).Fill(selMjj , w);
  }

  //////// plots filled after full selection
  if(status){
    NVtxFullSelection.at(t).Fill(Ntp->NVtx(),w);
    //std::cout << "ID after = " << id << std::endl;
  }

  //////// category specific plots, especially for background methods /////////

  ////// W+Jets Background estimation
  if(passedFullInclusiveSelNoMt && passed_ZeroJetLow){
	  Cat0JetLowMt.at(t).Fill(value.at(MT), w);
	  Cat0JetLowMtSideband.at(t).Fill(value.at(MT), w);
	  if(isWJetMC){
		  if(pass.at(MT))
			  Cat0JetLowMtExtrapolation.at(t).Fill(1, w); // fill with weight
		  if(value.at(MT) > 70.)
			  Cat0JetLowMtExtrapolation.at(t).Fill(2, w);
	  }
  }
  if(passedFullInclusiveSelNoMt && passed_ZeroJetHigh){
	  Cat0JetHighMt.at(t).Fill(value.at(MT), w);
	  Cat0JetHighMtSideband.at(t).Fill(value.at(MT), w);
	  if(isWJetMC){
		  if(pass.at(MT))
			  Cat0JetHighMtExtrapolation.at(t).Fill(1, w); // fill with weight
		  if(value.at(MT) > 70.)
			  Cat0JetHighMtExtrapolation.at(t).Fill(2, w); // fill with weight
	  }
  }
  if(passedFullInclusiveSelNoMt && passed_OneJetLow){
	  Cat1JetLowMt.at(t).Fill(value.at(MT), w);
	  Cat1JetLowMtSideband.at(t).Fill(value.at(MT), w);
	  if(isWJetMC){
		  if(pass.at(MT))
			  Cat1JetLowMtExtrapolation.at(t).Fill(1, w); // fill with weight
		  if(value.at(MT) > 70.)
			  Cat1JetLowMtExtrapolation.at(t).Fill(2, w); // fill with weight
	  }
  }
  if(passedFullInclusiveSelNoMt && passed_OneJetHigh){
	  Cat1JetHighMt.at(t).Fill(value.at(MT), w);
	  Cat1JetHighMtSideband.at(t).Fill(value.at(MT), w);
	  if(isWJetMC){
		  if(pass.at(MT))
			  Cat1JetHighMtExtrapolation.at(t).Fill(1, w); // fill with weight
		  if(value.at(MT) > 70.)
			  Cat1JetHighMtExtrapolation.at(t).Fill(2, w); // fill with weight
	  }
  }
  if(passedFullInclusiveSelNoMt && passed_OneJetBoost){
	  Cat1JetBoostMt.at(t).Fill(value.at(MT), w);
	  Cat1JetBoostMtSideband.at(t).Fill(value.at(MT), w);
	  if(isWJetMC){
		  if(pass.at(MT))
			  Cat1JetBoostMtExtrapolation.at(t).Fill(1, w); // fill with weight
		  if(value.at(MT) > 70.)
			  Cat1JetBoostMtExtrapolation.at(t).Fill(2, w); // fill with weight
	  }
  }
  if(passedFullInclusiveSelNoMt && passed_VBFLoose){
	  CatVBFLooseMt.at(t).Fill(value.at(MT), w);
	  CatVBFLooseMtSideband.at(t).Fill(value.at(MT), w);
  }
  if(passedFullInclusiveSelNoMtNoOS && passed_VBFLooseRelaxed){
	  CatVBFLooseRelaxMt.at(t).Fill(value.at(MT), w);
	  if(isWJetMC){
		  if(pass.at(MT))
			  CatVBFLooseRelaxMtExtrapolation.at(t).Fill(1, w); // fill with weight
		  if(value.at(MT) > 60. && value.at(MT) < 120.)
			  CatVBFLooseRelaxMtExtrapolation.at(t).Fill(2, w); // fill with weight
	  }
  }
  if(passedFullInclusiveSelNoMt && passed_VBFTight){
	  CatVBFTightMt.at(t).Fill(value.at(MT), w);
	  CatVBFTightMtSideband.at(t).Fill(value.at(MT), w);
  }
  if(passedFullInclusiveSelNoMtNoOS && passed_VBFTightRelaxed){
	  CatVBFTightRelaxMt.at(t).Fill(value.at(MT), w);
	  if(isWJetMC){
		  if(pass.at(MT))
			  CatVBFTightRelaxMtExtrapolation.at(t).Fill(1, w); // fill with weight
		  if(value.at(MT) > 60. && value.at(MT) < 120.)
			  CatVBFTightRelaxMtExtrapolation.at(t).Fill(2, w); // fill with weight
	  }
  }
  if(passedFullInclusiveSelNoMt){
	  CatInclusiveMt.at(t).Fill(value.at(MT), w);
	  CatInclusiveMtSideband.at(t).Fill(value.at(MT), w);
	  if(isWJetMC){
		  if(pass.at(MT))
			  CatInclusiveMtExtrapolation.at(t).Fill(1, w); // fill with weight
		  if(value.at(MT) > 70.)
			  CatInclusiveMtExtrapolation.at(t).Fill(2, w); // fill with weight
	  }
  }

  ////// QCD Background estimation
  //     OS/SS
  //       ^
  //    C  |  D
  //   ---------> relIso(mu)
  //    A  |  B
  if (passedFullInclusiveNoIsoNoCharge){
	  bool isSS = ( abs(value.at(OppCharge)) == 2 );
	  // (!passedObjects && hasRelaxedIsoTau) means that there is no signal muon (as relaxedIsoTaus fully contains signal taus)
	  // --> veto events with signal muon AND antiIsoMuon, as in these cases mT etc. are calculated using the signal muon
	  bool isA = pass.at(OppCharge) && passedObjects;
	  bool isB = pass.at(OppCharge) && !passedObjects && hasRelaxedIsoTau && hasAntiIsoMuon;
	  bool isC = isSS && passedObjects;
	  bool isD = isSS && !passedObjects && hasRelaxedIsoTau && hasAntiIsoMuon;
	  if (isA+isB+isC+isD > 1) printf("WARNING: Event %i enters more than 1 ABCD region! (A%i, B%i, C%i D%i)\n", Ntp->EventNumber(), isA, isB, isC, isD);
//	  if (isA+isB+isC+isD == 0) {
//		  printf("ATTENTION: Event %9d enters no ABCD region! Sum(q) = %d, passedMu = %d, passedTau = %d\n", Ntp->EventNumber(), value.at(OppCharge), passedMu, passedTau);
//		  printf("       		                                    hasRelTau = %d, hasAntiIsoMu = %d\n", hasRelaxedIsoTau, hasAntiIsoMuon);
//	  }
	  double iso(-10.), q(-10.);
	  if (isA || isB) q = -0.5;
	  if (isC || isD) q = +0.5;
	  if (isA || isC) iso = -0.5;
	  if (isB || isD) iso = +0.5;

	  CatInclusiveQcdAbcd.at(t).Fill(iso,q,w);
	  if (passed_ZeroJetLow) Cat0JetLowQcdAbcd.at(t).Fill(iso,q,w);
	  if (passed_ZeroJetHigh) Cat0JetHighQcdAbcd.at(t).Fill(iso,q,w);
	  if (passed_OneJetLow) Cat1JetLowQcdAbcd.at(t).Fill(iso,q,w);
	  if (passed_OneJetHigh) Cat1JetHighQcdAbcd.at(t).Fill(iso,q,w);
	  if (passed_OneJetBoost) Cat1JetBoostQcdAbcd.at(t).Fill(iso,q,w);
	  if (passed_VBFLoose) CatVBFLooseQcdAbcd.at(t).Fill(iso,q,w);
	  if (passed_VBFTight) CatVBFTightQcdAbcd.at(t).Fill(iso,q,w);
  }
  // todo: fill the ABCD plot
  // todo: create plots showing isolation for OS and SS events

}


void HToTaumuTauh::Finish() {
	if (verbose)
		std::cout << "HToTaumuTauh::Finish()" << std::endl;

	if (wJetsBGSource == "Data") {
		if (mode == RECONSTRUCT) { // only apply data-driven numbers on "combine" level
			std::cout << "WJet BG: Using data driven yield method." << std::endl;

			double sumSelEvts = 0;
			for (unsigned id = 20; id < 24; id++) {
				if (!HConfig.hasID(id))
					continue;
				int type = HConfig.GetType(id);
				// check that cross-section for WJet processes is set to -1 in Histo.txt
				double oldXSec = HConfig.GetCrossSection(id);
				if (oldXSec != -1) {
					// Histo.txt has WJet xsec unequal -1, so set it to -1 to avoid scaling by framework
					if (!HConfig.SetCrossSection(id, -1))
						std::cout << "WARNING: Could not change cross-section for id " << id << std::endl;
					printf("WJet process %i had xsec = %6.1f. Setting to %6.1f for data-driven WJet yield.\n", id, oldXSec, HConfig.GetCrossSection(id));
				}
				sumSelEvts += Npassed.at(type).GetBinContent(NCuts);
			}

			// second loop, now the total sum of all Wjets events is known, so we can scale
			for (unsigned id = 20; id < 24; id++) {
				if (!HConfig.hasID(id))
					continue;
				int type = HConfig.GetType(id);
				double rawSelEvts = Npassed.at(type).GetBinContent(NCuts);

				// scale all WJet histograms to data-driven yield
				ScaleAllHistOfType(type, wJetsYieldMap[categoryFlag] / sumSelEvts);
				printf("WJet process %i was scaled from yield %f to yield %f \n", id, rawSelEvts, Npassed.at(type).GetBinContent(NCuts));
			}
		}
		else
			std::cout << "WJet BG: Data driven will be used at Combine stage, but not in this individual set." << std::endl;
	}
	else if (wJetsBGSource == "MC")
		std::cout << "WJet BG: Using MC." << std::endl;
	else
		std::cout << "WJet BG: Please specify \"MC\" or \"Data\". Using MC for this run..." << std::endl;
	// call GetHistoInfo here (instead of in Configure function), otherwise the SetCrossSection calls are not reflected
	HConfig.GetHistoInfo(types, CrossSectionandAcceptance, legend, colour);
	Selection::Finish();

	// todo: implement cross check: Integral in WJet plot <-> yield from background method

}


/////////////////////////////////////////
// Definition of selection and helper functions
/////////////////////////////////////////

///////// Helper functions
double HToTaumuTauh::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return (-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt();
}

double HToTaumuTauh::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2);
}

double HToTaumuTauh::matchTrigger(unsigned int i_obj, std::vector<TString> trigger, std::string objectType){
	unsigned int id = 0;
	TLorentzVector particle(0.,0.,0.,0.);
	TLorentzVector triggerObj(0.,0.,0.,0.);
	if(objectType=="tau"){
		id = 84;
		particle = Ntp->PFTau_p4(i_obj);
	}
	if(objectType=="muon"){
		id = 83;
		particle = Ntp->Muon_p4(i_obj);
	}

	double minDR = 100.;
	for(unsigned i_trig = 0; i_trig < trigger.size(); i_trig++){
		for(unsigned i=0;i<Ntp->NHLTTrigger_objs();i++){
			if(Ntp->HLTTrigger_objs_trigger(i).find(trigger.at(i_trig)) != string::npos){
				for(unsigned j=0;j<Ntp->NHLTTrigger_objs(i);j++){
					if(Ntp->HLTTrigger_objs_Id(i,j)==id){
						triggerObj.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
								Ntp->HLTTrigger_objs_Eta(i,j),
								Ntp->HLTTrigger_objs_Phi(i,j),
								Ntp->HLTTrigger_objs_E(i,j));
					}
					if( triggerObj.Pt()>0. && particle.Pt()>0. ) {
						double dr = particle.DeltaR(triggerObj);
						if (dr < minDR) minDR = dr;
					}
				}
			}
		}
	}
	return minDR;
}

///////// Muons

bool HToTaumuTauh::selectMuon_Id(unsigned i, unsigned vertex){
	if(	Ntp->isSelectedMuon(i,vertex,cMu_dxy,cMu_dz) &&
		Ntp->Muon_RelIso(i) < cMu_relIso &&
		matchTrigger(i,cTriggerNames,"muon") < cMu_dRHltMatch
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_Kinematics(unsigned i){
	if(	Ntp->Muon_p4(i).Pt() >= cMu_pt &&
		fabs(Ntp->Muon_p4(i).Eta()) <= cMu_eta
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_antiIso(unsigned i, unsigned vertex) {
	if (Ntp->isSelectedMuon(i, vertex, cMu_dxy, cMu_dz) &&
		Ntp->Muon_RelIso(i) <= 0.5 &&
		Ntp->Muon_RelIso(i) >= 0.2 &&
		matchTrigger(i, cTriggerNames, "muon") < cMu_dRHltMatch &&
		selectMuon_Kinematics(i)) {
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_diMuonVeto(unsigned i, unsigned i_vtx){
	if(	Ntp->Muon_p4(i).Pt() > 15.0 &&
		fabs(Ntp->Muon_p4(i).Eta()) < 2.4 &&
		Ntp->Muon_isPFMuon(i) &&
		Ntp->Muon_isGlobalMuon(i) &&
		Ntp->Muon_isTrackerMuon(i) &&
		Ntp->Muon_RelIso(i) < 0.3 &&
		fabs(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx))) < 0.2
		) {
	  return true;
	}
	return false;
}

bool HToTaumuTauh::selectMuon_triLeptonVeto(unsigned i, int selectedMuon, unsigned i_vtx){
	if(	i != selectedMuon &&
		Ntp->isTightMuon(i,i_vtx) &&
		fabs(dxy(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx))) < cMu_dxy &&
		fabs(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(i_vtx))) < cMu_dz &&
		Ntp->Muon_RelIso(i) < 0.3 &&
		Ntp->Muon_p4(i).Pt() > cMuTriLep_pt &&
		fabs(Ntp->Muon_p4(i).Eta()) < cMuTriLep_eta
			){
			return true;
	}
	return false;
}


///////// Electrons
bool HToTaumuTauh::selectElectron_triLeptonVeto(unsigned i, unsigned i_vtx, std::vector<int> muonCollection){
	// check if elec is matched to a muon, if so this is not a good elec (should become obsolete when using top projections)
//	for(std::vector<int>::iterator it_mu = muonCollection.begin(); it_mu != muonCollection.end(); ++it_mu){
//	  if( Ntp->Electron_p4(i).DeltaR(Ntp->Muon_p4(*it_mu)) < cMuTau_dR ) {
//		  return false;
//	  }
//	}

	if ( 	Ntp->isSelectedElectron(i,i_vtx,0.045,0.2) &&
			//TODO: This electron isolation is using rho corrections, but should use deltaBeta corrections
			//documentation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaEARhoCorrection
			Ntp->Electron_RelIso03(i) < 0.3 && // "For electron the default cone size if 0.3" https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPFBasedIsolation#Alternate_code_to_calculate_PF_I
			Ntp->Electron_p4(i).Pt() > 10.0 &&
			fabs(Ntp->Electron_p4(i).Eta()) < 2.5
			){
		return true;
	}
	return false;
}

///////// Taus

bool HToTaumuTauh::selectPFTau_Id(unsigned i){
	if ( 	Ntp->PFTau_isHPSByDecayModeFinding(i) &&
			Ntp->PFTau_isHPSAgainstElectronsLoose(i) &&
			Ntp->PFTau_isHPSAgainstMuonTight(i)
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFTau_Id(unsigned i, std::vector<int> muonCollection){
	// check if tau is matched to a muon, if so this is not a good tau
	// https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013#Sync_Issues
	for(std::vector<int>::iterator it_mu = muonCollection.begin(); it_mu != muonCollection.end(); ++it_mu){
	  if( Ntp->PFTau_p4(i).DeltaR(Ntp->Muon_p4(*it_mu)) < cMuTau_dR ) {
		  return false;
	  }
	}
	// trigger matching
	if (matchTrigger(i,cTriggerNames,"tau") > cMu_dRHltMatch) {
		return false;
	}

	if ( 	selectPFTau_Id(i) ){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFTau_Iso(unsigned i){
	if ( 	Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(i) < cTau_rawIso
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFTau_Kinematics(unsigned i){
	if ( 	Ntp->PFTau_p4(i).Pt() >= cTau_pt &&
			fabs(Ntp->PFTau_p4(i).Eta()) <= cTau_eta
			){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFTau_relaxedIso(unsigned i, std::vector<int> muonCollection){
	if (selectPFTau_Id(i, muonCollection) &&
		Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(i) < 10. &&
		selectPFTau_Kinematics(i)){
		return true;
	}
	return false;
}

std::vector<int> HToTaumuTauh::sortPFjets(){
	// create vector of pairs to allow for sorting by jet pt
	std::vector< std::pair<int, double> > jetIdxPt;
	for (unsigned i = 0; i<Ntp->NPFJets(); i++ ){
		jetIdxPt.push_back( std::make_pair(i,Ntp->PFJet_p4(i).Pt()) );
	}
	// sort vector of pairs
	std::sort(jetIdxPt.begin(), jetIdxPt.end(), sortIdxByValue());
	// create vector of indices in correct order
	std::vector<int> sortedJets;
	for (unsigned i = 0; i<jetIdxPt.size(); i++){
		sortedJets.push_back(jetIdxPt.at(i).first);
	}
	return sortedJets;
}

bool HToTaumuTauh::selectPFJet_Cleaning(unsigned i, int selectedMuon, int selectedTau){
	// clean against selected muon and tau
	if (selectedMuon >= 0) {
		if (Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(selectedMuon)) < cJetClean_dR) return false;
	}
	if (selectedTau >= 0){
		if (Ntp->PFJet_p4(i).DeltaR(Ntp->PFTau_p4(selectedTau)) < cJetClean_dR) return false;
	}
	return true;
}


bool HToTaumuTauh::selectPFJet_Kinematics(unsigned i){
	if ( 	fabs(Ntp->PFJet_p4(i).Eta()) < cCat_jetEta &&
			Ntp->PFJet_p4(i).Pt() > cCat_jetPt){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectPFJet_Id(unsigned i){
	if (	Ntp->PFJet_PUJetID_looseWP(i) ){
		return true;
	}
	return false;
}

bool HToTaumuTauh::selectBJet(unsigned i, int selectedMuon, int selectedTau){
	// clean against selected muon and tau
	if (selectedMuon >= 0) {
		if (Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(selectedMuon)) < cJetClean_dR) return false;
	}
	if (selectedTau >= 0){
		if (Ntp->PFJet_p4(i).DeltaR(Ntp->PFTau_p4(selectedTau)) < cJetClean_dR) return false;
	}

	if (	fabs(Ntp->PFJet_p4(i).Eta()) < cCat_bjetEta &&
			Ntp->PFJet_p4(i).Pt() > cCat_bjetPt &&
			Ntp->PFJet_bDiscriminator(i) > cCat_btagDisc){
		return true;
	}
	return false;
}


//// *****functions defining the categories*****

void HToTaumuTauh::configure_VBFTight(){
	//todo: move to category classes
	// to be called only if VBFTight is chosen category

	// set cut values to be the cut values of this category
	cut.at(VbfTight_NJet) = cut_VBFTight.at(VbfTight_NJet);
	cut.at(VbfTight_DeltaEta) = cut_VBFTight.at(VbfTight_DeltaEta);
	cut.at(VbfTight_NJetRapGap) = cut_VBFTight.at(VbfTight_NJetRapGap);
	cut.at(VbfTight_JetInvM) = cut_VBFTight.at(VbfTight_JetInvM);
	cut.at(VbfTight_HiggsPt) = cut_VBFTight.at(VbfTight_HiggsPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(VbfTight_NJet)="Number VBF Jets $>=$";
	title.at(VbfTight_NJet)+=cut.at(VbfTight_NJet);
	htitle=title.at(VbfTight_NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jet_{VBF}";
	c="_Cut_";c+=VbfTight_NJet;
	Nminus1.at(VbfTight_NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfTight_NJet_",htitle,6,-0.5,5.5,hlabel,"Events");
	Nminus0.at(VbfTight_NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfTight_NJet_",htitle,6,-0.5,5.5,hlabel,"Events");

	title.at(VbfTight_DeltaEta)="$\\Delta\\eta(jj) >$";
	title.at(VbfTight_DeltaEta)+=cut.at(VbfTight_DeltaEta);
	htitle=title.at(VbfTight_DeltaEta);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="#Delta#eta(Jet_{VBF}^{1},Jet_{VBF}^{2})";
	c="_Cut_";c+=VbfTight_DeltaEta;
	Nminus1.at(VbfTight_DeltaEta) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfTight_DeltaEta_",htitle,32,-8.,8.,hlabel,"Events");
	Nminus0.at(VbfTight_DeltaEta) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfTight_DeltaEta_",htitle,32,-8.,8.,hlabel,"Events");

	title.at(VbfTight_NJetRapGap)="Number Jets in $\\eta$ gap $<=$";
	title.at(VbfTight_NJetRapGap)+=cut.at(VbfTight_NJetRapGap);
	htitle=title.at(VbfTight_NJetRapGap);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jet_{VBF} in rapidity gap";
	c="_Cut_";c+=VbfTight_NJetRapGap;
	Nminus1.at(VbfTight_NJetRapGap) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfTight_NJetRapGap_",htitle,6,-0.5,5.5,hlabel,"Events");
	Nminus0.at(VbfTight_NJetRapGap) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfTight_NJetRapGap_",htitle,6,-0.5,5.5,hlabel,"Events");

	title.at(VbfTight_JetInvM)="$m_{jj}(VBF) >$";
	title.at(VbfTight_JetInvM)+=cut.at(VbfTight_JetInvM);
	title.at(VbfTight_JetInvM)+=" GeV";
	htitle=title.at(VbfTight_JetInvM);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="m_{inv}(jj) of VBF-jets";
	c="_Cut_";c+=VbfTight_JetInvM;
	Nminus1.at(VbfTight_JetInvM) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfTight_JetInvM_",htitle,25,0.,2000.,hlabel,"Events");
	Nminus0.at(VbfTight_JetInvM) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfTight_JetInvM_",htitle,25,0.,2000.,hlabel,"Events");

	title.at(VbfTight_HiggsPt)="$p_{T}(H) >$";
	title.at(VbfTight_HiggsPt)+=cut.at(VbfTight_HiggsPt);
	title.at(VbfTight_HiggsPt)+=" GeV";
	htitle=title.at(VbfTight_HiggsPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T} of Higgs candidate";
	c="_Cut_";c+=VbfTight_HiggsPt;
	Nminus1.at(VbfTight_HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfTight_HiggsPt_",htitle,25,0.,250.,hlabel,"Events");
	Nminus0.at(VbfTight_HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfTight_HiggsPt_",htitle,25,0.,250.,hlabel,"Events");
}
bool HToTaumuTauh::category_VBFTight(unsigned NJets, double DEta, int NJetsInGap, double Mjj, double higgsPt){
	std::vector<float> value_VBFTight;
	std::vector<float> pass_VBFTight;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_VBFTight.push_back(-10.);
	pass_VBFTight.push_back(false);
	}

	value_VBFTight.at(VbfTight_NJet) = NJets;
	pass_VBFTight.at(VbfTight_NJet) = (value_VBFTight.at(VbfTight_NJet) >= cut_VBFTight.at(VbfTight_NJet));

	if(pass_VBFTight.at(VbfTight_NJet)){
		value_VBFTight.at(VbfTight_DeltaEta) = DEta;
		pass_VBFTight.at(VbfTight_DeltaEta) = (fabs(value_VBFTight.at(VbfTight_DeltaEta)) > cut_VBFTight.at(VbfTight_DeltaEta));

		value_VBFTight.at(VbfTight_NJetRapGap) = NJetsInGap;
		pass_VBFTight.at(VbfTight_NJetRapGap) = (value_VBFTight.at(VbfTight_NJetRapGap) <= cut_VBFTight.at(VbfTight_NJetRapGap));

		value_VBFTight.at(VbfTight_JetInvM) = Mjj;
		pass_VBFTight.at(VbfTight_JetInvM) = (value_VBFTight.at(VbfTight_JetInvM) > cut_VBFTight.at(VbfTight_JetInvM));
	}
	else{
		pass_VBFTight.at(VbfTight_DeltaEta) = true;
		pass_VBFTight.at(VbfTight_NJetRapGap) = true;
		pass_VBFTight.at(VbfTight_JetInvM) = true;
	}

	value_VBFTight.at(VbfTight_HiggsPt) = higgsPt;
	pass_VBFTight.at(VbfTight_HiggsPt) = (value_VBFTight.at(VbfTight_HiggsPt) > cut_VBFTight.at(VbfTight_HiggsPt));

	// migrate into main analysis if this is chosen category
	return migrateCategoryIntoMain("VBFTight",value_VBFTight, pass_VBFTight,VbfTight_NCuts);
}

void HToTaumuTauh::configure_VBFLoose(){
	// to be called only if VBFLoose is chosen category

	// set cut values to be the cut values of this category
	cut.at(VbfLoose_NJet) = cut_VBFLoose.at(VbfLoose_NJet);
	cut.at(VbfLoose_DeltaEta) = cut_VBFLoose.at(VbfLoose_DeltaEta);
	cut.at(VbfLoose_NJetRapGap) = cut_VBFLoose.at(VbfLoose_NJetRapGap);
	cut.at(VbfLoose_JetInvM) = cut_VBFLoose.at(VbfLoose_JetInvM);
	cut.at(VbfLoose_NotVbfTight) = cut_VBFLoose.at(VbfLoose_NotVbfTight);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(VbfLoose_NJet)="Number VBF Jets $>=$";
	title.at(VbfLoose_NJet)+=cut.at(VbfLoose_NJet);
	htitle=title.at(VbfLoose_NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jet_{VBF}";
	c="_Cut_";c+=VbfLoose_NJet;
	Nminus1.at(VbfLoose_NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfLoose_NJet_",htitle,6,-0.5,5.5,hlabel,"Events");
	Nminus0.at(VbfLoose_NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfLoose_NJet_",htitle,6,-0.5,5.5,hlabel,"Events");

	title.at(VbfLoose_DeltaEta)="$\\Delta\\eta(jj) >$";
	title.at(VbfLoose_DeltaEta)+=cut.at(VbfLoose_DeltaEta);
	htitle=title.at(VbfLoose_DeltaEta);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="#Delta#eta(Jet_{VBF}^{1},Jet_{VBF}^{2})";
	c="_Cut_";c+=VbfLoose_DeltaEta;
	Nminus1.at(VbfLoose_DeltaEta) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfLoose_DeltaEta_",htitle,32,-8.,8.,hlabel,"Events");
	Nminus0.at(VbfLoose_DeltaEta) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfLoose_DeltaEta_",htitle,32,-8.,8.,hlabel,"Events");

	title.at(VbfLoose_NJetRapGap)="Number Jets in $\\eta$ gap $<=$";
	title.at(VbfLoose_NJetRapGap)+=cut.at(VbfLoose_NJetRapGap);
	htitle=title.at(VbfLoose_NJetRapGap);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jet_{VBF} in rapidity gap";
	c="_Cut_";c+=VbfLoose_NJetRapGap;
	Nminus1.at(VbfLoose_NJetRapGap) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfLoose_NJetRapGap_",htitle,6,-0.5,5.5,hlabel,"Events");
	Nminus0.at(VbfLoose_NJetRapGap) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfLoose_NJetRapGap_",htitle,6,-0.5,5.5,hlabel,"Events");

	title.at(VbfLoose_JetInvM)="$m_{jj}(VBF) >$";
	title.at(VbfLoose_JetInvM)+=cut.at(VbfLoose_JetInvM);
	title.at(VbfLoose_JetInvM)+=" GeV";
	htitle=title.at(VbfLoose_JetInvM);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="m_{inv}(jj) of VBF-jets";
	c="_Cut_";c+=VbfLoose_JetInvM;
	Nminus1.at(VbfLoose_JetInvM) = HConfig.GetTH1D(Name+c+"_Nminus1_VbfLoose_JetInvM_",htitle,25,0.,2000.,hlabel,"Events");
	Nminus0.at(VbfLoose_JetInvM) = HConfig.GetTH1D(Name+c+"_Nminus0_VbfLoose_JetInvM_",htitle,25,0.,2000.,hlabel,"Events");

	title.at(VbfLoose_NotVbfTight)="Not VBFTight $==$";
	title.at(VbfLoose_NotVbfTight)+=cut.at(VbfLoose_NotVbfTight);
	htitle=title.at(VbfLoose_NotVbfTight);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Did not pass VBFTight cat.";
	c="_Cut_";c+=VbfLoose_NotVbfTight;
	Nminus1.at(VbfLoose_NotVbfTight) = HConfig.GetTH1D(Name+c+"_Nminus1_NotVbfTight_",htitle,2,-0.5,1.5,hlabel,"Events");
	Nminus0.at(VbfLoose_NotVbfTight) = HConfig.GetTH1D(Name+c+"_Nminus0_NotVbfTight_",htitle,2,-0.5,1.5,hlabel,"Events");
}
bool HToTaumuTauh::category_VBFLoose(unsigned NJets, double DEta, int NJetsInGap, double Mjj, bool passedVBFTight){
	std::vector<float> value_VBFLoose;
	std::vector<float> pass_VBFLoose;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_VBFLoose.push_back(-10.);
	pass_VBFLoose.push_back(false);
	}

	value_VBFLoose.at(VbfLoose_NJet) = NJets;
	pass_VBFLoose.at(VbfLoose_NJet) = (value_VBFLoose.at(VbfLoose_NJet) >= cut_VBFLoose.at(VbfLoose_NJet));

	if(pass_VBFLoose.at(VbfLoose_NJet)){
		value_VBFLoose.at(VbfLoose_DeltaEta) = DEta;
		pass_VBFLoose.at(VbfLoose_DeltaEta) = (fabs(value_VBFLoose.at(VbfLoose_DeltaEta)) > cut_VBFLoose.at(VbfLoose_DeltaEta));

		value_VBFLoose.at(VbfLoose_NJetRapGap) = NJetsInGap;
		pass_VBFLoose.at(VbfLoose_NJetRapGap) = (value_VBFLoose.at(VbfLoose_NJetRapGap) <= cut_VBFLoose.at(VbfLoose_NJetRapGap));

		value_VBFLoose.at(VbfLoose_JetInvM) = Mjj;
		pass_VBFLoose.at(VbfLoose_JetInvM) = (value_VBFLoose.at(VbfLoose_JetInvM) > cut_VBFLoose.at(VbfLoose_JetInvM));
	}
	else{
		pass_VBFLoose.at(VbfLoose_DeltaEta) = true;
		pass_VBFLoose.at(VbfLoose_NJetRapGap) = true;
		pass_VBFLoose.at(VbfLoose_JetInvM) = true;
	}

	value_VBFLoose.at(VbfLoose_NotVbfTight) = !passedVBFTight;
	pass_VBFLoose.at(VbfLoose_NotVbfTight) = ( value_VBFLoose.at(VbfLoose_NotVbfTight) == cut_VBFLoose.at(VbfLoose_NotVbfTight) );

	// migrate into main analysis if this is chosen category
	return migrateCategoryIntoMain("VBFLoose",value_VBFLoose, pass_VBFLoose,VbfLoose_NCuts);
}

void HToTaumuTauh::configure_OneJetLow(){
	// to be called only if OneJetLow is chosen category

	// set cut values to be the cut values of this category
	cut.at(OneJetLow_NJet) = cut_OneJetLow.at(OneJetLow_NJet);
	cut.at(OneJetLow_NotVbf) = cut_OneJetLow.at(OneJetLow_NotVbf);
	cut.at(OneJetLow_TauPt) = cut_OneJetLow.at(OneJetLow_TauPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(OneJetLow_NJet)="Number Jets $>=$";
	title.at(OneJetLow_NJet)+=cut.at(OneJetLow_NJet);
	htitle=title.at(OneJetLow_NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=OneJetLow_NJet;
	Nminus1.at(OneJetLow_NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetLow_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(OneJetLow_NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetLow_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(OneJetLow_NotVbf)="Not VBFTight or VBFLoose $==$";
	title.at(OneJetLow_NotVbf)+=cut.at(OneJetLow_NotVbf);
	htitle=title.at(OneJetLow_NotVbf);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Did not pass VBF cat.";
	c="_Cut_";c+=OneJetLow_NotVbf;
	Nminus1.at(OneJetLow_NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetLow_NotVBF_",htitle,2,-0.5,1.5,hlabel,"Events");
	Nminus0.at(OneJetLow_NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetLow_NotVBF_",htitle,2,-0.5,1.5,hlabel,"Events");

	title.at(OneJetLow_TauPt)="$p_{T}(\\tau_{h}) <$";
	title.at(OneJetLow_TauPt)+=cut.at(OneJetLow_TauPt);
	title.at(OneJetLow_TauPt)+=" GeV";
	htitle=title.at(OneJetLow_TauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=OneJetLow_TauPt;
	Nminus1.at(OneJetLow_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetLow_TauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(OneJetLow_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetLow_TauPt_",htitle,50,0.,200.,hlabel,"Events");
}
bool HToTaumuTauh::category_OneJetLow(unsigned NJets, double TauPt, bool passedVBF){
	bool categoryPass = true;
	std::vector<float> value_OneJetLow;
	std::vector<float> pass_OneJetLow;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_OneJetLow.push_back(-10.);
	pass_OneJetLow.push_back(false);
	}

	value_OneJetLow.at(OneJetLow_NJet) = NJets;
	pass_OneJetLow.at(OneJetLow_NJet) = ( value_OneJetLow.at(OneJetLow_NJet) >= cut_OneJetLow.at(OneJetLow_NJet) );

	value_OneJetLow.at(OneJetLow_NotVbf) = !passedVBF;
	pass_OneJetLow.at(OneJetLow_NotVbf) = ( value_OneJetLow.at(OneJetLow_NotVbf) == cut_OneJetLow.at(OneJetLow_NotVbf) );

	if (selTau == -1){
		// TauPt cut is set to true for nice N-0 and N-1 plots
		value_OneJetLow.at(OneJetLow_TauPt) = -10.;
		pass_OneJetLow.at(OneJetLow_TauPt) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		value_OneJetLow.at(OneJetLow_TauPt) = TauPt;
		pass_OneJetLow.at(OneJetLow_TauPt) = ( value_OneJetLow.at(OneJetLow_TauPt) < cut_OneJetLow.at(OneJetLow_TauPt) );
	}


	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("OneJetLow",value_OneJetLow, pass_OneJetLow,OneJetLow_NCuts) && categoryPass;
	return categoryPass;
}

void HToTaumuTauh::configure_OneJetHigh(){
	// to be called only if OneJetHigh is chosen category

	// set cut values to be the cut values of this category
	cut.at(OneJetHigh_NJet) = cut_OneJetHigh.at(OneJetHigh_NJet);
	cut.at(OneJetHigh_NotVbf) = cut_OneJetHigh.at(OneJetHigh_NotVbf);
	cut.at(OneJetHigh_TauPt) = cut_OneJetHigh.at(OneJetHigh_TauPt);
	cut.at(OneJetHigh_HiggsPt) = cut_OneJetHigh.at(OneJetHigh_HiggsPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(OneJetHigh_NJet)="Number Jets $>=$";
	title.at(OneJetHigh_NJet)+=cut.at(OneJetHigh_NJet);
	htitle=title.at(OneJetHigh_NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=OneJetHigh_NJet;
	Nminus1.at(OneJetHigh_NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetHigh_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(OneJetHigh_NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetHigh_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(OneJetHigh_NotVbf)="No VBF $==$";
	title.at(OneJetHigh_NotVbf)+=cut.at(OneJetHigh_NotVbf);
	htitle=title.at(OneJetHigh_NotVbf);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Did not pass VBF cat.";
	c="_Cut_";c+=OneJetHigh_NotVbf;
	Nminus1.at(OneJetHigh_NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetHigh_NotVbf_",htitle,2,-0.5,1.5,hlabel,"Events");
	Nminus0.at(OneJetHigh_NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetHigh_NotVbf_",htitle,2,-0.5,1.5,hlabel,"Events");

	title.at(OneJetHigh_TauPt)="$p_{T}(\\tau_{h}) >=$";
	title.at(OneJetHigh_TauPt)+=cut.at(OneJetHigh_TauPt);
	title.at(OneJetHigh_TauPt)+=" GeV";
	htitle=title.at(OneJetHigh_TauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=OneJetHigh_TauPt;
	Nminus1.at(OneJetHigh_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetHigh_TauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(OneJetHigh_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetHigh_TauPt_",htitle,50,0.,200.,hlabel,"Events");

	title.at(OneJetHigh_HiggsPt)="$p_{T}(H) <$";
	title.at(OneJetHigh_HiggsPt)+=cut.at(OneJetHigh_HiggsPt);
	title.at(OneJetHigh_HiggsPt)+=" GeV";
	htitle=title.at(OneJetHigh_HiggsPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T} of Higgs candidate";
	c="_Cut_";c+=OneJetHigh_HiggsPt;
	Nminus1.at(OneJetHigh_HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetHigh_HiggsPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(OneJetHigh_HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetHigh_HiggsPt_",htitle,50,0.,200.,hlabel,"Events");
}
bool HToTaumuTauh::category_OneJetHigh(unsigned NJets, double TauPt, double higgsPt, bool passedVBF){
	bool categoryPass = true;
	std::vector<float> value_OneJetHigh;
	std::vector<float> pass_OneJetHigh;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_OneJetHigh.push_back(-10.);
	pass_OneJetHigh.push_back(false);
	}

	value_OneJetHigh.at(OneJetHigh_NJet) = NJets;
	pass_OneJetHigh.at(OneJetHigh_NJet) = ( value_OneJetHigh.at(OneJetHigh_NJet) >= cut_OneJetHigh.at(OneJetHigh_NJet) );

	value_OneJetHigh.at(OneJetHigh_NotVbf) = !passedVBF;
	pass_OneJetHigh.at(OneJetHigh_NotVbf) = ( value_OneJetHigh.at(OneJetHigh_NotVbf) == cut_OneJetHigh.at(OneJetHigh_NotVbf) );

	if (selTau == -1){
		// TauPt cut is set to true for nice N-0 and N-1 plots
		value_OneJetHigh.at(OneJetHigh_TauPt) = -10.;
		pass_OneJetHigh.at(OneJetHigh_TauPt) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		value_OneJetHigh.at(OneJetHigh_TauPt) = TauPt;
		pass_OneJetHigh.at(OneJetHigh_TauPt) = ( value_OneJetHigh.at(OneJetHigh_TauPt) >= cut_OneJetHigh.at(OneJetHigh_TauPt) );
	}

	value_OneJetHigh.at(OneJetHigh_HiggsPt) = higgsPt;
	pass_OneJetHigh.at(OneJetHigh_HiggsPt) = (value_OneJetHigh.at(OneJetHigh_HiggsPt) < cut_OneJetHigh.at(OneJetHigh_HiggsPt));

	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("OneJetHigh",value_OneJetHigh, pass_OneJetHigh,OneJetHigh_NCuts) && categoryPass;
	return categoryPass;
}

void HToTaumuTauh::configure_OneJetBoost(){
	// to be called only if OneJetBoost is chosen category

	// set cut values to be the cut values of this category
	cut.at(OneJetBoost_NJet) = cut_OneJetBoost.at(OneJetBoost_NJet);
	cut.at(OneJetBoost_NotVbf) = cut_OneJetBoost.at(OneJetBoost_NotVbf);
	cut.at(OneJetBoost_TauPt) = cut_OneJetBoost.at(OneJetBoost_TauPt);
	cut.at(OneJetBoost_HiggsPt) = cut_OneJetBoost.at(OneJetBoost_HiggsPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(OneJetBoost_NJet)="Number Jets $>=$";
	title.at(OneJetBoost_NJet)+=cut.at(OneJetBoost_NJet);
	htitle=title.at(OneJetBoost_NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=OneJetBoost_NJet;
	Nminus1.at(OneJetBoost_NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetBoost_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(OneJetBoost_NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetBoost_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(OneJetBoost_NotVbf)="No VBF $==$";
	title.at(OneJetBoost_NotVbf)+=cut.at(OneJetBoost_NotVbf);
	htitle=title.at(OneJetBoost_NotVbf);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Did not pass VBF cat.";
	c="_Cut_";c+=OneJetBoost_NotVbf;
	Nminus1.at(OneJetBoost_NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetBoost_NotVbf_",htitle,2,-0.5,1.5,hlabel,"Events");
	Nminus0.at(OneJetBoost_NotVbf) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetBoost_NotVbf_",htitle,2,-0.5,1.5,hlabel,"Events");

	title.at(OneJetBoost_TauPt)="$p_{T}(\\tau_{h}) >=$";
	title.at(OneJetBoost_TauPt)+=cut.at(OneJetBoost_TauPt);
	title.at(OneJetBoost_TauPt)+=" GeV";
	htitle=title.at(OneJetBoost_TauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=OneJetBoost_TauPt;
	Nminus1.at(OneJetBoost_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetBoost_TauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(OneJetBoost_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetBoost_TauPt_",htitle,50,0.,200.,hlabel,"Events");

	title.at(OneJetBoost_HiggsPt)="$p_{T}(H) >=$";
	title.at(OneJetBoost_HiggsPt)+=cut.at(OneJetBoost_HiggsPt);
	title.at(OneJetBoost_HiggsPt)+=" GeV";
	htitle=title.at(OneJetBoost_HiggsPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T} of Higgs candidate";
	c="_Cut_";c+=OneJetBoost_HiggsPt;
	Nminus1.at(OneJetBoost_HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus1_OneJetBoost_HiggsPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(OneJetBoost_HiggsPt) = HConfig.GetTH1D(Name+c+"_Nminus0_OneJetBoost_HiggsPt_",htitle,50,0.,200.,hlabel,"Events");
}
bool HToTaumuTauh::category_OneJetBoost(unsigned NJets, double TauPt, double higgsPt, bool passedVBF){
	bool categoryPass = true;
	std::vector<float> value_OneJetBoost;
	std::vector<float> pass_OneJetBoost;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_OneJetBoost.push_back(-10.);
	pass_OneJetBoost.push_back(false);
	}

	value_OneJetBoost.at(OneJetBoost_NJet) = NJets;
	pass_OneJetBoost.at(OneJetBoost_NJet) = ( value_OneJetBoost.at(OneJetBoost_NJet) >= cut_OneJetBoost.at(OneJetBoost_NJet) );

	value_OneJetBoost.at(OneJetBoost_NotVbf) = !passedVBF;
	pass_OneJetBoost.at(OneJetBoost_NotVbf) = ( value_OneJetBoost.at(OneJetBoost_NotVbf) == cut_OneJetBoost.at(OneJetBoost_NotVbf) );

	if (selTau == -1){
		// TauPt cut is set to true for nice N-0 and N-1 plots
		value_OneJetBoost.at(OneJetBoost_TauPt) = -10.;
		pass_OneJetBoost.at(OneJetBoost_TauPt) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		value_OneJetBoost.at(OneJetBoost_TauPt) = TauPt;
		pass_OneJetBoost.at(OneJetBoost_TauPt) = ( value_OneJetBoost.at(OneJetBoost_TauPt) >= cut_OneJetBoost.at(OneJetBoost_TauPt) );
	}

	value_OneJetBoost.at(OneJetBoost_HiggsPt) = higgsPt;
	pass_OneJetBoost.at(OneJetBoost_HiggsPt) = (value_OneJetBoost.at(OneJetBoost_HiggsPt) >= cut_OneJetBoost.at(OneJetBoost_HiggsPt));

	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("OneJetBoost",value_OneJetBoost, pass_OneJetBoost,OneJetBoost_NCuts) && categoryPass;
	return categoryPass;
}

void HToTaumuTauh::configure_ZeroJetHigh(){
	// to be called only if ZeroJetHigh is chosen category

	// set cut values to be the cut values of this category
	cut.at(ZeroJetHigh_NJet) = cut_ZeroJetHigh.at(ZeroJetHigh_NJet);
	cut.at(ZeroJetHigh_TauPt) = cut_ZeroJetHigh.at(ZeroJetHigh_TauPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(ZeroJetHigh_NJet)="Number Jets $<=$";
	title.at(ZeroJetHigh_NJet)+=cut.at(ZeroJetHigh_NJet);
	htitle=title.at(ZeroJetHigh_NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=ZeroJetHigh_NJet;
	Nminus1.at(ZeroJetHigh_NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetHigh_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(ZeroJetHigh_NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetHigh_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(ZeroJetHigh_TauPt)="$p_{T}(\\tau_{h}) >=$";
	title.at(ZeroJetHigh_TauPt)+=cut.at(ZeroJetHigh_TauPt);
	title.at(ZeroJetHigh_TauPt)+=" GeV";
	htitle=title.at(ZeroJetHigh_TauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=ZeroJetHigh_TauPt;
	Nminus1.at(ZeroJetHigh_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetHigh_TauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(ZeroJetHigh_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetHigh_TauPt_",htitle,50,0.,200.,hlabel,"Events");
}
bool HToTaumuTauh::category_ZeroJetHigh(unsigned NJets, double TauPt){
	bool categoryPass = true;
	std::vector<float> value_ZeroJetHigh;
	std::vector<float> pass_ZeroJetHigh;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_ZeroJetHigh.push_back(-11.);
	pass_ZeroJetHigh.push_back(false);
	}

	value_ZeroJetHigh.at(ZeroJetHigh_NJet) = NJets;
	pass_ZeroJetHigh.at(ZeroJetHigh_NJet) = ( value_ZeroJetHigh.at(ZeroJetHigh_NJet) <= cut_ZeroJetHigh.at(ZeroJetHigh_NJet) );

	if (selTau == -1){
		// TauPt cut is set to true for nice N-0 and N-1 plots
		value_ZeroJetHigh.at(ZeroJetHigh_TauPt) = -10.;
		pass_ZeroJetHigh.at(ZeroJetHigh_TauPt) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		value_ZeroJetHigh.at(ZeroJetHigh_TauPt) = TauPt;
		pass_ZeroJetHigh.at(ZeroJetHigh_TauPt) = ( value_ZeroJetHigh.at(ZeroJetHigh_TauPt) >= cut_ZeroJetHigh.at(ZeroJetHigh_TauPt) );
	}


	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("ZeroJetHigh",value_ZeroJetHigh, pass_ZeroJetHigh,ZeroJetHigh_NCuts) && categoryPass;
	return categoryPass;
}

void HToTaumuTauh::configure_ZeroJetLow(){
	// to be called only if ZeroJetLow is chosen category

	// set cut values to be the cut values of this category
	cut.at(ZeroJetLow_NJet) = cut_ZeroJetLow.at(ZeroJetLow_NJet);
	cut.at(ZeroJetLow_TauPt) = cut_ZeroJetLow.at(ZeroJetLow_TauPt);

	// set histograms of category cuts
	TString hlabel;
	TString htitle;
	TString c;

	title.at(ZeroJetLow_NJet)="Number Jets $<=$";
	title.at(ZeroJetLow_NJet)+=cut.at(ZeroJetLow_NJet);
	htitle=title.at(ZeroJetLow_NJet);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="Number of Jets";
	c="_Cut_";c+=ZeroJetLow_NJet;
	Nminus1.at(ZeroJetLow_NJet) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetLow_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");
	Nminus0.at(ZeroJetLow_NJet) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetLow_NJet_",htitle,11,-0.5,10.5,hlabel,"Events");

	title.at(ZeroJetLow_TauPt)="$p_{T}(\\tau_{h}) <$";
	title.at(ZeroJetLow_TauPt)+=cut.at(ZeroJetLow_TauPt);
	title.at(ZeroJetLow_TauPt)+=" GeV";
	htitle=title.at(ZeroJetLow_TauPt);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="p_{T}(\\tau_{h})/GeV";
	c="_Cut_";c+=ZeroJetLow_TauPt;
	Nminus1.at(ZeroJetLow_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus1_ZeroJetLow_TauPt_",htitle,50,0.,200.,hlabel,"Events");
	Nminus0.at(ZeroJetLow_TauPt) = HConfig.GetTH1D(Name+c+"_Nminus0_ZeroJetLow_TauPt_",htitle,50,0.,200.,hlabel,"Events");
}

bool HToTaumuTauh::category_ZeroJetLow(unsigned NJets, double TauPt) {
	bool categoryPass = true;
	std::vector<float> value_ZeroJetLow(NCuts,-10);
	std::vector<float> pass_ZeroJetLow(NCuts,false);


	value_ZeroJetLow.at(ZeroJetLow_NJet) = NJets;
	pass_ZeroJetLow.at(ZeroJetLow_NJet) = ( value_ZeroJetLow.at(ZeroJetLow_NJet) <= cut_ZeroJetLow.at(ZeroJetLow_NJet) );

	if (selTau == -1){
		// TauPt cut is set to true for nice N-0 and N-1 plots
		value_ZeroJetLow.at(ZeroJetLow_TauPt) = -10.;
		pass_ZeroJetLow.at(ZeroJetLow_TauPt) = true;
		// whole category is failing selection, to avoid NCat > 1
		categoryPass = false;
	}
	else{
		value_ZeroJetLow.at(ZeroJetLow_TauPt) = TauPt;
		pass_ZeroJetLow.at(ZeroJetLow_TauPt) = ( value_ZeroJetLow.at(ZeroJetLow_TauPt) < cut_ZeroJetLow.at(ZeroJetLow_TauPt) );
	}


	// migrate into main analysis if this is chosen category
	categoryPass = migrateCategoryIntoMain("ZeroJetLow",value_ZeroJetLow, pass_ZeroJetLow,ZeroJetLow_NCuts) && categoryPass;
	return categoryPass;
}

void HToTaumuTauh::configure_NoCategory(){
	// to be called only if No Category shall be run

	// set cut values to be the cut values of this category
	// nothing to do

	// set histograms of category cuts
	// no histograms to be set
}
bool HToTaumuTauh::category_NoCategory(){
	std::vector<float> value_NoCategory;
	std::vector<float> pass_NoCategory;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_NoCategory.push_back(-10.);
	pass_NoCategory.push_back(false);
	}

	// no cuts to compute

	// migrate into main analysis if this is chosen category
	TString cat = (categoryFlag == "Inclusive") ? "Inclusive" : "NoCategory"; // make sure that Inclusive category is handled as NoCategory
	return migrateCategoryIntoMain(cat,value_NoCategory, pass_NoCategory,NoCategory_NCuts);
}

// migrate a category into main analysis if this is chosen category
// return value: if category passed
bool HToTaumuTauh::migrateCategoryIntoMain(TString thisCategory, std::vector<float> categoryValueVector, std::vector<float> categoryPassVector, int categoryNCuts) {
	bool catPassed = true;
	for (unsigned i_cut = CatCut1; i_cut < NCuts; i_cut++) {

		// migrate only if this category is the chosen one
		if (categoryFlag == thisCategory) {
			if (i_cut < categoryNCuts) {
				value.at(i_cut) = categoryValueVector.at(i_cut);
				pass.at(i_cut) = categoryPassVector.at(i_cut);
			} else {
				// cut not implemented in this category -> set to true
				value.at(i_cut) = -10.;
				pass.at(i_cut) = true;
			}
		}
		if (i_cut < categoryNCuts) {
			catPassed = catPassed && categoryPassVector.at(i_cut);
		}
	}

	return catPassed;
}

// helper category definitions for background methods
bool HToTaumuTauh::helperCategory_VBFLooseRelaxed(bool useRelaxedForPlots, unsigned NJets, double DEta, int NJetsInGap, double Mjj){
	std::vector<float> value_VBFLooseRelaxed;
	std::vector<float> pass_VBFLooseRelaxed;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_VBFLooseRelaxed.push_back(-10.);
	pass_VBFLooseRelaxed.push_back(false);
	}

	value_VBFLooseRelaxed.at(VbfLoose_NJet) = NJets;
	pass_VBFLooseRelaxed.at(VbfLoose_NJet) = (value_VBFLooseRelaxed.at(VbfLoose_NJet) >= cut_VBFLooseRelaxed.at(VbfLoose_NJet));

	if(pass_VBFLooseRelaxed.at(VbfLoose_NJet)){
		value_VBFLooseRelaxed.at(VbfLoose_DeltaEta) = DEta;
		pass_VBFLooseRelaxed.at(VbfLoose_DeltaEta) = (fabs(value_VBFLooseRelaxed.at(VbfLoose_DeltaEta)) > cut_VBFLooseRelaxed.at(VbfLoose_DeltaEta));

		value_VBFLooseRelaxed.at(VbfLoose_NJetRapGap) = NJetsInGap;
		pass_VBFLooseRelaxed.at(VbfLoose_NJetRapGap) = (value_VBFLooseRelaxed.at(VbfLoose_NJetRapGap) <= cut_VBFLooseRelaxed.at(VbfLoose_NJetRapGap));

		value_VBFLooseRelaxed.at(VbfLoose_JetInvM) = Mjj;
		pass_VBFLooseRelaxed.at(VbfLoose_JetInvM) = (value_VBFLooseRelaxed.at(VbfLoose_JetInvM) > cut_VBFLooseRelaxed.at(VbfLoose_JetInvM));
	}
	else{
		pass_VBFLooseRelaxed.at(VbfLoose_DeltaEta) = true;
		pass_VBFLooseRelaxed.at(VbfLoose_NJetRapGap) = true;
		pass_VBFLooseRelaxed.at(VbfLoose_JetInvM) = true;
	}

	value_VBFLooseRelaxed.at(VbfLoose_NotVbfTight) = true;
	pass_VBFLooseRelaxed.at(VbfLoose_NotVbfTight) = true; // disabled cut

	// migrate into main analysis if this is chosen category
	TString cat = useRelaxedForPlots ? "Loose" : "DoNotUseThisCategoryForPlotting";
	return migrateCategoryIntoMain(cat,value_VBFLooseRelaxed, pass_VBFLooseRelaxed,VbfLoose_NCuts);
}
bool HToTaumuTauh::helperCategory_VBFTightRelaxed(bool useRelaxedForPlots, unsigned NJets, double DEta, int NJetsInGap, double Mjj, double higgsPt){
	std::vector<float> value_VBFTightRelaxed;
	std::vector<float> pass_VBFTightRelaxed;

	// cut implementation
	for(int i=0; i<NCuts;i++){
	value_VBFTightRelaxed.push_back(-10.);
	pass_VBFTightRelaxed.push_back(false);
	}

	value_VBFTightRelaxed.at(VbfTight_NJet) = NJets;
	pass_VBFTightRelaxed.at(VbfTight_NJet) = (value_VBFTightRelaxed.at(VbfTight_NJet) >= cut_VBFTightRelaxed.at(VbfTight_NJet));

	if(pass_VBFTightRelaxed.at(VbfTight_NJet)){
		value_VBFTightRelaxed.at(VbfTight_DeltaEta) = DEta;
		pass_VBFTightRelaxed.at(VbfTight_DeltaEta) = (fabs(value_VBFTightRelaxed.at(VbfTight_DeltaEta)) > cut_VBFTightRelaxed.at(VbfTight_DeltaEta));

		value_VBFTightRelaxed.at(VbfTight_NJetRapGap) = NJetsInGap;
		pass_VBFTightRelaxed.at(VbfTight_NJetRapGap) = (value_VBFTightRelaxed.at(VbfTight_NJetRapGap) <= cut_VBFTightRelaxed.at(VbfTight_NJetRapGap));

		value_VBFTightRelaxed.at(VbfTight_JetInvM) = Mjj;
		pass_VBFTightRelaxed.at(VbfTight_JetInvM) = (value_VBFTightRelaxed.at(VbfTight_JetInvM) > cut_VBFTightRelaxed.at(VbfTight_JetInvM));
	}
	else{
		pass_VBFTightRelaxed.at(VbfTight_DeltaEta) = true;
		pass_VBFTightRelaxed.at(VbfTight_NJetRapGap) = true;
		pass_VBFTightRelaxed.at(VbfTight_JetInvM) = true;
	}

	value_VBFTightRelaxed.at(VbfTight_HiggsPt) = higgsPt;
	pass_VBFTightRelaxed.at(VbfTight_HiggsPt) = (value_VBFTightRelaxed.at(VbfTight_HiggsPt) > cut_VBFTightRelaxed.at(VbfTight_HiggsPt));

	// migrate into main analysis if this is chosen category
	TString cat = useRelaxedForPlots ? "VBFTight" : "DoNotUseThisCategoryForPlotting";
	return migrateCategoryIntoMain(cat,value_VBFTightRelaxed, pass_VBFTightRelaxed,VbfTight_NCuts);
}
