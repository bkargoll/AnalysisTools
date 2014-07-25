#include "TriggerStudyMC.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"

TriggerStudyMC::TriggerStudyMC(TString Name_, TString id_):
  Selection(Name_,id_)
{
  doKFTau=false;
  doHPSLooseTau=false;
  doHPSMediumTau=false;
  doHPSTightTau=false;
  verbose=true;
  if(Name_.Contains("KFTau"))        doKFTau=true;
  if(Name_.Contains("HPSLooseTau"))  doHPSLooseTau=true;
  if(Name_.Contains("HPSMediumTau")) doHPSMediumTau=true;
  if(Name_.Contains("HPSTightTau"))  doHPSTightTau=true;
  // doKFTau=true;
  doHPSMediumTau=true;
}

TriggerStudyMC::~TriggerStudyMC(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "TriggerStudyMC::~TriggerStudyMC Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "TriggerStudyMC::~TriggerStudyMC()" << std::endl;
}

void  TriggerStudyMC::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==isSignal)          cut.at(isSignal)=1;
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    Accumdist.push_back(std::vector<TH1D>());
    Nminus1dist.push_back(std::vector<TH1D>());

    TString c="_Cut_";
    if(i<10)c+="0";
    c+=i;
  
    if(i==isSignal){
      title.at(i)="isSignal";
      htitle=title.at(i);
      hlabel="isSignal";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_isSignal_",htitle,31,-0.5,30.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_isSignal_",htitle,31,-0.5,30.5,hlabel,"Events"));
    }
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  for(unsigned int h=0;h<4;h++){
    TString NP;
    if(h==0)NP="All";
    if(h==1)NP="1prong1prong";
    if(h==2)NP="1prong3prong";
    if(h==3)NP="3prong3prong";
    N_pt.push_back(HConfig.GetTH1D(Name+NP+"_N_pt","N_pt",40,0,200,"P_{T} (GeV)","Events"));
    n_pt.push_back(HConfig.GetTH1D(Name+NP+"_n_pt","n_pt",40,0,200,"P_{T} (GeV)","Events"));
    Eff_pt.push_back(HConfig.GetTH1D(Name+NP+"_Eff_pt","Eff_pt",40,0,200,"P_{T} (GeV)","Efficiency"));
    N_eta.push_back(HConfig.GetTH1D(Name+NP+"_N_eta","N_eta",40,-5.0,5.0,"#eta","Events"));
    n_eta.push_back(HConfig.GetTH1D(Name+NP+"_n_eta","n_eta",40,-5.0,5.0,"#eta","Events"));
    Eff_eta.push_back(HConfig.GetTH1D(Name+NP+"_Eff_eta","Eff_eta",40,-5.0,5.0,"#eta","Efficiency"));
    N_abs.push_back(HConfig.GetTH1D(Name+NP+"_N_abs","N_abs",4,-0.5,3.5,"0=Our 1=1prong 2=our had 3=1prong had","Events"));
    n_abs.push_back(HConfig.GetTH1D(Name+NP+"_n_abs","n_abs",4,-0.5,3.5,"0=Our 1=1prong 2=our had 3=1prong had","Events"));
    Eff_abs.push_back(HConfig.GetTH1D(Name+NP+"_Eff_abs","Eff_abs",4,-0.5,3.5,"0=Our 1=1prong 2=our had 3=1prong had","Efficiency"));

    N_pt_vis.push_back(HConfig.GetTH1D(Name+NP+"_N_pt_vis","N_pt_vis",40,0,200,"P_{T} (GeV)","Events"));
    n_pt_vis.push_back(HConfig.GetTH1D(Name+NP+"_n_pt_vis","n_pt_vis",40,0,200,"P_{T} (GeV)","Events"));
    Eff_pt_vis.push_back(HConfig.GetTH1D(Name+NP+"_Eff_pt_vis","Eff_pt_vis",40,0,200,"P_{T} (GeV)","Efficiency"));
    N_eta_vis.push_back(HConfig.GetTH1D(Name+NP+"_N_eta_vis","N_eta_vis",40,-5.0,5.0,"#eta","Events"));
    n_eta_vis.push_back(HConfig.GetTH1D(Name+NP+"_n_eta_vis","n_eta_vis",40,-5.0,5.0,"#eta","Events"));
    Eff_eta_vis.push_back(HConfig.GetTH1D(Name+NP+"_Eff_eta_vis","Eff_eta_vis",40,-5.0,5.0,"#eta","Efficiency"));

    N_pt_visreco.push_back(HConfig.GetTH1D(Name+NP+"_N_pt_visreco","N_pt_visreco",40,0,200,"P_{T} (GeV)","Events"));
    n_pt_visreco.push_back(HConfig.GetTH1D(Name+NP+"_n_pt_visreco","n_pt_visreco",40,0,200,"P_{T} (GeV)","Events"));
    Eff_pt_visreco.push_back(HConfig.GetTH1D(Name+NP+"_Eff_pt_visreco","Eff_pt_visreco",40,0,200,"P_{T} (GeV)","Efficiency"));
    N_eta_visreco.push_back(HConfig.GetTH1D(Name+NP+"_N_eta_visreco","N_eta_visreco",40,-5.0,5.0,"#eta","Events"));
    n_eta_visreco.push_back(HConfig.GetTH1D(Name+NP+"_n_eta_visreco","n_eta_visreco",40,-5.0,5.0,"#eta","Events"));
    Eff_eta_visreco.push_back(HConfig.GetTH1D(Name+NP+"_Eff_eta_visreco","Eff_eta_visreco",40,-5.0,5.0,"#eta","Efficiency"));

    N_2d_vis.push_back(HConfig.GetTH2D(Name+NP+"_N_2d_vis","N_2d_vis",20,0,200,20,-5,5,"P_{T} (GeV)","#eta"));
    n_2d_vis.push_back(HConfig.GetTH2D(Name+NP+"_n_2d_vis","n_2d_vis",20,0,200,20,-5,5,"P_{T} (GeV)","#eta"));
    Eff_2d_vis.push_back(HConfig.GetTH2D(Name+NP+"_Eff_2d_vis","Eff_2d_vis",20,0,200,20,-5,5,"P_{T} (GeV)","#eta"));
    n_2d_vis_Triggershift5GeV.push_back(HConfig.GetTH2D(Name+NP+"_n_2d_vis_Triggershift5GeV","Eff_2d_vis_Triggershift5GeV",20,0,200,20,-5,5,"P_{T} (GeV)","#eta"));
    n_2d_vis_Triggershift10GeV.push_back(HConfig.GetTH2D(Name+NP+"_n_2d_vis_Triggershift10GeV","Eff_2d_vis_Triggershift10GeV",20,0,200,20,-5,5,"P_{T} (GeV)","#eta"));
  }

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  TriggerStudyMC::Store_ExtraDist(){
  for(unsigned int h=0;h<Eff_pt.size();h++){
    Extradist1d.push_back(&N_pt.at(h));
    Extradist1d.push_back(&N_eta.at(h));
    Extradist1d.push_back(&n_pt.at(h));
    Extradist1d.push_back(&n_eta.at(h));
    Extradist1d.push_back(&Eff_pt.at(h));
    Extradist1d.push_back(&Eff_eta.at(h));
    Extradist1d.push_back(&N_abs.at(h));
    Extradist1d.push_back(&n_abs.at(h));
    Extradist1d.push_back(&Eff_abs.at(h));

    Extradist1d.push_back(&N_pt_vis.at(h));
    Extradist1d.push_back(&N_eta_vis.at(h));
    Extradist1d.push_back(&n_pt_vis.at(h));
    Extradist1d.push_back(&n_eta_vis.at(h));
    Extradist1d.push_back(&Eff_pt_vis.at(h));
    Extradist1d.push_back(&Eff_eta_vis.at(h));

    Extradist1d.push_back(&N_pt_visreco.at(h));
    Extradist1d.push_back(&N_eta_visreco.at(h));
    Extradist1d.push_back(&n_pt_visreco.at(h));
    Extradist1d.push_back(&n_eta_visreco.at(h));
    Extradist1d.push_back(&Eff_pt_visreco.at(h));
    Extradist1d.push_back(&Eff_eta_visreco.at(h));

    Extradist2d.push_back(&N_2d_vis.at(h));
    Extradist2d.push_back(&n_2d_vis.at(h));
    Extradist2d.push_back(&Eff_2d_vis.at(h));

    Extradist2d.push_back(&n_2d_vis_Triggershift5GeV.at(h));
    Extradist2d.push_back(&n_2d_vis_Triggershift10GeV.at(h));

  }
}

void  TriggerStudyMC::doEvent(){
  unsigned int t;
  int id(2);
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << id <<std::endl; return;}

  if(verbose)std::cout << "void  TriggerStudyMC::doEvent() A" << std::endl;

  // cout << "ntaus" <<Ntp->NMCTaus() << endl;
  value.at(isSignal)=0;
  if(Ntp->NMCTaus()==2){
    //cout<< Ntp->MCTau_JAK(0) <<  " " << Ntp->MCTau_JAK(1) << endl;
    if(Ntp->MCTau_JAK(0)!=TauDecay::JAK_ELECTRON && Ntp->MCTau_JAK(1)!=TauDecay::JAK_ELECTRON &&
       Ntp->MCTau_JAK(0)!=TauDecay::JAK_MUON && Ntp->MCTau_JAK(1)!=TauDecay::JAK_MUON){
      pass.at(isSignal)=true;
      value.at(isSignal)=1;
    }
  }
  /*  if(pass.at(isSignal)){
    for(unsigned p=0;p<Ntp->NHLTTriggers();p++){
      cout << Ntp->HTLTriggerName(p) <<  " Prescale: " << Ntp->HLTPrescale(p) << " " << Ntp->L1SEEDPrescale(p) <<endl;
    }
    int ntau=0;
    for(unsigned int i=0;i<Ntp->NKFTau() && doKFTau;i++){
      if(Ntp->isGoodKFTau(i)) ntau++;
    }
    for(unsigned int i=0;i<Ntp->NPFTaus() && doHPSLooseTau;i++){
      if(Ntp->PFTau_isLooseIsolation(i))ntau++;
    }
    for(unsigned int i=0;i<Ntp->NPFTaus() && doHPSMediumTau;i++){
      if(Ntp->PFTau_isMediumIsolation(i))ntau++;
    }
    for(unsigned int i=0;i<Ntp->NPFTaus() && doHPSTightTau;i++){
      if(Ntp->PFTau_isTightIsolation(i))ntau++;
    }
    if(ntau>=2 || (!doKFTau && !doHPSLooseTau && !doHPSMediumTau && !doHPSTightTau)){pass.at(isSignal)=true;value.at(isSignal)+=4;}
    else{pass.at(isSignal)=false;value.at(isSignal)+=2;}
    }*/
  //////////////////////////////////////////////////////////
  double wobs(1),w(1);
  if(!Ntp->isData()){
    if(verbose)std::cout << "void  TriggerStudyMC::doEvent() J" << std::endl;
    //w*=Ntp->PUWeightFineBins();
    if(verbose)std::cout << "void  TriggerStudyMC::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  if(verbose)std::cout << "void  TriggerStud::doEvent() L" << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  TauDecay TD;
  int idx=0;
  unsigned int taubit1(Ntp->MCTau_DecayBitMask(0)), taubit2(Ntp->MCTau_DecayBitMask(1));
  if( TD.nProng(taubit1)==1 && TD.nProng(taubit2)==1){
    idx=1;
    }
  else if(TD.nProng(taubit1)==1 && TD.nProng(taubit2)==3 ||
	  TD.nProng(taubit1)==3 && TD.nProng(taubit2)==1){
    idx=2;
  }
  else if( TD.nProng(taubit1)==3 && TD.nProng(taubit2)==3){
    idx=3;
  }
  if(TD.nProng(taubit1)==1 && TD.nProng(taubit2)==3 && Ntp->MCTau_JAK(1)==5 ||
     TD.nProng(taubit1)==3 && TD.nProng(taubit2)==1 && Ntp->MCTau_JAK(0)==5){
    idx=2;
  }

  bool flag1prong=Ntp->TriggerAccept("HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1_v2");
  bool flag=Ntp->TriggerAccept("HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_v2");
  N_abs.at(0).at(t).Fill(0.,w);
  N_abs.at(idx).at(t).Fill(0.,w);
  N_abs.at(0).at(t).Fill(1.,w);
  N_abs.at(idx).at(t).Fill(1.,w);
  if(flag){
    n_abs.at(0).at(t).Fill(0.,w);
    n_abs.at(idx).at(t).Fill(0.,w);
  }
  if(flag1prong){
    n_abs.at(0).at(t).Fill(1.,w);
    n_abs.at(idx).at(t).Fill(1.,w);
  }
  if(status){
    N_abs.at(0).at(t).Fill(2.,w);
    N_abs.at(idx).at(t).Fill(2.,w);
    N_abs.at(0).at(t).Fill(3.,w);
    N_abs.at(idx).at(t).Fill(3.,w);
    if(flag){
      n_abs.at(0).at(t).Fill(2.,w);
      n_abs.at(idx).at(t).Fill(2.,w);
    }
    if(flag1prong){
      n_abs.at(0).at(t).Fill(3.,w);
      n_abs.at(idx).at(t).Fill(3.,w);
    }

    for(unsigned int tau=0;tau<Ntp->NMCTaus();tau++){
      double tau_pt=Ntp->MCTau_p4(tau).Pt();
      double tau_eta=Ntp->MCTau_p4(tau).Eta();
      double tau_ptvis=0;
      double tau_etavis=0;
      TLorentzVector tautrig_p4(0,0,0,0);
      for(unsigned int p=0; p<Ntp->NMCTauDecayProducts(tau);p++){
	if(abs(Ntp->MCTauandProd_pdgid(tau,p))==abs(PDGInfo::pi_plus) || abs(Ntp->MCTauandProd_pdgid(tau,p))==abs(PDGInfo::pi0) || abs(Ntp->MCTauandProd_pdgid(tau,p))==abs(PDGInfo::K_plus)){
	  //cout << tau << " " << p << " PID "<< Ntp->MCTauandProd_pdgid(tau,p) << endl;
	  tautrig_p4+=Ntp->MCTauandProd_p4(tau,p);
	}
      }
      tau_ptvis=tautrig_p4.Pt();
      tau_etavis=tautrig_p4.Eta();
      ////////////////////////////////////////////////////

      ////////////////////////////////////////////////////
      N_pt.at(0).at(t).Fill(tau_pt,w);
      if(flag)n_pt.at(0).at(t).Fill(tau_pt,w);
      N_eta.at(0).at(t).Fill(tau_eta,w);
      if(flag)n_eta.at(0).at(t).Fill(tau_eta,w);
      if(idx!=0){
	N_pt.at(idx).at(t).Fill(tau_pt,w);
	if(flag)n_pt.at(idx).at(t).Fill(tau_pt,w);
	N_eta.at(idx).at(t).Fill(tau_eta,w);
	if(flag)n_eta.at(idx).at(t).Fill(tau_eta,w);
      }
      ///////////////////////////////////////////////////
      N_pt_vis.at(0).at(t).Fill(tau_ptvis,w);
      if(flag)n_pt_vis.at(0).at(t).Fill(tau_ptvis,w);
      N_eta_vis.at(0).at(t).Fill(tau_etavis,w);
      if(flag)n_eta_vis.at(0).at(t).Fill(tau_etavis,w);
      if(idx!=0){
        N_pt_vis.at(idx).at(t).Fill(tau_ptvis,w);
        if(flag)n_pt_vis.at(idx).at(t).Fill(tau_ptvis,w);
        N_eta_vis.at(idx).at(t).Fill(tau_etavis,w);
        if(flag)n_eta_vis.at(idx).at(t).Fill(tau_etavis,w);
      }


      N_2d_vis.at(0).at(t).Fill(tau_ptvis,tau_etavis,w);
      if(flag)n_2d_vis.at(0).at(t).Fill(tau_ptvis,tau_etavis,w);
      if(idx!=0){
        N_2d_vis.at(idx).at(t).Fill(tau_ptvis,tau_etavis,w);
        if(flag)n_2d_vis.at(idx).at(t).Fill(tau_ptvis,tau_etavis,w);
      }

      ///////////////////////////////////////////////////
      unsigned int trig=0;
      float tau_ptvisreco(0),tau_etavisreco(0);
      if(Ntp->GetTriggerIndex("HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_v2",trig) && tau==0){
	for(unsigned int jidx=0; jidx<Ntp->NHLTTriggerObject(trig);jidx++){
	  tau_ptvisreco=Ntp->HLTTriggerObject_p4(trig,jidx).Pt();
	  tau_etavisreco=Ntp->HLTTriggerObject_p4(trig,jidx).Eta();
	  N_pt_visreco.at(0).at(t).Fill(tau_ptvisreco,w);
	  if(flag)n_pt_visreco.at(0).at(t).Fill(tau_ptvisreco,w);
	  N_eta_visreco.at(0).at(t).Fill(tau_etavisreco,w);
	  if(flag)n_eta_visreco.at(0).at(t).Fill(tau_etavisreco,w);
	  if(idx!=0){
	    N_pt_visreco.at(idx).at(t).Fill(tau_ptvisreco,w);
	    if(flag)n_pt_visreco.at(idx).at(t).Fill(tau_ptvisreco,w);
	    N_eta_visreco.at(idx).at(t).Fill(tau_etavisreco,w);
	    if(flag)n_eta_visreco.at(idx).at(t).Fill(tau_etavisreco,w);
	  }
	}
      }
      ///////////////////////////////////////////////////
    }
  }
}


void TriggerStudyMC::Finish(){
  for(unsigned int i=0;i<N_pt.size();i++){
    for(unsigned int j=0;j<N_pt.at(i).size();j++){
      // n_pt.at(i).at(j).Rebin(5);N_pt.at(i).at(j).Rebin(5);
      // n_eta.at(i).at(j).Rebin(5);N_eta.at(i).at(j).Rebin(5);
      Eff_pt.at(i).at(j).Reset();//Eff_pt.at(i).at(j).Rebin(5);
      Eff_pt.at(i).at(j).Divide(&n_pt.at(i).at(j),&N_pt.at(i).at(j),1.0,1.0,"B");
      Eff_eta.at(i).at(j).Reset(); //Eff_eta.at(i).at(j).Rebin(5);
      Eff_eta.at(i).at(j).Divide(&n_eta.at(i).at(j),&N_eta.at(i).at(j),1.0,1.0,"B");
      Eff_abs.at(i).at(j).Reset(); //Eff_abs.at(i).at(j).Rebin(5);
      Eff_abs.at(i).at(j).Divide(&n_abs.at(i).at(j),&N_abs.at(i).at(j),1.0,1.0,"B");

      Eff_pt_vis.at(i).at(j).Reset();
      Eff_pt_vis.at(i).at(j).Divide(&n_pt_vis.at(i).at(j),&N_pt_vis.at(i).at(j),1.0,1.0,"B");
      Eff_eta_vis.at(i).at(j).Reset(); 
      Eff_eta_vis.at(i).at(j).Divide(&n_eta_vis.at(i).at(j),&N_eta_vis.at(i).at(j),1.0,1.0,"B");

      Eff_pt_visreco.at(i).at(j).Reset();
      Eff_pt_visreco.at(i).at(j).Divide(&n_pt_visreco.at(i).at(j),&N_pt_visreco.at(i).at(j),1.0,1.0,"B");
      Eff_eta_visreco.at(i).at(j).Reset();
      Eff_eta_visreco.at(i).at(j).Divide(&n_eta_visreco.at(i).at(j),&N_eta_visreco.at(i).at(j),1.0,1.0,"B");

      Eff_2d_vis.at(i).at(j).Reset();
      Eff_2d_vis.at(i).at(j).Divide(&n_2d_vis.at(i).at(j),&N_2d_vis.at(i).at(j),1.0,1.0,"B");

      n_2d_vis_Triggershift5GeV.at(i).at(j).Reset();
      n_2d_vis_Triggershift10GeV.at(i).at(j).Reset();
      for(unsigned int l=0; l<n_2d_vis_Triggershift5GeV.at(i).at(j).GetNbinsX();l++){
	for(unsigned int k=0; k<n_2d_vis_Triggershift5GeV.at(i).at(j).GetNbinsY();k++){
	  double value5GeV=N_2d_vis.at(i).at(j).GetBinContent(l,k);
	  value5GeV*=Eff_2d_vis.at(i).at(j).GetBinContent(l+1,k);
	  double value10GeV=N_2d_vis.at(i).at(j).GetBinContent(l,k);
          value10GeV*=Eff_2d_vis.at(i).at(j).GetBinContent(l+2,k);
	  if(k<=8 || k>=13){ value5GeV=0;  value10GeV=0;}
	  n_2d_vis_Triggershift5GeV.at(i).at(j).SetBinContent(l,k,value5GeV);
          n_2d_vis_Triggershift10GeV.at(i).at(j).SetBinContent(l,k,value10GeV);
	}
      }
    }
  }
  Selection::Finish();
}

