#include "TauLifeTime.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauSolver.h"
//#include "TCanvas.h"
//#include "NTuple_Controller.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "math.h"

TauLifeTime::TauLifeTime(TString Name_, TString id_):
  Selection(Name_,id_)
  ,muoneta(2.1)
  ,jeteta(2.3)
  ,TauTrackPtThreshold(5.0)
{
  verbose=true;//false;
}

TauLifeTime::~TauLifeTime(){
  for(int j=0; j<Npassed.size(); j++){
    //std::cout << "TauLifeTime::~TauLifeTime Selection Summary before: "
	// << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 //<< Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  //std::cout << "TauLifeTime::~TauLifeTime()" << std::endl;
}

void  TauLifeTime::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==hasTag)	      cut.at(hasTag)=1;
    if(i==TagPtMin)           cut.at(TagPtMin)=20;
    if(i==TagIso)             cut.at(TagIso)=0.2;
//    if(i==numIsoTags)         cut.at(numIsoTags)=1;
    if(i==TauPt)              cut.at(TauPt)=20;
    if(i==TauEta)             cut.at(TauEta)=2.0;
    if(i==TauIsIsolated)      cut.at(TauIsIsolated)=true;
    if(i==TauFit)             cut.at(TauFit)=1;
    if(i==deltaPhi)           cut.at(deltaPhi)=TMath::Pi()*3.0/4.0;
    if(i==Charge)             cut.at(Charge)=0.0;
    //if(i==ZMassmax)           cut.at(ZMassmax)=80;
  //std::cout << "Setting cut no. i=" << i << std::endl;
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
  
    if(i==TriggerOk){
      title.at(i)="$N_{\\mu} passing the Trigger$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu} passing the Trigger";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }        
    else if(i==hasTag){
      title.at(i)="$Number of good muons$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of good muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TagPtMin){
      title.at(i)="$N_{\\mu} with P_{T}>$";
      title.at(i)+=cut.at(TagPtMin);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPtMin_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPtMin_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TagIso){
      title.at(i)="$N_{\\mu} with rel. isolation <=$";
      title.at(i)+=cut.at(TagIso);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#mu}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagIso_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauPt){
      title.at(i)="$N_{\\tau} with P_{T} >$";
      title.at(i)+=cut.at(TauPt);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauEta){
      title.at(i)="$N_{\\tau} | \\eta <$";
      title.at(i)+=cut.at(TauEta);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==TauIsIsolated){
      title.at(i)="$N_{\\tau} | isolated$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauIsIsolated_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauIsIsolated_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
else if(i==TauFit){
      title.at(i)="$N_{\\tau} | fitted good$";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="N_{#tau}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauFit_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauFit_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==deltaPhi){
      title.at(i)="$\\Delta \\phi distribution >$";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(deltaPhi));
      title.at(i)+=buffer;
      title.at(i)+="(rad)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#Delta #phi (rad)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
    }
    else if(i==Charge){
      title.at(i)="$|Q_{\\tau}+Q_{\\mu}|-0.5<$";
      title.at(i)+=cut.at(Charge);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="|Q_{#tau}+Q_{#mu}|";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Charge_",htitle,44,0,2.2,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Charge_",htitle,44,0,2.2,hlabel,"Events"));
    }
    //-----------
   /*else if(i==ZMassmax){
     title.at(i)="$M_{Z}<$";
     title.at(i)+=cut.at(ZMassmax);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{Z} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
   }*/
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");


  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  TauFlightLength=HConfig.GetTH1D(Name+"_TauFlightLength","TauFlightLength",62,-0.55,2.55,"L (cm)","Events");
  TauFlightLengthTransverse=HConfig.GetTH1D(Name+"_TauFlightLengthTransverse","TauFlightLengthTransverse",62,-0.55,2.55,"L_{T} (cm)","Events");
  TauMomentum=HConfig.GetTH1D(Name+"_TauMomentum","TauMomentum",50,0,400,"|P| (GeV)","Events");
  TauMomentumTransverse=HConfig.GetTH1D(Name+"_TauMomentumTransverse","TauMomentum",50,0,400,"|P_{T}| (GeV)","Events");
  TauLife=HConfig.GetTH1D(Name+"_TauLife","TauLife",1000,-0.5,2,"#tau_{#tau} (ps)","Events");
  TauLifeTransverse=HConfig.GetTH1D(Name+"_TauLifeTransverse","TauLifeTransverse",1000,-0.5,2,"transverse #tau_{#tau} (ps)","Events");
  ResTauFlightLength=HConfig.GetTH1D(Name+"_ResTauFlightLength","ResTauFlightLength",100,-2,2,"L (cm)","Events");
  ResTauFlightLengthTransverse=HConfig.GetTH1D(Name+"_ResTauFlightLengthTransverse","ResTauFlightLengthTransverse",100,-2,2,"L_{T} (cm)","Events");
  ResTauMomentum=HConfig.GetTH1D(Name+"_ResTauMomentum","ResTauMomentum",80,-200,200,"|P| (GeV)","Events");
  ResTauMomentumTransverse=HConfig.GetTH1D(Name+"_ResTauMomentumTransverse","ResTauMomentum",80,-200,200,"|P_{T}| (GeV)","Events");
  ResTauLife=HConfig.GetTH1D(Name+"_ResTauLife","ResTauLife",1000,-1,1,"#tau_{#tau} (ps)","Events");
  ResTauLifeTransverse=HConfig.GetTH1D(Name+"_ResTauLifeTransverse","ResTauLifeTransverse",1000,-1,1,"transverse #tau_{#tau} (ps)","Events");

  //
  TauMomentumAvg=HConfig.GetTH1D(Name+"_TauMomentumAvg","TauMomentumAvg",50,0,400,"|P| (GeV)","Events");
  TauMomentumTransverseAvg=HConfig.GetTH1D(Name+"_TauMomentumTransverseAvg","TauMomentumAvg",50,0,400,"|P_{T}| (GeV)","Events");
  TauLifeAvg=HConfig.GetTH1D(Name+"_TauLifeAvg","TauLifeAvg",1000,-0.5,2,"#tau_{#tau} (ps)","Events");
  TauLifeTransverseAvg=HConfig.GetTH1D(Name+"_TauLifeTransverseAvg","TauLifeTransverseAvg",1000,-0.5,2,"transverse #tau_{#tau} (ps)","Events");
  ResTauMomentumAvg=HConfig.GetTH1D(Name+"_ResTauMomentumAvg","ResTauMomentumAvg",80,-200,200,"|P| (GeV)","Events");
  ResTauMomentumTransverseAvg=HConfig.GetTH1D(Name+"_ResTauMomentumTransverseAvg","ResTauMomentumAvg",80,-200,200,"|P_{T}| (GeV)","Events");
  ResTauLifeAvg=HConfig.GetTH1D(Name+"_ResTauLifeAvg","ResTauLifeAvg",1000,-1,1,"#tau_{#tau} (ps)","Events");
  ResTauLifeTransverseAvg=HConfig.GetTH1D(Name+"_ResTauLifeTransverseAvg","ResTauLifeTransverseAvg",1000,-1,1,"transverse #tau_{#tau} (ps)","Events");



  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    //std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}



void  TauLifeTime::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&TauFlightLength);
 Extradist1d.push_back(&TauFlightLengthTransverse);
 Extradist1d.push_back(&TauMomentum);
 Extradist1d.push_back(&TauMomentumTransverse);
 Extradist1d.push_back(&TauLife);
 Extradist1d.push_back(&TauLifeTransverse);
 Extradist1d.push_back(&ResTauFlightLength);
 Extradist1d.push_back(&ResTauFlightLengthTransverse);
 Extradist1d.push_back(&ResTauMomentum);
 Extradist1d.push_back(&ResTauMomentumTransverse);
 Extradist1d.push_back(&ResTauLife);
 Extradist1d.push_back(&ResTauLifeTransverse);

 Extradist1d.push_back(&TauMomentumAvg);
 Extradist1d.push_back(&TauMomentumTransverseAvg);
 Extradist1d.push_back(&ResTauMomentumAvg);
 Extradist1d.push_back(&ResTauMomentumTransverseAvg);
 Extradist1d.push_back(&TauLifeAvg);
 Extradist1d.push_back(&TauLifeTransverseAvg);
 Extradist1d.push_back(&ResTauLifeAvg);
 Extradist1d.push_back(&ResTauLifeTransverseAvg);
}

void  TauLifeTime::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  std::cout << "ID: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

  if(verbose)std::cout << "void  TauLifeTime::doEvent() Mu" << std::endl;

  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("IsoMu24"))value.at(TriggerOk)+=1;
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);

  // Apply Selection
  std::vector<unsigned int> mu_idx_good, mu_idx_pt, mu_idx_iso;
  unsigned mu_idx(999);
  //double mu_pt(0);
  for(unsigned int i=0;i<Ntp->NMuons();i++){
    //std::cout << "Checking if muon number = " << i << " is good..." << std::endl;
    if(Ntp->isGoodMuon(i) && fabs(Ntp->Muon_p4(i).Eta())<muoneta){
      mu_idx_good.push_back(i);
    }  
  }
  for(unsigned int i=0;i<mu_idx_good.size();i++){	  
    if(Ntp->Muon_p4(mu_idx_good.at(i)).Pt() > cut.at(TagPtMin)){
      mu_idx_pt.push_back(mu_idx_good.at(i));
    }
  }	
  for(unsigned int i=0;i<mu_idx_pt.size();i++){
    if((Ntp->Muon_emEt05(mu_idx_pt.at(i)) + Ntp->Muon_hadEt05(mu_idx_pt.at(i)) + Ntp->Muon_sumPt05(mu_idx_pt.at(i)))/Ntp->Muon_p4(mu_idx_pt.at(i)).Pt()<=cut.at(TagIso)){
      mu_idx=mu_idx_pt.at(i);//mu_idx is the idx of the effectively selected muon.
      mu_idx_iso.push_back(mu_idx);
    }
  }
  if(verbose)std::cout << "void  TauLifeTime::doEvent() MuA" << std::endl;
  //std::cout << "nmus = " << mu_idx_iso.size() << std::endl;
  value.at(hasTag)=mu_idx_good.size();
  pass.at(hasTag)=(value.at(hasTag)>=cut.at(hasTag));
  value.at(TagPtMin)=mu_idx_pt.size();
  pass.at(TagPtMin)=(value.at(TagPtMin)>0);
  value.at(TagIso)=mu_idx_iso.size();
  pass.at(TagIso)=(value.at(TagIso)==1);


  if(verbose)std::cout << "void  TauLifeTime::doEvent() MuEnd - mu_idx = " << mu_idx << " NMuons = " << Ntp->NMuons() << std::endl;
  if(verbose)std::cout << "void  TauLifeTime::doEvent() Tau" << std::endl;
  
  std::vector<unsigned int> tau_idx_pt, tau_idx_eta, tau_idx_iso, tau_idx_fit, tau_idx_phi, tau_idx_charge;
  unsigned int tau_idx(999);
  if(pass.at(TagIso)==true){
    for(unsigned int i=0;i<Ntp->NPFTaus();i++){
      if(Ntp->PFTau_p4(i).Pt()>cut.at(TauPt)){
        tau_idx_pt.push_back(i);
      }
    }
    for(unsigned int i=0;i<tau_idx_pt.size();i++){
      if(fabs(Ntp->PFTau_p4(tau_idx_pt.at(i)).Eta())<cut.at(TauEta)){
        tau_idx_eta.push_back(tau_idx_pt.at(i));
      }
    }
    for(unsigned int i=0;i<tau_idx_eta.size();i++){
      if(Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA(tau_idx_eta.at(i))){
        tau_idx_iso.push_back(tau_idx_eta.at(i));
      }
    }    
    for(unsigned int i=0;i<tau_idx_iso.size();i++){
      if(Ntp->PFTau_hpsDecayMode(tau_idx_iso.at(i))==10 && 
	 Ntp->PFTau_isHPSByDecayModeFinding(tau_idx_iso.at(i)) && 
	 Ntp->PFTau_TIP_hassecondaryVertex(tau_idx_iso.at(i)) &&
     Ntp->PFTau_TIP_hasA1Momentum(tau_idx_iso.at(i))){
        tau_idx=tau_idx_iso.at(i);
        tau_idx_fit.push_back(tau_idx);
      }
    }
    double delta_phi(0), charge(999);
    if(verbose)std::cout << "void  TauLifeTime::doEvent() Tau -a" << std::endl;
    //event
    //should have only one element, shouldn't it? for(int i=0;i<tau_idx_iso.size();i++){
    if(tau_idx_fit.size()>=1 && mu_idx_iso.size()>=1){
      //std::cout << "The index of the tau is " << tau_idx_fit.at(0) << " = " << tau_idx << std::endl;
      if(fabs(Tools::DeltaPhi(Ntp->Muon_p4(mu_idx),Ntp->PFTau_p4(tau_idx)))>=cut.at(deltaPhi)){
        delta_phi=fabs(Tools::DeltaPhi(Ntp->Muon_p4(mu_idx),Ntp->PFTau_p4(tau_idx)));
        tau_idx_phi.push_back(tau_idx_fit.at(0));
      }
      if(fabs(Ntp->PFTau_Charge(tau_idx) + Ntp->Muon_Charge(mu_idx))<0.5){ 
        charge=Ntp->PFTau_Charge(tau_idx) + Ntp->Muon_Charge(mu_idx);
        tau_idx_charge.push_back(tau_idx_fit.at(0));      
      }
    }

    if(verbose)std::cout << "void  TauLifeTime::doEvent() Tau -b" << std::endl;
  /*std::cout << "tau_idx_charge.size = " << tau_idx_charge.size() << std::endl;
  std::cout << "tau_idx_phi.size = " << tau_idx_phi.size() << std::endl;
  std::cout << "tau_idx_iso.size = " << tau_idx_iso.size() << std::endl;
  std::cout << "tau_idx_eta.size = " << tau_idx_eta.size() << std::endl;
  std::cout << "tau_idx_pt.size = " << tau_idx_pt.size() << std::endl;*/

  value.at(TauPt)=tau_idx_pt.size();
  value.at(TauEta)=tau_idx_eta.size();
  value.at(TauIsIsolated)=tau_idx_iso.size();
  value.at(TauFit)=tau_idx_fit.size();
  value.at(deltaPhi)=delta_phi;
  value.at(Charge)=charge;
  pass.at(TauPt)=(value.at(TauPt)>0);
  pass.at(TauEta)=(value.at(TauEta)>0);
  pass.at(TauIsIsolated)=(value.at(TauIsIsolated)>0);
  pass.at(TauFit)=(value.at(TauFit)==cut.at(TauFit));
  pass.at(deltaPhi)=(value.at(deltaPhi)>=cut.at(deltaPhi));
  pass.at(Charge)=(fabs(charge)<0.5);
  }
  //else{std::cout << "no Tau cuts made since no mu passed mu-cuts" << std::endl;}  
  if(verbose)std::cout << "void  TauLifeTime::doEvent() Tau -c " << std::endl;
  if(mu_idx==999 || tau_idx==999) { 
    //std::cout << "Either mu or tau did not pass" << std::endl;
    value.at(deltaPhi)=-1;
    pass.at(deltaPhi)=false;
    value.at(Charge)=-10.0;
    pass.at(Charge)=false;
  }
  //////////////////////////////////////////////////////////////////////////////////
  // QCD Control sample
  /*  if(!pass.at(Charge) && fabs(value.at(Charge))==1){
    if(Ntp->isData()){
      if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
      pass.at(Charge)=true;
    }
  }
  */

  double wobs(1),w(1);
  if(!Ntp->isData()){
    w*=Ntp->PUWeightFineBins();
  }
  else{w=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
 
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose)std::cout << "void  TauLifeTime::doEvent() -Tau done" << std::endl;
  if(status){
    if(verbose)std::cout << "MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;
    //std::cout << "tau_idx: " << tau_idx << "; mu_idx: " << mu_idx << std::endl;
    //std::cout << "tau_idx_iso.size: " << tau_idx_iso.size() << " mu_idx_iso.size: " << mu_idx_iso.size() << std::endl;
    if (tau_idx!=999 && mu_idx!=999){
      if(verbose)std::cout << "void  TauLifeTime::doEvent() A" << std::endl;
      std::vector<bool> Tau_FitOk; 
      std::vector<TLorentzVector> Tau_sol;
      for(unsigned int j=0;j< MultiProngTauSolver::NAmbiguity;j++){
	LorentzVectorParticle theTau;
	std::vector<LorentzVectorParticle> daughter;
	double LC_chi2(0);
	if(verbose)std::cout << "void  TauLifeTime::doEvent() A1" << std::endl;
	Tau_FitOk.push_back(Ntp->ThreeProngTauFit(tau_idx,j,theTau,daughter,LC_chi2));
	if(verbose)std::cout << "void  TauLifeTime::doEvent() A2" << std::endl;
	Tau_sol.push_back(theTau.LV());
      }
      if(verbose)std::cout << "void  TauLifeTime::doEvent() B" << std::endl;
      double   tauMass=PDGInfo::tau_mass();
      double   c=299792458*100;//units cm/s
      double   MCtauLifetime=2.906*1e-13*1e12;;
      TVector3 flightLength=(Ntp->PFTau_TIP_secondaryVertex_pos(tau_idx)-Ntp->PFTau_TIP_primaryVertex_pos(tau_idx));
      TVector3 MCFlightLength=Ntp->MCSignalParticle_Poca(tau_idx)-Ntp->MCTauandProd_Vertex(tau_idx,0);
      TVector3 a1Momentum = Ntp->PFTau_3PS_A1_LV(tau_idx).Vect();
      TVector3 flightLengthPt=MCFlightLength;flightLengthPt.SetZ(0);
      TVector3 MCFlightLengthPt=MCFlightLength;MCFlightLengthPt.SetZ(0);
      TVector3 a1MomentumPt=a1Momentum;a1MomentumPt.SetZ(0);
      double   signOfFlight=flightLength.Dot(a1Momentum)/fabs(flightLength.Dot(a1Momentum));
      double   signOfFlightPt=flightLengthPt.Dot(a1MomentumPt)/fabs(flightLengthPt.Dot(a1MomentumPt));
      double   tauFlightLength=signOfFlight*flightLength.Mag();
      double   tauFlightLengthPt=signOfFlightPt*flightLengthPt.Mag();

      if(verbose)std::cout << "void  TauLifeTime::doEvent() C" << std::endl;
      // Fill Plots
      TauFlightLength.at(t).Fill(tauFlightLength,w);
      ResTauFlightLength.at(t).Fill(tauFlightLength-MCFlightLength.Mag(),w);
      
      TauFlightLengthTransverse.at(t).Fill(tauFlightLengthPt,w);
      ResTauFlightLengthTransverse.at(t).Fill(tauFlightLengthPt-MCFlightLength.Pt(),w);
      if(verbose)std::cout << "void  TauLifeTime::doEvent() D" << std::endl;
      if(Tau_FitOk.at(MultiProngTauSolver::zero)){
	if(verbose)std::cout << "void  TauLifeTime::doEvent() E" << std::endl;
	TVector3 tauMomentum0 = Tau_sol.at(MultiProngTauSolver::zero).Vect();
	TVector3 resTauMomentum0 = tauMomentum0-Ntp->MCTau_p4(tau_idx).Vect();
        double tauLife = 1e12*tauMass*tauFlightLength/(c*tauMomentum0.Mag());
	double tauLifePt = 1e12*tauMass*tauFlightLengthPt/(c*tauMomentum0.Pt());
	// Fill Plots
	TauMomentum.at(t).Fill(tauMomentum0.Mag(),w);
	ResTauMomentum.at(t).Fill(tauMomentum0.Mag()-MCFlightLength.Mag(),w);
	ResTauLife.at(t).Fill(tauLife-MCtauLifetime,w);
	TauLife.at(t).Fill(tauLife,w);

	TauMomentumTransverse.at(t).Fill(tauMomentum0.Pt(),w);
	ResTauMomentumTransverse.at(t).Fill(tauMomentum0.Pt()-MCFlightLength.Pt(),w);
	TauLifeTransverse.at(t).Fill(tauLifePt,w);
	ResTauLifeTransverse.at(t).Fill(tauLifePt-MCtauLifetime,w);
	if(verbose)std::cout << "void  TauLifeTime::doEvent() F" << std::endl;
      }
      if(Tau_FitOk.at(MultiProngTauSolver::minus) && Tau_FitOk.at(MultiProngTauSolver::plus)){
	if(verbose)std::cout << "void  TauLifeTime::doEvent() G" << std::endl;
	TVector3  tauMomentum0((Tau_sol.at(MultiProngTauSolver::minus).Px()+Tau_sol.at(MultiProngTauSolver::plus).Px())/2,
			       (Tau_sol.at(MultiProngTauSolver::minus).Py()+Tau_sol.at(MultiProngTauSolver::plus).Py())/2,
			       (Tau_sol.at(MultiProngTauSolver::minus).Pz()+Tau_sol.at(MultiProngTauSolver::plus).Pz())/2);
        double tauLife = 1e12*tauMass*tauFlightLength/(c*tauMomentum0.Mag());
        double tauLifePt = 1e12*tauMass*tauFlightLengthPt/(c*tauMomentum0.Pt());
        // Fill Plots
        TauMomentumAvg.at(t).Fill(tauMomentum0.Mag(),w);
        ResTauMomentumAvg.at(t).Fill(tauMomentum0.Mag()-MCFlightLength.Mag(),w);
        ResTauLifeAvg.at(t).Fill(tauLife-MCtauLifetime,w);
        TauLifeAvg.at(t).Fill(tauLife,w);

        TauMomentumTransverseAvg.at(t).Fill(tauMomentum0.Pt(),w);
        ResTauMomentumTransverseAvg.at(t).Fill(tauMomentum0.Pt()-MCFlightLength.Pt(),w);
        TauLifeTransverseAvg.at(t).Fill(tauLifePt,w);
        ResTauLifeTransverseAvg.at(t).Fill(tauLifePt-MCtauLifetime,w);
	if(verbose)std::cout << "void  TauLifeTime::doEvent() H" << std::endl;
      }
    }
  }
}

