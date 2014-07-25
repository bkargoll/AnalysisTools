#include "ChargedHiggs_dilepontic.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"

ChargedHiggs_dilepontic::ChargedHiggs_dilepontic(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mutrigname("IsoMu24")
  ,tau_pt(20)
  ,tau_eta(2.0)
  ,jet_pt(30)
  ,jet_eta(2.4)
  ,muon_pt(26)
  ,muon_eta(2.4)
{
  //  verbose=true;
}

ChargedHiggs_dilepontic::~ChargedHiggs_dilepontic(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ChargedHiggs_dilepontic::~ChargedHiggs_dilepontic Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ChargedHiggs_dilepontic::~ChargedHiggs_dilepontic()" << std::endl;
}

void  ChargedHiggs_dilepontic::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==NMuons)             cut.at(NMuons)=1;
    if(i==MuonTrigMatch)      cut.at(MuonTrigMatch)=0.2;
    if(i==MuRelIso)           cut.at(MuRelIso)=0.2;
    if(i==N2Jets)             cut.at(N2Jets)=2;
    if(i==MaxJetPT)           cut.at(MaxJetPT)=60;
    if(i==ElectronVeto)       cut.at(ElectronVeto)=0;
    if(i==MuonVeto)           cut.at(MuonVeto)=0;
    if(i==NBJets)             cut.at(NBJets)=2;
    if(i==MET)                cut.at(MET)=40;
    if(i==MTplusMET)          cut.at(MTplusMET)=50;
    if(i==NTauKinFit)         cut.at(NTauKinFit)=1;
    if(i==NTauPt)             cut.at(NTauPt)=1;
    if(i==NTauEta)            cut.at(NTauEta)=1 ;
    if(i==MT)                 cut.at(MT)=30;
    if(i==MuMETJetdphi)       cut.at(MuMETJetdphi)=-1.0/sqrt(2.0);
    if(i==MTandMuMETdphi)     cut.at(MTandMuMETdphi)=20;
    if(i==Charge)             cut.at(Charge)=0.0;
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
  
    if(i==PrimeVtx){
      title.at(i)="Number of Prime Vertices $(N>$";
      title.at(i)+=cut.at(PrimeVtx);
      title.at(i)+=")";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Prime Vertices";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_PrimeVtx_",htitle,31,-0.5,30.5,hlabel,"Events"));
    }
    else if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,9,-0.5,8.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,9,-0.5,8.5,hlabel,"Events"));
    }
    else if(i==NMuons){
      title.at(i)="Number of Muons $>=$";
      title.at(i)+=cut.at(NMuons);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of NMuons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuons_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuons_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==MuRelIso){
      title.at(i)="Relative Muon Isolation $<$";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(MuRelIso));
      title.at(i)+=buffer;
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Relative Muon Isolation";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuRelIso_",htitle,20,0.0,1.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuRelIso_",htitle,20,0.0,1.0,hlabel,"Events"));
    }
    else if(i==MuonTrigMatch){
      title.at(i)="Muon Trigger Matching $<$";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(MuonTrigMatch));
      title.at(i)+=buffer;
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Muon Trigger Matching (dR)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonTrigMatch_",htitle,20,0.0,1.0,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonTrigMatch_",htitle,20,0.0,1.0,hlabel,"Events"));
    }
    else if(i==N2Jets){
      title.at(i)="Number of Jets $>=$";
      title.at(i)+=cut.at(N2Jets);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Jets";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_N2Jets_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_N2Jets_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==MaxJetPT){
      title.at(i)="$P_{T}^{Leading-Jet}>";
      title.at(i)+=cut.at(MaxJetPT);
      title.at(i)+="$ (GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="P_{T}^{Leading-Jet} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MaxJetPT_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MaxJetPT_",htitle,40,0,200,hlabel,"Events"));
    }
    else if(i==ElectronVeto){
      title.at(i)="Number of extra Electrons $>=$";
      htitle=title.at(i);
      title.at(i)+=cut.at(ElectronVeto);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of extra Electrons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ElectronVeto_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ElectronVeto_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==MuonVeto){
      title.at(i)="Number of extra Muons $>=$";
      htitle=title.at(i);
      title.at(i)+=cut.at(MuonVeto);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of extra Muons";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuonVeto_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuonVeto_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==NBJets){
      title.at(i)="Number of BJets $>=$";
      title.at(i)+=cut.at(NBJets);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of BJets";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NBJets_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NBJets_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==NTauKinFit){
      title.at(i)="Number of Kin. Fit $\\tau >=$";
      title.at(i)+=cut.at(NTauKinFit);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Kin. Fit #tau";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauKinFit_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauKinFit_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==NTauPt){
      title.at(i)="Number of $\\tau$ [$P_{T}^{\\tau}>$";
      title.at(i)+=tau_pt;
      title.at(i)+="(GeV)] $>=$";
      title.at(i)+=cut.at(NTauPt);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #tau [Kin. Fit+P^{T}]";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NTauPt_","P_{T,#tau} (N-1 Distribution)",20,0,200,"P_{T,#tau} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_NTauPt_","P_{T,#tau} (Accumulative Distribution)",20,0,200,"P_{T,#tau} (GeV)","Events");
    }
    else if(i==NTauEta){
      title.at(i)="Number of $\\tau$ [$|\\eta^{\\tau}|<$";
      title.at(i)+=tau_eta;
      title.at(i)+="] $>=$";
      title.at(i)+=cut.at(NTauPt);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #tau [Kin. Fit+P^{T}+#eta]";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NTauEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NTauEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NTauEta_","#eta{#tau} (N-1 Distribution)",20,-5,5,"#eta_{#tau} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accumdist_NTauEta_","#eta_{#tau} (Accumulative Distribution)",20,-5,5,"#eta_{#tau} (GeV)","Events");
    }
    else if(i==Charge){
      title.at(i)="$\\tau-\\mu$ charge =";
      title.at(i)+=cut.at(Charge);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#tau-#mu Charge";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_Charge_",htitle,7,-3.5,3.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_Charge_",htitle,7,-3.5,3.5,hlabel,"Events"));
    }
    else if(i==MET){
      title.at(i)="$E_{T}^{Miss} > $";
      title.at(i)+=cut.at(MET);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="E_{T}^{Miss} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MET_",htitle,15,0,300,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MET_",htitle,15,0,300,hlabel,"Events"));
    } 
   else if(i==MTplusMET){
     title.at(i)="$M_{T}+E_{T}^{Miss} > $";
     title.at(i)+=cut.at(MTplusMET);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{T}+E_{T}^{Miss} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MTplusMET_",htitle,25,0,750,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MTplusMET_",htitle,25,0,750,hlabel,"Events"));
   }
   else if(i==MT){
     title.at(i)="$M_{T} > $";
     title.at(i)+=cut.at(MT);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{T} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MT_",htitle,15,0,300,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MT_",htitle,15,0,300,hlabel,"Events"));
   }
   else if(i==MuMETJetdphi){
    title.at(i)="$cos(\\phi(\\mu+MET,Jet^{Lead})) < $";
    char buffer[50];
    sprintf(buffer,"%5.2f",cut.at(MuMETJetdphi));
    title.at(i)+=buffer;
    title.at(i)+="(rad)";
    htitle=title.at(i);
    htitle.ReplaceAll("$","");
    htitle.ReplaceAll("\\","#");
    hlabel="cos(#phi(#mu+MET,Jet^{Lead})) (rad)";
    Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuMETJetdphi_",htitle,10,-1,1,hlabel,"Events"));
    Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuMETJetdphi_",htitle,10,-1,1,hlabel,"Events"));
  }
   else if(i==MTandMuMETdphi){
      title.at(i)="$[\\alpha(1- cos(\\phi(\\mu,MET)))^2+M_{T}^2]^{1/2} < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(MTandMuMETdphi));
      title.at(i)+=buffer;
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="[#alpha(1- cos(#phi(#mu,MET)))^2+M_{T}^2]^{1/2}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MTandMuMETdphi_",htitle,40,0,800,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MTandMuMETdphi_",htitle,40,0,800,hlabel,"Events"));
   } 

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx_","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx_","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx_","NTracksperVtx",30,-0.5,150.5,"Number of Track per Vertex","Events");
  ChargedHiggsMT=HConfig.GetTH1D(Name+"_ChargedHiggsMT_","ChargedHiggsMT",100,0,250,"M_{T,#tau,MET}","Events");
  TagEtaPT=HConfig.GetTH2D(Name+"_TagEtaPT_","TagEtaPT",25,0,2.5,50,0,50,"#eta","P_{T}^{Tag}");
  METvsMT=HConfig.GetTH2D(Name+"_METvsMT_","METvsMT",30,0,300,30,0,300,"M_{T} (GeV)","E_{T}^{Miss} (GeV)");
  MTvsMuMETdphi=HConfig.GetTH2D(Name+"_MTvsMuMETdphi_","MTvsMuMETdphi",20,-1,1,30,0,300,"cos(\\phi(\\mu,E_{T}^{Miss}))","M_{T} (GeV)");
  METvsMuMETdphi=HConfig.GetTH2D(Name+"_METvsMuMETdphi_","METvsMuMETdphi",20,-1,1,30,0,300,"cos(\\phi(\\mu,E_{T}^{Miss}))","E_{T}^{Miss} (GeV)");
  METvsMTvsMuMETdphi=HConfig.GetTH3F(Name+"_METvsMTvsMuMETdphi_","METvsMTvsMuMETdphi",30,0,300,30,0,300,20,-1,1,"M_{T} (GeV)","E_{T}^{Miss} (GeV)","cos(\\phi(\\mu,E_{T}^{Miss}))");
  Jet1stEtaPT=HConfig.GetTH2D(Name+"_Jet1stEtaPT_","Jet1stEtaPT",24,-2.4,2.4,20,0,200,"#eta^{Jet,1st}","P_{T}^{Jet,1st} (GeV)");
  Jet2ndEtaPT=HConfig.GetTH2D(Name+"_Jet2ndEtaPT_","Jet2ndEtaPT",24,-2.4,2.4,20,0,200,"#eta^{Jet,2nd}","P_{T}^{Jet,2nd} (GeV)");
  METWHypvsMT=HConfig.GetTH2D(Name+"_METWHypvsMT_","METWHypvsMT",30,0,300,30,0,300,"M_{T} (GeV)","E_{T}^{Miss,W-Hypothesis} (GeV)");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ChargedHiggs_dilepontic::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist2d.push_back(&TagEtaPT);
 Extradist2d.push_back(&METvsMT);
 Extradist2d.push_back(&MTvsMuMETdphi);
 Extradist2d.push_back(&METvsMuMETdphi);
 Extradist2d.push_back(&Jet1stEtaPT);
 Extradist2d.push_back(&Jet2ndEtaPT);
 Extradist2d.push_back(&METWHypvsMT);
 Extradist3d.push_back(&METvsMTvsMuMETdphi);
}

void  ChargedHiggs_dilepontic::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() A" << std::endl;
  //////////////////////////////////////////////
  //
  // Setup Trigger
  //
  unsigned int mutrig=0;
  Ntp->GetTriggerIndex(mutrigname,mutrig);
  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept(mutrigname))value.at(TriggerOk)+=1;
  pass.at(TriggerOk)=value.at(TriggerOk)>=cut.at(TriggerOk);

  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() B" << std::endl;

  ///////////////////////////////////////////////
  //
  // Tau Cuts
  //
  std::vector<unsigned int> GoodTaus;
  for(unsigned i=0;i<Ntp->NKFTau();i++){
    if(Ntp->isGoodKFTau(i))GoodTaus.push_back(i);
  }
  value.at(NTauKinFit)=GoodTaus.size();
  pass.at(NTauKinFit)=(value.at(NTauKinFit)>=cut.at(NTauKinFit));

  for(unsigned i=0;i<GoodTaus.size();i++){
    dist.at(NTauPt).push_back(Ntp->KFTau_TauFit_p4(GoodTaus.at(i)).Pt());
    if(Ntp->KFTau_TauFit_p4(GoodTaus.at(i)).Pt()<tau_pt){
      GoodTaus.erase(GoodTaus.begin()+i);
      i--;
    }
  }
  value.at(NTauPt)=GoodTaus.size();
  pass.at(NTauPt)=(value.at(NTauPt)>=cut.at(NTauPt));
  for(unsigned i=0;i<GoodTaus.size();i++){
    dist.at(NTauEta).push_back(Ntp->KFTau_TauFit_p4(GoodTaus.at(i)).Eta());
    if(fabs(Ntp->KFTau_TauFit_p4(GoodTaus.at(i)).Eta())>tau_eta){
      GoodTaus.erase(GoodTaus.begin()+i);
      i--;
    }
  }
  value.at(NTauEta)=GoodTaus.size();
  pass.at(NTauEta)=(value.at(NTauEta)==cut.at(NTauEta));
  unsigned int tauidx(999);
  if(GoodTaus.size()>0)tauidx=GoodTaus.at(0);
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() C " << tauidx << std::endl;

  ////////////////////////////////////////////////
  //
  // Muons
  //
  std::vector<unsigned int> mu_idx;
  for(unsigned int i=0; i<Ntp->NMuons(); i++){
    if(Ntp->isGoodMuon(i) && Ntp->Muon_p4(i).Pt()>muon_pt && fabs(Ntp->Muon_p4(i).Eta())<muon_eta){mu_idx.push_back(i);}
  }
  value.at(NMuons)=mu_idx.size();
  pass.at(NMuons)=(value.at(NMuons)==cut.at(NMuons));
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() C-1 " << std::endl;
  float max_mu_pt=0;
  value.at(MuRelIso)=1.1;
  for(unsigned int i=0; i<mu_idx.size();i++){
    if(Ntp->Muon_p4(mu_idx.at(i)).Pt()>max_mu_pt){
      value.at(MuRelIso)=Ntp->Muon_RelIso(mu_idx.at(i));
    }
  }
  pass.at(MuRelIso)=(value.at(MuRelIso)<cut.at(MuRelIso));
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() C-2 " << std::endl;

  value.at(MuonTrigMatch)=998; 
  if(mu_idx.size()>0)value.at(MuonTrigMatch)=Ntp->MuonTriggerMatch(mutrig,mu_idx.at(0));
  pass.at(MuonTrigMatch)=(value.at(MuonTrigMatch)<cut.at(MuonTrigMatch));

  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() C-3 " << std::endl;
  value.at(ElectronVeto)=0;
  pass.at(ElectronVeto)=true;

  value.at(MuonVeto)=0;
  pass.at(MuonVeto)=true;
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() C-4 " << std::endl;
  ///////////////////////////////////////////////
  //
  // Jet Cuts
  //
  std::vector<unsigned int> GoodJets;
  for(int i=0;i<Ntp->NPFJets();i++){
    if(Ntp->isGoodJet(i) && Ntp->PFJet_p4(i).Pt()>jet_pt && fabs(Ntp->PFJet_p4(i).Eta())<jet_eta){
      bool overlap=false;
      for(unsigned j=0;j<GoodTaus.size();j++){
	if(Tools::dr(Ntp->KFTau_TauFit_p4(GoodTaus.at(j)),Ntp->PFJet_p4(i))<0.5) overlap=true;
      }
      if(!overlap)GoodJets.push_back(i);
    }
  }

  value.at(N2Jets)=GoodJets.size();
  pass.at(N2Jets)=(value.at(N2Jets)>=cut.at(N2Jets));
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() C-5 " << std::endl;
  value.at(MaxJetPT)=0;
  for(unsigned int i=0; i<GoodJets.size();i++){
    if(value.at(MaxJetPT)<Ntp->PFJet_p4(GoodJets.at(i)).Pt())value.at(MaxJetPT)=Ntp->PFJet_p4(GoodJets.at(i)).Pt();
  }
  pass.at(MaxJetPT)=(value.at(MaxJetPT)>cut.at(MaxJetPT));

  // use highest 2 pt jets for GoodBJets
  std::vector<unsigned int> GoodBJets(2,0);
  if(GoodJets.size()>0) GoodBJets.at(0)=GoodJets.at(0);
  if(GoodJets.size()>1) GoodBJets.at(1)=GoodJets.at(1);
  for(unsigned i=2;i<GoodJets.size();i++){ 
    if(Ntp->PFJet_p4(GoodJets.at(i)).Pt()>Ntp->PFJet_p4(GoodBJets.at(0)).Pt()){
      GoodBJets.at(1)=GoodBJets.at(0);
      GoodBJets.at(0)=GoodJets.at(i);
    }
    else if(Ntp->PFJet_p4(GoodJets.at(i)).Pt()>Ntp->PFJet_p4(GoodBJets.at(1)).Pt()){
      GoodBJets.at(1)=GoodJets.at(i);
    }
  }
  if(Ntp->PFJet_p4(GoodBJets.at(1)).Pt()>Ntp->PFJet_p4(GoodBJets.at(0)).Pt()){
    unsigned int Temp=GoodBJets.at(0);
    GoodBJets.at(0)=GoodBJets.at(1);
    GoodBJets.at(1)=Temp;
  }
  value.at(NBJets)=GoodBJets.size();
  pass.at(NBJets)=true;//(value.at(NBJets)>=cut.at(NBJets));
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() D " << tauidx << std::endl;
  ///////////////////////////////////////////////
  //
  // Event Shape and Energy Cuts
  // 1) compute variables
  // 2) run cuts

  double MET_Ex(Ntp->MET_CorrMVA_ex()),MET_Ey(Ntp->MET_CorrMVA_ey());
  // correct for neutrino from taus
  for(unsigned int i=0; i<GoodTaus.size();i++){
    MET_Ex-=Ntp->KFTau_Neutrino_p4(GoodTaus.at(i)).Px();
    MET_Ey-=Ntp->KFTau_Neutrino_p4(GoodTaus.at(i)).Py();
  }

  TLorentzVector Met(MET_Ex,MET_Ey,0,sqrt(MET_Ex*MET_Ex+MET_Ey*MET_Ey));

  TLorentzVector Tau(0,0,0,0);
  if(GoodTaus.size()>0){
    Tau=Ntp->KFTau_TauFit_p4(GoodTaus.at(0));
  }

  value.at(MT)=0;
  if(mu_idx.size()>0){
    value.at(MT)=sqrt(2*(Ntp->MET_CorrMVA_et())*Ntp->Muon_p4(mu_idx.at(0)).Pt()*fabs(1-cos(Ntp->Muon_p4(mu_idx.at(0)).Phi()-Ntp->MET_CorrMVA_phi())));
    pass.at(MT)=(value.at(MT)>cut.at(MT));
  }
  pass.at(MT)=true;

  value.at(MET)=sqrt(MET_Ex*MET_Ex+MET_Ey*MET_Ey);
  pass.at(MET)=(value.at(MET)>cut.at(MET));

  float dphiMuMET=0;
  if(mu_idx.size()>0){
    dphiMuMET=cos(Tools::DeltaPhi(Ntp->Muon_p4(mu_idx.at(0)),Met));
  }
  value.at(MTandMuMETdphi)=sqrt(pow((double)(1-dphiMuMET)*200,2.0)+pow((double)value.at(MT),2.0));
  pass.at(MTandMuMETdphi)=(value.at(MTandMuMETdphi)>=cut.at(MTandMuMETdphi));

  value.at(MTplusMET)=value.at(MET)+value.at(MT);
  pass.at(MTplusMET)=(value.at(MTplusMET)>=cut.at(MTplusMET));


  value.at(MuMETJetdphi)=0;
  if(mu_idx.size()>0 && GoodBJets.size()>0){
    TLorentzVector LVMuMET=Met;LVMuMET+=Ntp->Muon_p4(mu_idx.at(0));
    value.at(MuMETJetdphi)=cos(Tools::DeltaPhi(Ntp->PFJet_p4(GoodBJets.at(0)),LVMuMET));
  }
  pass.at(MuMETJetdphi)=(value.at(MuMETJetdphi)>cut.at(MuMETJetdphi));
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() H " << tauidx << std::endl;


  ///////////////////////////////////////////////////////////
  value.at(Charge)=0;
  pass.at(Charge)=false;
  if(tauidx!=999 && mu_idx.size()>0){
    value.at(Charge)=Ntp->KFTau_Fit_charge(tauidx)+Ntp->Muon_Charge(mu_idx.at(0));
    pass.at(Charge)=value.at(Charge)==cut.at(Charge);
  }
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() I" << std::endl;

  ///////////////////////////////////////////////////////////
  double wobs(1),w(1);
  if(!Ntp->isData()){
    if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() J" << std::endl;
    w*=Ntp->PUWeightFineBins();
    if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1;}
  if(verbose) std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  if(verbose)std::cout << "void  ChargedHiggs_dilepontic::doEvent() L" << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w*wobs);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w*wobs);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w*wobs);;
    if(GoodBJets.at(0)>=0)Jet1stEtaPT.at(t).Fill(Ntp->PFJet_p4(GoodBJets.at(0)).Eta(),Ntp->PFJet_p4(GoodBJets.at(0)).Pt(),w*wobs);
    if(GoodBJets.at(1)>=0)Jet2ndEtaPT.at(t).Fill(Ntp->PFJet_p4(GoodBJets.at(1)).Eta(),Ntp->PFJet_p4(GoodBJets.at(1)).Pt(),w*wobs);
  }
  if(NMinusL(MT,MET/*,MTplusMET*/)){
    METvsMT.at(t).Fill(value.at(MT),value.at(MET),w*wobs);
    MTvsMuMETdphi.at(t).Fill(dphiMuMET,value.at(MT),w*wobs);
    METvsMuMETdphi.at(t).Fill(dphiMuMET,value.at(MET),w*wobs);
    METvsMTvsMuMETdphi.at(t).Fill(value.at(MT),value.at(MET),dphiMuMET,w*wobs);
    METWHypvsMT.at(t).Fill(value.at(MT),Ntp->MET_CorrMVA_et(),w*wobs);
  }

}




