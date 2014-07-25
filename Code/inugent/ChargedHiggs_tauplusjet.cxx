#include "ChargedHiggs_tauplusjet.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"

ChargedHiggs_tauplusjet::ChargedHiggs_tauplusjet(TString Name_, TString id_):
  Selection(Name_,id_)
  ,tau_pt(20)
  ,tau_eta(2.0)
  ,jet_pt(45)
  ,jet_eta(2.4)
{

}

ChargedHiggs_tauplusjet::~ChargedHiggs_tauplusjet(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ChargedHiggs_tauplusjet::~ChargedHiggs_tauplusjet Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ChargedHiggs_tauplusjet::~ChargedHiggs_tauplusjet()" << std::endl;
}

void  ChargedHiggs_tauplusjet::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==N1Jets)             cut.at(N1Jets)=1;
    if(i==N2Jets)             cut.at(N2Jets)=2;
    if(i==N3Jets)             cut.at(N3Jets)=3;
    if(i==NJets)              cut.at(NJets)=4;
    if(i==NBJets)             cut.at(NBJets)=2;
    if(i==MET)                cut.at(MET)=50;
    if(i==HT)                 cut.at(HT)=200;
    if(i==NTauKinFit)         cut.at(NTauKinFit)=1;
    if(i==NTauPt)             cut.at(NTauPt)=1;
    if(i==NTauEta)            cut.at(NTauEta)=1 ;
    if(i==HadWMass)           cut.at(HadWMass)=30;
    if(i==HadTopMass)         cut.at(HadTopMass)=50;
    if(i==TauMETTopMT)        cut.at(TauMETTopMT)=50;
    if(i==TauMETdphi)         cut.at(TauMETdphi)=-1/sqrt(2);
    if(i==etaq)               cut.at(etaq)=0.1;
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
    else if(i==NJets){
      title.at(i)="Number of Jets $>=$";
      title.at(i)+=cut.at(NJets);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Jets";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NJets_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NJets_",htitle,11,-0.5,10.5,hlabel,"Events"));
    }
    else if(i==N1Jets){
      title.at(i)="Number of Jets $>=$";
      title.at(i)+=cut.at(N1Jets);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Jets";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_N1Jets_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_N1Jets_",htitle,11,-0.5,10.5,hlabel,"Events"));
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
    else if(i==N3Jets){
      title.at(i)="Number of Jets $>=$";
      htitle=title.at(i);
      title.at(i)+=cut.at(N3Jets);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Jets";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_N3Jets_",htitle,11,-0.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_N3Jets_",htitle,11,-0.5,10.5,hlabel,"Events"));
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
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NTauPt_","P_{T,#tau} (N-1 Distribution)",100,0,200,"P_{T,#tau} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_NTauPt_","P_{T,#tau} (Accumulative Distribution)",100,0,200,"P_{T,#tau} (GeV)","Events");
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
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NTauEta_","#eta{#tau} (N-1 Distribution)",100,0,200,"#eta_{#tau} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accumdist_NTauEta_","#eta_{#tau} (Accumulative Distribution)",100,0,200,"#eta_{#tau} (GeV)","Events");
    }
    else if(i==MET){
      title.at(i)="$E_{T}^{Miss} < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="E_{T}^{Miss} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MET_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MET_",htitle,40,0,200,hlabel,"Events"));
    } 
   else if(i==HT){
     title.at(i)="$H_{T} < $";
     title.at(i)+=cut.at(HT);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="H_{T} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HT_",htitle,40,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HT_",htitle,40,0,200,hlabel,"Events"));
   }
   else if(i==HadWMass){
     title.at(i)="$M_{W,had} < $";
     title.at(i)+=cut.at(HadWMass);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{W,had} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HadWMass_",htitle,75,0,250,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HadWMass_",htitle,75,0,250,hlabel,"Events"));
   }
   else if(i==HadTopMass){
    title.at(i)="$M_{Top,had} < $";
    title.at(i)+=cut.at(HadTopMass);
    title.at(i)+="(GeV)";
    htitle=title.at(i);
    htitle.ReplaceAll("$","");
    htitle.ReplaceAll("\\","#");
    hlabel="M_{Top,had} (GeV)";
    Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HadTopMass_",htitle,100,0,500,hlabel,"Events"));
    Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HadTopMass_",htitle,100,0,500,hlabel,"Events"));
  }
   else if(i==TauMETTopMT){
    title.at(i)="$M_{T,Top,\\tau-side} < $";
    title.at(i)+=cut.at(TauMETTopMT);
    title.at(i)+="(GeV)";
    htitle=title.at(i);
    htitle.ReplaceAll("$","");
    htitle.ReplaceAll("\\","#");
    hlabel="M_{T,Top,#tau-side} (GeV)";
    Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMETTopMT_",htitle,100,0,500,hlabel,"Events"));
    Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMETTopMT_",htitle,100,0,500,hlabel,"Events"));
   }
   else if(i==TauMETdphi){
      title.at(i)="$cos(\\phi(\\tau,MET)) < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(TauMETdphi));
      title.at(i)+=buffer;
      title.at(i)+="(rad)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="cos(#phi(#tau,MET)) (rad)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TauMETdphi_",htitle,32,-1,1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TauMETdphi_",htitle,32,-1,1,hlabel,"Events"));
   } 
   else if(i==etaq){
     title.at(i)="$q_{\\tau}|\\eta_{\\tau}|$";
     char buffer[50];
     sprintf(buffer,"%5.2f",cut.at(etaq));
     title.at(i)+=buffer;
     title.at(i)+="";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="q_{#tau}|#eta_{#tau}|";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_etaq_",htitle,40,-2.0,2.0,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_etaq_",htitle,40,-2.0,2.0,hlabel,"Events"));
   }

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  ChargedHiggsMT=HConfig.GetTH1D(Name+"_ChargedHiggsMT","ChargedHiggsMT",100,0,250,"M_{T,#tau,MET}","Events");
  TagEtaPT=HConfig.GetTH2D(Name+"_TagEtaPT","TagEtaPT",25,0,2.5,50,0,50,"#eta","P_{T}^{Tag}");
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ChargedHiggs_tauplusjet::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist2d.push_back(&TagEtaPT);
}

void  ChargedHiggs_tauplusjet::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() A" << std::endl;

  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("HLT_QuadJet40_IsoPFTau40"))value.at(TriggerOk)+=1;
  if(Ntp->TriggerAccept("HLT_QuadJet45_IsoPFTau45"))value.at(TriggerOk)+=2;
  if(Ntp->TriggerAccept("HLT_QuadJet50_IsoPFTau50"))value.at(TriggerOk)+=4;
  pass.at(TriggerOk)=value.at(TriggerOk)>=cut.at(TriggerOk);

  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() B" << std::endl;

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
  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() C " << tauidx << std::endl;
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
  value.at(N1Jets)=GoodJets.size();
  pass.at(N1Jets)=(value.at(N1Jets)>=cut.at(N1Jets));

  value.at(N2Jets)=GoodJets.size();
  pass.at(N2Jets)=(value.at(N2Jets)>=cut.at(N2Jets));

  value.at(N3Jets)=GoodJets.size();
  pass.at(N3Jets)=(value.at(N3Jets)>=cut.at(N3Jets));

  value.at(NJets)=GoodJets.size();
  pass.at(NJets)=(value.at(NJets)>=cut.at(NJets));

  // not implemented yet just use highest 2 pt jets for GoodBJets
  std::vector<unsigned int> GoodBJets(2,0);
  for(unsigned i=0;i<GoodJets.size();i++){ 
    if(Ntp->PFJet_p4(GoodJets.at(i)).Pt()>Ntp->PFJet_p4(GoodBJets.at(0)).Pt()){
      GoodBJets.at(0)=GoodJets.at(i);
      continue;
    }
    if(Ntp->PFJet_p4(GoodJets.at(i)).Pt()>Ntp->PFJet_p4(GoodBJets.at(1)).Pt())GoodBJets.at(1)=GoodJets.at(i);
  }
  value.at(NBJets)=GoodBJets.size();
  pass.at(NBJets)=(value.at(NBJets)>=cut.at(NBJets));
  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() D " << tauidx << std::endl;
  ///////////////////////////////////////////////
  //
  // Event Shape and Energy Cuts
  //
  double MET_Ex(Ntp->MET_CorrMVA_ex()),MET_Ey(Ntp->MET_CorrMVA_ey());
  // correct for neutrino from taus
  for(unsigned int i=0; i<GoodTaus.size();i++){
    MET_Ex+=Ntp->KFTau_Neutrino_p4(GoodTaus.at(i)).Px();
    MET_Ey+=Ntp->KFTau_Neutrino_p4(GoodTaus.at(i)).Py();
  }
  value.at(MET)=sqrt(MET_Ex*MET_Ex+MET_Ey*MET_Ey);
  pass.at(MET)=(value.at(MET)<cut.at(MET));

  value.at(HT)=0;
  for(unsigned int i=0; i<GoodTaus.size();i++){
    value.at(HT)+=Ntp->KFTau_TauFit_p4(GoodTaus.at(i)).Et();
  }
  for(unsigned i=0;i<GoodJets.size();i++){
    value.at(HT)+=Ntp->PFJet_p4(GoodJets.at(i)).Et();
  }
  pass.at(HT)=(value.at(HT)>=cut.at(HT));
  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() E " << tauidx << std::endl;
  ///////////////////////////////////////////////
  //
  // Top and W Masses
  //
  std::vector<unsigned int> WJetCand=GoodJets;
  for(unsigned int i=0;i<WJetCand.size();i++){
    for(unsigned int j=0;j<GoodBJets.size();j++){
      if(WJetCand.at(i)==GoodBJets.at(j)){
	WJetCand.erase(WJetCand.begin()+i);
	i--;
	break;
      }
    }
  }

  std::vector<unsigned int> WJets(2,0);
  value.at(HadWMass)=999;
  for(unsigned int i=0;i<WJetCand.size();i++){
    for(unsigned int j=i+1;j<WJetCand.size();j++){
      TLorentzVector W=Ntp->PFJet_p4(WJetCand.at(i));
      W+=Ntp->PFJet_p4(WJetCand.at(j));
      if(fabs(PDG_Var::W_mass()-W.M())<fabs(PDG_Var::W_mass()-value.at(HadWMass))){
	WJets.at(0)=WJetCand.at(i);
	WJets.at(1)=WJetCand.at(j);
	value.at(HadWMass)=W.M();
      }
    }
  }
  pass.at(HadWMass)=(fabs(PDG_Var::W_mass()-value.at(HadWMass))<=cut.at(HadWMass));

  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() F " << tauidx << std::endl;
  unsigned int HadTopB(0);
  value.at(HadTopMass)=999;
  for(unsigned int i=0;i<GoodBJets.size();i++){
    TLorentzVector Top=Ntp->PFJet_p4(GoodBJets.at(i));
    Top+=Ntp->PFJet_p4(WJets.at(0));
    Top+=Ntp->PFJet_p4(WJets.at(1));
    if(fabs(PDG_Var::Top_mass()-Top.M())<fabs(PDG_Var::Top_mass()-value.at(HadTopMass))){
      HadTopB=GoodBJets.at(i);
      value.at(HadTopMass)=Top.M();
    }
  }
  pass.at(HadWMass)=(fabs(PDG_Var::Top_mass()-value.at(HadWMass))<=cut.at(HadWMass));

  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() G " << tauidx << std::endl;

  TLorentzVector MET;
  TLorentzVector Tau;
  TLorentzVector BTau;
  value.at(TauMETTopMT)=999;
  for(unsigned int i=0;i<GoodBJets.size();i++){
    if(GoodBJets.at(i)!=HadTopB){
      BTau=Ntp->PFJet_p4(GoodBJets.at(i));
      MET_Ex=Ntp->MET_CorrMVA_ex();
      MET_Ey=Ntp->MET_CorrMVA_ey();
      // correct for neutrino from taus   
      if(GoodTaus.size()>0){
	Tau=Ntp->KFTau_TauFit_p4(GoodTaus.at(0));
	BTau+=Tau;
	MET_Ex+=Ntp->KFTau_Neutrino_p4(GoodTaus.at(0)).Px();
	MET_Ey+=Ntp->KFTau_Neutrino_p4(GoodTaus.at(0)).Py();
      }
      TLorentzVector LV(MET_Ex,MET_Ey,0,sqrt(MET_Ex*MET_Ex+MET_Ey*MET_Ey));
      MET=LV;
      value.at(TauMETTopMT)=sqrt(2*(MET.Pt())*BTau.Pt()*fabs(1-cos(BTau.Phi()-MET.Phi())));
      break;
    }
  }
  pass.at(TauMETTopMT)=true;//(fabs(PDG_Var::Top_mass()-value.at(HadWMass))<=cut.at(HadWMass));

  value.at(TauMETdphi)=cos(Tools::DeltaPhi(Tau,MET));
  pass.at(TauMETdphi)=(value.at(TauMETdphi)<=cut.at(TauMETdphi));
  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() H " << tauidx << std::endl;

  ///////////////////////////////////////////////////////////
  value.at(etaq)=0;
  pass.at(etaq)=false;
  if(tauidx!=999){
    value.at(etaq)=Ntp->KFTau_Fit_charge(tauidx)*fabs(Ntp->KFTau_TauFit_p4(tauidx).Eta());
    pass.at(etaq)=value.at(etaq)>cut.at(etaq);
  }
  pass.at(etaq)=true;
  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() I" << std::endl;

  pass.at(HT)=true;
  pass.at(etaq)=true;
  pass.at(HadWMass)=true;
  pass.at(HadTopMass)=true;
  pass.at(TauMETTopMT)=true;
  pass.at(TauMETdphi)=true;


  ///////////////////////////////////////////////////////////
  double wobs(1),w(1);
  if(!Ntp->isData()){
    if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() J" << std::endl;
    // w*=Ntp->PUWeightFineBins();
    if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1;}
  if(verbose) std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  if(verbose)std::cout << "void  ChargedHiggs_tauplusjet::doEvent() L" << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(status){
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;

  }
}




