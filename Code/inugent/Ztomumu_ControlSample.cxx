#include "Ztomumu_ControlSample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

Ztomumu_ControlSample::Ztomumu_ControlSample(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_pt(22)
  ,mu_eta(2.0)
  ,jet_pt(20)
  ,jet_eta(2.4)
{
  //  verbose=true;
}

Ztomumu_ControlSample::~Ztomumu_ControlSample(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "Ztomumu_ControlSample::~Ztomumu_ControlSample Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "Ztomumu_ControlSample::~Ztomumu_ControlSample()" << std::endl;
}

void  Ztomumu_ControlSample::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==NJets)              cut.at(NJets)=1;
    if(i==MET)                cut.at(MET)=40;
    if(i==NMu)                cut.at(NMu)=2;
    if(i==NMuPt)              cut.at(NMuPt)=2;
    if(i==NMuEta)             cut.at(NMuEta)=2 ;
    if(i==MuIso)              cut.at(MuIso)=0.4 ;
    if(i==MuMuVertex)         cut.at(MuMuVertex)=0.5;
    if(i==deltaPhi)           cut.at(deltaPhi)=-1/sqrt(2);
    if(i==charge)             cut.at(charge)=0;
    if(i==ZMassmax)           cut.at(ZMassmax)=PDG_Var::Z_mass()+30;
    if(i==ZMassmin)           cut.at(ZMassmin)=PDG_Var::Z_mass()-30;
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
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,17,-0.5,16.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,17,-0.5,16.5,hlabel,"Events"));
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
    else if(i==NMu){
      title.at(i)="Number $\\mu >=$";
      title.at(i)+=cut.at(NMu);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of Kin. Fit #mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==NMuPt){
      title.at(i)="Number of $\\mu$ [$P_{T}^{\\mu}>$";
      title.at(i)+=mu_pt;
      title.at(i)+="(GeV)] $>=$";
      title.at(i)+=cut.at(NMuPt);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu [Kin. Fit+P^{T}]";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuPt_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NMuPt_","P_{T,#mu} (N-1 Distribution)",100,0,200,"P_{T,#mu} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_NMuPt_","P_{T,#mu} (Accumulative Distribution)",100,0,200,"P_{T,#mu} (GeV)","Events");
    }
    else if(i==NMuEta){
      title.at(i)="Number of $\\mu$ [$|\\eta^{\\mu}|<$";
      title.at(i)+=mu_eta;
      title.at(i)+="] $>=$";
      title.at(i)+=cut.at(NMuPt);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu [Kin. Fit+P^{T}+#eta]";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMuEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMuEta_",htitle,6,-0.5,5.5,hlabel,"Events"));
      distindx.at(i)=true;
      Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_NMuEta_","#eta{#mu} (N-1 Distribution)",100,0,200,"#eta_{#mu} (GeV)","Events");
      Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accumdist_NMuEta_","#eta_{#mu} (Accumulative Distribution)",100,0,200,"#eta_{#mu} (GeV)","Events");
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
    else if(i==charge){
      title.at(i)="$\\mu-\\mu$ Charge = ";
      title.at(i)+=cut.at(charge);
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#mu-#mu Charge ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HT_",htitle,21,-5.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HT_",htitle,21,-5.5,5.5,hlabel,"Events"));
    }
    else if(i==deltaPhi){
      title.at(i)="$cos(\\phi_{\\mu-\\mu}) < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(deltaPhi));
      title.at(i)+=buffer;
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="cos(#phi_{#mu-#mu}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HadTopMass_",htitle,40,-1,1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HadTopMass_",htitle,40,-1,1,hlabel,"Events"));
    }
    else if(i==ZMassmax){
      title.at(i)="$M_{\\mu-\\mu} > $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(ZMassmax));
      title.at(i)+=buffer;
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{#mu,#mu} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuMETTopMT_",htitle,60,0,300,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuMETTopMT_",htitle,60,0,300,hlabel,"Events"));
    }
    else if(i==ZMassmin){
      title.at(i)="$M_{\\mu-\\mu} < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(ZMassmin));
      title.at(i)+=buffer;
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{#mu,#mu} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuMETdphi_",htitle,60,0,300,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuMETdphi_",htitle,60,0,300,hlabel,"Events"));
   } 
    else if(i==MuIso){
      title.at(i)="$dr(\\mu,Jet) > $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(MuIso));
      title.at(i)+=buffer;
      title.at(i)+="";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="dr(#mu,Jet) ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuMETdphi_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuMETdphi_",htitle,40,0,200,hlabel,"Events"));
    }
    else if(i==MuMuVertex){
      title.at(i)="has $\\mu-\\mu$ vertex";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="has #mu-#mu vertex ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MuMuVertex_",htitle,100,-5,5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MuMuVertex_",htitle,100,-5,5,hlabel,"Events"));
    }
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");


  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  Ztomumu_ControlSample::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);

}

void  Ztomumu_ControlSample::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() A" << std::endl;

  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("IsoMu24"))value.at(TriggerOk)+=1;
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);

  // Apply Selection
  unsigned int nGoodVtx=0;
  for(unsigned int i=0;i<Ntp->NVtx();i++){
    if(Ntp->isVtxGood(i))nGoodVtx++;
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() B" << std::endl;

  ///////////////////////////////////////////////
  //
  // Mu Cuts
  //
  std::vector<unsigned int> GoodMuons;
  for(unsigned i=0;i<Ntp->NMuons();i++){
    if(Ntp->isGoodMuon(i))GoodMuons.push_back(i);
  }
  value.at(NMu)=GoodMuons.size();
  pass.at(NMu)=(value.at(NMu)>=cut.at(NMu));
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() C" << std::endl;

  for(unsigned i=0;i<GoodMuons.size();i++){
    dist.at(NMuPt).push_back(Ntp->Muon_p4(GoodMuons.at(i)).Pt());
    if(Ntp->Muon_p4(GoodMuons.at(i)).Pt()<mu_pt){
      GoodMuons.erase(GoodMuons.begin()+i);
      i--;
    }
  }
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() D" << std::endl;
  value.at(NMuPt)=GoodMuons.size();
  pass.at(NMuPt)=(value.at(NMuPt)>=cut.at(NMuPt));
  for(unsigned i=0;i<GoodMuons.size();i++){
    dist.at(NMuEta).push_back(Ntp->Muon_p4(GoodMuons.at(i)).Eta());
    if(fabs(Ntp->Muon_p4(GoodMuons.at(i)).Eta())>mu_eta){
      GoodMuons.erase(GoodMuons.begin()+i);
      i--;
    }
  }
  value.at(NMuEta)=GoodMuons.size();
  pass.at(NMuEta)=(value.at(NMuEta)==cut.at(NMuEta));
  unsigned int muidx1(999),muidx2(999),muInfoidx1(999),muInfoidx2(999);
  if(GoodMuons.size()>=1){muidx1=GoodMuons.at(0);muInfoidx1=muidx1;}
  if(GoodMuons.size()>=2){muidx2=GoodMuons.at(1);muInfoidx2=muidx2;}
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() E " << muidx1 <<" "<< muInfoidx1 << " " << muidx2 << " "  << muInfoidx2 << std::endl;

  // mu-mu Charge
  value.at(charge)=-5;
  if(muInfoidx1!=999 && muInfoidx2!=999){
    value.at(charge)=Ntp->Muon_Charge(muInfoidx1)+Ntp->Muon_Charge(muInfoidx2);
  }
  pass.at(charge)=(value.at(charge)==0);
 
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() F" << std::endl;
  value.at(deltaPhi)=1;
  if(muidx1!=999 && muidx2!=999){
    TLorentzVector Mu1=Ntp->Muon_p4(muidx1);
    TLorentzVector Mu2=Ntp->Muon_p4(muidx2);
    value.at(deltaPhi)=cos(Tools::DeltaPhi(Mu1.Phi(),Mu2.Phi()));
    TLorentzVector Z=Mu1+Mu2;
    value.at(ZMassmax)=Z.M();
    value.at(ZMassmin)=Z.M();
  }
  pass.at(deltaPhi)=(value.at(deltaPhi)<cut.at(deltaPhi));
  pass.at(ZMassmax)=(value.at(ZMassmax)<cut.at(ZMassmax));
  pass.at(ZMassmin)=(value.at(ZMassmin)>cut.at(ZMassmin));
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() F0" << std::endl;
  value.at(MuIso)=5;
  int overlap1(0), overlap2(0);
  for(int i=0;i<Ntp->NPFJets();i++){
    double dR=5;
    if(muidx1!=999) dR=Tools::dr(Ntp->Muon_p4(muidx1),Ntp->PFJet_p4(i));
    if(dR<cut.at(MuIso))overlap1++;
    if(muidx2!=999) dR=Tools::dr(Ntp->Muon_p4(muidx2),Ntp->PFJet_p4(i));
    if(dR<cut.at(MuIso))overlap2++;
  }
  pass.at(MuIso)=(2==overlap1+overlap2);
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() F1" << std::endl;
  pass.at(MuMuVertex)=false;
  value.at(MuMuVertex)=0;
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() F2" << std::endl;
  
  if(muInfoidx1!=999 && muInfoidx2!=999){
    value.at(MuMuVertex)=Ntp->Muon_Poca(muInfoidx1).Z()-Ntp->Muon_Poca(muInfoidx2).Z();
    pass.at(MuMuVertex)=fabs(value.at(MuMuVertex))<cut.at(MuMuVertex);
  }

  ///////////////////////////////////////////////
  //
  // Jet Cuts
  //
  std::vector<unsigned int> GoodJets;
  for(int i=0;i<Ntp->NPFJets();i++){
    if(/*Ntp->isGoodJet(i) &&*/ Ntp->PFJet_p4(i).Pt()>jet_pt && fabs(Ntp->PFJet_p4(i).Eta())<jet_eta){
      bool overlap=false;
      for(unsigned j=0;j<GoodMuons.size();j++){
	//	if(Tools::dr(Ntp->Muon_p4(GoodMuons.at(j)),Ntp->PFJet_p4(i))<0.5) overlap=true;
      }
      if(!overlap)GoodJets.push_back(i);
    }
  }

  value.at(NJets)=GoodJets.size();
  //pass.at(NJets)=(value.at(NJets)>=cut.at(NJets));
  pass.at(NJets)=true;
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() G" << std::endl;
  ///////////////////////////////////////////////
  //
  // Event Shape and Energy Cuts
  //
  double MET_Ex(Ntp->MET_CorrMVA_ex()),MET_Ey(Ntp->MET_CorrMVA_ey());
  // correct for neutrino from mus
  value.at(MET)=sqrt(MET_Ex*MET_Ex+MET_Ey*MET_Ey);
  pass.at(MET)=(value.at(MET)<cut.at(MET));

  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() H" << std::endl;
  ///////////////////////////////////////////////////////////
  //Do QCD bkg
  /*if(!pass.at(charge)){
    if(Ntp->isData()){
      if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
      pass.at(charge)=true;
    }
    }*/
  //////////////////////////////////////////////////////////
  double wobs(1),w(1);
  if(!Ntp->isData()){
    if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() J" << std::endl;
    //w*=Ntp->PUWeightFineBins();
    if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
  if(verbose)std::cout << "void  Ztomumu_ControlSample::doEvent() L" << std::endl;
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




