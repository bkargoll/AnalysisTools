#include "ZtoEMu_Skim.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include <stdlib.h>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

#include "Parameters.h"
#include "TMath.h"

#include <TFile.h>
#include <sstream>

ZtoEMu_Skim::ZtoEMu_Skim(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_ptlow(10)
  ,mu_pthigh(20)
  ,e_ptlow(10)
  ,e_pthigh(20)
  ,mu_eta(2.4)
  ,e_eta(2.5)
  ,mmin(20)
  ,jet_pt(18)
  ,jet_eta(5.2) // 2.5 in dileptonic top selection
  ,singlejet(40)
  ,zmin(88)
  ,zmax(94)
  ,mtmu(60)
  ,ptbalance(10)
  ,csvl(0.244)
  ,csvm(0.679)
  ,csvt(0.898)
{
    //verbose=true;

	doHiggsObjects = false;
	doWWObjects = true;
	useMadgraphZ = true;
	if(useMadgraphZ) mmin = 50;
	if(doHiggsObjects){
		mu_eta = 2.1;
		e_eta = 2.3;
	}

	// corrections to be applied to particle candidates
	mucorr = "roch";
	ecorr = "";
	jetcorr = "runJER";
}

ZtoEMu_Skim::~ZtoEMu_Skim(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZtoEMu::~ZtoEMu Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZtoEMu::~ZtoEMu()" << std::endl;
}

void  ZtoEMu_Skim::Configure(){

  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
    if(i==NMu)                cut.at(NMu)=1;
    if(i==NE)                 cut.at(NE)=1;
    if(i==ptthreshold)        cut.at(ptthreshold)=1;
    if(i==mll)                cut.at(mll)=mmin;
    if(i==charge)             cut.at(charge)=0;
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
      title.at(i)="Number of Prime Vertices $(N>=$";
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
    else if(i==NMu){
      title.at(i)="Number $\\mu >=$";
      title.at(i)+=cut.at(NMu);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of #mu";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NMu_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==NE){
      title.at(i)="Number $e >=$";
      title.at(i)+=cut.at(NE);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Number of e";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NE_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NE_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
    else if(i==ptthreshold){
	  title.at(i)="ptthreshold ";
	  hlabel="ptthreshold ";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ptthreshold_",htitle,17,-0.5,16.5,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ptthreshold_",htitle,17,-0.5,16.5,hlabel,"Events"));
    }
    else if(i==mll){
	  title.at(i)="$M_{e,\\mu} > $";
	  char buffer[50];
	  sprintf(buffer,"%5.2f",cut.at(mll));
	  title.at(i)+=buffer;
	  title.at(i)+="(GeV)";
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="m_{e#mu} / GeV";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_mll_",htitle,41,19,142,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_mll_",htitle,41,19,142,hlabel,"Events"));
    }
    else if(i==charge){
	  title.at(i)="$e-\\mu$ Charge = ";
	  title.at(i)+=cut.at(charge);
	  title.at(i)+="";
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="e-#mu Charge";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,21,-5.5,5.5,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,21,-5.5,5.5,hlabel,"Events"));
    }

    // calling external files (e.g. root files for efficiencies)
    TString baseCode = "";
    TString baseEff = "";
	if(runtype==GRID){
		baseEff = (TString)std::getenv("PWD")+"/Code/"+"/nehrkorn/";
	}
	else if(runtype==Local){
		baseEff = (TString)std::getenv("workdir")+"/Code/"+"/nehrkorn/";
	}
	RSF = new ReferenceScaleFactors(runtype);
	FRFile = new TFile(baseEff+"FakeRates_2012_19ifb_rereco.root");
	ZptCorrFile = new TFile(baseEff+"zpt_correction_2012_roch.root");

	ZptCorrection = (TH1D*)(ZptCorrFile->Get("zptratio"));

	ElectronFakeRate35 = (TH2D*)(FRFile->Get("ElectronFakeRateHist_35"));
	ElectronFakeRate20 = (TH2D*)(FRFile->Get("ElectronFakeRateHist_20"));
	ElectronFakeRate50 = (TH2D*)(FRFile->Get("ElectronFakeRateHist_50"));
	MuonFakeRate15 = (TH2D*)(FRFile->Get("MuonFakeRateHist_15"));
	MuonFakeRate5 = (TH2D*)(FRFile->Get("MuonFakeRateHist_5"));
	MuonFakeRate30 = (TH2D*)(FRFile->Get("MuonFakeRateHist_30"));


    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu_Skim::Store_ExtraDist(){

}

void  ZtoEMu_Skim::doEvent(){
  if(verbose)std::cout << "ZtoEMu::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}

  ///////////////////////////////////////////////
  //
  // Trigger passed?
  //
  if(verbose)std::cout << "Trigger" << std::endl;
  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") || Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")){
	  value.at(TriggerOk)=1;
  }
  if(Ntp->GetMCID()==DataMCType::DY_emu_embedded)value.at(TriggerOk)=1;
  pass.at(TriggerOk)=(value.at(TriggerOk)==cut.at(TriggerOk));

  // Apply Selection

  ///////////////////////////////////////////////
  //
  // Vertex selection
  //
  if(verbose)std::cout << "Vertex selection" << std::endl;
  unsigned int nGoodVtx=0;
  int vertex = -1;
  for(unsigned i=0;i<Ntp->NVtx();i++){
	  if(Ntp->isGoodVtx(i)){
		  if(vertex==-1)vertex=i;
		  nGoodVtx++;
	  }
  }
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));

  ///////////////////////////////////////////////
  //
  // Mu Cuts
  //
  if(verbose) std::cout << "Muon cuts" << std::endl;
  std::vector<unsigned int> GoodMuons;
  std::vector<unsigned int> Fakemuons;
  
  // muon ID cuts
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->Muon_p4(i,mucorr).Pt()>mu_ptlow
			  && fabs(Ntp->Muon_p4(i,mucorr).Eta())<mu_eta
			  && vertex>=0
			  && (Ntp->matchTrigger(Ntp->Muon_p4(i,mucorr),0.2,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","muon") || Ntp->matchTrigger(Ntp->Muon_p4(i,mucorr),0.2,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","muon") || Ntp->GetMCID()==DataMCType::DY_emu_embedded)
			  ){
		  if(doHiggsObjects){
			  if(Ntp->isSelectedMuon(i,vertex,0.02,0.1,mucorr)
					  && ((fabs(Ntp->Muon_p4(i,mucorr).Eta())<1.479 && Ntp->Muon_RelIso(i,mucorr)<0.15) || (fabs(Ntp->Muon_p4(i,mucorr).Eta())>=1.479 && Ntp->Muon_RelIso(i,mucorr)<0.10))
							  ){
				  GoodMuons.push_back(i);
			  }else if(isFakeMuon(i,vertex,mucorr)
					  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
					  ){
				  Fakemuons.push_back(i);
				  GoodMuons.push_back(i);
			  }
		  }else{
			  if(Ntp->isTightMuon(i,vertex,mucorr)
					  && Ntp->Muon_RelIso(i,mucorr)<0.12
					  ){
				  GoodMuons.push_back(i);
			  }else if(doWWObjects
					  && isFakeMuon(i,vertex,mucorr)
					  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
					  ){
				  Fakemuons.push_back(i);
				  GoodMuons.push_back(i);
			  }
		  }
	  }
  }
  
  value.at(NMu)=GoodMuons.size();
  pass.at(NMu)=(value.at(NMu)>=cut.at(NMu));
  
  unsigned int muidx(999);
  double hardestmu(0);
  if(GoodMuons.size()>1){
	  for(unsigned i=0;i<GoodMuons.size();i++){
		  if(Ntp->Muon_p4(GoodMuons.at(i),mucorr).Pt()>hardestmu){
			  hardestmu = Ntp->Muon_p4(GoodMuons.at(i),mucorr).Pt();
			  muidx = GoodMuons.at(i);
		  }
	  }
  }
  if(GoodMuons.size()==1){muidx=GoodMuons.at(0);}

  ///////////////////////////////////////////////
  //
  // E Cuts
  //
  if(verbose) std::cout << "electrons cuts" << std::endl;
  std::vector<unsigned int> GoodElectrons;
  std::vector<unsigned int> Fakeelectrons;
  bool hasMuonTrack = false;
  bool matchRecoMuon = false;
  
  // electron ID cuts (eta-dependent MVA or simple cut-based)
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  matchRecoMuon = false;
	  if(Ntp->Electron_p4(i,ecorr).Et()>e_ptlow
			  && fabs(Ntp->Electron_supercluster_eta(i))<e_eta
			  && vertex>=0
			  && (Ntp->matchTrigger(Ntp->Electron_p4(i,ecorr),0.2,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","electron") || Ntp->matchTrigger(Ntp->Electron_p4(i,ecorr),0.2,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","electron") || Ntp->GetMCID()==DataMCType::DY_emu_embedded)
			  ){
		  // no overlapping reco muons
		  for(unsigned j=0;j<Ntp->NMuons();j++){
			  if(Ntp->Muon_p4(j,mucorr).Pt()<3) continue;
			  if(fabs(Ntp->Muon_p4(j,mucorr).Eta())>2.4) continue;
			  if(Ntp->Electron_p4(i,ecorr).DeltaR(Ntp->Muon_p4(j,mucorr))<0.3) matchRecoMuon = true;
		  }
		  if(matchRecoMuon) continue;
		  if(doHiggsObjects){
			  if(Ntp->isSelectedElectron(i,vertex,0.02,0.1,ecorr)
					  && ((fabs(Ntp->Electron_supercluster_eta(i))<1.479 && Ntp->Electron_RelIso04(i,ecorr)<0.15) || (fabs(Ntp->Electron_supercluster_eta(i))>=1.479 && Ntp->Electron_RelIso04(i,ecorr)<0.10))
					  ){
				  GoodElectrons.push_back(i);
			  }else if(isFakeElectron(i,vertex,ecorr)
					  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
					  ){
				  Fakeelectrons.push_back(i);
				  GoodElectrons.push_back(i);
			  }
		  }else{
			  if(isWWElectron(i,vertex,ecorr)
					  && Ntp->Electron_RelIso04(i,ecorr)<0.15
					  ){
				  GoodElectrons.push_back(i);
			  }else if(doWWObjects
					  && isFakeElectron(i,vertex,ecorr)
					  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
					  ){
				  Fakeelectrons.push_back(i);
				  GoodElectrons.push_back(i);
			  }
		  }
	  }
  }
  
  value.at(NE)=GoodElectrons.size();
  pass.at(NE)=(value.at(NE)>=cut.at(NE));
  
  unsigned int eidx(999);
  double hardeste(0);
  if(GoodElectrons.size()>1){
	  for(unsigned i=0;i<GoodElectrons.size();i++){
		  if(Ntp->Electron_p4(GoodElectrons.at(i),ecorr).Et()>hardeste){
			  hardeste = Ntp->Electron_p4(GoodElectrons.at(i),ecorr).Et();
			  eidx = GoodElectrons.at(i);
		  }
	  }
  }
  if(GoodElectrons.size()==1){eidx=GoodElectrons.at(0);}
  
  ///////////////////////////////////////////////
  //
  // pt thresholds
  //
  if(verbose) std::cout << "Setting pt thresholds" << std::endl;
  bool passembed = false;
  bool leadingmu = false;
  value.at(ptthreshold)=0;
  if(muidx!=999 && eidx!=999){
	  value.at(ptthreshold)=1;
	  if(Ntp->Muon_p4(muidx,mucorr).Pt()<=mu_ptlow || Ntp->Electron_p4(eidx,ecorr).Et()<=e_ptlow) value.at(ptthreshold)=0;
	  if(Ntp->Muon_p4(muidx,mucorr).Pt()<mu_pthigh && Ntp->Electron_p4(eidx,ecorr).Et()<e_pthigh) value.at(ptthreshold)=0;
	  if(value.at(ptthreshold)==1 && Ntp->GetMCID()==DataMCType::DY_emu_embedded) passembed = true;
	  if(Ntp->Muon_p4(muidx,mucorr).Pt()<mu_pthigh){
		  if(!Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) value.at(ptthreshold)=0;
	  }
	  else if(Ntp->Electron_p4(eidx,ecorr).Et()<e_pthigh){
		  if(!Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) value.at(ptthreshold)=0;
		  else leadingmu = true;
	  }
	  else if(Ntp->Muon_p4(muidx,mucorr).Pt()>mu_pthigh && Ntp->Electron_p4(eidx,ecorr).Et()>e_pthigh){
		  if(Ntp->Muon_p4(muidx,mucorr).Pt()>Ntp->Electron_p4(eidx,ecorr).Et()) leadingmu = true;
	  }
	  if(passembed) value.at(ptthreshold)=1;
  }
  pass.at(ptthreshold)=(value.at(ptthreshold)==cut.at(ptthreshold));

  ///////////////////////////////////////////////
  //
  // m(ll)
  //
  if(verbose) std::cout << "m(ll)" << std::endl;
  //value.at(mll)=mmin+1;
  value.at(mll)=61;
  if(muidx!=999 && eidx!=999){
	  value.at(mll)=(Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).M();
  }
  //pass.at(mll)=(value.at(mll)>cut.at(mll));
  pass.at(mll)=(value.at(mll)>=60. && value.at(mll)<120.);

  ///////////////////////////////////////////////
  //
  // charge cut
  //
  if(verbose) std::cout << "Charge cut" << std::endl;
  value.at(charge)=-1;
  if(eidx!=999 && muidx!=999){
	  value.at(charge)=Ntp->Electron_Charge(eidx)*Ntp->Muon_Charge(muidx);
	  if(Ntp->Electron_Charge(eidx)==-999 || Ntp->Muon_Charge(muidx)==-999) value.at(charge)=1.;
  }
  pass.at(charge)=(value.at(charge)<cut.at(charge));

  ////////////////////////////////////////////////
  //
  // QCD
  //
  if(verbose) std::cout << "QCD" << std::endl;
  fakeRate = 1.;
  bool fakemu = false;
  bool fakee = false;
  if(doHiggsObjects || doWWObjects){
	  for(unsigned i=0;i<Fakemuons.size();i++){
		  if(Fakemuons.at(i)==muidx){
			  fakemu=true;
			  if(doHiggsObjects || doWWObjects) fakeRateMu = Fakerate(Ntp->Muon_p4(muidx,mucorr).Pt(),Ntp->Muon_p4(muidx,mucorr).Eta(),MuonFakeRate15);
			  break;
		  }
	  }
	  for(unsigned i=0;i<Fakeelectrons.size();i++){
		  if(Fakeelectrons.at(i)==eidx){
			  fakee=true;
			  if(doHiggsObjects || doWWObjects) fakeRateE = Fakerate(Ntp->Electron_p4(eidx,ecorr).Et(),Ntp->Electron_supercluster_eta(eidx),ElectronFakeRate35);
			  break;
		  }
	  }
	  if(pass.at(charge)
			  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
			  ){
		  if(fakemu || fakee) fakeRate = 0.;
	  }
	  if(!pass.at(charge)
			  && Ntp->isData()
			  ){
		  if(fakemu && !fakee){
			  fakeRate = fakeRateMu;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge)=true;
		  }else if(fakee && !fakemu){
			  fakeRate = fakeRateE;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge)=true;
		  }else if(fakemu && fakee){
			  fakeRate = fakeRateMu*fakeRateE;
			  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
			  pass.at(charge)=true;
		  }
	  }
	  if(!pass.at(charge)
			  && (Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
			  ){
		  if(fakemu && !fakee){
			  fakeRate = -fakeRateMu;
			  pass.at(charge)=true;
		  }else if(fakee && !fakemu){
			  fakeRate = -fakeRateE;
			  pass.at(charge)=true;
		  }else if(fakemu && fakee){
			  fakeRate = -fakeRateMu*fakeRateE;
			  pass.at(charge)=true;
		  }
	  }
	  if(fabs(fakeRate)>0 && fabs(fakeRate)<1) fakeRate*=0.83;
  }

  ///////////////////////////////////////////////
  //
  // Weights
  //
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1);
  if(!Ntp->isData() && Ntp->GetMCID()!=DataMCType::DY_emu_embedded){
    w*=Ntp->PUWeight()*fakeRate;
    if(pass.at(NE)){
    	if(doHiggsObjects){
    		w*=RSF->HiggsTauTau_EMu_Id_E(Ntp->Electron_p4(eidx,ecorr).Et(),Ntp->Electron_supercluster_eta(eidx));
    		w*=RSF->ElectronReconstruction2012(Ntp->Electron_p4(eidx,ecorr).Et(),Ntp->Electron_supercluster_eta(eidx));
    		w*=RSF->HiggsTauTau_EMu_Trigger_E(Ntp->Electron_p4(eidx,ecorr).Et(),Ntp->Electron_supercluster_eta(eidx));
    	}else{
    		w*=RSF->ElectronIdTrig2012(Ntp->Electron_p4(eidx,ecorr).Et(),Ntp->Electron_supercluster_eta(eidx));
    		w*=RSF->ElectronReconstruction2012(Ntp->Electron_p4(eidx,ecorr).Et(),Ntp->Electron_supercluster_eta(eidx));
    	}
    }
    if(pass.at(NMu)){ // for systematics: add systematic uncertainty of 1.5%(0.5%) when pt<20(>20) for Id and 0.2% for isolation when pt>20 to statistical in quadrature.
    	if(doHiggsObjects){
    		w*=RSF->HiggsEMuId_Mu(Ntp->Muon_p4(muidx,mucorr));
    		w*=RSF->HiggsTauTau_EMu_Trigger_Mu(Ntp->Muon_p4(muidx,mucorr));
    	}else{
    		w*=RSF->MuonIdTight2012(Ntp->Muon_p4(muidx,mucorr));
    		w*=RSF->MuonIsoTight2012(Ntp->Muon_p4(muidx,mucorr));
    		w*=RSF->TrackingEfficiency2012(Ntp->Muon_p4(muidx,mucorr));
    	}
    }
    if(pass.at(TriggerOk)
    		&& pass.at(NMu)
    		&& pass.at(NE)
    		&& !doHiggsObjects
    		){
    	if(leadingmu) w*=RSF->HiggsWW_EMu_Trigger(Ntp->Muon_p4(muidx,mucorr),Ntp->Electron_p4(eidx,ecorr).Et(),Ntp->Electron_supercluster_eta(eidx),"Mu17_Ele8");
    	else w*=RSF->HiggsWW_EMu_Trigger(Ntp->Muon_p4(muidx,mucorr),Ntp->Electron_p4(eidx,ecorr).Et(),Ntp->Electron_supercluster_eta(eidx),"Mu8_Ele17");
    }
    if(pass.at(NMu)
    		&& pass.at(NE)
    		&& !useMadgraphZ
    		){
    	//if(Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau)w*=PowhegReweight((Ntp->Muon_p4(muidx)+Ntp->Electron_p4(eidx,ecorr)).Pt());
    }
    if(verbose)std::cout << "void  ZtoEMu::doEvent() k" << w << " " << wobs << std::endl;
  }
  else if(!Ntp->isData() && Ntp->GetMCID()==DataMCType::DY_emu_embedded){
	  w*=Ntp->EmbeddedWeight();
	  if(pass.at(NE)) w*=RSF->ElectronEmbedding2012(Ntp->Electron_p4(eidx,ecorr).Et(),Ntp->Electron_supercluster_eta(eidx));
  }
  else{w=1*fakeRate;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;

  ///////////////////////////////////////////////
  //
  // Add plots
  //
  if(verbose) std::cout << "add plots" << std::endl;

  if(verbose)std::cout << "ZtoEMu::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

double ZtoEMu_Skim::ZPtReweight(double zpt){
	double weight = 1;
	if(zpt<2.5)	weight/=1.03445;
	else if(zpt>=2.5 && zpt<5.) weight/=1.02996;
	else if(zpt>=5. && zpt<7.5)	weight/=1.01111;
	else if(zpt>=7.5 && zpt<10.) weight/=1.02506;
	else if(zpt>=10. && zpt<12.5) weight/=1.03364;
	else if(zpt>=12.5 && zpt<15.) weight/=0.999692;
	else if(zpt>=15. && zpt<17.5) weight/=0.980893;
	else if(zpt>=17.5 && zpt<20.) weight/=1.03963;
	else if(zpt>=20. && zpt<22.5) weight/=0.889207;
	else if(zpt>=22.5 && zpt<25.) weight/=0.886043;
	else if(zpt>=25. && zpt<27.5) weight/=0.851487;
	else if(zpt>=27.5 && zpt<30.) weight/=0.951472;
	return weight;
}

double ZtoEMu_Skim::PowhegReweight(double zpt){
	double weight = 1.;
	weight = ZptCorrection->GetBinContent(ZptCorrection->FindFixBin(zpt));
	return weight;
}

// relative uncertainties obtained from table 16 in AN-2012-225
double ZtoEMu_Skim::ZPtRelUnc(double zpt){
	double unc = 0.;
	if(zpt<2.5) unc = 0.027;
	else if(zpt>=2.5 && zpt<5.0) unc = 0.014;
	else if(zpt>=5.0 && zpt<7.5) unc = 0.0050;
	else if(zpt>=7.5 && zpt<10.) unc = 0.013;
	else if(zpt>=10. && zpt<12.5) unc = 0.013;
	else if(zpt>=12.5 && zpt<15.0) unc = 0.023;
	else if(zpt>=15.0 && zpt<17.5) unc = 0.013;
	else if(zpt>=17.5 && zpt<20.0) unc = 0.016;
	else if(zpt>=20.0 && zpt<30.0) unc = 0.0042;
	else if(zpt>=30.0 && zpt<40.0) unc = 0.0068;
	else if(zpt>=40.0 && zpt<50.0) unc = 0.013;
	else if(zpt>=50.0 && zpt<70.0) unc = 0.016;
	else if(zpt>=70.0 && zpt<90.0) unc = 0.0027;
	else if(zpt>=90.0 && zpt<110.) unc = 0.0037;
	else if(zpt>=110. && zpt<150.) unc = 0.037;
	else if(zpt>=150. && zpt<190.) unc = 0.0020;
	else if(zpt>=190. && zpt<250.) unc = 0.12;
	else if(zpt>=250. && zpt<600.) unc = 0.020;
	return unc;
}

double ZtoEMu_Skim::ZPtMadgraphRelUnc(double zpt){
	double unc = 0.;
	if(zpt<2.5) unc = 0.027;
	else if(zpt>=2.5 && zpt<5.0) unc = 0.013;
	else if(zpt>=5.0 && zpt<7.5) unc = 0.0028;
	else if(zpt>=7.5 && zpt<10.) unc = 0.013;
	else if(zpt>=10. && zpt<12.5) unc = 0.014;
	else if(zpt>=12.5 && zpt<15.0) unc = 0.023;
	else if(zpt>=15.0 && zpt<17.5) unc = 0.013;
	else if(zpt>=17.5 && zpt<20.0) unc = 0.016;
	else if(zpt>=20.0 && zpt<30.0) unc = 0.0041;
	else if(zpt>=30.0 && zpt<40.0) unc = 0.0056;
	else if(zpt>=40.0 && zpt<50.0) unc = 0.010;
	else if(zpt>=50.0 && zpt<70.0) unc = 0.0026;
	else if(zpt>=70.0 && zpt<90.0) unc = 0.0037;
	else if(zpt>=90.0 && zpt<110.) unc = 0.0067;
	else if(zpt>=110. && zpt<150.) unc = 0.011;
	else if(zpt>=150. && zpt<190.) unc = 0.0014;
	else if(zpt>=190. && zpt<250.) unc = 0.099;
	else if(zpt>=250. && zpt<600.) unc = 0.0045;
	return unc;
}

double ZtoEMu_Skim::calculatePzeta(int muiterator, int eiterator){
  double pex=Ntp->Electron_p4(eiterator).Px();
  double pey=Ntp->Electron_p4(eiterator).Py();
  double pmux=Ntp->Muon_p4(muiterator).Px();
  double pmuy=Ntp->Muon_p4(muiterator).Py();
  double phie=Ntp->Electron_p4(eiterator).Phi();
  double phimu=Ntp->Muon_p4(muiterator).Phi();
  double combpt=TMath::Sqrt(pow(pex+pmux,2)+pow(pey+pmuy,2));
  double aemu=TMath::ACos(pmux*pex+pmuy*pey/(Ntp->Muon_p4(muiterator).Pt()*Ntp->Electron_p4(eiterator).Pt()));
  double phismall = 0.;
  if(phie<phimu && fabs(phie-phimu)<TMath::Pi()) phismall=phie;
  else if(phimu<phie && fabs(phie-phimu)>TMath::Pi())phismall=phie;
  else if(phie<phimu && fabs(phie-phimu)>TMath::Pi())phismall=phimu;
  else if(phimu<phie && fabs(phie-phimu)<TMath::Pi())phismall=phimu;
  double beta=TMath::ACos(((pex+pmux)*TMath::Cos(phismall+0.5*aemu)+(pey+pmuy)*TMath::Sin(phismall+0.5*aemu))/combpt);
  double gamma=TMath::ACos((Ntp->MET_CorrT0pcT1_ex()*TMath::Cos(phismall+0.5*aemu)+Ntp->MET_CorrT0pcT1_ey()*TMath::Sin(phismall+0.5*aemu))/Ntp->MET_CorrT0pcT1_et());
  if(Ntp->MET_CorrT0pcT1_phi()>(phismall+0.5*aemu+0.5*TMath::Pi()) && Ntp->MET_CorrT0pcT1_phi()<(phismall+0.5*aemu+1.5*TMath::Pi()))gamma*=-1;
  double pvis=TMath::Sin(beta)*combpt;
  double pmiss=TMath::Sin(gamma)*Ntp->MET_CorrT0pcT1_et();
  return pmiss-pvis;
}

double ZtoEMu_Skim::calculatePzetaDQM(int muiterator, int eiterator){
	double cosPhi1 = TMath::Cos(Ntp->Electron_p4(eiterator).Phi());
	double sinPhi1 = TMath::Sin(Ntp->Electron_p4(eiterator).Phi());
	double cosPhi2 = TMath::Cos(Ntp->Muon_p4(muiterator).Phi());
	double sinPhi2 = TMath::Sin(Ntp->Muon_p4(muiterator).Phi());
	double zetaX = cosPhi1 + cosPhi2;
	double zetaY = sinPhi1 + sinPhi2;
	double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
	if(zetaR>0.){
		zetaX/=zetaR;
		zetaY/=zetaR;
	}
	double pxVis=Ntp->Electron_p4(eiterator).Px()+Ntp->Muon_p4(muiterator).Px();
	double pyVis=Ntp->Electron_p4(eiterator).Py()+Ntp->Muon_p4(muiterator).Py();
	double pZetaVis=pxVis*zetaX+pyVis*zetaY;
	double px=pxVis+Ntp->MET_CorrT0pcT1_ex();
	double py=pyVis+Ntp->MET_CorrT0pcT1_ey();
	double pZeta=px*zetaX+py*zetaY;
	return pZeta-1.5*pZetaVis;
}

double ZtoEMu_Skim::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double ZtoEMu_Skim::cosphi3d(TVector3 vec1, TVector3 vec2){
	return (vec1.Dot(vec2))/vec1.Mag()/vec2.Mag();
}

int ZtoEMu_Skim::findBin(TGraphAsymmErrors* graph, double xval){
	int bin = -1;
	for(unsigned i=0;i<graph->GetN()-1;i++){
		if(fabs(graph->GetX()[i]-xval)<fabs(graph->GetX()[i+1]-xval)) bin = i;
		else if(fabs(graph->GetX()[i]-xval)>=fabs(graph->GetX()[i+1]-xval)) bin = i+1;
	}
	return bin;
}

//////////////////////////////
//
// Muon related functions
//

bool ZtoEMu_Skim::isFakeMuon(unsigned int idx, TString corr){
	if(!Ntp->Muon_isGlobalMuon(idx)) return false;
	if(Ntp->Muon_p4(idx,corr).Pt()<=10) return false;
	if(fabs(Ntp->Muon_p4(idx,corr).Eta())>2.4) return false;
	if(Ntp->Muon_p4(idx,corr).Pt()<=20){
		if(Ntp->Muon_sumPt03(idx)>=8.) return false;
		if(Ntp->Muon_emEt03(idx)>=8.) return false;
		if(Ntp->Muon_hadEt03(idx)>=8.) return false;
	}
	if(Ntp->Muon_p4(idx,corr).Pt()>20){
		if(Ntp->Muon_sumPt03(idx)/Ntp->Muon_p4(idx,corr).Pt()>=0.4) return false;
		if(Ntp->Muon_emEt03(idx)/Ntp->Muon_p4(idx,corr).Pt()>=0.4) return false;
		if(Ntp->Muon_hadEt03(idx)/Ntp->Muon_p4(idx,corr).Pt()>=0.4) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isFakeMuon(unsigned int idx, unsigned int vtx, TString corr){
	if(vtx<0 || vtx>=Ntp->NVtx()) return false;
	if(!isFakeMuon(idx,corr)) return false;
	if(Ntp->dxy(Ntp->Muon_p4(idx,corr),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.2) return false;
	if(Ntp->dz(Ntp->Muon_p4(idx,corr),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.1) return false;
	return true;
}

//////////////////////////////
//
// Electron related functions
//

bool ZtoEMu_Skim::isWWElectron(unsigned int idx, unsigned int vtx, TString corr){
	double mvapt = Ntp->Electron_p4(idx,corr).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(idx));
	if(mvapt<10.) return false;
	if(mvaeta>2.5) return false;
	if(!isFakeElectron(idx,vtx,corr)) return false;
	if(mvapt>10. && mvapt<20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_Trig_discriminator(idx)<=0.00) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_Trig_discriminator(idx)<=0.10) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_Trig_discriminator(idx)<=0.62) return false;
	}
	if(mvapt>=20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_Trig_discriminator(idx)<=0.94) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_Trig_discriminator(idx)<=0.85) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_Trig_discriminator(idx)<=0.92) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isFakeElectron(unsigned int idx, TString corr){
	if(Ntp->Electron_p4(idx,corr).Pt()<10) return false;
	if(fabs(Ntp->Electron_supercluster_eta(idx))>2.5) return false;
	if(Ntp->Electron_HasMatchedConversions(idx)) return false;
	if(Ntp->Electron_numberOfMissedHits(idx)>0) return false;
	if(Ntp->Electron_tkSumPt03(idx)/Ntp->Electron_p4(idx,corr).Pt()>0.2) return false;
	if(std::max(Ntp->Electron_ecalRecHitSumEt03(idx)-1.,0.)/Ntp->Electron_p4(idx,corr).Pt()>0.2) return false;
	if((Ntp->Electron_hcalDepth1TowerSumEt03(idx)+Ntp->Electron_hcalDepth2TowerSumEt03(idx))/Ntp->Electron_p4(idx,corr).Pt()>0.2) return false;
	if(fabs(Ntp->Electron_supercluster_eta(idx))<1.479){
		if(Ntp->Electron_sigmaIetaIeta(idx)>0.01) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(idx)>0.15) return false;
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(idx)>0.007) return false;
		if(Ntp->Electron_hadronicOverEm(idx)>0.12) return false;
	}
	if(fabs(Ntp->Electron_supercluster_eta(idx))>=1.479){
		if(Ntp->Electron_sigmaIetaIeta(idx)>0.03) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(idx)>0.10) return false;
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(idx)>0.009) return false;
		if(Ntp->Electron_hadronicOverEm(idx)>0.10) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isFakeElectron(unsigned int idx, unsigned int vtx, TString corr){
	if(vtx<0 || vtx>=Ntp->NVtx()) return false;
	if(!isFakeElectron(idx,corr)) return false;
	if(Ntp->dz(Ntp->Electron_p4(idx,corr),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>0.1) return false;
	if(Ntp->dxy(Ntp->Electron_p4(idx,corr),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>0.02) return false;
	return true;
}

//////////////////////////////
//
// Calculate fakerate
//

double ZtoEMu_Skim::Fakerate(double pt, double eta, TH2D *fakeRateHist){
	double fakerate = 0.;
	if(pt>=35.) pt = 34.99;
	fakerate = fakeRateHist->GetBinContent(fakeRateHist->FindFixBin(pt,eta));
	return fakerate/(1-fakerate);
}

double ZtoEMu_Skim::FakerateError(double pt, double eta, TH2D *fakeRateHist){
	double error = 0.;
	if(pt>=35.) pt = 34.99;
	error = fakeRateHist->GetBinError(fakeRateHist->FindFixBin(pt,eta));
	return error/(1-error);
}

//////////////////////////////
//
// Finish function
//

void ZtoEMu_Skim::Finish(){
	Selection::Finish();
	Selection::ScaleAllHistOfType(DataMCType::DY_emu_embedded,10152445./9760728.);
}
