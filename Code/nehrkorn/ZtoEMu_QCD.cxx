#include "ZtoEMu_QCD.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

#include "Parameters.h"
#include "TMath.h"

#include <TFile.h>

ZtoEMu_QCD::ZtoEMu_QCD(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_ptlow(10)
  ,mu_pthigh(20)
  ,e_ptlow(10)
  ,e_pthigh(20)
  ,mu_eta(2.4)
  ,e_eta(2.5)
  ,jet_pt(18)
  ,jet_eta(4.7)
  ,jet_sum(70)
  ,zmin(88)
  ,zmax(94)
{
    //verbose=true;
}

ZtoEMu_QCD::~ZtoEMu_QCD(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZtoEMu_QCD::~ZtoEMu_QCD Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZtoEMu_QCD::~ZtoEMu_QCD()" << std::endl;
}

void  ZtoEMu_QCD::Configure(){

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
    if(i==drEMu)              cut.at(drEMu)=0.3;
    if(i==diMuonVeto)         cut.at(diMuonVeto)=0;
    if(i==triLeptonVeto)      cut.at(triLeptonVeto)=0;
    if(i==charge)             cut.at(charge)=0;
    if(i==jetVeto)            cut.at(jetVeto)=jet_sum;
    if(i==mtmu)               cut.at(mtmu)=50;
    if(i==ptbal)              cut.at(ptbal)=20;
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
    else if(i==drEMu){
      title.at(i)="$\\Delta R(e,\\mu) < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(drEMu));
      title.at(i)+=buffer;
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="dR(e,mu)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_drEMu_",htitle,50,0,5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_drEMu_",htitle,50,0,5,hlabel,"Events"));
    }
    else if(i==diMuonVeto){
  	  title.at(i)="di-muon veto";
  	  htitle=title.at(i);
  	  htitle.ReplaceAll("$","");
  	  htitle.ReplaceAll("\\","#");
  	  hlabel="dimuon veto";
  	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_diMuonVeto_",htitle,40,0,20,hlabel,"Events"));
  	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_diMuonVeto_",htitle,40,0,20,hlabel,"Events"));
  	}
    else if(i==triLeptonVeto){
	  title.at(i)="tri-lepton veto";
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="trilepton veto";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_triLeptonVeto_",htitle,40,0,20,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_triLeptonVeto_",htitle,40,0,20,hlabel,"Events"));
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
    else if(i==jetVeto){
        title.at(i)="jet veto";
        htitle=title.at(i);
        htitle.ReplaceAll("$","");
        htitle.ReplaceAll("\\","#");
        hlabel="jet P_{T} / GeV";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_jetVeto_",htitle,33,35,200,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_jetVeto_",htitle,33,35,200,hlabel,"Events"));
    }
    else if(i==mtmu){
        title.at(i)="$m_{T}^{\\mu,Miss} < $";
        title.at(i)+=cut.at(mtmu);
        title.at(i)+="(GeV)";
        htitle=title.at(i);
        htitle.ReplaceAll("$","");
        htitle.ReplaceAll("\\","#");
        hlabel="m_{T}^{#mu,Miss} / GeV";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_mtmu_",htitle,40,0,200,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_mtmu_",htitle,40,0,200,hlabel,"Events"));
    }
    else if(i==ptbal){
        title.at(i)="$p_{t,e+\\mu} < $";
        title.at(i)+=cut.at(ptbal);
        title.at(i)+="(GeV)";
        htitle=title.at(i);
        htitle.ReplaceAll("$","");
        htitle.ReplaceAll("\\","#");
        hlabel="p_t balance / GeV";
        Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ptbal_",htitle,40,0,200,hlabel,"Events"));
        Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ptbal_",htitle,40,0,200,hlabel,"Events"));
    }

    // calling external files (e.g. root files for efficiencies)
    TString base = "";
	if(runtype==GRID){
		base+=std::getenv("PWD");
		base+="/Code/nehrkorn/";
	}
	else if(runtype==Local){
		//base+=Selection::splitString(std::getenv("PWD"),'/',"workdir");
		//base+="/Code/nehrkorn/";
		base+="/net/scratch_cms/institut_3b/nehrkorn/";
	}
	FRFile = new TFile(base+"FakeRates_2012_19ifb_rereco.root");
	EmbEffFile = new TFile(base+"RecHitElectronEfficiencies.root");
	MuIdEffFile = new TFile(base+"MuonEfficiencies_Run2012ReReco_53X.root");
	MuIsoEffFile = new TFile(base+"MuonEfficiencies_ISO_Run_2012ReReco_53X.root");
	ETrigIdEffFile = new TFile(base+"ElectronEfficiencies_Run2012ReReco_53X_Trig.root");
	ENonTrigIdEffFile = new TFile(base+"ElectronEfficiencies_Run2012ReReco_53X_NonTrig.root");
	TriggerEfficiencies = new TFile(base+"TriggerEfficienciesWW.root");

	ElectronFakeRate = (TH2D*)(FRFile->Get("ElectronFakeRateHist"));
	MuonFakeRate = (TH2D*)(FRFile->Get("MuonFakeRateHist"));
	EmbEff = (TH2D*)(EmbEffFile->Get("hPtEtaSFL"));

	ElectronTrigEff = (TH2D*)(ETrigIdEffFile->Get("electronsDATAMCratio_FO_ID_ISO"));
	ElectronNonTrigEff = (TH2D*)(ENonTrigIdEffFile->Get("h_electronScaleFactor_IdIsoSip"));
	MuIdEff09 = (TGraphAsymmErrors*)(MuIdEffFile->Get("DATA_over_MC_Tight_pt_abseta<0.9"));
	MuIdEff12 = (TGraphAsymmErrors*)(MuIdEffFile->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2"));
	MuIdEff21 = (TGraphAsymmErrors*)(MuIdEffFile->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1"));
	MuIdEff24 = (TGraphAsymmErrors*)(MuIdEffFile->Get("DATA_over_MC_Tight_pt_abseta2.1-2.4"));
	MuIsoEff09 = (TGraphAsymmErrors*)(MuIsoEffFile->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta<0.9"));
	MuIsoEff12 = (TGraphAsymmErrors*)(MuIsoEffFile->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta0.9-1.2"));
	MuIsoEff21 = (TGraphAsymmErrors*)(MuIsoEffFile->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta1.2-2.1"));
	MuIsoEff24 = (TGraphAsymmErrors*)(MuIsoEffFile->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta2.1-2.4"));

	SingleEle15 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("SingleEle15"));
	SingleEle25 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("SingleEle25"));
	DoubleEleLead15 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleEleLead15"));
	DoubleEleLead25 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleEleLead25"));
	DoubleEleTrail15 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleEleTrail15"));
	DoubleEleTrail25 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleEleTrail25"));
	SingleMu08 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("SingleMu08"));
	SingleMu12 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("SingleMu12"));
	SingleMu21 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("SingleMu21"));
	SingleMu25 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("SingleMu25"));
	DoubleMuLead12 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleMuLead12"));
	DoubleMuLead21 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleMuLead21"));
	DoubleMuLead25 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleMuLead25"));
	DoubleMuTrail12 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleMuTrail12"));
	DoubleMuTrail21 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleMuTrail21"));
	DoubleMuTrail25 = (TGraphAsymmErrors*)(TriggerEfficiencies->Get("DoubleMuTrail25"));

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  qcd_a=HConfig.GetTH1D(Name+"_qcd_a","qcd_a",2,0.,2.,"number of entries");
  qcd_b=HConfig.GetTH1D(Name+"_qcd_b","qcd_b",2,0.,2.,"number of entries");
  qcd_c=HConfig.GetTH1D(Name+"_qcd_c","qcd_c",2,0.,2.,"number of entries");
  qcd_d=HConfig.GetTH1D(Name+"_qcd_d","qcd_d",2,0.,2.,"number of entries");
  wa=HConfig.GetTH1D(Name+"_wa","wa",2,0.,2.,"number of entries");
  wb=HConfig.GetTH1D(Name+"_wb","wb",2,0.,2.,"number of entries");
  wc=HConfig.GetTH1D(Name+"_wc","wc",2,0.,2.,"number of entries");
  wd=HConfig.GetTH1D(Name+"_wd","wd",2,0.,2.,"number of entries");
  qcd_control_a=HConfig.GetTH1D(Name+"_qcd_control_a","qcd_control_a",20,60.,120.,"m_{e,#mu} / GeV");
  qcd_control_b=HConfig.GetTH1D(Name+"_qcd_control_b","qcd_control_b",20,60.,120.,"m_{e,#mu} / GeV");
  qcd_control_c=HConfig.GetTH1D(Name+"_qcd_control_c","qcd_control_c",20,60.,120.,"m_{e,#mu} / GeV");
  qcd_control_d=HConfig.GetTH1D(Name+"_qcd_control_d","qcd_control_d",20,60.,120.,"m_{e,#mu} / GeV");
  qcd=HConfig.GetTH2D(Name+"_qcd","qcd",2,1,2,2,1,2,"opposite sign | same sign","#mu isolation < 0.12 | #mu isolation > 0.12");
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu_QCD::Store_ExtraDist(){
	Extradist1d.push_back(&qcd_a);
	Extradist1d.push_back(&qcd_b);
	Extradist1d.push_back(&qcd_c);
	Extradist1d.push_back(&qcd_d);
	Extradist1d.push_back(&wa);
	Extradist1d.push_back(&wb);
	Extradist1d.push_back(&wc);
	Extradist1d.push_back(&wd);
	Extradist1d.push_back(&qcd_control_a);
	Extradist1d.push_back(&qcd_control_b);
	Extradist1d.push_back(&qcd_control_c);
	Extradist1d.push_back(&qcd_control_d);

	Extradist2d.push_back(&qcd);
}

void  ZtoEMu_QCD::doEvent(){
  if(verbose)std::cout << "ZtoEMu_QCD::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}

  value.at(TriggerOk)=0;
  /*std::cout << "Number of HLT Triggers:" << Ntp->NHLTTriggers() << std::endl;
  for(unsigned int i=0; i<Ntp->NHLTTriggers(); i++){
	  std::cout << "######################################" << std::endl;
	  std::cout << "Trigger: " << Ntp->HTLTriggerName(i) << std::endl;
	  std::cout << "Acceptet (1) or not (0): " << Ntp->TriggerAccept(Ntp->HTLTriggerName(i)) << std::endl;
	  std::cout << "######################################" << std::endl;
  }*/
  if(Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") || Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")){
	  value.at(TriggerOk)=1;
  }
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
	  if(isGoodVtx(i)){
		  if(vertex==-1)vertex=i;
		  nGoodVtx++;
	  }
  }
  if(verbose)std::cout << "selected vertex: " << vertex << std::endl;
  value.at(PrimeVtx)=nGoodVtx;
  pass.at(PrimeVtx)=(value.at(PrimeVtx)>=cut.at(PrimeVtx));

  ///////////////////////////////////////////////
  //
  // Mu Cuts
  //
  if(verbose) std::cout << "Muon cuts" << std::endl;
  std::vector<unsigned int> GoodMuons;

  // muon ID cuts
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->Muon_p4(i).Pt()>mu_ptlow
			  && fabs(Ntp->Muon_p4(i).Eta())<mu_eta
			  && (matchTrigger(i,0.2,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","muon") || matchTrigger(i,0.2,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","muon"))
			  && vertex>=0
			  ){
		  if((Ntp->GetMCID()==31 || Ntp->GetMCID()==32) && !matchTruth(Ntp->Muon_p4(i),13,0.2)) continue;
		  if(Ntp->GetMCID()==33 && !matchTruth(Ntp->Muon_p4(i),13,0.2) && !matchTruth(Ntp->Muon_p4(i),15,0.2)) continue;
		  if(isTightMuon(i,vertex)){
			  GoodMuons.push_back(i);
		  }
	  }
  }

  value.at(NMu)=GoodMuons.size();
  pass.at(NMu)=(value.at(NMu)>=cut.at(NMu));

  unsigned int muidx(999);
  double hardestmu(0);
  if(GoodMuons.size()>1){
	  for(unsigned i=0;i<GoodMuons.size();i++){
		  if(Ntp->Muon_p4(GoodMuons.at(i)).Pt()>hardestmu){
			  hardestmu = Ntp->Muon_p4(GoodMuons.at(i)).Pt();
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
  bool hasMuonTrack = false;
  bool notThisOne = false;
  bool matchRecoMuon = false;

  // electron ID cuts (eta-dependent MVA or simple cut-based)
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(Ntp->Electron_p4(i).Et()>e_ptlow
			  && fabs(Ntp->Electron_supercluster_eta(i))<e_eta
			  && (matchTrigger(i,0.2,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","electron") || matchTrigger(i,0.2,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","electron"))
			  && vertex>=0
			  ){
		  for(unsigned j=0;j<GoodMuons.size();j++){
			  if(Ntp->Electron_Track_idx(i)==Ntp->Muon_Track_idx(GoodMuons.at(j))) hasMuonTrack = true;
		  }
		  if(hasMuonTrack) continue;
		  for(unsigned j=0;j<Ntp->NMuons();j++){
			  for(unsigned k=0;k<GoodMuons.size();k++){
				  notThisOne = false;
				  if(j==k){
					  notThisOne = true;
				  }
				  if(notThisOne)break;
			  }
			  if(!notThisOne && Ntp->Electron_p4(i).DeltaR(Ntp->Muon_p4(j))<0.3) matchRecoMuon = true;
		  }
		  if(matchRecoMuon) continue;
		  //if(isMVATrigElectron(i) && Electron_RelIso(i)<0.15) GoodElectrons.push_back(i);
		  if(Ntp->GetMCID()==31 && !matchTruth(Ntp->Electron_p4(i),11,0.2) && !matchTruth(Ntp->Electron_p4(i),22,0.2)) continue;
		  if(Ntp->GetMCID()==32 && !matchTruth(Ntp->Electron_p4(i),11,0.2) && !matchTruth(Ntp->Electron_p4(i),22,0.2) && !matchTruth(Ntp->Electron_p4(i),13,0.2)) continue;
		  if(Ntp->GetMCID()==33 && !matchTruth(Ntp->Electron_p4(i),11,0.2) && !matchTruth(Ntp->Electron_p4(i),13,0.2) && !matchTruth(Ntp->Electron_p4(i),22,0.2) && !matchTruth(Ntp->Electron_p4(i),15,0.2)) continue;
		  if(isWWElectron(i,vertex) && Electron_RelIso(i)<0.15) GoodElectrons.push_back(i);
	  }
  }

  value.at(NE)=GoodElectrons.size();
  pass.at(NE)=(value.at(NE)>=cut.at(NE));

  unsigned int eidx(999);
  double hardeste(0);
  if(GoodElectrons.size()>1){
	  for(unsigned i=0;i<GoodElectrons.size();i++){
		  if(Ntp->Electron_p4(GoodElectrons.at(i)).Et()>hardeste){
			  hardeste = Ntp->Electron_p4(GoodElectrons.at(i)).Et();
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
  bool leadingmu = false;
  value.at(ptthreshold)=0;
  if(muidx!=999 && eidx!=999){
	  value.at(ptthreshold)=1;
	  if(Ntp->Muon_p4(muidx).Pt()<=mu_ptlow || Ntp->Electron_p4(eidx).Et()<=e_ptlow) value.at(ptthreshold)=0;
	  if(Ntp->Muon_p4(muidx).Pt()<mu_pthigh && Ntp->Electron_p4(eidx).Et()<e_pthigh) value.at(ptthreshold)=0;
	  if(Ntp->Muon_p4(muidx).Pt()<mu_pthigh){
		  if(!Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) value.at(ptthreshold)=0;
	  }else if(Ntp->Electron_p4(eidx).Et()<e_pthigh){
		  if(!Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) value.at(ptthreshold)=0;
		  else leadingmu = true;
	  }
  }
  pass.at(ptthreshold)=(value.at(ptthreshold)==cut.at(ptthreshold));

  ///////////////////////////////////////////////
  //
  // dR(e,mu)
  //
  if(verbose) std::cout << "dR(e,mu)" << std::endl;
  value.at(drEMu)=0;
  if(muidx!=999 && eidx!=999){
	  value.at(drEMu)=Ntp->Muon_p4(muidx).DeltaR(Ntp->Electron_p4(eidx));
  }
  pass.at(drEMu)=(value.at(drEMu)>cut.at(drEMu));

  ///////////////////////////////////////////////
  //
  // Di Muon Veto
  //
  if(verbose)std::cout << "dimuon veto" << std::endl;
  unsigned int dimu(0);
  if(eidx!=999 && muidx!=999){
	  for(unsigned i=0;i<Ntp->NMuons();i++){
		  if(i==muidx) continue;
		  if(Ntp->Muon_p4(i).Pt()<3) continue;
		  if(fabs(Ntp->Muon_p4(i).Eta())>2.4) continue;
		  if(Ntp->Electron_p4(eidx).DeltaR(Ntp->Muon_p4(i))>0.3) continue;
		  dimu++;
	  }
  }
  value.at(diMuonVeto)=dimu;
  pass.at(diMuonVeto)=(value.at(diMuonVeto)==cut.at(diMuonVeto));

  ///////////////////////////////////////////////
  //
  // Tri Lepton Veto
  //
  if(verbose)std::cout << "trilepton veto" << std::endl;
  unsigned int trilep(0);
  if(eidx!=999 && muidx!=999){
	  for(unsigned i=0;i<Ntp->NMuons();i++){
		  if(i==muidx) continue;
		  if(vertex<0) continue;
		  if(Ntp->Muon_p4(i).Pt()<10) continue;
		  if(fabs(Ntp->Muon_p4(i).Eta())>2.4) continue;
		  if(!isTightMuon(i,vertex)) continue;
		  if(Muon_RelIso(i)>0.3) continue;
		  trilep++;
	  }
	  for(unsigned i=0;i<Ntp->NElectrons();i++){
		  if(i==eidx) continue;
		  if(vertex<0) continue;
		  if(Ntp->Electron_p4(i).Et()<10) continue;
		  if(fabs(Ntp->Electron_supercluster_eta(i))>2.5) continue;
		  //if(!isMVATrigElectron(i)) continue;
		  if(!isWWElectron(i,vertex)) continue;
		  if(Electron_RelIso(i)>0.3) continue;
		  trilep++;
	  }
  }
  value.at(triLeptonVeto)=trilep;
  pass.at(triLeptonVeto)=(value.at(triLeptonVeto)==cut.at(triLeptonVeto));

  ///////////////////////////////////////////////
  //
  // charge cut
  //
  if(verbose) std::cout << "Charge cut" << std::endl;
  value.at(charge)=999;
  if(eidx!=999 && muidx!=999){
	  value.at(charge)=Ntp->Electron_Charge(eidx)*Ntp->Muon_Charge(muidx);
  }
  pass.at(charge)=(value.at(charge)<cut.at(charge));

  ///////////////////////////////////////////////
  //
  // jet veto
  //
  if(verbose)std::cout << "jet veto" << std::endl;
  if(verbose)std::cout << "Cleaning jets" << std::endl;
  bool etrackjet(false);
  bool mutrackjet(false);
  std::vector<int> jetsfromvtx;

  /*if(doOldJetVeto){
	  // loop over all jets & only save the ones that do not overlap with the selected muon/electron
	  std::vector<int> jetidx;
	  if(eidx!=999 && muidx!=999){
		  for(unsigned i=0;i<Ntp->NPFJets();i++){
			  etrackjet=false;
			  mutrackjet=false;
			  for(unsigned j=0;j<Ntp->PFJet_nTrk(i);j++){
				  if(Ntp->PFJet_TracksP4(i,j).DeltaR(Ntp->Muon_p4(muidx))<0.001){
					  mutrackjet=true;
				  }
				  if(Ntp->PFJet_TracksP4(i,j).DeltaR(Ntp->Electron_p4(eidx))<0.1){
					  etrackjet=true;
				  }
			  }
			  if(Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(muidx))<0.3) mutrackjet = true;
			  if(Ntp->PFJet_p4(i).DeltaR(Ntp->Electron_p4(eidx))<0.3) etrackjet = true;
			  if(!mutrackjet && !etrackjet){
				  jetidx.push_back(i);
			  }
		  }
	  }

	  if(verbose)std::cout << "Looking for jets from Vtx" << std::endl;
	  int leadingtrack;
	  bool jetfromvtx;
	  std::vector<int> leadingtracks;

	  // loop over selected jets, find each leading track & check if it's assigned to selected vertex. save only jets from vertex
	  for(unsigned i=0;i<jetidx.size();i++){
		  leadingtrack = 0;
		  jetfromvtx = false;
		  int counter = 0;
		  for(unsigned j=0;j<Ntp->PFJet_nTrk(jetidx.at(i));j++){
			  if(Ntp->PFJet_TracksP4(jetidx.at(i),j).Pt()>Ntp->PFJet_TracksP4(jetidx.at(i),leadingtrack).Pt()){
				  leadingtrack=j;
			  }
		  }
		  if(Ntp->PFJet_nTrk(jetidx.at(i))==0)leadingtrack = -1;
		  if(vertex!=-1){
			  for(unsigned j=0;j<Ntp->Vtx_nTrk(vertex);j++){
				  if(leadingtrack>=0 && Ntp->PFJet_TracksP4(jetidx.at(i),leadingtrack).DeltaR(Ntp->Vtx_TracksP4(vertex,j))<0.00001){
					  jetfromvtx=true;
					  counter++;
				  }
			  }
		  }
		  if(counter>1)std::cout << "More than one vtx track associated to leading track from jet" << std::endl;
		  if(jetfromvtx){
			  jetsfromvtx.push_back(jetidx.at(i));
			  leadingtracks.push_back(leadingtrack);
		  }
	  }
  }else{
	  //using pileup jet id (loose wp)

	  // clean jets
	  for(unsigned int i=0;i<Ntp->NPFJets();i++){
		  if(eidx!=999 && muidx!=999){
			  etrackjet = false;
			  mutrackjet = false;
			  for(unsigned int j=0;j<Ntp->PFJet_Track_idx(i).size();j++){
				  if(Ntp->PFJet_Track_idx(i).at(j)==Ntp->Electron_Track_idx(eidx)) etrackjet = true;
				  if(Ntp->PFJet_Track_idx(i).at(j)==Ntp->Muon_Track_idx(muidx)) mutrackjet = true;
			  }
			  if(Ntp->PFJet_p4(i).DeltaR(Ntp->Electron_p4(eidx))<0.3) etrackjet = true;
			  if(Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(muidx))<0.3) mutrackjet = true;
		  }
		  if(etrackjet || mutrackjet) continue;
		  //if(Ntp->PFJet_PUJetID_looseWP(i)>0.5){ // loose
		  if(Ntp->PFJet_PUJetID_tightWP(i)>0.5){ // tight
			  jetsfromvtx.push_back(i);
		  }
	  }
  }*/

  if(verbose)std::cout << "Finding jets from vtx" << std::endl;
  for(unsigned i=0;i<Ntp->NPFJets();i++){
	  // clean jets against signal objects
	  if(Ntp->PFJet_p4(i).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(i).Eta(),Ntp->RunNumber())<20) continue;
	  if(fabs(Ntp->PFJet_p4(i).Eta())>jet_eta) continue;
	  if(muidx!=999){
		  mutrackjet = false;
		  if(Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(muidx))<0.3) mutrackjet = true;
	  }
	  if(eidx!=999){
		  etrackjet = false;
		  if(Ntp->PFJet_p4(i).DeltaR(Ntp->Electron_p4(eidx))<0.3) etrackjet = true;
	  }
	  if(etrackjet || mutrackjet) continue;
	  // find jets from vertex: use pileup jet id for jets with pt>20 GeV
	  if(Ntp->PFJet_PUJetID_tightWP(i)>0.5) jetsfromvtx.push_back(i);
  }

  if(verbose)std::cout<< "Find two highest pt jets" << std::endl;
  int firstjet_idx=-1;
  int secondjet_idx=-1;
  double initialpt=0.;

  // loop over jets from selected vertex & find the two jets with the highest pt
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(jetsfromvtx.at(i)).Eta(),Ntp->RunNumber())>initialpt){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(jetsfromvtx.at(i)).Eta(),Ntp->RunNumber());
		  firstjet_idx=jetsfromvtx.at(i);
	  }
  }
  initialpt=0.;
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(jetsfromvtx.size()>1 && firstjet_idx!=-1
			  && Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(jetsfromvtx.at(i)).Eta(),Ntp->RunNumber())>initialpt
			  && Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(jetsfromvtx.at(i)).Eta(),Ntp->RunNumber())<Ntp->PFJet_p4(firstjet_idx).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(firstjet_idx).Eta(),Ntp->RunNumber())
			  ){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(jetsfromvtx.at(i)).Eta(),Ntp->RunNumber());
		  secondjet_idx=jetsfromvtx.at(i);
	  }
  }

  if(verbose)std::cout << "applying veto" << std::endl;
  value.at(jetVeto)=0;
  if(jetsfromvtx.size()>1 && firstjet_idx!=-1 && secondjet_idx!=-1){
	  value.at(jetVeto)=Ntp->PFJet_p4(firstjet_idx).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(firstjet_idx).Eta(),Ntp->RunNumber())+Ntp->PFJet_p4(secondjet_idx).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(secondjet_idx).Eta(),Ntp->RunNumber());
  }else if(jetsfromvtx.size()==1 && firstjet_idx!=-1){
	  value.at(jetVeto)=Ntp->PFJet_p4(firstjet_idx).Pt()*rundependentJetPtCorrection(Ntp->PFJet_p4(firstjet_idx).Eta(),Ntp->RunNumber());
	  cut.at(jetVeto)=40;
  }
  /*int nb(0),nj(0);
  for(unsigned i=0;i<Ntp->NPFJets();i++){
	  if(Ntp->PFJet_bDiscriminator(i)>csvl) nb++;
	  if(Ntp->PFJet_PUJetID_looseWP(i)>0.5) nj++;
  }
  pass.at(jetVeto)=(nb==0 && nj<4);*/
  pass.at(jetVeto)=(value.at(jetVeto)<cut.at(jetVeto));

  ///////////////////////////////////////////////
  //
  // mtmu cut
  //
  if(verbose) std::cout << "Mt Mu cut" << std::endl;
  value.at(mtmu)=0;
  if(muidx!=999){
	  value.at(mtmu)=sqrt(2*Ntp->Muon_p4(muidx).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(muidx).Px(),Ntp->Muon_p4(muidx).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey())));
  }
  pass.at(mtmu)=(value.at(mtmu)<cut.at(mtmu));

  ///////////////////////////////////////////////
  //
  // Pt balance cut
  //
  if(verbose) std::cout << "pt balance cut" << std::endl;
  value.at(ptbal)=0.;
  if(muidx!=999
		  && eidx!=999
		  ){
	  value.at(ptbal) = (Ntp->Muon_p4(muidx)+Ntp->Electron_p4(eidx)).Pt();
  }
  pass.at(ptbal)=(value.at(ptbal)<cut.at(ptbal));

  ///////////////////////////////////////////////
  //
  // W+jets
  //
  if(verbose) std::cout << "w+jets" << std::endl;
  double wrate = 1.;
  /*if(Ntp->isData()){
	  if(!pass.at(mtmu)){
		  wrate = 3.15;
		  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::W_lnu,t)){ std::cout << "failed to find id "<< DataMCType::W_lnu <<std::endl; return;}
		  pass.at(mtmu)=true;
	  }
  }*/

  //////////////////////////////////////////////////////////
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1);
  if(!Ntp->isData() && Ntp->GetMCID()!=34){
    w*=Ntp->PUWeight();
    if(verbose)std::cout << "void  ZtoEMu_QCD::doEvent() k" << w << " " << wobs << std::endl;
    if(pass.at(NMu)){
    	w*=MuonIDeff(muidx);
    }
    if(pass.at(NE)){
    	w*=ElectronIDeff(eidx,"Trig");
    }
  }
  if(pass.at(TriggerOk)
  		&& pass.at(NMu)
  		&& pass.at(NE)
  		){
  	if(leadingmu) w*=TriggerEff(muidx,eidx,"Mu17_Ele8");
  	else w*=TriggerEff(muidx,eidx,"Mu8_Ele17");
  }
  else{w=1*wrate;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;

  // qcd abcd plots
  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(NMu)
		  && pass.at(NE)
		  && pass.at(ptthreshold)
		  && pass.at(diMuonVeto)
		  && pass.at(triLeptonVeto)
		  //&& pass.at(jetVeto)
		  //&& pass.at(mtmu)
		  //&& pass.at(ptbal)
		  ){
	  double m = (Ntp->Muon_p4(muidx)+Ntp->Electron_p4(eidx)).M();
	  if(m>60. && m<120.){
		  if(pass.at(charge) && Muon_RelIso(muidx)<0.12){
			  qcd_a.at(t).AddBinContent(1);
			  qcd_control_a.at(t).Fill(m,w);
		  }
		  if(!pass.at(charge) && Muon_RelIso(muidx)<0.12){
			  qcd_b.at(t).AddBinContent(1);
			  qcd_control_b.at(t).Fill(m,w);
		  }
		  if(pass.at(charge) && Muon_RelIso(muidx)>0.2 && Muon_RelIso(muidx)<0.5){
			  qcd_c.at(t).AddBinContent(1);
			  qcd_control_c.at(t).Fill(m,w);
		  }
		  if(!pass.at(charge) && Muon_RelIso(muidx)>0.2 && Muon_RelIso(muidx)<0.5){
			  qcd_d.at(t).AddBinContent(1);
			  qcd_control_d.at(t).Fill(m,w);
		  }
	  }
  }

  if(verbose)std::cout << "ZtoEMu_QCD::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

bool ZtoEMu_QCD::isGoodVtx(unsigned int i){
	if(fabs(Ntp->Vtx(i).z())>=24) return false;
	if(Ntp->Vtx(i).Perp()>=2) return false;
	if(Ntp->Vtx_ndof(i)<=4) return false;
	if(Ntp->Vtx_isFake(i)!=0) return false;
	return true;
}

bool ZtoEMu_QCD::jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx){
	for(unsigned i=0;i<vtx_track_idx.size();i++){
		if(vtx_track_idx.at(i)==leadingtrack_idx)return true;
	}
	return false;
}

double ZtoEMu_QCD::rundependentJetPtCorrection(double jeteta, int runnumber){
	if(!Ntp->isData()) return 1.;
	const double corrs[5] = {0.0, -0.454e-6, -0.952e-6, 1.378e-6, 0.0};
	const int run0 = 201000;
	double eta = fabs(jeteta);
	double corr = 0.;
	if(eta<1.3) corr = corrs[0];
	else if(eta<2.0) corr = corrs[1];
	else if(eta<2.5) corr = corrs[2];
	else if(eta<3.) corr = corrs[3];
	else if(eta<5.) corr = corrs[4];
	return (1.+corr*(runnumber-run0));
}

double ZtoEMu_QCD::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double ZtoEMu_QCD::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

double ZtoEMu_QCD::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}

double ZtoEMu_QCD::vertexSignificance(TVector3 vec, unsigned int vertex){
	if(vertex>=0 && vertex<Ntp->NVtx()){
		const float elm[3] = {(vec.X()-Ntp->Vtx(vertex).X()),(vec.Y()-Ntp->Vtx(vertex).Y()),(vec.Z()-Ntp->Vtx(vertex).Z())};
		TVectorF diff(3,elm);
		TMatrixF M(Ntp->Vtx_Cov(vertex));
		if(M.IsValid()){
			double mag = diff.Norm2Sqr();
			double sim = M.Similarity(diff);
			return mag/sqrt(sim);
		}
	}
	return 999;
}

bool ZtoEMu_QCD::matchTrigger(unsigned int i, double dr, std::string trigger, std::string object){
	unsigned int id = 0;
	TLorentzVector particle(0.,0.,0.,0.);
	TLorentzVector triggerObj(0.,0.,0.,0.);
	if(object=="electron"){
		id = 82;
		particle = Ntp->Electron_p4(i);
	}
	if(object=="muon"){
		id = 83;
		particle = Ntp->Muon_p4(i);
	}
	for(unsigned i=0;i<Ntp->NHLTTrigger_objs();i++){
		if(Ntp->HLTTrigger_objs_trigger(i).find(trigger) != string::npos){
			for(unsigned j=0;j<Ntp->NHLTTrigger_objs(i);j++){
				if(Ntp->HLTTrigger_objs_Id(i,j)==id){
					triggerObj.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(i,j),
							Ntp->HLTTrigger_objs_Eta(i,j),
							Ntp->HLTTrigger_objs_Phi(i,j),
							Ntp->HLTTrigger_objs_E(i,j));
				}
				if(triggerObj.Pt()>0.
						&& particle.Pt()>0.
						&& particle.DeltaR(triggerObj)<dr) return true;
			}
		}
	}
	return false;
}

int ZtoEMu_QCD::matchTruth(TLorentzVector tvector){
	double testdr=1.;
	int pdgid = 0;
	for(unsigned i=0;i<Ntp->NMCParticles();i++){
		if(Ntp->MCParticle_p4(i).Pt()>0.){
			if(tvector.DeltaR(Ntp->MCParticle_p4(i))<testdr){
				testdr = tvector.DeltaR(Ntp->MCParticle_p4(i));
				pdgid = Ntp->MCParticle_pdgid(i);
			}
		}
	}
	return pdgid;
}

bool ZtoEMu_QCD::matchTruth(TLorentzVector tvector, int pid, double dr){
	for(unsigned i=0;i<Ntp->NMCParticles();i++){
		if(Ntp->MCParticle_p4(i).Pt()>0.){
			if(fabs(Ntp->MCParticle_pdgid(i))==pid){
				if(tvector.DeltaR(Ntp->MCParticle_p4(i))<dr) return true;
			}
		}
	}
	return false;
}

//////////////////////////////
//
// Muon related functions
//

bool ZtoEMu_QCD::isTightMuon(unsigned int i){
	if(!Ntp->Muon_isGlobalMuon(i)) return false;
	if(!Ntp->Muon_isPFMuon(i)) return false;
	if(Ntp->Muon_normChi2(i)>=10.) return false;
	if(Ntp->Muon_hitPattern_numberOfValidMuonHits(i)<=0) return false;
	if(Ntp->Muon_numberOfMatchedStations(i)<=1) return false;
	if(Ntp->Muon_numberofValidPixelHits(i)<=0) return false;
	if(Ntp->Muon_trackerLayersWithMeasurement(i)<=5) return false;
	return true;
}

bool ZtoEMu_QCD::isTightMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.2) return false;
	if(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.5) return false;
	return true;
}

bool ZtoEMu_QCD::isHiggsMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	if(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
}

bool ZtoEMu_QCD::isLooseMuon(unsigned int i){
	if(!Ntp->Muon_isPFMuon(i)) return false;
	if(!(Ntp->Muon_isGlobalMuon(i) || Ntp->Muon_isTrackerMuon(i))) return false;
	return true;
}

bool ZtoEMu_QCD::isFakeMuon(unsigned int i){
	if(!Ntp->Muon_isGlobalMuon(i)) return false;
	if(Ntp->Muon_p4(i).Pt()<=10) return false;
	if(fabs(Ntp->Muon_p4(i).Eta())>2.4) return false;
	if(Ntp->Muon_p4(i).Pt()<=20){
		if(Ntp->Muon_sumPt03(i)>=8.) return false;
		if(Ntp->Muon_emEt03(i)>=8.) return false;
		if(Ntp->Muon_hadEt03(i)>=8.) return false;
	}
	if(Ntp->Muon_p4(i).Pt()>20){
		if(Ntp->Muon_sumPt03(i)/Ntp->Muon_p4(i).Pt()>=0.4) return false;
		if(Ntp->Muon_emEt03(i)/Ntp->Muon_p4(i).Pt()>=0.4) return false;
		if(Ntp->Muon_hadEt03(i)/Ntp->Muon_p4(i).Pt()>=0.4) return false;
	}
	return true;
}

bool ZtoEMu_QCD::isFakeMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isFakeMuon(i)) return false;
	if(dxy(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.2) return false;
	if(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
}

double ZtoEMu_QCD::Muon_RelIso(unsigned int i){
	return (Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i)))/Ntp->Muon_p4(i).Pt();
}

double ZtoEMu_QCD::Muon_AbsIso(unsigned int i){
	return Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i));
}

//////////////////////////////
//
// Electron related functions
//

bool ZtoEMu_QCD::isTrigPreselElectron(unsigned int i){
	if(fabs(Ntp->Electron_supercluster_eta(i))>2.5) return false;
	if(fabs(Ntp->Electron_supercluster_eta(i))<1.479){
		if(Ntp->Electron_sigmaIetaIeta(i)>0.014) return false;
		if(Ntp->Electron_hadronicOverEm(i)>0.15) return false;
		if(Ntp->Electron_Gsf_dr03TkSumPt(i)/Ntp->Electron_p4(i).Pt()>0.2) return false;
		if(Ntp->Electron_Gsf_dr03HcalTowerSumEt(i)/Ntp->Electron_p4(i).Pt()>0.2) return false;
		if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	}else{
		if(Ntp->Electron_sigmaIetaIeta(i)>0.035) return false;
		if(Ntp->Electron_hadronicOverEm(i)>0.1) return false;
		if(Ntp->Electron_Gsf_dr03TkSumPt(i)/Ntp->Electron_p4(i).Pt()>0.2) return false;
		if(Ntp->Electron_Gsf_dr03HcalTowerSumEt(i)/Ntp->Electron_p4(i).Pt()>0.2) return false;
		if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	}
	return true;
}

bool ZtoEMu_QCD::isTrigNoIPPreselElectron(unsigned int i){
	if(fabs(Ntp->Electron_supercluster_eta(i))>2.5) return false;
	if(fabs(Ntp->Electron_supercluster_eta(i))<1.479){
		if(Ntp->Electron_sigmaIetaIeta(i)>0.01) return false;
		if(Ntp->Electron_hadronicOverEm(i)>0.12) return false;
		if(fabs(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i))>0.007) return false;
		if(fabs(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i))>0.15) return false;
		if(Ntp->Electron_Gsf_dr03TkSumPt(i)/Ntp->Electron_p4(i).Pt()>0.2) return false;
		if(Ntp->Electron_Gsf_dr03HcalTowerSumEt(i)/Ntp->Electron_p4(i).Pt()>0.2) return false;
		if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	}else{
		if(Ntp->Electron_sigmaIetaIeta(i)>0.03) return false;
		if(Ntp->Electron_hadronicOverEm(i)>0.1) return false;
		if(fabs(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i))>0.009) return false;
		if(fabs(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i))>0.1) return false;
		if(Ntp->Electron_Gsf_dr03TkSumPt(i)/Ntp->Electron_p4(i).Pt()>0.2) return false;
		if(Ntp->Electron_Gsf_dr03HcalTowerSumEt(i)/Ntp->Electron_p4(i).Pt()>0.2) return false;
		if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	}
	return true;
}

bool ZtoEMu_QCD::isMVATrigNoIPElectron(unsigned int i){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	if(!isTrigNoIPPreselElectron(i)) return false;
	if(mvapt<20){
		if(mvaeta<0.8){
			if(Electron_RelIso(i)>=0.15) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=-0.5375) return false;
		}
		if(mvaeta>=0.8 && mvaeta<1.479){
			if(Electron_RelIso(i)>=0.15) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=-0.375) return false;
		}
		if(mvaeta>=1.479 && mvaeta<2.5){
			if(Electron_RelIso(i)>=0.10) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=-0.025) return false;
		}
	}
	if(mvapt>=20){
		if(mvaeta<0.8){
			if(Electron_RelIso(i)>=0.15) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=0.325) return false;
		}
		if(mvaeta>=0.8 && mvaeta<1.479){
			if(Electron_RelIso(i)>=0.15) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
		}
		if(mvaeta>=1.479 && mvaeta<2.5){
			if(Electron_RelIso(i)>=0.10) return false;
			if(Ntp->Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
		}
	}
	return true;
}

bool ZtoEMu_QCD::isMVANonTrigElectron(unsigned int i, unsigned int j){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_numberOfMissedHits(i)>1) return false;
	if(vertexSignificance(Ntp->Electron_Poca(i),j)>=4) return false;
	if(mvapt>7. && mvapt<10.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.47) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.004) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.295) return false;
	}
	if(mvapt>=10.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=-0.34) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=-0.65) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.6) return false;
	}
	return true;
}

bool ZtoEMu_QCD::isHiggsElectron(unsigned int i, unsigned int j){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	if(dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	if(mvapt<20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.925) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.915) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.965) return false;
	}
	if(mvapt>=20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.905) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.955) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.975) return false;
	}
	return true;
}

bool ZtoEMu_QCD::isWWElectron(unsigned int i, unsigned int j){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(!isFakeElectron(i,j)) return false;
	if(mvapt>10. && mvapt<20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.00) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.10) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.62) return false;
	}else if(mvapt>=20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.94) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.85) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.92) return false;
	}
	return true;
}

bool ZtoEMu_QCD::isMVATrigElectron(unsigned int i){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(!isTrigPreselElectron(i)) return false;
	if(mvapt>10. && mvapt<20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.00) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.10) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.62) return false;
	}else if(mvapt>=20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.94) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.85) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_Trig_discriminator(i)<=0.92) return false;
	}
	return true;
}

bool ZtoEMu_QCD::isTightElectron(unsigned int i){
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_RelIso(i)>=0.1) return false;
	if(fabs(1/Ntp->Electron_ecalEnergy(i)-1/Ntp->Electron_trackMomentumAtVtx(i))>=0.05) return false;
	if(fabs(Ntp->Electron_supercluster_eta(i))<=1.479){
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.004) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.03) return false;
		if(Ntp->Electron_sigmaIetaIeta(i)>=0.01) return false;
		if(Ntp->Electron_hadronicOverEm(i)>=0.12) return false;
	}
	if(fabs(Ntp->Electron_supercluster_eta(i))>1.479 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.005) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.02) return false;
		if(Ntp->Electron_sigmaIetaIeta(i)>=0.03) return false;
		if(Ntp->Electron_hadronicOverEm(i)>=0.10) return false;
		if(Ntp->Electron_p4(i).Pt()<20 && Electron_RelIso(i)>=0.07) return false;
	}
	return true;
}

bool ZtoEMu_QCD::isTightElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightElectron(i)) return false;
	if(dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	if(dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
}

bool ZtoEMu_QCD::isFakeElectron(unsigned int i){
	if(Ntp->Electron_p4(i).Pt()<10) return false;
	if(fabs(Ntp->Electron_supercluster_eta(i))>2.5) return false;
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	if(Ntp->Electron_tkSumPt03(i)/Ntp->Electron_p4(i).Pt()>0.2) return false;
	if(std::max(Ntp->Electron_ecalRecHitSumEt03(i)-1.,0.)/Ntp->Electron_p4(i).Pt()>0.2) return false;
	if((Ntp->Electron_hcalDepth1TowerSumEt03(i)+Ntp->Electron_hcalDepth2TowerSumEt03(i))/Ntp->Electron_p4(i).Pt()>0.2) return false;
	if(fabs(Ntp->Electron_supercluster_eta(i))<1.479){
		if(Ntp->Electron_sigmaIetaIeta(i)>0.01) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>0.15) return false;
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>0.007) return false;
		if(Ntp->Electron_hadronicOverEm(i)>0.12) return false;
	}
	if(fabs(Ntp->Electron_supercluster_eta(i))>=1.479 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){
		if(Ntp->Electron_sigmaIetaIeta(i)>0.03) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>0.10) return false;
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>0.009) return false;
		if(Ntp->Electron_hadronicOverEm(i)>0.10) return false;
	}
	return true;
}

bool ZtoEMu_QCD::isFakeElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isFakeElectron(i)) return false;
	if(dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>0.1) return false;
	if(dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>0.02) return false;
	return true;
}

double ZtoEMu_QCD::Electron_RelIso(unsigned int i){
	return (Ntp->Electron_chargedHadronIso(i)+std::max((double)0.,Ntp->Electron_neutralHadronIso(i)+Ntp->Electron_photonIso(i)-Ntp->RhoIsolationAllInputTags()*Electron_Aeff_R04(Ntp->Electron_supercluster_eta(i))))/Ntp->Electron_p4(i).Pt();
}

double ZtoEMu_QCD::Electron_Aeff_R04(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.208;
	else if(eta>=1. && eta<1.479) return 0.209;
	else if(eta>=1.479 && eta<2.) return 0.115;
	else if(eta>=2. && eta<2.2) return 0.143;
	else if(eta>=2.2 && eta<2.3) return 0.183;
	else if(eta>=2.3 && eta<2.3) return 0.194;
	else if(eta>=2.4) return 0.261;
}

double ZtoEMu_QCD::Electron_Aeff_R03(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.13;
	else if(eta>=1. && eta<1.479) return 0.14;
	else if(eta>=1.479 && eta<2.) return 0.07;
	else if(eta>=2. && eta<2.2) return 0.09;
	else if(eta>=2.2 && eta<2.3) return 0.11;
	else if(eta>=2.3 && eta<2.3) return 0.11;
	else if(eta>=2.4) return 0.14;
}

//////////////////////////////
//
// Trigger & ID efficiencies
//

double ZtoEMu_QCD::MuonIDeff(unsigned int i){
	double pt = Ntp->Muon_p4(i).Pt();
	double eta = fabs(Ntp->Muon_p4(i).Eta());
	double eff = 1.;
	if(pt<20) pt=20;
	if(pt>100) pt=100;
	if(eta<0.9) eff = MuIdEff09->Eval(pt)*MuIsoEff09->Eval(pt);
	if(eta>=0.9 && eta<1.2) eff = MuIdEff12->Eval(pt)*MuIsoEff12->Eval(pt);
	if(eta>=1.2 && eta<2.1) eff = MuIdEff21->Eval(pt)*MuIsoEff21->Eval(pt);
	if(eta>=2.1 && eta<2.4) eff = MuIdEff24->Eval(pt)*MuIsoEff24->Eval(pt);
	return eff;
}

double ZtoEMu_QCD::MuonHiggsIDeff(unsigned int i){
	double pt = Ntp->Muon_p4(i).Pt();
	double eta = fabs(Ntp->Muon_p4(i).Eta());
	double eff = 1.;
	if(pt>10 && pt<=15){
		if(eta>=0 && eta<0.8){
			eff = 0.9771;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9746;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9644;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9891;
		}
	}else if(pt>15 && pt<=20){
		if(eta>=0 && eta<0.8){
			eff = 0.9548;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9701;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9766;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9892;
		}
	}else if(pt>20 && pt<=25){
		if(eta>=0 && eta<0.8){
			eff = 0.9648;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9836;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9820;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9909;
		}
	}else if(pt>25 && pt<=30){
		if(eta>=0 && eta<0.8){
			eff = 0.9676;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9817;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9886;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9883;
		}
	}else if(pt>30 && pt<=35){
		if(eta>=0 && eta<0.8){
			eff = 0.9883;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9833;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9910;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9900;
		}
	}else if(pt>35){
		if(eta>=0 && eta<0.8){
			eff = 0.9826;
		}else if(eta>=0.8 && eta<1.2){
			eff = 0.9841;
		}else if(eta>=1.2 && eta<1.6){
			eff = 0.9900;
		}else if(eta>=1.6 && eta<2.1){
			eff = 0.9886;
		}
	}
	return eff;
}

double ZtoEMu_QCD::ElectronIDeff(unsigned int i, std::string id){
	if(id=="Trig") return ElectronTrigIDeff(i);
	if(id=="NonTrig") return ElectronNonTrigIDeff(i);
	if(id=="Higgs") return ElectronHiggsIDeff(i);
	return 1.;
}

double ZtoEMu_QCD::ElectronTrigIDeff(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	double eff = 1.;
	if(eta<2.5){
		eff = ElectronTrigEff->GetBinContent(ElectronTrigEff->FindFixBin(eta,pt));
	}
	return eff;
}

double ZtoEMu_QCD::ElectronNonTrigIDeff(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = Ntp->Electron_supercluster_eta(i);
	double eff = 1.;
	if(fabs(eta)<2.5){
		eff = ElectronNonTrigEff->GetBinContent(ElectronNonTrigEff->FindFixBin(pt,eta));
	}
	return eff;
}

double ZtoEMu_QCD::ElectronHiggsIDeff(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = Ntp->Electron_supercluster_eta(i);
	double eff = 1.;
	if(fabs(eta)<2.5){
		if(pt>10 && pt<=15){
			if(eta>=0 && eta<0.8){
				eff = 0.7654;
			}else if(eta>=0.8 && eta<1.5){
				eff = 0.7693;
			}else if(eta>=1.5 && eta<2.3){
				eff = 0.5719;
			}
		}else if(pt>15 && pt<=20){
			if(eta>=0 && eta<0.8){
				eff = 0.8394;
			}else if(eta>=0.8 && eta<1.5){
				eff = 0.8457;
			}else if(eta>=1.5 && eta<2.3){
				eff = 0.7024;
			}
		}else if(pt>20 && pt<=25){
			if(eta>=0 && eta<0.8){
				eff = 0.8772;
			}else if(eta>=0.8 && eta<1.5){
				eff = 0.8530;
			}else if(eta>=1.5 && eta<2.3){
				eff = 0.7631;
			}
		}else if(pt>25 && pt<=30){
			if(eta>=0 && eta<0.8){
				eff = 0.9006;
			}else if(eta>=0.8 && eta<1.5){
				eff = 0.8874;
			}else if(eta>=1.5 && eta<2.3){
				eff = 0.8092;
			}
		}else if(pt>30 && pt<=35){
			if(eta>=0 && eta<0.8){
				eff = 0.9261;
			}else if(eta>=0.8 && eta<1.5){
				eff = 0.9199;
			}else if(eta>=1.5 && eta<2.3){
				eff = 0.8469;
			}
		}else if(pt>35){
			if(eta>=0 && eta<0.8){
				eff = 0.9514;
			}else if(eta>=0.8 && eta<1.5){
				eff = 0.9445;
			}else if(eta>=1.5 && eta<2.3){
				eff = 0.9078;
			}
		}
	}
	return eff;
}

double ZtoEMu_QCD::TriggerEff(unsigned int muid, unsigned int eid, TString path){
	double eff = 1.;
	if(path.Contains("Mu17_Ele8")){
		eff = SingleMu(muid) + (1-SingleMu(muid))*SingleEle(eid)
				+ (1-SingleMu(muid))*(1-SingleEle(eid))*
				(DoubleMuLeading(muid)*DoubleEleTrailing(eid) + (1-DoubleMuLeading(muid)*DoubleEleTrailing(eid))*DoubleEleLeading(eid)*DoubleMuTrailing(muid));
	}
	if(path.Contains("Mu8_Ele17")){
		eff = SingleEle(eid) + (1-SingleEle(eid))*SingleMu(muid)
				+ (1-SingleEle(eid))*(1-SingleMu(muid))*
				(DoubleEleLeading(eid)*DoubleMuTrailing(muid) + (1-DoubleEleLeading(eid)*DoubleMuTrailing(muid))*DoubleMuLeading(muid)*DoubleEleTrailing(eid));
	}
	return eff;
}

double ZtoEMu_QCD::SingleEle(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	double eff = 1.;
	if(eta<1.5) eff = SingleEle15->Eval(pt);
	else if(eta<2.5) eff = SingleEle25->Eval(pt);
	return eff;
}

double ZtoEMu_QCD::DoubleEleLeading(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	double eff = 1.;
	if(eta<1.5) eff = DoubleEleLead15->Eval(pt);
	else if(eta<2.5) eff = DoubleEleLead25->Eval(pt);
	return eff;
}

double ZtoEMu_QCD::DoubleEleTrailing(unsigned int i){
	double pt = Ntp->Electron_p4(i).Pt();
	double eta = fabs(Ntp->Electron_supercluster_eta(i));
	double eff = 1.;
	if(eta<1.5) eff = DoubleEleTrail15->Eval(pt);
	else if(eta<2.5) eff = DoubleEleTrail25->Eval(pt);

	return eff;
}

double ZtoEMu_QCD::SingleMu(unsigned int i){
	double pt = Ntp->Muon_p4(i).Pt();
	double eta = fabs(Ntp->Muon_p4(i).Eta());
	double eff = 1.;
	if(eta<0.8) eff = SingleMu08->Eval(pt);
	else if(eta<1.2) eff = SingleMu12->Eval(pt);
	else if(eta<2.1) eff = SingleMu21->Eval(pt);
	else if(eta<2.5) eff = SingleMu25->Eval(pt);
	return eff;
}

double ZtoEMu_QCD::DoubleMuLeading(unsigned int i){
	double pt = Ntp->Muon_p4(i).Pt();
	double eta = fabs(Ntp->Muon_p4(i).Eta());
	double eff = 1.;
	if(eta<1.2) eff = DoubleMuLead12->Eval(pt);
	else if(eta<2.1) eff = DoubleMuLead21->Eval(pt);
	else if(eta<2.5) eff = DoubleMuLead25->Eval(pt);
	return eff;
}

double ZtoEMu_QCD::DoubleMuTrailing(unsigned int i){
	double pt = Ntp->Muon_p4(i).Pt();
	double eta = fabs(Ntp->Muon_p4(i).Eta());
	double eff = 1.;
	if(eta<1.2)	eff = DoubleMuTrail12->Eval(pt);
	else if(eta<2.1) eff = DoubleMuTrail21->Eval(pt);
	else if(eta<2.5) eff = DoubleMuTrail25->Eval(pt);
	return eff;
}

//////////////////////////////
//
// Calculate fakerate
//

double ZtoEMu_QCD::Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type){

	double eta1, eta2;
	int ptbin = 0;
	int etabin = 0;
	double FakePt = vec.Pt();
	double FakeEta = vec.Eta();
	double fakerate = 0.;

	if(type=="muon"){
		eta1 = 2.1;
		eta2 = 1.2;
	}else if(type=="electron"){
		eta1 = 2.;
		eta2 = 1.479;
	}

	if(FakePt<15.)ptbin = 1;
	else if(FakePt>=15. && FakePt<20.)ptbin = 2;
	else if(FakePt>=20. && FakePt<25.)ptbin = 3;
	else if(FakePt>=25. && FakePt<30.)ptbin = 4;
	else if(FakePt>=30.)ptbin = 5;

	if(FakeEta<-eta1)etabin = 1;
	else if(FakeEta>=-eta1 && FakeEta<-eta2)etabin = 2;
	else if(FakeEta>=-eta2 && FakeEta<-0.8)etabin = 3;
	else if(FakeEta>=-0.8 && FakeEta<0.)etabin = 4;
	else if(FakeEta>=0. && FakeEta<0.8)etabin = 5;
	else if(FakeEta>=0.8 && FakeEta<eta2)etabin = 6;
	else if(FakeEta>=eta2 && FakeEta<eta1)etabin = 7;
	else if(FakeEta>=eta1)etabin = 8;

	if(ptbin==0 || etabin==0){
		fakerate = 0;
	}else{
		fakerate = fakeRateHist->GetBinContent(ptbin,etabin);
	}

	return fakerate;//(1-fakerate);
}

//////////////////////////////
//
// Finish function
//

void ZtoEMu_QCD::Finish(){
	unsigned int nhistos = 14;

	double lumi = 19712.;
	double xdytautau = 1966.7;
	double xdyee = 1966.7;
	double xdymumu = 1966.7;
	double xww = 5.824;
	double xtt = 245.8;
	double xtw = 11.1;
	double xtbarw = 11.1;
	double xwz3lnu = 1.058;
	double xwz2l2q = 2.207;
	double xzz4l = 0.181;
	double xzz2l2q = 2.502;
	double xzz2l2nu = 0.716;

	int ndytautau = 48790562;
	int ndyee = 3297045;
	int ndymumu = 3293740;
	int nww = 1933235;
	int ntt = 21675970;
	int ntw = 497658;
	int ntbarw = 493460;
	int nwz3lnu = 2017254;
	int nwz2l2q = 3212348;
	int nzz4l = 3854021;
	int nzz2l2q = 1577042;
	int nzz2l2nu = 777964;

	double nwjets = 3.15;

	std::vector<double> scale;
	scale.push_back(1.);
	scale.push_back(nwjets);
	scale.push_back(lumi*xww/nww);
	scale.push_back(lumi*xwz2l2q/nwz2l2q);
	scale.push_back(lumi*xwz3lnu/nwz3lnu);
	scale.push_back(lumi*xzz4l/nzz4l);
	scale.push_back(lumi*xzz2l2nu/nzz2l2nu);
	scale.push_back(lumi*xzz2l2q/nzz2l2q);
	scale.push_back(lumi*xtt/ntt);
	scale.push_back(lumi*xtw/ntw);
	scale.push_back(lumi*xtbarw/ntbarw);
	scale.push_back(lumi*xdytautau/ndytautau);
	scale.push_back(lumi*xdymumu/ndymumu);
	scale.push_back(lumi*xdyee/ndyee);

	wa.at(0).SetBinContent(1,qcd_a.at(0).GetBinContent(1));
	wb.at(0).SetBinContent(1,qcd_b.at(0).GetBinContent(1));
	wc.at(0).SetBinContent(1,qcd_c.at(0).GetBinContent(1));
	wd.at(0).SetBinContent(1,qcd_d.at(0).GetBinContent(1));
	for(unsigned int i=1;i<nhistos;i++){
		wa.at(0).SetBinContent(1,wa.at(0).GetBinContent(1)-qcd_a.at(i).GetBinContent(1)*scale.at(i));
		wb.at(0).SetBinContent(1,wb.at(0).GetBinContent(1)-qcd_b.at(i).GetBinContent(1)*scale.at(i));
		wc.at(0).SetBinContent(1,wc.at(0).GetBinContent(1)-qcd_c.at(i).GetBinContent(1)*scale.at(i));
		wd.at(0).SetBinContent(1,wd.at(0).GetBinContent(1)-qcd_d.at(i).GetBinContent(1)*scale.at(i));
	}
	qcd.at(0).SetBinContent(1,1,wa.at(0).Integral());
	qcd.at(0).SetBinContent(2,1,wb.at(0).Integral());
	qcd.at(0).SetBinContent(1,2,wc.at(0).Integral());
	qcd.at(0).SetBinContent(2,2,wd.at(0).Integral());

	std::cout << "events in region a: " << qcd.at(0).GetBinContent(1,1) << std::endl;
	std::cout << "events in region b: " << qcd.at(0).GetBinContent(2,1) << std::endl;
	std::cout << "events in region c: " << qcd.at(0).GetBinContent(1,2) << std::endl;
	std::cout << "events in region d: " << qcd.at(0).GetBinContent(2,2) << std::endl;
	std::cout << "Use ratio of " << (double)qcd.at(0).GetBinContent(1,2)/(double)qcd.at(0).GetBinContent(2,2) << " to extrapolate shape of region b to region a" << std::endl;
	std::cout << "Uncertainty on ratio: " << TMath::Sqrt(qcd.at(0).GetBinContent(1,2)/TMath::Power(qcd.at(0).GetBinContent(2,2),2)+TMath::Power(qcd.at(0).GetBinContent(1,2),2)/TMath::Power(qcd.at(0).GetBinContent(2,2),3)) << std::endl;

	Selection::Finish();
}
