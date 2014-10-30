#include "ZtoEMu.h"
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

ZtoEMu::ZtoEMu(TString Name_, TString id_):
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
  ,normunc_dy(0.033)
  ,normunc_tt(0.047)
  ,normunc_tw(0.09)
  ,normunc_diboson(0.15)
  ,normunc_qcd(0.387)
{
    //verbose=true;

	doHiggsObjects = false;
	doWWObjects = true;
	useMadgraphZ = true;
	doPDFuncertainty = false;
	if(useMadgraphZ) mmin = 50;
	if(doHiggsObjects){
		mu_eta = 2.1;
		e_eta = 2.3;
	}

	// corrections to be applied to particle candidates
	mucorr = "roch";
	ecorr = "";
	jetcorr = "runJER";

	// initialize pdf reweighting for systematics
	if(doPDFuncertainty){
		pdfname1 = "cteq6ll.LHgrid";
		pdfname2 = "MSTW2008nnlo68cl.LHgrid"; // CT10nnlo.LHgrid works only locally for some reason
		//pdfname2 = "NNPDF23_nnlo_as_0119.LHgrid";
		if(Ntp->GetMCID()==DataMCType::ttbar){ // todo: correct pdfs?
			pdfname1 = "CT10.LHgrid";
			pdfname2 = "MSTW2008nnlo68cl.LHgrid";
			//pdfname2 = "NNPDF23_nnlo_as_0119.LHgrid";
		}
		pdf = new PDFweights(pdfname1,pdfname2);
		nPDFmembers = pdf->numberOfMembers();
	}
}

ZtoEMu::~ZtoEMu(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZtoEMu::~ZtoEMu Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZtoEMu::~ZtoEMu()" << std::endl;
}

void  ZtoEMu::Configure(){

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
    if(i==triLeptonVeto)      cut.at(triLeptonVeto)=0;
    if(i==charge)             cut.at(charge)=0;
    if(i==oneJet)             cut.at(oneJet)=singlejet;
    if(i==MtMu)               cut.at(MtMu)=mtmu;
    if(i==ptBalance)          cut.at(ptBalance)=ptbalance;
    if(i==ZMassmax)           cut.at(ZMassmax)=zmax;
    if(i==ZMassmin)           cut.at(ZMassmin)=zmin;
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
    else if(i==oneJet){
  	  title.at(i)="jet pt < ";
  	  char buffer[50];
  	  sprintf(buffer,"%5.2f",cut.at(oneJet));
  	  title.at(i)+=buffer;
  	  htitle=title.at(i);
  	  htitle.ReplaceAll("$","");
  	  htitle.ReplaceAll("\\","#");
  	  hlabel="jet P_{T} / GeV";
  	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_oneJet_",htitle,40,0,200,hlabel,"Events"));
  	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_oneJet_",htitle,40,0,200,hlabel,"Events"));
	}
    else if(i==MtMu){
	  title.at(i)="$m_{T}^{\\mu,Miss} < $";
	  title.at(i)+=cut.at(MtMu);
	  title.at(i)+="(GeV)";
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="m_{T}^{#mu,Miss} / GeV";
	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MtMu_",htitle,40,0,200,hlabel,"Events"));
	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MtMu_",htitle,40,0,200,hlabel,"Events"));
    }
    else if(i==ptBalance){
      title.at(i)="$p_{t,e+\\mu} < $";
      title.at(i)+=cut.at(ptBalance);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="p_t balance / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ptBalance_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ptBalance_",htitle,40,0,200,hlabel,"Events"));
    }
    else if(i==ZMassmax){
      title.at(i)="$M_{e,\\mu} < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(ZMassmax));
      title.at(i)+=buffer;
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="m_{e#mu} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmax_",htitle,41,19,142,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmax_",htitle,41,19,142,hlabel,"Events"));
    }
    else if(i==ZMassmin){
      title.at(i)="$M_{e,\\mu} > $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(ZMassmin));
      title.at(i)+=buffer;
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="m_{e#mu} / GeV";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmin_",htitle,41,19,142,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmin_",htitle,41,19,142,hlabel,"Events"));
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
  RelIsoE=HConfig.GetTH1D(Name+"_RelIsoE","RelIsoE",20,0.,1.,"Relative Isolation of Electron");
  RelIsoMu=HConfig.GetTH1D(Name+"_RelIsoMu","RelIsoMu",20,0.,1.,"Relative Isolation of Muon");
  EPt=HConfig.GetTH1D(Name+"_PtE","PtE",40,0.,100.,"p_{T}^{e} / GeV");
  EEt=HConfig.GetTH1D(Name+"_EtE","EtE",40,0.,100.,"E_{T}^{e} / GeV");
  MuPt=HConfig.GetTH1D(Name+"_PtMu","PtMu",40,0.,100.,"p_{T}^{#mu} / GeV");
  mtMu=HConfig.GetTH1D(Name+"_mtMu","mtMu",40,0.,200.,"m_{t}^{#mu} / GeV");
  mtE=HConfig.GetTH1D(Name+"_mtE","mtE",40,0.,200.,"m_{t}^{e} / GeV");
  NJets=HConfig.GetTH1D(Name+"_NJets","NJets",20,0,20,"number of jets");
  NJetsLoose=HConfig.GetTH1D(Name+"_NJetsLoose","NJetsLoose",20,0,20,"number of jets (loose wp)");
  NJetsMedium=HConfig.GetTH1D(Name+"_NJetsMedium","NJetsMedium",20,0,20,"number of jets (medium wp)");
  NJetsTight=HConfig.GetTH1D(Name+"_NJetsTight","NJetsTight",20,0,20,"number of jets (tight wp)");
  PUJetId=HConfig.GetTH1D(Name+"_PUJetId","PUJetId",50,-1.,1.,"pu jet id discr.");
  etaMu=HConfig.GetTH1D(Name+"_etaMu","etaMu",20,-2.5,2.5,"#eta_{#mu}");
  etaE=HConfig.GetTH1D(Name+"_etaE","etaE",20,-2.5,2.5,"#eta_{e}");
  jetsum=HConfig.GetTH1D(Name+"_jetsum","jetsum",80,40.,440,"p_{T}^{leading jet}+p_{T}^{subleading jet} / GeV");
  chargesum=HConfig.GetTH1D(Name+"_chargesum","chargesum",21,-5.5,5.5,"charge sum");
  drmue=HConfig.GetTH1D(Name+"_drmue","drmue",50,0.,5.,"dR(e,#mu)");
  deltaphi=HConfig.GetTH1D(Name+"_deltaphi","deltaphi",40,0.,2,"#phi_{e#mu}");
  ptbal=HConfig.GetTH1D(Name+"_ptbal","ptbal",40,0.,200.,"p_{T}^{e#mu} / GeV");
  chargesumsigned=HConfig.GetTH1D(Name+"_chargesumsigned","chargesumsigned",21,-5.5,5.5,"charge sum");
  FirstJetPt=HConfig.GetTH1D(Name+"_FirstJetPt","FirstJetPt",40,0.,200.,"p_{T}^{hardest jet} / GeV");
  SecondJetPt=HConfig.GetTH1D(Name+"_SecondJetPt","SecondJetPt",40,0.,200.,"p_{T}^{2nd jet} / GeV");
  
  invmass_zmass=HConfig.GetTH1D(Name+"_invmass_zmass","invmass_zmass",41,19,142,"m_{e#mu} / GeV");
  invmass_ptbalance=HConfig.GetTH1D(Name+"_invmass_ptbalance","invmass_ptbalance",41,19,142,"m_{e#mu} / GeV");
  invmass_mtmu=HConfig.GetTH1D(Name+"_invmass_mtmu","invmass_mtmu",41,19,142,"m_{e#mu} / GeV");
  invmass_jetveto=HConfig.GetTH1D(Name+"_invmass_jetveto","invmass_jetveto",41,19,142,"m_{e#mu} / GeV");
  invmass_vetos=HConfig.GetTH1D(Name+"_invmass_vetos","invmass_vetos",41,19,142,"m_{e#mu} / GeV");
  invmass_only_object_id=HConfig.GetTH1D(Name+"_invmass_only_object_id","invmass_only_object_id",41,19,142,"m_{e#mu} / GeV");

  invmass_zmass_m=HConfig.GetTH1D(Name+"_invmass_zmass_m","invmass_zmass_m_m",20,61,121,"m_{e#mu} / GeV");
  invmass_ptbalance_m=HConfig.GetTH1D(Name+"_invmass_ptbalance_m","invmass_ptbalance_m",20,61,121,"m_{e#mu} / GeV");
  invmass_mtmu_m=HConfig.GetTH1D(Name+"_invmass_mtmu_m","invmass_mtmu_m",20,61,121,"m_{e#mu} / GeV");
  invmass_jetveto_m=HConfig.GetTH1D(Name+"_invmass_jetveto_m","invmass_jetveto_m",20,61,121,"m_{e#mu} / GeV");
  invmass_vetos_m=HConfig.GetTH1D(Name+"_invmass_vetos_m","invmass_vetos_m",20,61,121,"m_{e#mu} / GeV");
  invmass_only_object_id_m=HConfig.GetTH1D(Name+"_invmass_only_object_id_m","invmass_only_object_id_m",20,61,121,"m_{e#mu} / GeV");

  invmass_trilepton_only=HConfig.GetTH1D(Name+"_invmass_trilepton_only","invmass_trilepton_only",41,19,142,"m_{e#mu} / GeV");
  invmass_charge_only=HConfig.GetTH1D(Name+"_invmass_charge_only","invmass_charge_only",41,19,142,"m_{e#mu} / GeV");
  invmass_jetveto_only=HConfig.GetTH1D(Name+"_invmass_jetveto_only","invmass_jetveto_only",41,19,142,"m_{e#mu} / GeV");
  invmass_mtmu_only=HConfig.GetTH1D(Name+"_invmass_mtmu_only","invmass_mtmu_only",41,19,142,"m_{e#mu} / GeV");
  invmass_ptbal_only=HConfig.GetTH1D(Name+"_invmass_ptbal_only","invmass_ptbal_only",41,19,142,"m_{e#mu} / GeV");

  nm0_met=HConfig.GetTH1D(Name+"_nm0_met","nm0_met",30,0.,150.,"MET / GeV");
  nm0_mvamet=HConfig.GetTH1D(Name+"_nm0_mvamet","nm0_mvamet",30,0.,150.,"MET / GeV");
  nm0_onejet=HConfig.GetTH1D(Name+"_nm0_onejet","nm0_onejet",40,0.,200.,"p_{T}^{jet} / GeV");
  nm0_mtmu=HConfig.GetTH1D(Name+"_nm0_mtmu","nm0_mtmu",40,0.,160.,"m_{t}^{#mu} / GeV");
  nm0_ptbalance=HConfig.GetTH1D(Name+"_nm0_ptbalance","nm0_ptbalance",40,0.,200.,"p_{T}^{e#mu} / GeV");
  nmm_met=HConfig.GetTH1D(Name+"_nmm_met","nmm_met",30,0.,150.,"MET / GeV");
  nmm_mvamet=HConfig.GetTH1D(Name+"_nmm_mvamet","nmm_mvamet",30,0.,150.,"MET / GeV");
  nmm_onejet=HConfig.GetTH1D(Name+"_nmm_onejet","nmm_onejet",40,0.,200.,"p_{T}^{jet} / GeV");
  nmm_mtmu=HConfig.GetTH1D(Name+"_nmm_mtmu","nmm_mtmu",40,0.,160.,"m_{t}^{#mu} / GeV");
  nmm_ptbalance=HConfig.GetTH1D(Name+"_nmm_ptbalance","nmm_ptbalance",40,0.,200.,"p_{T}^{e#mu} / GeV");
  
  NPV=HConfig.GetTH1D(Name+"_NPV","NPV",60,0.,60.,"number of vertices");
  NPV3d=HConfig.GetTH1D(Name+"_NPV3d","NPV3d",60,0.,60.,"number of vertices");
  NPVfine=HConfig.GetTH1D(Name+"_NPVfine","NPVfine",60,0.,60.,"number of vertices");
  evtweight=HConfig.GetTH1D(Name+"_evtweight","evtweight",100,-0.1,2.1,"puweight");
  
  met=HConfig.GetTH1D(Name+"_met","met",30,0.,150.,"MET / GeV");
  met_uncorr=HConfig.GetTH1D(Name+"_met_uncorr","met_uncorr",30,0.,150.,"uncorrected MET / GeV");
  onejet=HConfig.GetTH1D(Name+"_onejet","onejet",40,0.,200.,"p_{T}^{jet} / GeV");
  onejet_eta=HConfig.GetTH1D(Name+"_onejet_eta","onejet_eta",40,-5.,5.,"#eta_{jet}");
  NbJets=HConfig.GetTH1D(Name+"_NbJets","NbJets",20,0,20,"number of b jets");
  NbJetsVtxL=HConfig.GetTH1D(Name+"_NbJetsVtxL","NbJetsVtxL",20,0,20,"number of b jets from vtx loose");
  NbJetsVtxM=HConfig.GetTH1D(Name+"_NbJetsVtxM","NbJetsVtxM",20,0,20,"number of b jets from vtx medium");
  NbJetsVtxT=HConfig.GetTH1D(Name+"_NbJetsVtxT","NbJetsVtxT",20,0,20,"number of b jets from vtx tight");
  
  zpt=HConfig.GetTH1D(Name+"_zpt","zpt",12,0.,30.,"p_{t}^{Z}");
  zeta=HConfig.GetTH1D(Name+"_zeta","zeta",40,-5.,5.,"#eta_{Z}");
  zmass=HConfig.GetTH1D(Name+"_zmass","zmass",20,60.,120.,"m_{Z}");

  eta_mu_e=HConfig.GetTH2D(Name+"_eta_mu_e","eta_mu_e",20,-2.5,2.5,20,-2.5,2.5,"#eta_{#mu}","#eta_{e}");
  pt_vs_eta_mu=HConfig.GetTH2D(Name+"_pt_vs_eta_mu","pt_vs_eta_mu",40,0.,100.,20,-2.5,2.5,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  pt_vs_eta_e=HConfig.GetTH2D(Name+"_pt_vs_eta_e","pt_vs_eta_e",40,0.,100.,20,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  nfakes=HConfig.GetTH1D(Name+"_nfakes","nfakes",2,0.5,2.5,"number of fake leptons");
  pt_vs_eta_mu_gen=HConfig.GetTH2D(Name+"_pt_vs_eta_mu_gen","pt_vs_eta_mu_gen",40,0.,100.,20,-2.5,2.5,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  pt_vs_eta_e_gen=HConfig.GetTH2D(Name+"_pt_vs_eta_e_gen","pt_vs_eta_e_gen",40,0.,100.,20,-2.5,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  higgs_mass=HConfig.GetTH1D(Name+"_higgs_mass","higgs_mass",10,121,141,"m_{e,#mu}");
  ht_pseudo=HConfig.GetTH1D(Name+"_ht_pseudo","ht_pseudo",200,0.,1000.,"#Sigma p_{T}^{jet} / GeV");
  zmass_zoom=HConfig.GetTH1D(Name+"_zmass_zoom","zmass_zoom",5,82,97,"m_{e#mu} / GeV");
  mtmu_vs_ptbal=HConfig.GetTH2D(Name+"_mtmu_vs_ptbal","mtmu_vs_ptbal",40,0.,200.,40,0.,200.,"m_{T}^{#mu} / GeV","p_{T}^{e#mu} / GeV");
  met_vs_ptbal=HConfig.GetTH2D(Name+"_met_vs_ptbal","met_vs_ptbal",30,0.,150.,40,0.,200.,"MET / GeV","p_{T}^{e#mu} / GeV");
  met_vs_mtmu=HConfig.GetTH2D(Name+"_met_vs_mtmu","met_vs_mtmu",30,0.,150.,40,0.,200.,"MET / GeV","m_{T}^{#mu} / GeV");
  invmass_high=HConfig.GetTH1D(Name+"_invmass_high","invmass_high",50,60.,160.,"m_{e#mu}");
  ptsum=HConfig.GetTH1D(Name+"_ptsum","ptsum",40,0.,200.,"p_{T}^{e}+p_{T}^{#mu} / GeV");
  ptsum_nm0=HConfig.GetTH1D(Name+"_ptsum_nm0","ptsum_nm0",40,0.,200.,"p_{T}^{e}+p_{T}^{#mu} / GeV");
  mvamet=HConfig.GetTH1D(Name+"_mvamet","mvamet",30,0.,150.,"MET / GeV");
  mva_mtmu=HConfig.GetTH1D(Name+"_mva_mtmu","mva_mtmu",40,0.,200.,"m_{T}^{#mu} / GeV");
  mtmu_phicorr=HConfig.GetTH1D(Name+"_mtmu_phicorr","mtmu_phicorr",40,0.,200.,"m_{T}^{#mu} / GeV");
  mte_mtmu=HConfig.GetTH1D(Name+"_mte_mtmu","mte_mtmu",40,0.,200.,"m_{T}^{e} / GeV");
  mtmu_mufake=HConfig.GetTH1D(Name+"_mtmu_mufake","mtmu_mufake",40,0.,200.,"m_{T}^{#mu} / GeV");
  mtmu_efake=HConfig.GetTH1D(Name+"_mtmu_efake","mtmu_efake",40,0.,200.,"m_{T}^{#mu} / GeV");
  mtmu_onefake=HConfig.GetTH1D(Name+"_mtmu_onefake","mtmu_onefake",40,0.,200.,"m_{T}^{#mu} / GeV");
  mtmu_twofakes=HConfig.GetTH1D(Name+"_mtmu_twofakes","mtmu_twofakes",40,0.,200.,"m_{T}^{#mu} / GeV");
  mtmu_nmu=HConfig.GetTH1D(Name+"_mtmu_nmu","mtmu_nmu",40,0.,200.,"m_{T}^{#mu} / GeV");

  if(doPDFuncertainty){
	  pdf_w0=HConfig.GetTH1D(Name+"_pdf_w0","pdf_w0",nPDFmembers,0,nPDFmembers,"pdf member");
	  pdf_w1=HConfig.GetTH1D(Name+"_pdf_w1","pdf_w1",nPDFmembers,0,nPDFmembers,"pdf member");
  }

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu::Store_ExtraDist(){
 Extradist1d.push_back(&RelIsoE);
 Extradist1d.push_back(&RelIsoMu);
 Extradist1d.push_back(&EPt);
 Extradist1d.push_back(&EEt);
 Extradist1d.push_back(&MuPt);
 Extradist1d.push_back(&mtMu);
 Extradist1d.push_back(&mtE);
 Extradist1d.push_back(&NJets);
 Extradist1d.push_back(&NJetsLoose);
 Extradist1d.push_back(&NJetsMedium);
 Extradist1d.push_back(&NJetsTight);
 Extradist1d.push_back(&PUJetId);
 Extradist1d.push_back(&etaE);
 Extradist1d.push_back(&etaMu);
 Extradist1d.push_back(&jetsum);
 Extradist1d.push_back(&chargesum);
 Extradist1d.push_back(&drmue);
 Extradist1d.push_back(&deltaphi);
 Extradist1d.push_back(&ptbal);
 Extradist1d.push_back(&chargesumsigned);
 Extradist1d.push_back(&FirstJetPt);
 Extradist1d.push_back(&SecondJetPt);
 
 Extradist1d.push_back(&invmass_zmass);
 Extradist1d.push_back(&invmass_ptbalance);
 Extradist1d.push_back(&invmass_mtmu);
 Extradist1d.push_back(&invmass_jetveto);
 Extradist1d.push_back(&invmass_vetos);
 Extradist1d.push_back(&invmass_only_object_id);

 Extradist1d.push_back(&invmass_zmass_m);
 Extradist1d.push_back(&invmass_ptbalance_m);
 Extradist1d.push_back(&invmass_mtmu_m);
 Extradist1d.push_back(&invmass_jetveto_m);
 Extradist1d.push_back(&invmass_vetos_m);
 Extradist1d.push_back(&invmass_only_object_id_m);

 Extradist1d.push_back(&invmass_trilepton_only);
 Extradist1d.push_back(&invmass_charge_only);
 Extradist1d.push_back(&invmass_jetveto_only);
 Extradist1d.push_back(&invmass_mtmu_only);
 Extradist1d.push_back(&invmass_ptbal_only);

 Extradist1d.push_back(&nm0_met);
 Extradist1d.push_back(&nm0_mvamet);
 Extradist1d.push_back(&nm0_onejet);
 Extradist1d.push_back(&nm0_mtmu);
 Extradist1d.push_back(&nm0_ptbalance);
 Extradist1d.push_back(&nmm_met);
 Extradist1d.push_back(&nmm_mvamet);
 Extradist1d.push_back(&nmm_onejet);
 Extradist1d.push_back(&nmm_mtmu);
 Extradist1d.push_back(&nmm_ptbalance);
 
 Extradist1d.push_back(&NPV);
 Extradist1d.push_back(&NPV3d);
 Extradist1d.push_back(&NPVfine);
 Extradist1d.push_back(&evtweight);
 
 Extradist1d.push_back(&met);
 Extradist1d.push_back(&met_uncorr);
 Extradist1d.push_back(&onejet);
 Extradist1d.push_back(&onejet_eta);
 Extradist1d.push_back(&NbJets);
 Extradist1d.push_back(&NbJetsVtxL);
 Extradist1d.push_back(&NbJetsVtxM);
 Extradist1d.push_back(&NbJetsVtxT);

 Extradist1d.push_back(&zpt);
 Extradist1d.push_back(&zeta);
 Extradist1d.push_back(&zmass);

 Extradist2d.push_back(&eta_mu_e);
 Extradist2d.push_back(&pt_vs_eta_mu);
 Extradist2d.push_back(&pt_vs_eta_e);
 Extradist1d.push_back(&nfakes);
 Extradist2d.push_back(&pt_vs_eta_mu_gen);
 Extradist2d.push_back(&pt_vs_eta_e_gen);
 Extradist1d.push_back(&higgs_mass);
 Extradist1d.push_back(&ht_pseudo);
 Extradist1d.push_back(&zmass_zoom);
 Extradist2d.push_back(&mtmu_vs_ptbal);
 Extradist2d.push_back(&met_vs_ptbal);
 Extradist2d.push_back(&met_vs_mtmu);
 Extradist1d.push_back(&invmass_high);
 Extradist1d.push_back(&ptsum);
 Extradist1d.push_back(&ptsum_nm0);
 Extradist1d.push_back(&mvamet);
 Extradist1d.push_back(&mva_mtmu);
 Extradist1d.push_back(&mtmu_phicorr);
 Extradist1d.push_back(&mte_mtmu);
 Extradist1d.push_back(&mtmu_mufake);
 Extradist1d.push_back(&mtmu_efake);
 Extradist1d.push_back(&mtmu_onefake);
 Extradist1d.push_back(&mtmu_twofakes);
 Extradist1d.push_back(&mtmu_nmu);

 if(doPDFuncertainty){
	 Extradist1d.push_back(&pdf_w0);
	 Extradist1d.push_back(&pdf_w1);
 }

}

void  ZtoEMu::doEvent(){
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
  // Tri Lepton Veto
  //
  if(verbose)std::cout << "trilepton veto" << std::endl;
  unsigned int trilep(0);
  if(eidx!=999 && muidx!=999){
	  for(unsigned i=0;i<Ntp->NMuons();i++){
		  if(i==muidx) continue;
		  if(vertex<0) continue;
		  if(Ntp->Muon_p4(i,mucorr).Pt()<10) continue;
		  if(fabs(Ntp->Muon_p4(i,mucorr).Eta())>2.4) continue;
		  if(!Ntp->isTightMuon(i,vertex,mucorr)) continue;
		  if(Ntp->Muon_RelIso(i,mucorr)>0.3) continue;
		  if(doHiggsObjects){
			  if(Ntp->dxy(Ntp->Muon_p4(i,mucorr),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))<0.045
					  && Ntp->dz(Ntp->Muon_p4(i,mucorr),Ntp->Muon_Poca(i),Ntp->Vtx(vertex))<0.2
					  ){
				  trilep++;
			  }
		  }else{
			  trilep++;
		  }
	  }
	  for(unsigned i=0;i<Ntp->NElectrons();i++){
		  if(i==eidx) continue;
		  if(vertex<0) continue;
		  if(Ntp->Electron_p4(i,ecorr).Et()<10) continue;
		  if(fabs(Ntp->Electron_supercluster_eta(i))>2.5) continue;
		  if(doHiggsObjects){
			  if(Ntp->isSelectedElectron(i,vertex,0.045,0.2)
					  && Ntp->Electron_RelIso04(i,ecorr)<0.3
					  ){
				  trilep++;
			  }
		  }else{
			  if(isWWElectron(i,vertex)
					  && Ntp->Electron_RelIso04(i,ecorr)<0.3
					  ){
				  trilep++;
			  }
		  }
	  }
  }
  value.at(triLeptonVeto)=trilep;
  pass.at(triLeptonVeto)=(value.at(triLeptonVeto)==cut.at(triLeptonVeto));
  
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
  // jet veto
  //
  if(verbose)std::cout << "jet veto" << std::endl;
  if(verbose)std::cout << "Cleaning jets" << std::endl;
  std::vector<int> jetsfromvtx;

  if(verbose)std::cout << "Finding jets from vtx" << std::endl;
  for(unsigned i=0;i<Ntp->NPFJets();i++){
	  // clean jets against signal objects
	  if(Ntp->PFJet_p4(i,jetcorr).Pt()<20) continue;
	  if(fabs(Ntp->PFJet_p4(i,jetcorr).Eta())>jet_eta) continue;
	  if(!Ntp->isJetID(i,jetcorr)) continue;
	  if(muidx!=999){
		  if(Ntp->PFJet_p4(i,jetcorr).DeltaR(Ntp->Muon_p4(muidx,mucorr))<0.3) continue;
	  }
	  if(eidx!=999){
		  if(Ntp->PFJet_p4(i,jetcorr).DeltaR(Ntp->Electron_p4(eidx,ecorr))<0.3) continue;
	  }
	  // find jets from vertex: use pileup jet id for jets with pt>20 GeV
	  if(Ntp->PFJet_PUJetID_tightWP(i)>0.5) jetsfromvtx.push_back(i);
  }

  if(verbose)std::cout<< "Find highest pt jet" << std::endl;
  int firstjet_idx=-1;
  int secondjet_idx=-1;
  double initialpt=0.;

  // loop over jets from selected vertex & find the two jets with the highest pt
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(Ntp->PFJet_p4(jetsfromvtx.at(i),jetcorr).Pt()>initialpt){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i),jetcorr).Pt();
		  firstjet_idx=jetsfromvtx.at(i);
	  }
  }
  initialpt=0.;
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(jetsfromvtx.size()>1 && firstjet_idx!=-1
			  && Ntp->PFJet_p4(jetsfromvtx.at(i),jetcorr).Pt()>initialpt
			  && Ntp->PFJet_p4(jetsfromvtx.at(i),jetcorr).Pt()<Ntp->PFJet_p4(firstjet_idx,jetcorr).Pt()
			  ){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i),jetcorr).Pt();
		  secondjet_idx=jetsfromvtx.at(i);
	  }
  }

  if(verbose)std::cout << "applying veto" << std::endl;

  value.at(oneJet)=0;
  if(jetsfromvtx.size()>0 && firstjet_idx!=-1){
	  value.at(oneJet)=Ntp->PFJet_p4(firstjet_idx,jetcorr).Pt();
  }
  pass.at(oneJet)=(value.at(oneJet)<cut.at(oneJet));
  
  ///////////////////////////////////////////////
  //
  // Mt Mu cut
  //
  if(verbose) std::cout << "Mt Mu cut" << std::endl;
  value.at(MtMu)=0.;
  if(muidx!=999){
	  value.at(MtMu)=sqrt(2*Ntp->Muon_p4(muidx,mucorr).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(muidx,mucorr).Px(),Ntp->Muon_p4(muidx,mucorr).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey())));
  }
  pass.at(MtMu)=(value.at(MtMu)<cut.at(MtMu));
  
  ///////////////////////////////////////////////
  //
  // Pt balance cut
  //
  if(verbose) std::cout << "pt balance cut" << std::endl;
  value.at(ptBalance)=0.;
  if(muidx!=999 && eidx!=999){
	  value.at(ptBalance) = (Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).Pt();
	  //if(Ntp->GetMCID()==DataMCType::DY_emu)value.at(ptBalance)*=ZPtReweight(value.at(ptBalance));
  }
  pass.at(ptBalance)=(value.at(ptBalance)<cut.at(ptBalance));

  ///////////////////////////////////////////////
  //
  // Invariant mass cut
  //
  if(verbose) std::cout << "Invariant mass cut" << std::endl;
  value.at(ZMassmax)=zmax+1.;
  value.at(ZMassmin)=zmin-1.;
  if(eidx!=999 && muidx!=999){
	  value.at(ZMassmax)=(Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).M();
	  value.at(ZMassmin)=(Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).M();
  }
  pass.at(ZMassmax)=(value.at(ZMassmax)<cut.at(ZMassmax));
  pass.at(ZMassmin)=(value.at(ZMassmin)>=cut.at(ZMassmin));

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
  // Get PDF weights
  //
  if(doPDFuncertainty){
	  if(verbose) std::cout << "Calculating PDF weights" << std::endl;
	  for(unsigned int member=0;member<nPDFmembers;++member){
		  if(Ntp->isData() || Ntp->GetMCID()==DataMCType::QCD) break;
		  double pdfWeight = w*pdf->weight(Ntp->GenEventInfoProduct_id1(),Ntp->GenEventInfoProduct_id2(),Ntp->GenEventInfoProduct_x1(),Ntp->GenEventInfoProduct_x2(),Ntp->GenEventInfoProduct_scalePDF(),member);
		  pdf_w0.at(t).AddBinContent(member+1,pdfWeight);
		  if(status) pdf_w1.at(t).AddBinContent(member+1,pdfWeight);
	  }
  }

  ///////////////////////////////////////////////
  //
  // Add plots
  //
  if(verbose) std::cout << "add plots" << std::endl;
  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(NMu)
		  && pass.at(NE)
		  && pass.at(ptthreshold)
		  && pass.at(charge)
		  && pass.at(triLeptonVeto)
		  && pass.at(oneJet)
		  && pass.at(MtMu)
		  && pass.at(ptBalance)
                  ){
		  	if((Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).M()>60.) invmass_high.at(t).Fill((Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).M(),w);
		  }
  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(NMu)
		  ){
	  mtmu_nmu.at(t).Fill(Ntp->transverseMass(Ntp->Muon_p4(muidx,mucorr).Pt(),Ntp->Muon_p4(muidx,mucorr).Phi(),Ntp->MET_CorrT0pcT1_et(),Ntp->MET_CorrT0pcT1_phi()),w);
  }
  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(NMu)
		  && pass.at(NE)
		  && pass.at(ptthreshold)
		  && pass.at(mll)
		  //&& pass.at(charge)
		  ){
	  // often needed variables
	  double m = (Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).M();
	  double dp = Ntp->Muon_p4(muidx,mucorr).DeltaPhi(Ntp->Electron_p4(eidx,ecorr))/TMath::Pi();
	  if(dp<0)dp+=2;

	  // electron related histograms
	  EPt.at(t).Fill(Ntp->Electron_p4(eidx,ecorr).Pt(),w);
	  EEt.at(t).Fill(Ntp->Electron_p4(eidx,ecorr).Et(),w);
	  etaE.at(t).Fill(Ntp->Electron_supercluster_eta(eidx),w);
	  mtE.at(t).Fill(sqrt(2*Ntp->Electron_p4(eidx,ecorr).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Electron_p4(eidx,ecorr).Px(),Ntp->Electron_p4(eidx,ecorr).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
	  RelIsoE.at(t).Fill(Ntp->Electron_RelIso04(eidx),w);
	  pt_vs_eta_e.at(t).Fill(Ntp->Electron_p4(eidx,ecorr).Pt(),Ntp->Electron_supercluster_eta(eidx),w);

	  // muon related histograms
	  MuPt.at(t).Fill(Ntp->Muon_p4(muidx,mucorr).Pt(),w);
	  etaMu.at(t).Fill(Ntp->Muon_p4(muidx,mucorr).Eta(),w);
	  mtMu.at(t).Fill(sqrt(2*Ntp->Muon_p4(muidx,mucorr).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(muidx,mucorr).Px(),Ntp->Muon_p4(muidx,mucorr).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
	  RelIsoMu.at(t).Fill(Ntp->Muon_RelIso(muidx,mucorr),w);
	  pt_vs_eta_mu.at(t).Fill(Ntp->Muon_p4(muidx,mucorr).Pt(),Ntp->Muon_p4(muidx,mucorr).Eta(),w);

	  // histograms related to combination
	  eta_mu_e.at(t).Fill(Ntp->Muon_p4(muidx,mucorr).Eta(),Ntp->Electron_supercluster_eta(eidx),w);
	  drmue.at(t).Fill(Ntp->Muon_p4(muidx,mucorr).DeltaR(Ntp->Electron_p4(eidx,ecorr)),w);
	  met.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),w);
	  met_uncorr.at(t).Fill(Ntp->MET_Uncorr_et(),w);
	  mvamet.at(t).Fill(Ntp->MET_CorrMVA_et(),w);
	  mva_mtmu.at(t).Fill(sqrt(2*Ntp->Muon_p4(muidx,mucorr).Pt()*Ntp->MET_CorrMVA_et()*(1-cosphi2d(Ntp->Muon_p4(muidx,mucorr).Px(),Ntp->Muon_p4(muidx,mucorr).Py(),Ntp->MET_CorrMVA_ex(),Ntp->MET_CorrMVA_ey()))),w);
	  deltaphi.at(t).Fill(dp,w);
	  chargesum.at(t).Fill(fabs(Ntp->Muon_Charge(muidx)+Ntp->Electron_Charge(eidx)),w);
	  chargesumsigned.at(t).Fill(Ntp->Muon_Charge(muidx)+Ntp->Electron_Charge(eidx),w);
	  ptbal.at(t).Fill((Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).Pt(),w);
	  mtmu_vs_ptbal.at(t).Fill(sqrt(2*Ntp->Muon_p4(muidx,mucorr).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(muidx,mucorr).Px(),Ntp->Muon_p4(muidx,mucorr).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),(Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).Pt(),w);
	  ptsum.at(t).Fill(Ntp->Muon_p4(muidx,mucorr).Pt()+Ntp->Electron_p4(eidx,ecorr).Et(),w);

	  //mtmu cross checks
	  mtmu_phicorr.at(t).Fill(sqrt(2*Ntp->Muon_p4(muidx,mucorr).Pt()*Ntp->MET_CorrT0pcT1Txy_et()*(1-cosphi2d(Ntp->Muon_p4(muidx,mucorr).Px(),Ntp->Muon_p4(muidx,mucorr).Py(),Ntp->MET_CorrT0pcT1Txy_ex(),Ntp->MET_CorrT0pcT1Txy_ey()))),w);
	  double mt = sqrt(2*Ntp->Muon_p4(muidx,mucorr).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(muidx,mucorr).Px(),Ntp->Muon_p4(muidx,mucorr).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey())));
	  if(fakemu && !fakee && Ntp->isData()) mtmu_mufake.at(t).Fill(mt,w);
	  if(!fakemu && fakee && Ntp->isData()) mtmu_efake.at(t).Fill(mt,w);
	  if((fakemu || fakee) && !(fakemu && fakee) && Ntp->isData()) mtmu_onefake.at(t).Fill(mt,w);
	  if(fakemu && fakee && Ntp->isData()) mtmu_twofakes.at(t).Fill(mt,w);

	  if(fakemu && !fakee) nfakes.at(t).AddBinContent(1,w);
	  else if(!fakemu && fakee) nfakes.at(t).AddBinContent(1,w);
	  else if(fakemu && fakee) nfakes.at(t).AddBinContent(2,w);


	  if(jetsfromvtx.size()==1 && firstjet_idx!=-1){
		  onejet.at(t).Fill(Ntp->PFJet_p4(firstjet_idx,jetcorr).Pt(),w);
		  onejet_eta.at(t).Fill(Ntp->PFJet_p4(firstjet_idx,jetcorr).Eta(),w);
	  }
	  if(jetsfromvtx.size()>1 && firstjet_idx!=-1 && secondjet_idx!=-1){
		  jetsum.at(t).Fill(Ntp->PFJet_p4(firstjet_idx,jetcorr).Pt()+Ntp->PFJet_p4(secondjet_idx,jetcorr).Pt(),w);
	  }
	  if(pass.at(MtMu))mte_mtmu.at(t).Fill(sqrt(2*Ntp->Electron_p4(eidx,ecorr).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Electron_p4(eidx,ecorr).Px(),Ntp->Electron_p4(eidx,ecorr).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);

	  double sumht(0);
	  for(unsigned i=0;i<jetsfromvtx.size();i++){
		  sumht += Ntp->PFJet_p4(jetsfromvtx.at(i),jetcorr).Et();
	  }
	  ht_pseudo.at(t).Fill(sumht,w);

	  std::vector<unsigned int> jets;
	  std::vector<unsigned int> loosejets;
	  std::vector<unsigned int> mediumjets;
	  std::vector<unsigned int> tightjets;
	  int nbjets = 0;
	  int nbjetsvtxl = 0;
	  int nbjetsvtxm = 0;
	  int nbjetsvtxt = 0;

	  for(unsigned int i=0;i<Ntp->NPFJets();i++){
		  if(Ntp->PFJet_p4(i,jetcorr).Pt()<20) continue;
		  if(fabs(Ntp->PFJet_p4(i,jetcorr).Eta())>jet_eta) continue;
		  if(!Ntp->isJetID(i,jetcorr)) continue;
		  if(muidx!=999){
			  if(Ntp->PFJet_p4(i,jetcorr).DeltaR(Ntp->Muon_p4(muidx,mucorr))<0.3) continue;
		  }
		  if(eidx!=999){
			  if(Ntp->PFJet_p4(i,jetcorr).DeltaR(Ntp->Electron_p4(eidx,ecorr))<0.3) continue;
		  }
		  PUJetId.at(t).Fill(Ntp->PFJet_PUJetID_discr(i),w);
		  jets.push_back(i);
		  if(Ntp->PFJet_PUJetID_looseWP(i)>0.5) loosejets.push_back(i);
		  if(Ntp->PFJet_PUJetID_mediumWP(i)>0.5) mediumjets.push_back(i);
		  if(Ntp->PFJet_PUJetID_tightWP(i)>0.5) tightjets.push_back(i);
		  if(Ntp->PFJet_bDiscriminator(i)>csvl){
			  nbjets++;
			  if(Ntp->PFJet_PUJetID_looseWP(i)>0.5){
				  nbjetsvtxl++;
				  if(Ntp->PFJet_bDiscriminator(i)>csvm){
					  nbjetsvtxm++;
					  if(Ntp->PFJet_bDiscriminator(i)>csvt){
						  nbjetsvtxt++;
					  }
				  }
			  }
		  }
	  }

	  NJets.at(t).Fill(jets.size(),w);
	  NJetsLoose.at(t).Fill(loosejets.size(),w);
	  NJetsMedium.at(t).Fill(mediumjets.size(),w);
	  NJetsTight.at(t).Fill(tightjets.size(),w);
	  NbJets.at(t).Fill(nbjets,w);
	  NbJetsVtxL.at(t).Fill(nbjetsvtxl,w);
	  NbJetsVtxM.at(t).Fill(nbjetsvtxm,w);
	  NbJetsVtxT.at(t).Fill(nbjetsvtxt,w);

	  if(jetsfromvtx.size()>0 && firstjet_idx!=-1){
		  FirstJetPt.at(t).Fill(Ntp->PFJet_p4(firstjet_idx,jetcorr).Pt(),w);
	  }
	  if(jetsfromvtx.size()>1 && firstjet_idx!=-1 && secondjet_idx!=-1){
		  SecondJetPt.at(t).Fill(Ntp->PFJet_p4(secondjet_idx,jetcorr).Pt(),w);
	  }

	  NPV.at(t).Fill(Ntp->NVtx(),w);
	  if(!Ntp->isData()){
		  evtweight.at(t).Fill(Ntp->PUWeight());
		  if(Ntp->PUWeight()>0.)NPV3d.at(t).Fill(Ntp->NVtx(),w*Ntp->PUWeight3D()/Ntp->PUWeight());
		  if(Ntp->PUWeight()>0.)NPVfine.at(t).Fill(Ntp->NVtx(),w*Ntp->PUWeightFineBins()/Ntp->PUWeight());
	  }else{
		  NPV3d.at(t).Fill(Ntp->NVtx(),w);
		  NPVfine.at(t).Fill(Ntp->NVtx(),w);
	  }

	  if(Ntp->GetMCID()==30 || Ntp->GetMCID()==31 || Ntp->GetMCID()==32 || Ntp->GetMCID()==33 || Ntp->GetMCID()==40){
		  for(unsigned i=0;i<Ntp->NMCParticles();i++){
			  if(Ntp->MCParticle_pdgid(i)==23){
				  if(!useMadgraphZ && (Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau)){
					  zpt.at(t).Fill(Ntp->MCParticle_p4(i).Pt())*PowhegReweight(Ntp->MCParticle_p4(i).Pt());
				  }else{
					  zpt.at(t).Fill(Ntp->MCParticle_p4(i).Pt());
				  }
				  zeta.at(t).Fill(Ntp->MCParticle_p4(i).Eta());
				  zmass.at(t).Fill(Ntp->MCParticle_p4(i).M());
			  }
		  }
	  }

	  if(verbose)std::cout << "invariant mass with different cuts applied" << std::endl;

	  invmass_only_object_id.at(t).Fill(m,w);
	  invmass_only_object_id_m.at(t).Fill(m,w);
	  if(pass.at(triLeptonVeto)
			  && pass.at(charge)
			  ){
		  invmass_vetos.at(t).Fill(m,w);
		  invmass_vetos_m.at(t).Fill(m,w);
		  if(pass.at(oneJet)){
			  invmass_jetveto.at(t).Fill(m,w);
			  invmass_jetveto_m.at(t).Fill(m,w);
			  if(pass.at(MtMu)){
				  invmass_mtmu.at(t).Fill(m,w);
				  invmass_mtmu_m.at(t).Fill(m,w);
				  if(pass.at(ptBalance)){
					  invmass_ptbalance.at(t).Fill(m,w);
					  invmass_ptbalance_m.at(t).Fill(m,w);
					  if(jetsfromvtx.size()>0 && firstjet_idx!=-1) nm0_onejet.at(t).Fill(value.at(oneJet),w);
					  nm0_met.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),w);
					  nm0_mvamet.at(t).Fill(Ntp->MET_CorrMVA_et(),w);
					  nm0_mtmu.at(t).Fill(value.at(MtMu),w);
					  nm0_ptbalance.at(t).Fill(value.at(ptBalance),w);
					  zmass_zoom.at(t).Fill(m,w);
					  ptsum_nm0.at(t).Fill(Ntp->Muon_p4(muidx,mucorr).Pt()+Ntp->Electron_p4(eidx,ecorr).Pt(),w);
					  if(pass.at(ZMassmax)
							  && pass.at(ZMassmin)){
						  invmass_zmass.at(t).Fill(m,w);
						  invmass_zmass_m.at(t).Fill(m,w);
						  if(jetsfromvtx.size()>0 && firstjet_idx!=-1) nmm_onejet.at(t).Fill(Ntp->PFJet_p4(firstjet_idx,jetcorr).Pt(),w);
						  nmm_met.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),w);
						  nmm_mvamet.at(t).Fill(Ntp->MET_CorrMVA_et(),w);
						  nmm_mtmu.at(t).Fill(sqrt(2*Ntp->Muon_p4(muidx,mucorr).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(muidx,mucorr).Px(),Ntp->Muon_p4(muidx,mucorr).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
						  nmm_ptbalance.at(t).Fill(value.at(ptBalance),w);
					  }
					  if(m>=123 && m<129) higgs_mass.at(t).Fill(m,w);
				  }
			  	  met_vs_ptbal.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),(Ntp->Muon_p4(muidx,mucorr)+Ntp->Electron_p4(eidx,ecorr)).Pt(),w);
			  }
			  if(pass.at(ptBalance)) met_vs_mtmu.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),sqrt(2*Ntp->Muon_p4(muidx,mucorr).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(muidx,mucorr).Px(),Ntp->Muon_p4(muidx,mucorr).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
		  }
	  }
	  if(pass.at(triLeptonVeto)) invmass_trilepton_only.at(t).Fill(m,w);
	  if(pass.at(charge)) invmass_charge_only.at(t).Fill(m,w);
	  if(pass.at(oneJet)) invmass_jetveto_only.at(t).Fill(m,w);
	  if(pass.at(MtMu)) invmass_mtmu_only.at(t).Fill(m,w);
	  if(pass.at(ptBalance)) invmass_ptbal_only.at(t).Fill(m,w);
  }

  TLorentzVector egen(0.,0.,0.,0.);
  TLorentzVector mugen(0.,0.,0.,0.);
  for(unsigned i=0;i<Ntp->NMCParticles();i++){
	  if(fabs(Ntp->MCParticle_pdgid(i))==13 && (Ntp->MCParticle_midx(i)==23 || fabs(Ntp->MCParticle_midx(i))==24 || fabs(Ntp->MCParticle_midx(i))==15)){
		  pt_vs_eta_mu_gen.at(t).Fill(Ntp->MCParticle_p4(i).Pt(),Ntp->MCParticle_p4(i).Eta(),w);
	  }
	  if(fabs(Ntp->MCParticle_pdgid(i))==11 && (Ntp->MCParticle_midx(i)==23 || fabs(Ntp->MCParticle_midx(i))==24 || fabs(Ntp->MCParticle_midx(i))==15)){
		  pt_vs_eta_e_gen.at(t).Fill(Ntp->MCParticle_p4(i).Pt(),Ntp->MCParticle_p4(i).Eta(),w);
	  }
	  if(fabs(Ntp->MCParticle_pdgid(i))==11 && fabs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(i)))==15){
		  egen = Ntp->MCParticle_p4(i);
	  }
	  if(fabs(Ntp->MCParticle_pdgid(i))==13 && fabs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(i)))==15){
		  mugen = Ntp->MCParticle_p4(i);
	  }
  }

  if(verbose)std::cout << "ZtoEMu::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

double ZtoEMu::ZPtReweight(double zpt){
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

double ZtoEMu::PowhegReweight(double zpt){
	double weight = 1.;
	weight = ZptCorrection->GetBinContent(ZptCorrection->FindFixBin(zpt));
	return weight;
}

// relative uncertainties obtained from table 16 in AN-2012-225
double ZtoEMu::ZPtRelUnc(double zpt){
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

double ZtoEMu::ZPtMadgraphRelUnc(double zpt){
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

double ZtoEMu::calculatePzeta(int muiterator, int eiterator){
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

double ZtoEMu::calculatePzetaDQM(int muiterator, int eiterator){
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

double ZtoEMu::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double ZtoEMu::cosphi3d(TVector3 vec1, TVector3 vec2){
	return (vec1.Dot(vec2))/vec1.Mag()/vec2.Mag();
}

int ZtoEMu::findBin(TGraphAsymmErrors* graph, double xval){
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

bool ZtoEMu::isFakeMuon(unsigned int idx, TString corr){
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

bool ZtoEMu::isFakeMuon(unsigned int idx, unsigned int vtx, TString corr){
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

bool ZtoEMu::isWWElectron(unsigned int idx, unsigned int vtx, TString corr){
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

bool ZtoEMu::isFakeElectron(unsigned int idx, TString corr){
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

bool ZtoEMu::isFakeElectron(unsigned int idx, unsigned int vtx, TString corr){
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

double ZtoEMu::Fakerate(double pt, double eta, TH2D *fakeRateHist){
	double fakerate = 0.;
	if(pt>=35.) pt = 34.99;
	fakerate = fakeRateHist->GetBinContent(fakeRateHist->FindFixBin(pt,eta));
	return fakerate/(1-fakerate);
}

double ZtoEMu::FakerateError(double pt, double eta, TH2D *fakeRateHist){
	double error = 0.;
	if(pt>=35.) pt = 34.99;
	error = fakeRateHist->GetBinError(fakeRateHist->FindFixBin(pt,eta));
	return error/(1-error);
}

//////////////////////////////
//
// Finish function
//

void ZtoEMu::Finish(){
	Selection::Finish();
	Selection::ScaleAllHistOfType(DataMCType::DY_emu_embedded,10152445./9760728.);
	double sumdata(0), sumbkg(0), sumsignal(0);
	double stddata(87),stdbkg(80.629),stdsignal(74.680);
	double stdqcd(14.371),stdww(16.689),stdwz2l2q(0.000),stdwz3l1nu(0.299),stdzz4l(0.026),stdzz2l2q(0.000),stdzz2l2nu(0.009),stdtt(1.539),stdtw(0.592),stdtbarw(0.000),stddyll(7.141),stddytt(39.963);
	double std(0),stdunc(0);
	double bkgunc(0),individualunc(0);
	for(unsigned i=0;i<HConfig.GetNHisto();i++){
		if(HConfig.GetID(i)==DataMCType::Data) sumdata = Npassed.at(i).GetBinContent(NCuts+1);
		else if(HConfig.GetID(i)==DataMCType::DY_emu) sumsignal = Npassed.at(i).GetBinContent(NCuts+1);
		else if(HConfig.GetID(i)>0) sumbkg += Npassed.at(i).GetBinContent(NCuts+1);

		if(HConfig.GetID(i)==DataMCType::DY_tautau || HConfig.GetID(i)==DataMCType::DY_ll) bkgunc += pow(Npassed.at(i).GetBinContent(NCuts+1) * normunc_dy,2);
		if(HConfig.GetID(i)==DataMCType::ttbar) bkgunc += pow(Npassed.at(i).GetBinContent(NCuts+1) * normunc_tt,2);
		if(HConfig.GetID(i)==DataMCType::tw || HConfig.GetID(i)==DataMCType::tbarw) bkgunc += pow(Npassed.at(i).GetBinContent(NCuts+1) * normunc_tw,2);
		if(HConfig.GetID(i)==DataMCType::WW_2l2nu
				|| HConfig.GetID(i)==DataMCType::WZ_3l1nu
				|| HConfig.GetID(i)==DataMCType::WZ_2l2q
				|| HConfig.GetID(i)==DataMCType::ZZ_4l
				|| HConfig.GetID(i)==DataMCType::ZZ_2l2q
				|| HConfig.GetID(i)==DataMCType::ZZ_2l2nu
				){
			bkgunc += pow(Npassed.at(i).GetBinContent(NCuts+1) * normunc_diboson,2);
		}
		if(HConfig.GetID(i)==DataMCType::QCD) bkgunc += pow(Npassed.at(i).GetBinContent(NCuts+1) * normunc_qcd,2);

		if(HConfig.GetID(i)==DataMCType::Data) std = stddata;
		if(HConfig.GetID(i)==DataMCType::QCD) std = stdqcd;
		if(HConfig.GetID(i)==DataMCType::WW_2l2nu) std = stdww;
		if(HConfig.GetID(i)==DataMCType::WZ_2l2q) std = stdwz2l2q;
		if(HConfig.GetID(i)==DataMCType::WZ_3l1nu) std = stdwz3l1nu;
		if(HConfig.GetID(i)==DataMCType::ZZ_4l) std = stdzz4l;
		if(HConfig.GetID(i)==DataMCType::ZZ_2l2q) std = stdzz2l2q;
		if(HConfig.GetID(i)==DataMCType::ZZ_2l2nu) std = stdzz2l2nu;
		if(HConfig.GetID(i)==DataMCType::ttbar) std = stdtt;
		if(HConfig.GetID(i)==DataMCType::tw) std = stdtw;
		if(HConfig.GetID(i)==DataMCType::tbarw) std = stdtbarw;
		if(HConfig.GetID(i)==DataMCType::DY_ll) std = stddyll;
		if(HConfig.GetID(i)==DataMCType::DY_tautau) std = stddytt;
		if(HConfig.GetID(i)==DataMCType::DY_emu) std = stdsignal;
		if(HConfig.GetID(i)!=DataMCType::Data && HConfig.GetID(i)!=DataMCType::DY_emu) individualunc += pow(abs(std-Npassed.at(i).GetBinContent(NCuts+1)),2);
		std::cout << "Number of events for sample " << HConfig.GetName(i) << " ";
		printf("Standard: %.3f, Now: %.3f Resulting uncertainty: %.3f\n",std,Npassed.at(i).GetBinContent(NCuts+1),abs(std-Npassed.at(i).GetBinContent(NCuts+1)));

		if(zpt.at(i).Integral()>0)zpt.at(i).Scale(1./zpt.at(i).Integral());
		if(zeta.at(i).Integral()>0)zeta.at(i).Scale(1./zeta.at(i).Integral());
		if(zmass.at(i).Integral()>0)zmass.at(i).Scale(1./zmass.at(i).Integral());

	}

	std::cout << "################## Result ##################" << std::endl;
	printf("Total events in data: %.0f, background: %.3f and signal %.3f\n",sumdata,sumbkg,sumsignal);
	printf("Total background normalization uncertainty: %.3f%%\n",sqrt(bkgunc)/sumbkg*100);
	printf("Total background uncertainty not from cross section: %.3f%%\n",sqrt(individualunc)/sumbkg*100);
	printf("Cases with one fake lepton: %f. Cases with two fake leptons: %f.\n",nfakes.at(HConfig.GetType(DataMCType::QCD)).GetBinContent(1),nfakes.at(HConfig.GetType(DataMCType::QCD)).GetBinContent(2));

	printf("Difference in background events: %.3f%% and in signal events: %.3f%%\n",fabs(sumbkg/stdbkg-1)*100,fabs(sumsignal/stdsignal-1)*100);

	if(doPDFuncertainty){
		TString pdfname1_lower = pdfname1;
		TString pdfname2_lower = pdfname2;
		pdfname1_lower.ToLower();
		pdfname2_lower.ToLower();
		for(unsigned i=0;i<HConfig.GetNHisto();i++){
			if(HConfig.GetID(i)==DataMCType::Data || HConfig.GetID(i)==DataMCType::QCD) continue;
			double acc = pdf_w1.at(i).GetBinContent(1)/pdf_w0.at(i).GetBinContent(1);
			double yield = pdf_w1.at(i).GetBinContent(1);
			double acc2 = 0;
			double yield2 = 0;
			if(pdfname2_lower.Contains("nnpdf")){
				for(unsigned j=1;j<nPDFmembers;++j){
					acc2 += pow(pdf_w1.at(i).GetBinContent(j+1)/pdf_w0.at(i).GetBinContent(j+1)-acc,2);
					yield2 += pow(pdf_w1.at(i).GetBinContent(j+1)-yield,2);
				}
			}else{
				for(unsigned j=1;j<nPDFmembers;j+=2){
					acc2 += pow(pdf_w1.at(i).GetBinContent(j+1)/pdf_w0.at(i).GetBinContent(j+1)-pdf_w1.at(i).GetBinContent(j+2)/pdf_w0.at(i).GetBinContent(j+2),2);
					yield2 += pow(pdf_w1.at(i).GetBinContent(j+1)-pdf_w1.at(i).GetBinContent(j+2),2);
				}
			}
			double eacc = sqrt(acc2/(nPDFmembers-1));
			double eyield = sqrt(yield2/(nPDFmembers-1));
			if(pdfname1_lower.Contains("cteq") || pdfname2_lower.Contains("cteq")){ // CTEQ must be corrected from 90%CL to 68%CL
				eacc /= 1.645;
				eyield /= 1.645;
			}
			printf("====================\n");
			std::cout << "PDF uncertainties for sample " << HConfig.GetName(i) << std::endl;
			printf("Yield = %.3e +- %.3e (PDFs), i.e. %.2f%% relative uncertainty\n",yield,eyield,eyield/yield*100);
			printf("Acceptance = %.3e +- %.3e (PDFs), i.e. %.2f%% relative uncertainty\n",acc,eacc,eacc/acc*100);
		}
	}

}
