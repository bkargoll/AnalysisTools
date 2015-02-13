#include "Tvariable_Base.h"
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

Tvariable_Base::Tvariable_Base(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_ptlow(10)
  ,e_ptlow(10)
  ,mu_eta(2.4)
  ,e_eta(2.5)
  ,lead_pt(30)
  ,trail_pt(20)
  ,lead_eta(2.1)
  ,trail_eta(2.4)
  ,mmin(20)
  ,jet_pt(18)
  ,jet_eta(5.2) // 2.5 in dileptonic top selection
  ,singlejet(40)
  ,met(30)
  ,ptbalance(10)
  ,zmin(88)
  ,zmax(94)
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

	doWWObjects = true;
	useMadgraphZ = true;
	if(useMadgraphZ) mmin = 50;

	// Set category flag to standard value
	categoryFlag = "noCategory";

	// corrections to be applied to particle candidates
	mucorr = "roch";
	ecorr = "";
	jetcorr = "runJER";

	// decide if and which systematics should be done

	upwardUncertainty = false; // important for some uncertainties

	doPDFuncertainty = false;
	doTriggerUncertainty = false;
	doPileupUncertainty = false;
	doElectronIdUncertainty = false;
	doElectronScaleUncertainty = false;
	doElectronResUncertainty = false;
	doMuonIdUncertainty = false;
	doMuonScaleUncertainty = false;
	doMuonResUncertainty = false;
	doJECUncertainty = false;
	doJERUncertainty = false;
	doFakeRateUncertainty = false;
	doMetUncertainty = false;

	systValid = true;
	if((doPDFuncertainty+doTriggerUncertainty+doPileupUncertainty+doElectronIdUncertainty+doElectronScaleUncertainty+doElectronResUncertainty+doMuonIdUncertainty+doMuonScaleUncertainty+doMuonResUncertainty+doJECUncertainty+doJERUncertainty+doFakeRateUncertainty)>1){
		systValid = false;
	}

	// set correction strings for systematics
	if(doElectronResUncertainty) ecorr += "res";
	else if(doElectronScaleUncertainty){
		ecorr += "scale";
		if(!upwardUncertainty) ecorr += "down";
	}
	else if(doMuonResUncertainty) mucorr += "res";
	else if(doMuonScaleUncertainty){
		mucorr += "scale";
		if(!upwardUncertainty) mucorr += "down";
	}
	else if(doJECUncertainty){
		jetcorr += "JEC";
		if(upwardUncertainty) jetcorr += "plus";
	}
	else if(doJERUncertainty){
		if(upwardUncertainty) jetcorr += "up";
		else jetcorr += "down";
	}

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

Tvariable_Base::~Tvariable_Base(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    std::cout << "Tvariable_Base::~Tvariable_Base Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts+1) << std::endl;
  }
  std::cout << "Tvariable_Base::~Tvariable_Base()" << std::endl;
}

void  Tvariable_Base::Configure(){

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
    if(i==MET)                cut.at(MET)=met;
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
	  hlabel="m_{e#mu} (GeV)";
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
  	  hlabel="jet P_{T} (GeV)";
  	  Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_oneJet_",htitle,40,0,200,hlabel,"Events"));
  	  Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_oneJet_",htitle,40,0,200,hlabel,"Events"));
	}
    else if(i==MET){
	  title.at(i)="$m_{T}^{\\mu,Miss} < $";
	  title.at(i)+=cut.at(MET);
	  title.at(i)+="(GeV)";
	  htitle=title.at(i);
	  htitle.ReplaceAll("$","");
	  htitle.ReplaceAll("\\","#");
	  hlabel="m_{T}^{#mu,Miss} (GeV)";
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
      hlabel="p_t balance (GeV)";
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
      hlabel="m_{e#mu} (GeV)";
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
      hlabel="m_{e#mu} (GeV)";
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
  LeadPt=HConfig.GetTH1D(Name+"_LeadPt","LeadPt",40,0.,100.,"p_{T}^{lead} (GeV)");
  TrailPt=HConfig.GetTH1D(Name+"_TrailPt","TrailPt",40,0.,100.,"p_{T}^{trail} (GeV)");
  NJets=HConfig.GetTH1D(Name+"_NJets","NJets",20,0,20,"number of jets");
  NJetsLoose=HConfig.GetTH1D(Name+"_NJetsLoose","NJetsLoose",20,0,20,"number of jets (loose wp)");
  NJetsMedium=HConfig.GetTH1D(Name+"_NJetsMedium","NJetsMedium",20,0,20,"number of jets (medium wp)");
  NJetsTight=HConfig.GetTH1D(Name+"_NJetsTight","NJetsTight",20,0,20,"number of jets (tight wp)");
  PUJetId=HConfig.GetTH1D(Name+"_PUJetId","PUJetId",50,-1.,1.,"pu jet id discr.");
  etaLead=HConfig.GetTH1D(Name+"_etaLead","etaLead",20,-2.5,2.5,"#eta_{lead}");
  etaTrail=HConfig.GetTH1D(Name+"_etaTrail","etaTrail",20,-2.5,2.5,"#eta_{trail}");
  drleadtrail=HConfig.GetTH1D(Name+"_drleadtrail","drleadtrail",50,0.,5.,"dR(lead,trail)");
  deltaphi=HConfig.GetTH1D(Name+"_deltaphi","deltaphi",40,0.,2,"#phi_{e#mu}");
  ptbal=HConfig.GetTH1D(Name+"_ptbal","ptbal",40,0.,200.,"p_{T}^{e#mu} (GeV)");

  invmass_zmass=HConfig.GetTH1D(Name+"_invmass_zmass","invmass_zmass",41,19,142,"m_{e#mu} (GeV)");
  invmass_ptbalance=HConfig.GetTH1D(Name+"_invmass_ptbalance","invmass_ptbalance",41,19,142,"m_{e#mu} (GeV)");
  invmass_met=HConfig.GetTH1D(Name+"_invmass_met","invmass_met",41,19,142,"m_{e#mu} (GeV)");
  invmass_jetveto=HConfig.GetTH1D(Name+"_invmass_jetveto","invmass_jetveto",41,19,142,"m_{e#mu} (GeV)");
  invmass_vetos=HConfig.GetTH1D(Name+"_invmass_vetos","invmass_vetos",41,19,142,"m_{e#mu} (GeV)");
  invmass_only_object_id=HConfig.GetTH1D(Name+"_invmass_only_object_id","invmass_only_object_id",41,19,142,"m_{e#mu} (GeV)");

  invmass_zmass_m=HConfig.GetTH1D(Name+"_invmass_zmass_m","invmass_zmass_m_m",20,61,121,"m_{e#mu} (GeV)");
  invmass_ptbalance_m=HConfig.GetTH1D(Name+"_invmass_ptbalance_m","invmass_ptbalance_m",20,61,121,"m_{e#mu} (GeV)");
  invmass_met_m=HConfig.GetTH1D(Name+"_invmass_met_m","invmass_met_m",20,61,121,"m_{e#mu} (GeV)");
  invmass_jetveto_m=HConfig.GetTH1D(Name+"_invmass_jetveto_m","invmass_jetveto_m",20,61,121,"m_{e#mu} (GeV)");
  invmass_vetos_m=HConfig.GetTH1D(Name+"_invmass_vetos_m","invmass_vetos_m",20,61,121,"m_{e#mu} (GeV)");
  invmass_only_object_id_m=HConfig.GetTH1D(Name+"_invmass_only_object_id_m","invmass_only_object_id_m",20,61,121,"m_{e#mu} (GeV)");

  invmass_trilepton_only=HConfig.GetTH1D(Name+"_invmass_trilepton_only","invmass_trilepton_only",41,19,142,"m_{e#mu} (GeV)");
  invmass_charge_only=HConfig.GetTH1D(Name+"_invmass_charge_only","invmass_charge_only",41,19,142,"m_{e#mu} (GeV)");
  invmass_jetveto_only=HConfig.GetTH1D(Name+"_invmass_jetveto_only","invmass_jetveto_only",41,19,142,"m_{e#mu} (GeV)");
  invmass_met_only=HConfig.GetTH1D(Name+"_invmass_met_only","invmass_met_only",41,19,142,"m_{e#mu} (GeV)");
  invmass_ptbal_only=HConfig.GetTH1D(Name+"_invmass_ptbal_only","invmass_ptbal_only",41,19,142,"m_{e#mu} (GeV)");

  nm0_met=HConfig.GetTH1D(Name+"_nm0_met","nm0_met",30,0.,150.,"MET (GeV)");
  nm0_mvamet=HConfig.GetTH1D(Name+"_nm0_mvamet","nm0_mvamet",30,0.,150.,"MET (GeV)");
  nm0_onejet=HConfig.GetTH1D(Name+"_nm0_onejet","nm0_onejet",40,0.,200.,"p_{T}^{jet} (GeV)");
  nm0_ptbalance=HConfig.GetTH1D(Name+"_nm0_ptbalance","nm0_ptbalance",40,0.,200.,"p_{T}^{e#mu} (GeV)");
  nmm_met=HConfig.GetTH1D(Name+"_nmm_met","nmm_met",30,0.,150.,"MET (GeV)");
  nmm_mvamet=HConfig.GetTH1D(Name+"_nmm_mvamet","nmm_mvamet",30,0.,150.,"MET (GeV)");
  nmm_onejet=HConfig.GetTH1D(Name+"_nmm_onejet","nmm_onejet",40,0.,200.,"p_{T}^{jet} (GeV)");
  nmm_ptbalance=HConfig.GetTH1D(Name+"_nmm_ptbalance","nmm_ptbalance",40,0.,200.,"p_{T}^{e#mu} (GeV)");

  NPV=HConfig.GetTH1D(Name+"_NPV","NPV",60,0.,60.,"number of vertices");
  NPV3d=HConfig.GetTH1D(Name+"_NPV3d","NPV3d",60,0.,60.,"number of vertices");
  NPVfine=HConfig.GetTH1D(Name+"_NPVfine","NPVfine",60,0.,60.,"number of vertices");

  met_corr=HConfig.GetTH1D(Name+"_met_corr","met_corr",30,0.,150.,"MET (GeV)");
  met_uncorr=HConfig.GetTH1D(Name+"_met_uncorr","met_uncorr",30,0.,150.,"uncorrected MET (GeV)");
  onejet=HConfig.GetTH1D(Name+"_onejet","onejet",40,0.,200.,"p_{T}^{jet} (GeV)");
  onejet_eta=HConfig.GetTH1D(Name+"_onejet_eta","onejet_eta",40,-5.,5.,"#eta_{jet}");
  NbJets=HConfig.GetTH1D(Name+"_NbJets","NbJets",20,0,20,"number of b jets");
  NbJetsVtxL=HConfig.GetTH1D(Name+"_NbJetsVtxL","NbJetsVtxL",20,0,20,"number of b jets from vtx loose");
  NbJetsVtxM=HConfig.GetTH1D(Name+"_NbJetsVtxM","NbJetsVtxM",20,0,20,"number of b jets from vtx medium");
  NbJetsVtxT=HConfig.GetTH1D(Name+"_NbJetsVtxT","NbJetsVtxT",20,0,20,"number of b jets from vtx tight");

  ptbal_zoom=HConfig.GetTH1D(Name+"_ptbal_zoom","ptbal_zoom",15,0.,30.,"p_{T}^{e#mu} (GeV)");
  ht_pseudo=HConfig.GetTH1D(Name+"_ht_pseudo","ht_pseudo",200,0.,1000.,"#Sigma p_{T}^{jet} (GeV)");
  zmass_zoom=HConfig.GetTH1D(Name+"_zmass_zoom","zmass_zoom",5,82,97,"m_{e#mu} (GeV)");
  invmass_high=HConfig.GetTH1D(Name+"_invmass_high","invmass_high",50,60.,160.,"m_{e#mu}");
  ptsum=HConfig.GetTH1D(Name+"_ptsum","ptsum",40,0.,200.,"p_{T}^{e}+p_{T}^{#mu} (GeV)");
  ptsum_nm0=HConfig.GetTH1D(Name+"_ptsum_nm0","ptsum_nm0",40,0.,200.,"p_{T}^{e}+p_{T}^{#mu} (GeV)");
  mvamet=HConfig.GetTH1D(Name+"_mvamet","mvamet",30,0.,150.,"MET (GeV)");
  invmass_ptbalance_widerange=HConfig.GetTH1D(Name+"_invmass_ptbalance_widerange","invmass_ptbalance_widerange",41,19,142,"m_{e#mu} (GeV)");

  if(doPDFuncertainty){
	  pdf_w0=HConfig.GetTH1D(Name+"_pdf_w0","pdf_w0",nPDFmembers,0,nPDFmembers,"pdf member");
	  pdf_w1=HConfig.GetTH1D(Name+"_pdf_w1","pdf_w1",nPDFmembers,0,nPDFmembers,"pdf member");
  }

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(unsigned int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }

  // set category specific cuts
  if(categoryFlag == "ZtoEE"){
	  std::cout << "Using category: ZtoEE" << std::endl;
	  cut.at(NE) = 2;
	  cut.at(NMu) = 0;
  }else if(categoryFlag == "ZtoMuMu"){
	  std::cout << "Using category: ZtoMuMu" << std::endl;
	  cut.at(NE) = 0;
	  cut.at(NMu) = 2;
  }else if(categoryFlag == "ZtoEMu"){
	  std::cout << "Using category: ZtoEMu" << std::endl;
	  cut.at(NE) = 1;
	  cut.at(NMu) = 1;
  }
}




void  Tvariable_Base::Store_ExtraDist(){
 Extradist1d.push_back(&LeadPt);
 Extradist1d.push_back(&TrailPt);
 Extradist1d.push_back(&NJets);
 Extradist1d.push_back(&NJetsLoose);
 Extradist1d.push_back(&NJetsMedium);
 Extradist1d.push_back(&NJetsTight);
 Extradist1d.push_back(&PUJetId);
 Extradist1d.push_back(&etaLead);
 Extradist1d.push_back(&etaTrail);
 Extradist1d.push_back(&drleadtrail);
 Extradist1d.push_back(&deltaphi);
 Extradist1d.push_back(&ptbal);

 Extradist1d.push_back(&invmass_zmass);
 Extradist1d.push_back(&invmass_ptbalance);
 Extradist1d.push_back(&invmass_met);
 Extradist1d.push_back(&invmass_jetveto);
 Extradist1d.push_back(&invmass_vetos);
 Extradist1d.push_back(&invmass_only_object_id);

 Extradist1d.push_back(&invmass_zmass_m);
 Extradist1d.push_back(&invmass_ptbalance_m);
 Extradist1d.push_back(&invmass_met_m);
 Extradist1d.push_back(&invmass_jetveto_m);
 Extradist1d.push_back(&invmass_vetos_m);
 Extradist1d.push_back(&invmass_only_object_id_m);

 Extradist1d.push_back(&invmass_trilepton_only);
 Extradist1d.push_back(&invmass_charge_only);
 Extradist1d.push_back(&invmass_jetveto_only);
 Extradist1d.push_back(&invmass_met_only);
 Extradist1d.push_back(&invmass_ptbal_only);

 Extradist1d.push_back(&nm0_met);
 Extradist1d.push_back(&nm0_mvamet);
 Extradist1d.push_back(&nm0_onejet);
 Extradist1d.push_back(&nm0_ptbalance);
 Extradist1d.push_back(&nmm_met);
 Extradist1d.push_back(&nmm_mvamet);
 Extradist1d.push_back(&nmm_onejet);
 Extradist1d.push_back(&nmm_ptbalance);

 Extradist1d.push_back(&NPV);
 Extradist1d.push_back(&NPV3d);
 Extradist1d.push_back(&NPVfine);

 Extradist1d.push_back(&met_corr);
 Extradist1d.push_back(&met_uncorr);
 Extradist1d.push_back(&onejet);
 Extradist1d.push_back(&onejet_eta);
 Extradist1d.push_back(&NbJets);
 Extradist1d.push_back(&NbJetsVtxL);
 Extradist1d.push_back(&NbJetsVtxM);
 Extradist1d.push_back(&NbJetsVtxT);

 Extradist1d.push_back(&ptbal_zoom);
 Extradist1d.push_back(&ht_pseudo);
 Extradist1d.push_back(&zmass_zoom);
 Extradist1d.push_back(&invmass_high);
 Extradist1d.push_back(&ptsum);
 Extradist1d.push_back(&ptsum_nm0);
 Extradist1d.push_back(&mvamet);
 Extradist1d.push_back(&invmass_ptbalance_widerange);

 if(doPDFuncertainty){
	 Extradist1d.push_back(&pdf_w0);
	 Extradist1d.push_back(&pdf_w1);
 }

}

void  Tvariable_Base::doEvent(){
  if(verbose)std::cout << "Tvariable_Base::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}

  if(categoryFlag != "ZtoEE" && categoryFlag != "ZtoMuMu" && categoryFlag != "ZtoEMu"){
	  std::cout << "No category specified. Aborting event ..." << std::endl;
	  return;
  }

  ///////////////////////////////////////////////
  //
  // Set corrections
  //
  if(verbose)std::cout << "Set corrections" << std::endl;
  Ntp->SetMuonCorrections(mucorr);
  Ntp->SetElecCorrections(ecorr);
  Ntp->SetJetCorrections(jetcorr);

  ///////////////////////////////////////////////
  //
  // Trigger passed?
  //
  if(verbose)std::cout << "Trigger" << std::endl;
  value.at(TriggerOk)=0;
  if(categoryFlag=="ZtoEE" && Ntp->TriggerAccept("HLT_Ele27_WP80_v")) value.at(TriggerOk)=1;
  else if(categoryFlag=="ZtoMuMu" && Ntp->TriggerAccept("HLT_IsoMu24_eta2p1_v")) value.at(TriggerOk)=1;
  else if(categoryFlag=="ZtoEMu" && (Ntp->TriggerAccept("HLT_IsoMu24_eta2p1_v") || Ntp->TriggerAccept("HLT_Ele27_WP80_v"))) value.at(TriggerOk)=1;
  if(Ntp->GetMCID()==DataMCType::DY_emu_embedded)value.at(TriggerOk)=1;
  pass.at(TriggerOk)=(value.at(TriggerOk)==cut.at(TriggerOk));

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
	  if(Ntp->Muon_p4(i).Pt()>mu_ptlow
			  && fabs(Ntp->Muon_p4(i).Eta())<mu_eta
			  && vertex>=0
			  ){
		  if(Ntp->isTightMuon(i,vertex)
				  && Ntp->Muon_RelIso(i)<0.12
				  ){
			  GoodMuons.push_back(i);
		  }else if(doWWObjects
				  && isFakeMuon(i,vertex)
				  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll || Ntp->GetMCID()==DataMCType::DY_emu_embedded)
				  ){
			  Fakemuons.push_back(i);
			  GoodMuons.push_back(i);
		  }
	  }
  }

  value.at(NMu)=GoodMuons.size();
  pass.at(NMu)=(value.at(NMu)>=cut.at(NMu));

  unsigned int muidx1(999), muidx2(999);
  double hardestmu(0);
  if(GoodMuons.size()>1){
	  for(unsigned i=0;i<GoodMuons.size();i++){
		  if(Ntp->Muon_p4(GoodMuons.at(i)).Pt()>hardestmu){
			  hardestmu = Ntp->Muon_p4(GoodMuons.at(i)).Pt();
			  muidx1 = GoodMuons.at(i);
		  }
	  }
	  hardestmu = 0;
	  for(unsigned i=0;i<GoodMuons.size();i++){
		  if(GoodMuons.at(i)==muidx1) continue;
		  if(Ntp->Muon_p4(GoodMuons.at(i)).Pt()>hardestmu && Ntp->Muon_p4(GoodMuons.at(i)).Pt()<Ntp->Muon_p4(muidx1).Pt()){
			  hardestmu = Ntp->Muon_p4(GoodMuons.at(i)).Pt();
			  muidx2 = GoodMuons.at(i);
		  }
	  }
  }
  if(GoodMuons.size()==1){muidx1=GoodMuons.at(0);}

  ///////////////////////////////////////////////
  //
  // E Cuts
  //
  if(verbose) std::cout << "electrons cuts" << std::endl;
  std::vector<unsigned int> GoodElectrons;
  std::vector<unsigned int> Fakeelectrons;
  bool matchRecoMuon = false;

  // electron ID cuts (eta-dependent MVA or simple cut-based)
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  matchRecoMuon = false;
	  if(Ntp->Electron_p4(i).Pt()>e_ptlow
			  && fabs(Ntp->Electron_supercluster_eta(i))<e_eta
			  && vertex>=0
			  ){
		  // no overlapping reco muons
		  for(unsigned j=0;j<Ntp->NMuons();j++){
			  if(Ntp->Muon_p4(j).Pt()<3) continue;
			  if(fabs(Ntp->Muon_p4(j).Eta())>2.4) continue;
			  if(Ntp->Electron_p4(i).DeltaR(Ntp->Muon_p4(j))<0.3) matchRecoMuon = true;
		  }
		  if(matchRecoMuon) continue;
		  if(isWWElectron(i,vertex)
				  && Ntp->Electron_RelIsoDep04(i)<0.15
				  ){
			  GoodElectrons.push_back(i);
		  }else if(doWWObjects
				  && isFakeElectron(i,vertex)
				  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll || Ntp->GetMCID()==DataMCType::DY_emu_embedded)
				  ){
			  Fakeelectrons.push_back(i);
			  GoodElectrons.push_back(i);
		  }
	  }
  }

  value.at(NE)=GoodElectrons.size();
  pass.at(NE)=(value.at(NE)>=cut.at(NE));

  unsigned int eidx1(999), eidx2(999);
  double hardeste(0);
  if(GoodElectrons.size()>1){
	  for(unsigned i=0;i<GoodElectrons.size();i++){
		  if(Ntp->Electron_p4(GoodElectrons.at(i)).Pt()>hardeste){
			  hardeste = Ntp->Electron_p4(GoodElectrons.at(i)).Pt();
			  eidx1 = GoodElectrons.at(i);
		  }
	  }
	  hardeste = 0;
	  for(unsigned i=0;i<GoodElectrons.size();i++){
		  if(GoodElectrons.at(i)==eidx1) continue;
		  if(Ntp->Electron_p4(GoodElectrons.at(i)).Pt()>hardeste && Ntp->Electron_p4(GoodElectrons.at(i)).Pt()<Ntp->Electron_p4(eidx1).Pt()){
			  hardeste = Ntp->Electron_p4(GoodElectrons.at(i)).Pt();
			  eidx2 = GoodElectrons.at(i);
		  }
	  }
  }
  if(GoodElectrons.size()==1){eidx1=GoodElectrons.at(0);}

  ///////////////////////////////////////////////
  //
  // pt thresholds
  //
  if(verbose) std::cout << "Setting pt thresholds" << std::endl;
  bool passembed = false;
  bool mutrig = false;
  value.at(ptthreshold)=0;
  if(categoryFlag == "ZtoEE"){
	  if(eidx1 != 999 && eidx2 != 999){
		  value.at(ptthreshold) = 1;
		  if(Ntp->Electron_p4(eidx1).Pt()<lead_pt || fabs(Ntp->Electron_supercluster_eta(eidx1))>lead_eta) value.at(ptthreshold) = 0;
		  if(Ntp->Electron_p4(eidx2).Pt()<trail_pt || fabs(Ntp->Electron_supercluster_eta(eidx2))>trail_eta) value.at(ptthreshold) = 0;
		  if(value.at(ptthreshold) == 1 && Ntp->GetMCID() == DataMCType::DY_emu_embedded) passembed = true;
		  if(!Ntp->matchTrigger(Ntp->Electron_p4(eidx1),0.2,"HLT_Ele27_WP80_v","electron") && !Ntp->matchTrigger(Ntp->Electron_p4(eidx2),0.2,"HLT_Ele27_WP80_v","electron")) value.at(ptthreshold) = 0;
		  if(passembed) value.at(ptthreshold) = 1;
	  }
  }else if(categoryFlag == "ZtoMuMu"){
	  if(muidx1 != 999 && muidx2 != 999){
		  value.at(ptthreshold) = 1;
		  if(Ntp->Muon_p4(muidx1).Pt()<lead_pt || fabs(Ntp->Muon_p4(muidx1).Eta())>lead_eta) value.at(ptthreshold) = 0;
		  if(Ntp->Muon_p4(muidx2).Pt()<trail_pt || fabs(Ntp->Muon_p4(muidx2).Eta())>trail_eta) value.at(ptthreshold) = 0;
		  if(value.at(ptthreshold) == 1 && Ntp->GetMCID() == DataMCType::DY_emu_embedded) passembed = true;
		  if(!Ntp->matchTrigger(Ntp->Muon_p4(muidx1),0.2,"HLT_IsoMu24_eta2p1_v","muon") && !Ntp->matchTrigger(Ntp->Muon_p4(muidx2),0.2,"HLT_IsoMu24_eta2p1_v","muon")) value.at(ptthreshold) = 0;
		  if(passembed) value.at(ptthreshold) = 1;
	  }
  }else if(categoryFlag == "ZtoEMu"){
	  if(muidx1 != 999 && eidx1 != 999){
		  value.at(ptthreshold) = 1;
		  if(Ntp->Muon_p4(muidx1).Pt() > Ntp->Electron_p4(eidx1).Pt()){
			  if(Ntp->Muon_p4(muidx1).Pt()<lead_pt || fabs(Ntp->Muon_p4(muidx1).Eta())>lead_eta) value.at(ptthreshold) = 0;
			  if(Ntp->Electron_p4(eidx1).Pt()<trail_pt || fabs(Ntp->Electron_supercluster_eta(eidx1))>trail_eta) value.at(ptthreshold) = 0;
		  }else{
			  if(Ntp->Electron_p4(eidx1).Pt()<lead_pt || fabs(Ntp->Electron_supercluster_eta(eidx1))>lead_eta) value.at(ptthreshold) = 0;
			  if(Ntp->Muon_p4(muidx1).Pt()<trail_pt || fabs(Ntp->Muon_p4(muidx1).Eta())>trail_eta) value.at(ptthreshold) = 0;
		  }
		  if(value.at(ptthreshold) == 1 && Ntp->GetMCID() == DataMCType::DY_emu_embedded) passembed = true;
		  if(!Ntp->matchTrigger(Ntp->Muon_p4(muidx1),0.2,"HLT_IsoMu24_eta2p1_v","muon") && !Ntp->matchTrigger(Ntp->Electron_p4(eidx1),0.2,"HLT_Ele27_WP80_v","electron") && Ntp->GetMCID()!=DataMCType::DY_emu_embedded) value.at(ptthreshold)=0;
		  if(Ntp->matchTrigger(Ntp->Muon_p4(muidx1),0.2,"HLT_IsoMu24_eta2p1_v","muon") && Ntp->matchTrigger(Ntp->Electron_p4(eidx1),0.2,"HLT_Ele27_WP80_v","electron")) mutrig = true;
		  if(passembed) value.at(ptthreshold) = 1;
	  }
  }
  pass.at(ptthreshold)=(value.at(ptthreshold)==cut.at(ptthreshold));

  ///////////////////////////////////////////////
  //
  // m(ll)
  //
  if(verbose) std::cout << "m(ll)" << std::endl;
  //value.at(mll)=mmin+1;
  value.at(mll)=61;
  if(categoryFlag == "ZtoEE"){
	  if(eidx1 != 999 && eidx2 != 999){
		  value.at(mll) = (Ntp->Electron_p4(eidx1) + Ntp->Electron_p4(eidx2)).M();
	  }
  }else if(categoryFlag == "ZtoMuMu"){
	  if(muidx1 != 999 && muidx2 != 999){
		  value.at(mll) = (Ntp->Muon_p4(muidx1) + Ntp->Muon_p4(muidx2)).M();
	  }
  }else if(categoryFlag == "ZtoEMu"){
	  if(muidx1 !=999 && eidx1 != 999){
		  value.at(mll) = (Ntp->Muon_p4(muidx1) + Ntp->Electron_p4(eidx1)).M();
	  }
  }
  //pass.at(mll)=(value.at(mll)>cut.at(mll));
  pass.at(mll)=(value.at(mll)>=60. && value.at(mll)<120.);

  ///////////////////////////////////////////////
  //
  // Tri Lepton Veto
  //
  if(verbose)std::cout << "trilepton veto" << std::endl;
  unsigned int trilep(0);
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(i == muidx1 || i == muidx2) continue;
	  if(vertex<0) continue;
	  if(Ntp->Muon_p4(i).Pt()<10) continue;
	  if(fabs(Ntp->Muon_p4(i).Eta())>2.4) continue;
	  if(!Ntp->isTightMuon(i,vertex)) continue;
	  if(Ntp->Muon_RelIso(i)>0.3) continue;
	  trilep++;
  }
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  if(i == eidx1 || i == eidx2) continue;
	  if(vertex<0) continue;
	  if(Ntp->Electron_p4(i).Pt()<10) continue;
	  if(fabs(Ntp->Electron_supercluster_eta(i))>2.5) continue;
	  if(!isWWElectron(i,vertex)) continue;
	  if(Ntp->Electron_RelIsoDep04(i)>0.3) continue;
	  trilep++;
  }
  if(categoryFlag == "ZtoEE"){
	  if(muidx1 != 999 || muidx2 !=999) trilep++;
  }else if(categoryFlag == "ZtoMuMu"){
	  if(eidx1 != 999 || eidx2 !=999) trilep++;
  }else if(categoryFlag == "ZtoEMu"){
	  if(muidx2 != 999 || eidx2 != 999) trilep++;
  }
  value.at(triLeptonVeto)=trilep;
  pass.at(triLeptonVeto)=(value.at(triLeptonVeto)==cut.at(triLeptonVeto));

  ///////////////////////////////////////////////
  //
  // charge cut
  //
  if(verbose) std::cout << "Charge cut" << std::endl;
  value.at(charge)=-1;
  if(categoryFlag == "ZtoEE"){
	  if(eidx1 != 999 && eidx2 != 999){
		  value.at(charge) = Ntp->Electron_Charge(eidx1) * Ntp->Electron_Charge(eidx2);
		  if(Ntp->Electron_Charge(eidx1) == -999 || Ntp->Electron_Charge(eidx2) == -999) value.at(charge) = 1;
	  }
  }else if(categoryFlag == "ZtoMuMu"){
	  if(muidx1 != 999 && muidx2 != 999){
		  value.at(charge) = Ntp->Muon_Charge(muidx1) * Ntp->Muon_Charge(muidx2);
		  if(Ntp->Muon_Charge(muidx1) == -999 || Ntp->Muon_Charge(muidx2) == -999) value.at(charge) = 1;
	  }
  }else if(categoryFlag == "ZtoEMu"){
	  if(muidx1 != 999 && eidx1 != 999){
		  value.at(charge) = Ntp->Muon_Charge(muidx1) * Ntp->Electron_Charge(eidx1);
		  if(Ntp->Muon_Charge(muidx1) == -999 || Ntp->Electron_Charge(eidx1) == -999) value.at(charge) = 1;
	  }
  }
  pass.at(charge)=(value.at(charge)<cut.at(charge));

  ////////////////////////////////////////////////
  //
  // QCD
  //
  if(verbose) std::cout << "QCD" << std::endl;
  fakeRate = 1.;
  bool fakemu1(false), fakemu2(false);
  bool fakee1(false), fakee2(false);
  if(doWWObjects){
	  if(categoryFlag == "ZtoMuMu" || categoryFlag == "ZtoEMu"){
		  for(unsigned i=0;i<Fakemuons.size();i++){
			  if(Fakemuons.at(i) == muidx1){
				  fakemu1 = true;
				  fakeRateMu1 = Fakerate(Ntp->Muon_p4(muidx1).Pt(),Ntp->Muon_p4(muidx1).Eta(),MuonFakeRate15);
				  if(doFakeRateUncertainty){
					  if(upwardUncertainty) fakeRateMu1 = Fakerate(Ntp->Muon_p4(muidx1).Pt(),Ntp->Muon_p4(muidx1).Eta(),MuonFakeRate30);
					  else fakeRateMu1 = Fakerate(Ntp->Muon_p4(muidx1).Pt(),Ntp->Muon_p4(muidx1).Eta(),MuonFakeRate5);
				  }
			  }else if(Fakemuons.at(i) == muidx2){
				  fakemu2 = true;
				  fakeRateMu2 = Fakerate(Ntp->Muon_p4(muidx2).Pt(),Ntp->Muon_p4(muidx2).Eta(),MuonFakeRate15);
				  if(doFakeRateUncertainty){
					  if(upwardUncertainty) fakeRateMu2 = Fakerate(Ntp->Muon_p4(muidx2).Pt(),Ntp->Muon_p4(muidx2).Eta(),MuonFakeRate30);
					  else fakeRateMu2 = Fakerate(Ntp->Muon_p4(muidx2).Pt(),Ntp->Muon_p4(muidx2).Eta(),MuonFakeRate5);
				  }
			  }
			  if((fakemu1 && categoryFlag != "ZtoMuMu") || (fakemu1 && fakemu2)) break;
		  }
	  }
	  if(categoryFlag == "ZtoEE" || categoryFlag == "ZtoEMu"){
		  for(unsigned i=0;i<Fakeelectrons.size();i++){
			  if(Fakeelectrons.at(i) == eidx1){
				  fakee1 = true;
				  fakeRateE1 = Fakerate(Ntp->Electron_p4(eidx1).Pt(),Ntp->Electron_supercluster_eta(eidx1),ElectronFakeRate35);
				  if(doFakeRateUncertainty){
					  if(upwardUncertainty) fakeRateE1 = Fakerate(Ntp->Electron_p4(eidx1).Pt(),Ntp->Electron_supercluster_eta(eidx1),ElectronFakeRate50);
					  else fakeRateE1 = Fakerate(Ntp->Electron_p4(eidx1).Pt(),Ntp->Electron_supercluster_eta(eidx1),ElectronFakeRate20);
				  }
			  }else if(Fakeelectrons.at(i) == eidx2){
				  fakee2 = true;
				  fakeRateE2 = Fakerate(Ntp->Electron_p4(eidx2).Pt(),Ntp->Electron_supercluster_eta(eidx2),ElectronFakeRate35);
				  if(doFakeRateUncertainty){
					  if(upwardUncertainty) fakeRateE2 = Fakerate(Ntp->Electron_p4(eidx2).Pt(),Ntp->Electron_supercluster_eta(eidx2),ElectronFakeRate50);
					  else fakeRateE2 = Fakerate(Ntp->Electron_p4(eidx2).Pt(),Ntp->Electron_supercluster_eta(eidx2),ElectronFakeRate20);
				  }
			  }
			  if((fakee1 && categoryFlag != "ZtoEE") || (fakee1 && fakee2)) break;
		  }
	  }
	  if(pass.at(charge)
			  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll || Ntp->GetMCID()==DataMCType::DY_emu_embedded)
			  ){
		  if(fakemu1 || fakemu2 || fakee1 || fakee2) fakeRate = 0.;
	  }
	  if(!pass.at(charge)
			  && Ntp->isData()
			  ){
		  if(categoryFlag == "ZtoEE"){
			  if(fakee1 && !fakee2){
				  fakeRate = fakeRateE1;
				  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
				  pass.at(charge)=true;
			  }else if(fakee2 && !fakee1){
				  fakeRate = fakeRateE2;
				  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
				  pass.at(charge)=true;
			  }else if(fakee1 && fakee2){
				  fakeRate = fakeRateE1 * fakeRateE2;
				  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
				  pass.at(charge)=true;
			  }
		  }else if(categoryFlag == "ZtoMuMu"){
			  if(fakemu1 && !fakemu2){
				  fakeRate = fakeRateMu1;
				  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
				  pass.at(charge)=true;
			  }else if(fakemu2 && !fakemu1){
				  fakeRate = fakeRateMu2;
				  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
				  pass.at(charge)=true;
			  }else if(fakemu1 && fakemu2){
				  fakeRate = fakeRateMu1 * fakeRateMu2;
				  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
				  pass.at(charge)=true;
			  }
		  }else if(categoryFlag == "ZtoEMu"){
			  if(fakemu1 && !fakee1){
				  fakeRate = fakeRateMu1;
				  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
				  pass.at(charge)=true;
			  }else if(fakee1 && !fakemu1){
				  fakeRate = fakeRateE1;
				  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
				  pass.at(charge)=true;
			  }else if(fakemu1 && fakee1){
				  fakeRate = fakeRateMu1 * fakeRateE1;
				  if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
				  pass.at(charge)=true;
			  }
		  }
	  }
	  if(!pass.at(charge)
			  && (Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll || Ntp->GetMCID()==DataMCType::DY_emu_embedded)
			  ){
		  if(categoryFlag == "ZtoEE"){
			  if(fakee1 && !fakee2){
				  fakeRate = -fakeRateE1;
				  pass.at(charge)=true;
			  }else if(fakee2 && !fakee1){
				  fakeRate = -fakeRateE2;
				  pass.at(charge)=true;
			  }else if(fakee1 && fakee2){
				  fakeRate = -fakeRateE1 * fakeRateE2;
				  pass.at(charge)=true;
			  }
		  }else if(categoryFlag == "ZtoMuMu"){
			  if(fakemu1 && !fakemu2){
				  fakeRate = -fakeRateMu1;
				  pass.at(charge)=true;
			  }else if(fakemu2 && !fakemu1){
				  fakeRate = -fakeRateMu2;
				  pass.at(charge)=true;
			  }else if(fakemu1 && fakemu2){
				  fakeRate = -fakeRateMu1 * fakeRateMu2;
				  pass.at(charge)=true;
			  }
		  }else if(categoryFlag == "ZtoEMu"){
			  if(fakemu1 && !fakee1){
				  fakeRate = -fakeRateMu1;
				  pass.at(charge)=true;
			  }else if(fakee1 && !fakemu1){
				  fakeRate = -fakeRateE1;
				  pass.at(charge)=true;
			  }else if(fakemu1 && fakee1){
				  fakeRate = -fakeRateMu1 * fakeRateE1;
				  pass.at(charge)=true;
			  }
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
	  if(Ntp->PFJet_p4(i).Pt()<20) continue;
	  if(fabs(Ntp->PFJet_p4(i).Eta())>jet_eta) continue;
	  if(!Ntp->isJetID(i)) continue;
	  if(muidx1!=999 && Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(muidx1))<0.3) continue;
	  if(muidx2!=999 && Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(muidx2))<0.3) continue;
	  if(eidx1!=999 && Ntp->PFJet_p4(i).DeltaR(Ntp->Electron_p4(eidx1))<0.3) continue;
	  if(eidx2!=999 && Ntp->PFJet_p4(i).DeltaR(Ntp->Electron_p4(eidx2))<0.3) continue;
	  // find jets from vertex: use pileup jet id for jets with pt>20 GeV
	  if(Ntp->PFJet_PUJetID_tightWP(i)>0.5) jetsfromvtx.push_back(i);
  }

  if(verbose)std::cout<< "Find highest pt jet" << std::endl;
  int firstjet_idx=-1;
  int secondjet_idx=-1;
  double initialpt=0.;

  // loop over jets from selected vertex & find the two jets with the highest pt
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()>initialpt){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt();
		  firstjet_idx=jetsfromvtx.at(i);
	  }
  }
  initialpt=0.;
  for(unsigned i=0;i<jetsfromvtx.size();i++){
	  if(jetsfromvtx.size()>1 && firstjet_idx!=-1
			  && Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()>initialpt
			  && Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt()<Ntp->PFJet_p4(firstjet_idx).Pt()
			  ){
		  initialpt=Ntp->PFJet_p4(jetsfromvtx.at(i)).Pt();
		  secondjet_idx=jetsfromvtx.at(i);
	  }
  }

  if(verbose)std::cout << "applying veto" << std::endl;

  value.at(oneJet)=0;
  if(jetsfromvtx.size()>0 && firstjet_idx!=-1){
	  value.at(oneJet)=Ntp->PFJet_p4(firstjet_idx).Pt();
  }
  pass.at(oneJet)=(value.at(oneJet)<cut.at(oneJet));

  ///////////////////////////////////////////////
  //
  // Mt Mu cut
  //
  if(verbose) std::cout << "Mt Mu cut" << std::endl;
  value.at(MET)=0.;
  value.at(MET) = Ntp->MET_CorrT0pcT1Txy_et();
  if(doMetUncertainty && !Ntp->isData() && Ntp->GetMCID() != DataMCType::DY_emu_embedded) value.at(MET) = Ntp->MET_Type1CorrUnClusteredUp_et();
  pass.at(MET)=(value.at(MET)<cut.at(MET));

  ///////////////////////////////////////////////
  //
  // Pt balance cut
  //
  if(verbose) std::cout << "pt balance cut" << std::endl;
  value.at(ptBalance)=0.;
  if(categoryFlag == "ZtoEE"){
	  if(eidx1 != 999 && eidx2 != 999){
		  value.at(ptBalance) = (Ntp->Electron_p4(eidx1) + Ntp->Electron_p4(eidx2)).Pt();
	  }
  }else if(categoryFlag == "ZtoMuMu"){
	  if(muidx1 != 999 && muidx2 != 999){
		  value.at(ptBalance) = (Ntp->Muon_p4(muidx1) + Ntp->Muon_p4(muidx2)).Pt();
	  }
  }else if(categoryFlag == "ZtoEMu"){
	  if(muidx1 != 999 && eidx1 != 999){
		  value.at(ptBalance) = (Ntp->Muon_p4(muidx1) + Ntp->Electron_p4(eidx1)).Pt();
	  }
  }
  pass.at(ptBalance)=(value.at(ptBalance)<cut.at(ptBalance));

  ///////////////////////////////////////////////
  //
  // Invariant mass cut
  //
  if(verbose) std::cout << "Invariant mass cut" << std::endl;
  value.at(ZMassmax)=zmax-1.;
  value.at(ZMassmin)=zmin+1.;
  if(categoryFlag == "ZtoEE"){
	  if(eidx1 != 999 && eidx2 != 999){
		  value.at(ZMassmax) = (Ntp->Electron_p4(eidx1) + Ntp->Electron_p4(eidx2)).M();
		  value.at(ZMassmin) = (Ntp->Electron_p4(eidx1) + Ntp->Electron_p4(eidx2)).M();
	  }
  }else if(categoryFlag == "ZtoMuMu"){
	  if(muidx1 != 999 && muidx2 != 999){
		  value.at(ZMassmax) = (Ntp->Muon_p4(muidx1) + Ntp->Muon_p4(muidx2)).M();
		  value.at(ZMassmin) = (Ntp->Muon_p4(muidx1) + Ntp->Muon_p4(muidx2)).M();
	  }
  }else if(categoryFlag == "ZtoEMu"){
	  if(muidx1 != 999 && eidx1 != 999){
		  value.at(ZMassmax) = (Ntp->Muon_p4(muidx1) + Ntp->Electron_p4(eidx1)).M();
		  value.at(ZMassmin) = (Ntp->Muon_p4(muidx1) + Ntp->Electron_p4(eidx1)).M();
	  }
  }
  pass.at(ZMassmax)=(value.at(ZMassmax)<cut.at(ZMassmax));
  pass.at(ZMassmin)=(value.at(ZMassmin)>=cut.at(ZMassmin));

  ///////////////////////////////////////////////
  //
  // Weights
  //
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1);
  if(!Ntp->isData()){
	// pileup weights only for MC
	if(Ntp->GetMCID()!=DataMCType::DY_emu_embedded){
		if(!doPileupUncertainty) w*=Ntp->PUWeight();
	    else{
	    	if(upwardUncertainty) w*=Ntp->PUWeight_p5();
	    	else w*= Ntp->PUWeight3D_m5();
	    }
	}
	// weights just for embedded samples
    if(Ntp->GetMCID()==DataMCType::DY_emu_embedded && (categoryFlag == "ZtoEE" || categoryFlag == "ZtoEMu")){
    	w*=Ntp->EmbeddedWeight();
    	if(pass.at(NE)){
    		w*=RSF->ElectronEmbedding2012(Ntp->Electron_p4(eidx1).Pt(),Ntp->Electron_supercluster_eta(eidx1));
    		if(categoryFlag == "ZtoEE") w*=RSF->ElectronEmbedding2012(Ntp->Electron_p4(eidx2).Pt(),Ntp->Electron_supercluster_eta(eidx2));
    	}
    }
    // weights for MC & embedded samples
    w*=fakeRate;
    if(categoryFlag == "ZtoEE" || categoryFlag == "ZtoEMu"){
    	if(pass.at(NE)){
    		double eidunc1(0), eidunc2(0), erecunc1(0), erecunc2(0);
        	if(doElectronIdUncertainty){
        		eidunc1 = RSF->ElectronIdTrigUnc2012(Ntp->Electron_p4(eidx1).Pt(),Ntp->Electron_supercluster_eta(eidx1));
        		erecunc1 = RSF->ElectronReconstructionUnc2012(Ntp->Electron_p4(eidx1).Pt(),Ntp->Electron_supercluster_eta(eidx1));
        		if(categoryFlag == "ZtoEE"){
        			eidunc2 = RSF->ElectronIdTrigUnc2012(Ntp->Electron_p4(eidx2).Pt(),Ntp->Electron_supercluster_eta(eidx2));
            		erecunc2 = RSF->ElectronReconstructionUnc2012(Ntp->Electron_p4(eidx2).Pt(),Ntp->Electron_supercluster_eta(eidx2));
        		}
        		if(!upwardUncertainty){
        			eidunc1 *= -1;
        			eidunc2 *= -1;
        			erecunc1 *= -1;
        			erecunc2 *= -1;
        		}
        	}
        	w*=(RSF->ElectronIdTrig2012(Ntp->Electron_p4(eidx1).Pt(),Ntp->Electron_supercluster_eta(eidx1))+eidunc1);
    		w*=(RSF->ElectronReconstruction2012(Ntp->Electron_p4(eidx1).Pt(),Ntp->Electron_supercluster_eta(eidx1))+erecunc1);
    		if(categoryFlag == "ZtoEE"){
    			w*=(RSF->ElectronIdTrig2012(Ntp->Electron_p4(eidx2).Pt(),Ntp->Electron_supercluster_eta(eidx2))+eidunc2);
    			w*=(RSF->ElectronReconstruction2012(Ntp->Electron_p4(eidx2).Pt(),Ntp->Electron_supercluster_eta(eidx2))+erecunc2);
    		}
    	}
    }
    if(categoryFlag == "ZtoMuMu" || categoryFlag == "ZtoEMu"){
    	if(pass.at(NMu)){
    		// for systematics: add systematic uncertainty of 0.5%(1.5%) when pt>20(<20) for Id and 0.2% for isolation when pt>20 to statistical in quadrature.
    		double muidunc1(0), muidunc2(0), muisounc1(0), muisounc2(0), mutrkunc1(0), mutrkunc2(0);
    		if(doMuonIdUncertainty){
    			if(Ntp->Muon_p4(muidx1).Pt()>20){
    				muidunc1 = sqrt(pow(RSF->MuonIdUncTight2012(Ntp->Muon_p4(muidx1)),2)+pow(RSF->MuonIdTight2012(Ntp->Muon_p4(muidx1))*0.005,2));
    				muisounc1 = sqrt(pow(RSF->MuonIsoUncTight2012(Ntp->Muon_p4(muidx1)),2)+pow(0.002,2));
    			}else{
    				muidunc1 = sqrt(pow(RSF->MuonIdUncTight2012(Ntp->Muon_p4(muidx1)),2)+pow(RSF->MuonIdTight2012(Ntp->Muon_p4(muidx1))*0.015,2));
    				muisounc1 = RSF->MuonIsoUncTight2012(Ntp->Muon_p4(muidx1));
    			}
    			if(categoryFlag == "ZtoMuMu"){
    				if(Ntp->Muon_p4(muidx2).Pt()>20){
        				muidunc2 = sqrt(pow(RSF->MuonIdUncTight2012(Ntp->Muon_p4(muidx2)),2)+pow(RSF->MuonIdTight2012(Ntp->Muon_p4(muidx2))*0.005,2));
        				muisounc2 = sqrt(pow(RSF->MuonIsoUncTight2012(Ntp->Muon_p4(muidx2)),2)+pow(0.002,2));
        			}else{
        				muidunc2 = sqrt(pow(RSF->MuonIdUncTight2012(Ntp->Muon_p4(muidx2)),2)+pow(RSF->MuonIdTight2012(Ntp->Muon_p4(muidx2))*0.015,2));
        				muisounc2 = RSF->MuonIsoUncTight2012(Ntp->Muon_p4(muidx2));
        			}
    			}
    			mutrkunc1 = RSF->TrackingEfficiencyUnc2012(Ntp->Muon_p4(muidx1));
    			if(categoryFlag == "ZtoMuMu") mutrkunc2 = RSF->TrackingEfficiencyUnc2012(Ntp->Muon_p4(muidx2));
    			if(!upwardUncertainty){
    				muidunc1 *= -1;
    				muidunc2 *= -1;
    				muisounc1 *= -1;
    				muisounc2 *= -1;
    				mutrkunc1 *= -1;
    				mutrkunc2 *= -1;
    			}
    		}
    		w*=(RSF->MuonIdTight2012(Ntp->Muon_p4(muidx1))+muidunc1);
    		w*=(RSF->MuonIsoTight2012(Ntp->Muon_p4(muidx1))+muisounc1);
    		w*=(RSF->TrackingEfficiency2012(Ntp->Muon_p4(muidx1))+mutrkunc1);
    		if(categoryFlag == "ZtoMuMu"){
    			w*=(RSF->MuonIdTight2012(Ntp->Muon_p4(muidx2))+muidunc2);
        		w*=(RSF->MuonIsoTight2012(Ntp->Muon_p4(muidx2))+muisounc2);
        		w*=(RSF->TrackingEfficiency2012(Ntp->Muon_p4(muidx2))+mutrkunc2);
    		}
    	}
    }
    if(pass.at(TriggerOk)
    		&& pass.at(NMu)
    		&& pass.at(NE)
    		){
    	// for systematics: double mu uncertainty 0.5%, double e uncertainty 0.2%
    	double trigunc(0);
    	if(categoryFlag=="ZtoEE"){
    		//todo: implement trigger efficiencies
    		if(Ntp->matchTrigger(Ntp->Electron_p4(eidx1),0.2,"HLT_Ele27_WP80_v","electron")){
    			if(doTriggerUncertainty){
    				trigunc = RSF->HiggsWW_EMu_SingleEle(Ntp->Electron_p4(eidx1).Et(),Ntp->Electron_supercluster_eta(eidx1))*0.002;
    				if(!upwardUncertainty) trigunc *= -1;
    			}
    			w*=(RSF->HiggsWW_EMu_SingleEle(Ntp->Electron_p4(eidx1).Et(),Ntp->Electron_supercluster_eta(eidx1))+trigunc);
    		}else if(Ntp->matchTrigger(Ntp->Electron_p4(eidx2),0.2,"HLT_Ele27_WP80_v","electron")){
    			if(doTriggerUncertainty){
    				trigunc = RSF->HiggsWW_EMu_SingleEle(Ntp->Electron_p4(eidx2).Et(),Ntp->Electron_supercluster_eta(eidx2))*0.002;
    				if(!upwardUncertainty) trigunc *= -1;
    			}
    			w*=(RSF->HiggsWW_EMu_SingleEle(Ntp->Electron_p4(eidx2).Et(),Ntp->Electron_supercluster_eta(eidx2))+trigunc);
    		}
    	}else if(categoryFlag=="ZtoMuMu"){
    		//todo: implement trigger efficiencies
    		if(Ntp->matchTrigger(Ntp->Muon_p4(muidx1),0.2,"HLT_IsoMu24_eta2p1_v","muon")){
    			if(doTriggerUncertainty){
    				trigunc = RSF->IsoMu24_eta2p1_unc(Ntp->Muon_p4(muidx1));
    				if(!upwardUncertainty) trigunc *= -1;
    			}
    			w*=(RSF->IsoMu24_eta2p1(Ntp->Muon_p4(muidx1))+trigunc);
    		}else if(Ntp->matchTrigger(Ntp->Muon_p4(muidx2),0.2,"HLT_IsoMu24_eta2p1_v","muon")){
    			if(doTriggerUncertainty){
    				trigunc = RSF->IsoMu24_eta2p1_unc(Ntp->Muon_p4(muidx2));
    				if(!upwardUncertainty) trigunc *= -1;
    			}
    			w*=(RSF->IsoMu24_eta2p1(Ntp->Muon_p4(muidx2))+trigunc);
    		}
    	}else if(categoryFlag=="ZtoEMu"){
    		if(mutrig){
    			if(doTriggerUncertainty){
    				if(doTriggerUncertainty){
        				trigunc = RSF->IsoMu24_eta2p1_unc(Ntp->Muon_p4(muidx1));
        				if(!upwardUncertainty) trigunc *= -1;
        			}
        			w*=(RSF->IsoMu24_eta2p1(Ntp->Muon_p4(muidx1))+trigunc);
    			}
    		}else{
    			if(doTriggerUncertainty){
    				if(doTriggerUncertainty){
        				trigunc = RSF->HiggsWW_EMu_SingleEle(Ntp->Electron_p4(eidx1).Et(),Ntp->Electron_supercluster_eta(eidx1))*0.002;
        				if(!upwardUncertainty) trigunc *= -1;
        			}
        			w*=(RSF->HiggsWW_EMu_SingleEle(Ntp->Electron_p4(eidx1).Et(),Ntp->Electron_supercluster_eta(eidx1))+trigunc);
    			}
    		}
    	}
    }
    if(pass.at(NMu)
    		&& pass.at(NE)
    		&& !useMadgraphZ
    		){
    	//if(Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau)w*=PowhegReweight((Ntp->Muon_p4(muidx)+Ntp->Electron_p4(eidx)).Pt());
    }
    if(verbose)std::cout << "void  Tvariable_Base::doEvent() k" << w << " " << wobs << std::endl;
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
	  for(int member=0;member<nPDFmembers;++member){
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

  TLorentzVector leadlep(0.,0.,0.,0.);
  TLorentzVector traillep(0.,0.,0.,0.);
  bool leadele(false), trailele(false);

  if(categoryFlag=="ZtoEE"){
	  if(eidx1!=999 && eidx2!=999){
		  leadlep = Ntp->Electron_p4(eidx1);
		  traillep = Ntp->Electron_p4(eidx2);
		  leadele = true;
		  trailele = true;
	  }
  }else if(categoryFlag=="ZtoMuMu"){
	  if(muidx1!=999 && muidx2!=999){
		  leadlep = Ntp->Muon_p4(muidx1);
		  traillep = Ntp->Muon_p4(muidx2);
	  }
  }else if(categoryFlag=="ZtoEMu"){
	  if(eidx1!=999 && muidx1!=999){
		  if(Ntp->Electron_p4(eidx1).Pt()>Ntp->Muon_p4(muidx1).Pt()){
			  leadlep = Ntp->Electron_p4(eidx1);
			  traillep = Ntp->Muon_p4(muidx1);
			  leadele = true;
		  }else{
			  leadlep = Ntp->Muon_p4(muidx1);
			  traillep = Ntp->Electron_p4(eidx1);
			  trailele = true;
		  }
	  }
  }

  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(NMu)
		  && pass.at(NE)
		  && pass.at(ptthreshold)
		  && pass.at(mll)
		  && pass.at(triLeptonVeto)
		  && pass.at(charge)
		  	  	  ){
	  if((Ntp->GetMCID()==30 || Ntp->GetMCID()==33) && !(fakemu1 || fakemu2 || fakee1 || fakee2)){
		  printf("Eventnumber: %i, lead pt: %f, lead eta, %f, trail pt: %f, trail eta: %f\n",Ntp->EventNumber(),leadlep.Pt(),leadlep.Eta(),traillep.Pt(),traillep.Eta());
	  }
  }

  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(NMu)
		  && pass.at(NE)
		  && pass.at(ptthreshold)
		  && pass.at(charge)
		  && pass.at(triLeptonVeto)
		  && pass.at(oneJet)
		  && pass.at(MET)
		  && pass.at(ptBalance)
                  ){
		  	if((leadlep+traillep).M()>60.) invmass_high.at(t).Fill((leadlep+traillep).M(),w);
		  	invmass_ptbalance_widerange.at(t).Fill((leadlep+traillep).M(),w);
		  }
  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(NMu)
		  && pass.at(NE)
		  && pass.at(ptthreshold)
		  && pass.at(mll)
		  ){
	  // often needed variables
	  double m = (leadlep+traillep).M();
	  double dp = leadlep.DeltaPhi(traillep)/TMath::Pi();
	  if(dp<0)dp+=2;

	  // electron related histograms
	  LeadPt.at(t).Fill(leadlep.Pt(),w);
	  if(leadele) etaLead.at(t).Fill(Ntp->Electron_supercluster_eta(eidx1),w);
	  else etaLead.at(t).Fill(leadlep.Eta(),w);

	  // muon related histograms
	  TrailPt.at(t).Fill(traillep.Pt(),w);
	  if(trailele){
		  if(categoryFlag == "ZtoEE") etaTrail.at(t).Fill(Ntp->Electron_supercluster_eta(eidx2),w);
		  else etaTrail.at(t).Fill(Ntp->Electron_supercluster_eta(eidx1),w);
	  }
	  else etaTrail.at(t).Fill(traillep.Eta(),w);

	  // histograms related to combination
	  drleadtrail.at(t).Fill(leadlep.DeltaR(traillep),w);
	  met_corr.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),w);
	  met_uncorr.at(t).Fill(Ntp->MET_Uncorr_et(),w);
	  mvamet.at(t).Fill(Ntp->MET_CorrMVA_et(),w);
	  deltaphi.at(t).Fill(dp,w);
	  ptbal.at(t).Fill((leadlep+traillep).Pt(),w);
	  ptbal_zoom.at(t).Fill((leadlep+traillep).Pt(),w);
	  ptsum.at(t).Fill(leadlep.Pt()+traillep.Pt(),w);

	  if(jetsfromvtx.size()==1 && firstjet_idx!=-1){
		  onejet.at(t).Fill(Ntp->PFJet_p4(firstjet_idx).Pt(),w);
		  onejet_eta.at(t).Fill(Ntp->PFJet_p4(firstjet_idx).Eta(),w);
	  }

	  double sumht(0);
	  for(unsigned i=0;i<jetsfromvtx.size();i++){
		  sumht += Ntp->PFJet_p4(jetsfromvtx.at(i)).Et();
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
		  if(Ntp->PFJet_p4(i).Pt()<20) continue;
		  if(fabs(Ntp->PFJet_p4(i).Eta())>jet_eta) continue;
		  if(!Ntp->isJetID(i)) continue;
		  if(muidx1!=999 && Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(muidx1))<0.3) continue;
		  if(muidx2!=999 && Ntp->PFJet_p4(i).DeltaR(Ntp->Muon_p4(muidx2))<0.3) continue;
		  if(eidx1!=999 && Ntp->PFJet_p4(i).DeltaR(Ntp->Electron_p4(eidx1))<0.3) continue;
		  if(eidx2!=999 && Ntp->PFJet_p4(i).DeltaR(Ntp->Electron_p4(eidx2))<0.3) continue;
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

	  NPV.at(t).Fill(Ntp->NVtx(),w);
	  if(!Ntp->isData()){
		  if(Ntp->PUWeight()>0.)NPV3d.at(t).Fill(Ntp->NVtx(),w*Ntp->PUWeight3D()/Ntp->PUWeight());
		  if(Ntp->PUWeight()>0.)NPVfine.at(t).Fill(Ntp->NVtx(),w*Ntp->PUWeightFineBins()/Ntp->PUWeight());
	  }else{
		  NPV3d.at(t).Fill(Ntp->NVtx(),w);
		  NPVfine.at(t).Fill(Ntp->NVtx(),w);
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
			  if(pass.at(MET)){
				  invmass_met.at(t).Fill(m,w);
				  invmass_met_m.at(t).Fill(m,w);
				  if(pass.at(ptBalance)){
					  invmass_ptbalance.at(t).Fill(m,w);
					  invmass_ptbalance_m.at(t).Fill(m,w);
					  if(jetsfromvtx.size()>0 && firstjet_idx!=-1) nm0_onejet.at(t).Fill(value.at(oneJet),w);
					  nm0_met.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),w);
					  nm0_mvamet.at(t).Fill(Ntp->MET_CorrMVA_et(),w);
					  nm0_ptbalance.at(t).Fill(value.at(ptBalance),w);
					  zmass_zoom.at(t).Fill(m,w);
					  ptsum_nm0.at(t).Fill(leadlep.Pt()+traillep.Pt(),w);
					  if(pass.at(ZMassmax)
							  && pass.at(ZMassmin)){
						  invmass_zmass.at(t).Fill(m,w);
						  invmass_zmass_m.at(t).Fill(m,w);
						  if(jetsfromvtx.size()>0 && firstjet_idx!=-1) nmm_onejet.at(t).Fill(Ntp->PFJet_p4(firstjet_idx).Pt(),w);
						  nmm_met.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),w);
						  nmm_mvamet.at(t).Fill(Ntp->MET_CorrMVA_et(),w);
						  nmm_ptbalance.at(t).Fill(value.at(ptBalance),w);
					  }
				  }
			  }
		  }
	  }
	  if(pass.at(triLeptonVeto)) invmass_trilepton_only.at(t).Fill(m,w);
	  if(pass.at(charge)) invmass_charge_only.at(t).Fill(m,w);
	  if(pass.at(oneJet)) invmass_jetveto_only.at(t).Fill(m,w);
	  if(pass.at(MET)) invmass_met_only.at(t).Fill(m,w);
	  if(pass.at(ptBalance)) invmass_ptbal_only.at(t).Fill(m,w);
  }

  if(verbose)std::cout << "Tvariable_Base::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

double Tvariable_Base::ZPtReweight(double zpt){
	double weight = 1;
	/*if(zpt<2.5)	weight/=1.03445;
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
	else if(zpt>=27.5 && zpt<30.) weight/=0.951472;*/
	if(Ntp->GetMCID()==40){
		if(zpt<5.) weight/=1.84321;
		else if(zpt>=5. && zpt<10.) weight/=0.743062;
		else if(zpt>=10. && zpt<15.) weight/=0.899661;
		else if(zpt>=15. && zpt<20.) weight/=0.971033;
		else if(zpt>=20. && zpt<25.) weight/=0.970502;
		else if(zpt>=25. && zpt<30.) weight/=1.00624;
		else if(zpt>=30. && zpt<35.) weight/=1.0576;
		else if(zpt>=35. && zpt<40.) weight/=0.995763;
		else if(zpt>=40. && zpt<45.) weight/=1.07514;
		else if(zpt>=45. && zpt<50.) weight/=0.95596;
		else if(zpt>=50. && zpt<55.) weight/=1.00109;
		else if(zpt>=55. && zpt<60.) weight/=1.00435;
		else if(zpt>=60. && zpt<65.) weight/=1.05995;
		else if(zpt>=65. && zpt<70.) weight/=1.048;
		else if(zpt>=70. && zpt<75.) weight/=1.04944;
		else if(zpt>=75. && zpt<80.) weight/=0.981557;
		else if(zpt>=80. && zpt<85.) weight/=0.845718;
		else if(zpt>=85. && zpt<90.) weight/=0.851668;
		else if(zpt>=90. && zpt<95.) weight/=1.20121;
		else if(zpt>=95. && zpt<100.) weight/=0.898355;
	}else if(Ntp->GetMCID()==30 || Ntp->GetMCID()==33){
		if(zpt<5.) weight/=1.73327;
		else if(zpt>=5. && zpt<10.) weight/=0.709086;
		else if(zpt>=10. && zpt<15.) weight/=0.855728;
		else if(zpt>=15. && zpt<20.) weight/=0.932681;
		else if(zpt>=20. && zpt<25.) weight/=1.06339;
		else if(zpt>=25. && zpt<30.) weight/=1.08427;
		else if(zpt>=30. && zpt<35.) weight/=1.15234;
		else if(zpt>=35. && zpt<40.) weight/=1.17128;
		else if(zpt>=40. && zpt<45.) weight/=1.06381;
		else if(zpt>=45. && zpt<50.) weight/=1.27088;
		else if(zpt>=50. && zpt<55.) weight/=1.09183;
		else if(zpt>=55. && zpt<60.) weight/=1.18817;
		else if(zpt>=60. && zpt<65.) weight/=1.23817;
		else if(zpt>=65. && zpt<70.) weight/=1.24276;
		else if(zpt>=70. && zpt<75.) weight/=1.09174;
		else if(zpt>=75. && zpt<80.) weight/=1.13458;
		else if(zpt>=80. && zpt<85.) weight/=1.51522;
		else if(zpt>=85. && zpt<90.) weight/=1.12188;
		else if(zpt>=90. && zpt<95.) weight/=1.30456;
		else if(zpt>=95. && zpt<100.) weight/=1.18801;
	}
	return weight;
}

double Tvariable_Base::PowhegReweight(double zpt){
	double weight = 1.;
	weight = ZptCorrection->GetBinContent(ZptCorrection->FindFixBin(zpt));
	return weight;
}

// relative uncertainties obtained from table 16 in AN-2012-225
double Tvariable_Base::ZPtRelUnc(double zpt){
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

double Tvariable_Base::ZPtMadgraphRelUnc(double zpt){
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

double Tvariable_Base::calculatePzeta(int muiterator, int eiterator){
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

double Tvariable_Base::calculatePzetaDQM(int muiterator, int eiterator){
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

double Tvariable_Base::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double Tvariable_Base::cosphi3d(TVector3 vec1, TVector3 vec2){
	return (vec1.Dot(vec2))/vec1.Mag()/vec2.Mag();
}

int Tvariable_Base::findBin(TGraphAsymmErrors* graph, double xval){
	int bin = -1;
	for(int i=0;i<graph->GetN()-1;i++){
		if(fabs(graph->GetX()[i]-xval)<fabs(graph->GetX()[i+1]-xval)) bin = i;
		else if(fabs(graph->GetX()[i]-xval)>=fabs(graph->GetX()[i+1]-xval)) bin = i+1;
	}
	return bin;
}

//////////////////////////////
//
// Muon related functions
//

bool Tvariable_Base::isFakeMuon(unsigned int idx){
	if(!Ntp->Muon_isGlobalMuon(idx)) return false;
	if(Ntp->Muon_p4(idx).Pt()<=10) return false;
	if(fabs(Ntp->Muon_p4(idx).Eta())>2.4) return false;
	if(Ntp->Muon_p4(idx).Pt()<=20){
		if(Ntp->Muon_sumPt03(idx)>=8.) return false;
		if(Ntp->Muon_emEt03(idx)>=8.) return false;
		if(Ntp->Muon_hadEt03(idx)>=8.) return false;
	}
	if(Ntp->Muon_p4(idx).Pt()>20){
		if(Ntp->Muon_sumPt03(idx)/Ntp->Muon_p4(idx).Pt()>=0.4) return false;
		if(Ntp->Muon_emEt03(idx)/Ntp->Muon_p4(idx).Pt()>=0.4) return false;
		if(Ntp->Muon_hadEt03(idx)/Ntp->Muon_p4(idx).Pt()>=0.4) return false;
	}
	return true;
}

bool Tvariable_Base::isFakeMuon(unsigned int idx, unsigned int vtx){
	if(vtx<0 || vtx>=Ntp->NVtx()) return false;
	if(!isFakeMuon(idx)) return false;
	if(Ntp->dxy(Ntp->Muon_p4(idx),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.2) return false;
	if(Ntp->dz(Ntp->Muon_p4(idx),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.1) return false;
	return true;
}

//////////////////////////////
//
// Electron related functions
//

bool Tvariable_Base::isWWElectron(unsigned int idx, unsigned int vtx){
	double mvapt = Ntp->Electron_p4(idx).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(idx));
	if(mvapt<10.) return false;
	if(mvaeta>2.5) return false;
	if(!isFakeElectron(idx,vtx)) return false;
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

bool Tvariable_Base::isFakeElectron(unsigned int idx){
	if(Ntp->Electron_p4(idx).Pt()<10) return false;
	if(fabs(Ntp->Electron_supercluster_eta(idx))>2.5) return false;
	if(Ntp->Electron_HasMatchedConversions(idx)) return false;
	if(Ntp->Electron_numberOfMissedHits(idx)>0) return false;
	if(Ntp->Electron_tkSumPt03(idx)/Ntp->Electron_p4(idx).Pt()>0.2) return false;
	if(std::max(Ntp->Electron_ecalRecHitSumEt03(idx)-1.,0.)/Ntp->Electron_p4(idx).Pt()>0.2) return false;
	if((Ntp->Electron_hcalDepth1TowerSumEt03(idx)+Ntp->Electron_hcalDepth2TowerSumEt03(idx))/Ntp->Electron_p4(idx).Pt()>0.2) return false;
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

bool Tvariable_Base::isFakeElectron(unsigned int idx, unsigned int vtx){
	if(vtx<0 || vtx>=Ntp->NVtx()) return false;
	if(!isFakeElectron(idx)) return false;
	if(Ntp->dz(Ntp->Electron_p4(idx),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>0.1) return false;
	if(Ntp->dxy(Ntp->Electron_p4(idx),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>0.02) return false;
	return true;
}

//////////////////////////////
//
// Calculate fakerate
//

double Tvariable_Base::Fakerate(double pt, double eta, TH2D *fakeRateHist){
	double fakerate = 0.;
	if(pt>=35.) pt = 34.99;
	fakerate = fakeRateHist->GetBinContent(fakeRateHist->FindFixBin(pt,eta));
	return fakerate/(1-fakerate);
}

double Tvariable_Base::FakerateError(double pt, double eta, TH2D *fakeRateHist){
	double error = 0.;
	if(pt>=35.) pt = 34.99;
	error = fakeRateHist->GetBinError(fakeRateHist->FindFixBin(pt,eta));
	return error/(1-error);
}

//////////////////////////////
//
// Finish function
//

void Tvariable_Base::Finish(){
	/*if(mode==RECONSTRUCT){
		double ndymc(0),ndyemb(0);
		for(unsigned i=0;i<HConfig.GetNHisto();i++){
			if(HConfig.GetID(i)==DataMCType::DY_tautau){
				printf("DY_tautau after ptbalance cut: %.3f\n",Npassed.at(i).GetBinContent(ptBalance+2));
				ndymc = Npassed.at(i).GetBinContent(ptBalance+2);
			}
			if(HConfig.GetID(i)==DataMCType::DY_emu_embedded){
				printf("DY_emu_embedded after ptbalance cut: %.3f\n",Npassed.at(i).GetBinContent(ptBalance+2));
				ndyemb = Npassed.at(i).GetBinContent(ptBalance+2);
			}
		}
		if(ndymc>0 && ndyemb>0){
			Selection::ScaleAllHistOfType(HConfig.GetType(DataMCType::DY_emu_embedded),ndymc/ndyemb*10393634./10152445.);
			std::cout << "Scaled embedded histograms" << std::endl;
		}else{
			Selection::ScaleAllHistOfType(DataMCType::DY_emu_embedded,0);
			std::cout << "Embedded histograms not scaled properly" << std::endl;
		}
	}*/
	Selection::Finish();
	/*double sumdata(0), sumbkg(0), sumsignal(0);
	double stddata(87),stdbkg(80.629),stdsignal(74.177);
	double stdqcd(13.877),stdww(17.123),stdwz2l2q(0.000),stdwz3l1nu(0.249),stdzz4l(0.026),stdzz2l2q(0.000),stdzz2l2nu(0.009),stdtt(1.095),stdtw(0.592),stdtbarw(0.000),stddyll(7.213),stddytt(39.958);
	double std(0);
	double bkgunc(0),individualunc(0),statunc(0);
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
		if(HConfig.GetID(i)!=DataMCType::Data && HConfig.GetID(i)!=DataMCType::DY_emu){
			individualunc += pow(abs(std-Npassed.at(i).GetBinContent(NCuts+1)),2);
			statunc += pow(Npassed.at(i).GetBinError(NCuts+1),2);
		}
		//std::cout << "Number of events for sample " << HConfig.GetName(i) << " ";
		printf("Standard: %8.3f, Now: %8.3f, Resulting uncertainty: %8.3f. Sample: %s\n",std,Npassed.at(i).GetBinContent(NCuts+1),abs(std-Npassed.at(i).GetBinContent(NCuts+1)),HConfig.GetName(i).Data());
	}

	std::cout << "################## Result ##################" << std::endl;
	printf("Total events in data: %.0f, background: %.3f and signal %.3f\n",sumdata,sumbkg,sumsignal);
	printf("Total background normalization uncertainty: %.3f%%\n",sqrt(bkgunc)/sumbkg*100);
	printf("Total statistical background uncertainty: %.3f (abs), %.3f%%\n",sqrt(statunc),sqrt(statunc)/sumbkg*100);
	printf("Total background uncertainty not from cross section: %.3f%%\n",sqrt(individualunc)/sumbkg*100);

	printf("Difference in background events: %.3f%% and in signal events: %.3f%%\n",fabs(sumbkg/stdbkg-1)*100,fabs(sumsignal/stdsignal-1)*100);
	if(!systValid) std::cout << "!!! WARNING: MORE THAN ONE SYSTEMATIC SET TO TRUE !!!" << std::endl;
*/
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
				for(int j=1;j<nPDFmembers;++j){
					acc2 += pow(pdf_w1.at(i).GetBinContent(j+1)/pdf_w0.at(i).GetBinContent(j+1)-acc,2);
					yield2 += pow(pdf_w1.at(i).GetBinContent(j+1)-yield,2);
				}
			}else{
				for(int j=1;j<nPDFmembers;j+=2){
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
			printf("Yield = %8.3e +- %8.3e (PDFs), i.e. %8.2f%% relative uncertainty\n",yield,eyield,eyield/yield*100);
			printf("Acceptance = %8.3e +- %8.3e (PDFs), i.e. %8.2f%% relative uncertainty\n",acc,eacc,eacc/acc*100);
		}
	}

}
