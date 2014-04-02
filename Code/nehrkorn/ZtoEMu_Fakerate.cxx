#include "ZtoEMu_Fakerate.h"
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

ZtoEMu_Fakerate::ZtoEMu_Fakerate(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_pt(3)
  ,e_pt(3)
  ,mu_eta(2.4)
  ,e_eta(2.5)
  ,jet_pt(25)
{
    //verbose=true;
	//outfile = new TFile("ownfakerates.root","RECREATE");
	doHiggsObjects = false;
}

ZtoEMu_Fakerate::~ZtoEMu_Fakerate(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZtoEMu_Fakerate::~ZtoEMu_Fakerate Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZtoEMu_Fakerate::~ZtoEMu_Fakerate()" << std::endl;
}

void  ZtoEMu_Fakerate::Configure(){

  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==PrimeVtx)           cut.at(PrimeVtx)=1;
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
    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  tightmu=HConfig.GetTH2D(Name+"_tightmu","tightmu",100,0.,100.,24,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  tightmu_rebin=HConfig.GetTH2D(Name+"_tightmu_rebin","tightmu_rebin",9,10.,55.,4,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  tighte=HConfig.GetTH2D(Name+"_tighte","tighte",100,0.,100.,25,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  tighte_rebin=HConfig.GetTH2D(Name+"_tighte_rebin","tighte_rebin",9,10.,55.,5,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  fakemu=HConfig.GetTH2D(Name+"_fakemu","fakemu",100,0.,100.,24,0.,2.44,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  fakemu_rebin=HConfig.GetTH2D(Name+"_fakemu_rebin","fakemu_rebin",9,10.,55.,4,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  fakee=HConfig.GetTH2D(Name+"_fakee","fakee",100,0.,100.,25,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  fakee_rebin=HConfig.GetTH2D(Name+"_fakee_rebin","fakee_rebin",9,10.,55.,5,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  mueff=HConfig.GetTH2D(Name+"_mueff","mueff",9,10.,55.,4,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eeff=HConfig.GetTH2D(Name+"_eeff","eeff",9,10.,55.,5,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");

  muleg_numerator=HConfig.GetTH2D(Name+"_muleg_numerator","muleg_numerator",100,0.,100.,100,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  muleg_denominator=HConfig.GetTH2D(Name+"_muleg_denominator","muleg_denominator",100,0.,100.,100,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eleg_numerator=HConfig.GetTH2D(Name+"_eleg_numerator","eleg_numerator",100,0.,100.,100,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  eleg_denominator=HConfig.GetTH2D(Name+"_eleg_denominator","eleg_denominator",100,0.,100.,100,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");

  muleg_numerator_rebin=HConfig.GetTH2D(Name+"_muleg_numerator_rebin","muleg_numerator_rebin",10,5.,55.,4,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  muleg_denominator_rebin=HConfig.GetTH2D(Name+"_muleg_denominator_rebin","muleg_denominator_rebin",10,5.,55.,4,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eleg_numerator_rebin=HConfig.GetTH2D(Name+"_eleg_numerator_rebin","eleg_numerator_rebin",10,5.,55.,5,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  eleg_denominator_rebin=HConfig.GetTH2D(Name+"_eleg_denominator_rebin","eleg_denominator_rebin",10,5.,55.,5,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");

  muleg_eff=HConfig.GetTH2D(Name+"_muleg_eff","muleg_eff",10,5.,55.,4,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eleg_eff=HConfig.GetTH2D(Name+"_eleg_eff","eleg_eff",10,5.,55.,5,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");
  muleg_scale=HConfig.GetTH2D(Name+"_muleg_scale","muleg_scale",10,5.,55.,4,0.,2.4,"p_{T}^{#mu} / GeV","#eta_{#mu}");
  eleg_scale=HConfig.GetTH2D(Name+"_eleg_scale","eleg_scale",10,5.,55.,5,0.,2.5,"p_{T}^{e} / GeV","#eta_{e}");

  muleg_denom_09=HConfig.GetTH1D(Name+"_muleg_denom_09","muleg_denom_09",10,5.,55.,"p_{T}^{#mu} / GeV");
  muleg_denom_12=HConfig.GetTH1D(Name+"_muleg_denom_12","muleg_denom_12",10,5.,55.,"p_{T}^{#mu} / GeV");
  muleg_denom_21=HConfig.GetTH1D(Name+"_muleg_denom_21","muleg_denom_21",10,5.,55.,"p_{T}^{#mu} / GeV");
  muleg_denom_24=HConfig.GetTH1D(Name+"_muleg_denom_24","muleg_denom_24",10,5.,55.,"p_{T}^{#mu} / GeV");
  muleg_denom_all=HConfig.GetTH1D(Name+"_muleg_denom_all","muleg_denom_all",10,5.,55.,"p_{T}^{#mu} / GeV");
  muleg_num_09=HConfig.GetTH1D(Name+"_muleg_num_09","muleg_num_09",10,5.,55.,"p_{T}^{#mu} / GeV");
  muleg_num_12=HConfig.GetTH1D(Name+"_muleg_num_12","muleg_num_12",10,5.,55.,"p_{T}^{#mu} / GeV");
  muleg_num_21=HConfig.GetTH1D(Name+"_muleg_num_21","muleg_num_21",10,5.,55.,"p_{T}^{#mu} / GeV");
  muleg_num_24=HConfig.GetTH1D(Name+"_muleg_num_24","muleg_num_24",10,5.,55.,"p_{T}^{#mu} / GeV");
  muleg_num_all=HConfig.GetTH1D(Name+"_muleg_num_all","muleg_num_all",10,5.,55.,"p_{T}^{#mu} / GeV");

  eleg_denom_08=HConfig.GetTH1D(Name+"_eleg_denom_08","eleg_denom_10",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_denom_14=HConfig.GetTH1D(Name+"_eleg_denom_14","eleg_denom_14",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_denom_15=HConfig.GetTH1D(Name+"_eleg_denom_15","eleg_denom_15",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_denom_20=HConfig.GetTH1D(Name+"_eleg_denom_20","eleg_denom_20",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_denom_25=HConfig.GetTH1D(Name+"_eleg_denom_25","eleg_denom_25",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_denom_all=HConfig.GetTH1D(Name+"_eleg_denom_all","eleg_denom_all",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_num_08=HConfig.GetTH1D(Name+"_eleg_num_08","eleg_num_08",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_num_14=HConfig.GetTH1D(Name+"_eleg_num_14","eleg_num_14",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_num_15=HConfig.GetTH1D(Name+"_eleg_num_15","eleg_num_15",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_num_20=HConfig.GetTH1D(Name+"_eleg_num_20","eleg_num_20",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_num_25=HConfig.GetTH1D(Name+"_eleg_num_25","eleg_num_25",10,5.,55.,"p_{T}^{#mu} / GeV");
  eleg_num_all=HConfig.GetTH1D(Name+"_eleg_num_all","eleg_num_all",10,5.,55.,"p_{T}^{#mu} / GeV");

  tagmupt=HConfig.GetTH1D(Name+"_tagmupt","tagmupt",40,0.,100.,"p_{T}^{#mu} / GeV");
  tagmueta=HConfig.GetTH1D(Name+"_tagmueta","tagmueta",25,-2.5,2.5,"#eta_{#mu}");
  tagept=HConfig.GetTH1D(Name+"_tagept","tagept",40,0.,100.,"p_{T}^{e} / GeV");
  tageeta=HConfig.GetTH1D(Name+"_tageeta","tageeta",25,-2.5,2.5,"#eta_{e}");
  probemupt=HConfig.GetTH1D(Name+"_probemupt","probemupt",40,0.,100.,"p_{T}^{#mu} / GeV");
  probemueta=HConfig.GetTH1D(Name+"_probemueta","probemueta",25,-2.5,2.5,"#eta_{#mu}");
  probeept=HConfig.GetTH1D(Name+"_probeept","probeept",40,0.,100.,"p_{T}^{e} / GeV");
  probeeeta=HConfig.GetTH1D(Name+"_probeeeta","probeeeta",25,-2.5,2.5,"#eta_{e}");
  drmumu=HConfig.GetTH1D(Name+"_drmumu","drmumu",60,0.,3.,"#DeltaR(#mu,#mu)");
  dree=HConfig.GetTH1D(Name+"_dree","dree",60,0.,3.,"#DeltaR(e,e)");
  ptbalmumu=HConfig.GetTH1D(Name+"_ptbalmumu","ptbalmumu",40,0.,200.,"p_{T}^{#mu#mu} / GeV");
  ptbalee=HConfig.GetTH1D(Name+"_ptbalee","ptbalee",40,0.,200.,"p_{T}^{ee} / GeV");
  mttagmu=HConfig.GetTH1D(Name+"_mttagmu","mttagmu",40,0.,160.,"m_{t}^{#mu} / GeV");
  mttage=HConfig.GetTH1D(Name+"_mttage","mttage",40,0.,160.,"m_{t}^{e} / GeV");
  mtprobemu=HConfig.GetTH1D(Name+"_mtprobemu","mtprobemu",40,0.,160.,"m_{t}^{#mu} / GeV");
  mtprobee=HConfig.GetTH1D(Name+"_mtprobee","mtprobee",40,0.,160.,"m_{t}^{e} / GeV");
  mmumu=HConfig.GetTH1D(Name+"_mmumu","mmumu",30,60.,120.,"m_{#mu,#mu} / GeV");
  mmumupass=HConfig.GetTH1D(Name+"_mmumupass","mmumupass",30,60.,120.,"m_{#mu,#mu} / GeV");
  mmumufail=HConfig.GetTH1D(Name+"_mmumufail","mmumufail",30,60.,120.,"m_{#mu,#mu} / GeV");
  mee=HConfig.GetTH1D(Name+"_mee","mee",30,60.,120.,"m_{e,e} / GeV");
  meepass=HConfig.GetTH1D(Name+"_meepass","meepass",30,60.,120.,"m_{e,e} / GeV");
  meefail=HConfig.GetTH1D(Name+"_meefail","meefail",30,60.,120.,"m_{e,e} / GeV");
  nprobes=HConfig.GetTH1D(Name+"_nprobes","nprobes",10,0.,10.,"number of probes");
  ntags=HConfig.GetTH1D(Name+"_ntags","ntags",10,0.,10.,"number of tags");
  met=HConfig.GetTH1D(Name+"_met","met",20,0.,100.,"MET / GeV");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu_Fakerate::Store_ExtraDist(){
	Extradist2d.push_back(&tightmu);
	Extradist2d.push_back(&tightmu_rebin);
	Extradist2d.push_back(&tighte);
	Extradist2d.push_back(&tighte_rebin);
	Extradist2d.push_back(&fakemu);
	Extradist2d.push_back(&fakemu_rebin);
	Extradist2d.push_back(&fakee);
	Extradist2d.push_back(&fakee_rebin);
	Extradist2d.push_back(&mueff);
	Extradist2d.push_back(&eeff);

	Extradist2d.push_back(&muleg_numerator);
	Extradist2d.push_back(&muleg_denominator);
	Extradist2d.push_back(&eleg_numerator);
	Extradist2d.push_back(&eleg_denominator);

	Extradist2d.push_back(&muleg_numerator_rebin);
	Extradist2d.push_back(&muleg_denominator_rebin);
	Extradist2d.push_back(&eleg_numerator_rebin);
	Extradist2d.push_back(&eleg_denominator_rebin);

	Extradist2d.push_back(&muleg_eff);
	Extradist2d.push_back(&eleg_eff);
	Extradist2d.push_back(&muleg_scale);
	Extradist2d.push_back(&eleg_scale);

	Extradist1d.push_back(&muleg_denom_09);
	Extradist1d.push_back(&muleg_denom_12);
	Extradist1d.push_back(&muleg_denom_21);
	Extradist1d.push_back(&muleg_denom_24);
	Extradist1d.push_back(&muleg_denom_all);
	Extradist1d.push_back(&muleg_num_09);
	Extradist1d.push_back(&muleg_num_12);
	Extradist1d.push_back(&muleg_num_21);
	Extradist1d.push_back(&muleg_num_24);
	Extradist1d.push_back(&muleg_num_all);

	Extradist1d.push_back(&eleg_denom_08);
	Extradist1d.push_back(&eleg_denom_14);
	Extradist1d.push_back(&eleg_denom_15);
	Extradist1d.push_back(&eleg_denom_20);
	Extradist1d.push_back(&eleg_denom_25);
	Extradist1d.push_back(&eleg_denom_all);
	Extradist1d.push_back(&eleg_num_08);
	Extradist1d.push_back(&eleg_num_14);
	Extradist1d.push_back(&eleg_num_15);
	Extradist1d.push_back(&eleg_num_20);
	Extradist1d.push_back(&eleg_num_25);
	Extradist1d.push_back(&eleg_num_all);

	Extradist1d.push_back(&tagmupt);
	Extradist1d.push_back(&tagmueta);
	Extradist1d.push_back(&tagept);
	Extradist1d.push_back(&tageeta);
	Extradist1d.push_back(&probemupt);
	Extradist1d.push_back(&probemueta);
	Extradist1d.push_back(&probeept);
	Extradist1d.push_back(&probeeeta);
	Extradist1d.push_back(&drmumu);
	Extradist1d.push_back(&dree);
	Extradist1d.push_back(&ptbalmumu);
	Extradist1d.push_back(&ptbalee);
	Extradist1d.push_back(&mttagmu);
	Extradist1d.push_back(&mttage);
	Extradist1d.push_back(&mtprobemu);
	Extradist1d.push_back(&mtprobee);
	Extradist1d.push_back(&mmumu);
	Extradist1d.push_back(&mmumupass);
	Extradist1d.push_back(&mmumufail);
	Extradist1d.push_back(&mee);
	Extradist1d.push_back(&meepass);
	Extradist1d.push_back(&meefail);
	Extradist1d.push_back(&nprobes);
	Extradist1d.push_back(&ntags);
	Extradist1d.push_back(&met);
}

void  ZtoEMu_Fakerate::doEvent(){
  if(verbose)std::cout << "ZtoEMu_Fakerate::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}
  
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
  
  //////////////////////////////////////////////////////////
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1),ww(1);
  if(!Ntp->isData() && Ntp->GetMCID()!=34){
    w*=Ntp->PUWeight();
    ww*=Ntp->PUWeight();
    if(verbose)std::cout << "void  ZtoEMu_Fakerate::doEvent() k" << w << " " << wobs << std::endl;
  }
  else{w=1;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;

  ///////////////////////////////////////////////////////////
  //
  // Fakerate calculation
  //

  // muon fakerate

  double leadingJetPt = -1;
  int nfake = 0;

  if(Ntp->isData() && Ntp->MET_CorrT0pcT1_et()<=20.){
	  for(unsigned i=0;i<Ntp->NMuons();i++){
		  if(nfake>0) break;
		  if(Ntp->Muon_p4(i).Pt()>mu_pt
				  && fabs(Ntp->Muon_p4(i).Eta())<mu_eta
				  && vertex>=0
				  && (Ntp->TriggerAccept("HLT_Mu8_v")
						  || Ntp->TriggerAccept("HLT_Mu12_v")
						  || Ntp->TriggerAccept("HLT_Mu17_v"))
				  && passZVeto(vertex,"muon")
				  && sqrt(2*Ntp->Muon_p4(i).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(i).Px(),Ntp->Muon_p4(i).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey())))<=20.
				  ){
			  for(unsigned j=0;j<Ntp->NPFJets();j++){
				  if(Ntp->PFJet_p4(j).Pt()>leadingJetPt
						  && Ntp->PFJet_p4(j).DeltaR(Ntp->Muon_p4(i))>1.
						  ){
					  leadingJetPt = Ntp->PFJet_p4(j).Pt();
				  }
			  }
			  if(!leadingJetPt>jet_pt) continue;
			  if(!matchTrigger(i,0.2,"HLT_Mu8_v","muon")
					  && !matchTrigger(i,0.2,"HLT_Mu12_v","muon")
					  && !matchTrigger(i,0.2,"HLT_Mu17_v","muon")){
				  continue;
			  }
			  if(isFakeMuon(i,vertex)){
				  fakemu.at(t).Fill(Ntp->Muon_p4(i).Pt(),fabs(Ntp->Muon_p4(i).Eta()));
				  fakemu_rebin.at(t).Fill(Ntp->Muon_p4(i).Pt(),fabs(Ntp->Muon_p4(i).Eta()));
				  nfake++;
				  if(doHiggsObjects){
					  if(isHiggsMuon(i,vertex)
							  && ((fabs(Ntp->Muon_p4(i).Eta())<1.479 && Muon_RelIso(i)<0.15) || (fabs(Ntp->Muon_p4(i).Eta())>=1.479 && Muon_RelIso(i)<0.10))
							  ){
						  tightmu.at(t).Fill(Ntp->Muon_p4(i).Pt(),fabs(Ntp->Muon_p4(i).Eta()));
						  tightmu_rebin.at(t).Fill(Ntp->Muon_p4(i).Pt(),fabs(Ntp->Muon_p4(i).Eta()));
					  }
				  }else{
					  if(isTightMuon(i,vertex)
							  && Muon_RelIso(i)<0.12){
						  tightmu.at(t).Fill(Ntp->Muon_p4(i).Pt(),fabs(Ntp->Muon_p4(i).Eta()));
						  tightmu_rebin.at(t).Fill(Ntp->Muon_p4(i).Pt(),fabs(Ntp->Muon_p4(i).Eta()));
					  }
				  }
			  }
		  }
	  }

	  leadingJetPt = -1;
	  nfake = 0;

	  // electron fakerate
	  for(unsigned i=0;i<Ntp->NElectrons();i++){
		  if(nfake>0) break;
		  if(Ntp->Electron_p4(i).Et()>e_pt
				  && fabs(Ntp->Electron_supercluster_eta(i))<e_eta
				  && vertex>=0
				  && (Ntp->TriggerAccept("HLT_Ele8_CaloIdT_TrkIdVL_v")
						  || Ntp->TriggerAccept("HLT_Ele8_CaloIdL_CaloIsoVL_v")
						  || Ntp->TriggerAccept("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v"))
				  //&& passZVeto(vertex,"electron")
				  ){
			  for(unsigned j=0;j<Ntp->NPFJets();j++){
				  if(Ntp->PFJet_p4(j).Pt()>leadingJetPt
						  && Ntp->PFJet_p4(j).DeltaR(Ntp->Electron_p4(i))>1.
						  ){
					  leadingJetPt = Ntp->PFJet_p4(j).Pt();
				  }
			  }
			  if(!leadingJetPt>jet_pt) continue;
			  if(!matchTrigger(i,02.,"HLT_Ele8_CaloIdT_TrkIdVL_v","electron")
					  && !matchTrigger(i,0.2,"HLT_Ele8_CaloIdL_CaloIsoVL_v","electron")
					  && !matchTrigger(i,0.2,"HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","electron")
					  ){
				  continue;
			  }
			  if(isFakeElectron(i,vertex)){
				  fakee.at(t).Fill(Ntp->Electron_p4(i).Pt(),fabs(Ntp->Electron_supercluster_eta(i)));
				  fakee_rebin.at(t).Fill(Ntp->Electron_p4(i).Pt(),fabs(Ntp->Electron_supercluster_eta(i)));
				  nfake++;
				  if(doHiggsObjects){
					  if(isHiggsElectron(i,vertex)
							  && ((fabs(Ntp->Electron_supercluster_eta(i))<1.479 && Electron_RelIso(i)<0.15) || (fabs(Ntp->Electron_supercluster_eta(i))>=1.479 && Electron_RelIso(i)<0.10))
							  ){
						  tighte.at(t).Fill(Ntp->Electron_p4(i).Pt(),fabs(Ntp->Electron_supercluster_eta(i)));
						  tighte_rebin.at(t).Fill(Ntp->Electron_p4(i).Pt(),fabs(Ntp->Electron_supercluster_eta(i)));
					  }
				  }else{
					  if(isMVATrigElectron(i)
							  && Electron_RelIso(i)<0.15
					  ){
				  tighte.at(t).Fill(Ntp->Electron_p4(i).Pt(),fabs(Ntp->Electron_supercluster_eta(i)));
				  tighte_rebin.at(t).Fill(Ntp->Electron_p4(i).Pt(),fabs(Ntp->Electron_supercluster_eta(i)));
					  }
				  }
			  }
		  }
	  }
  }

  ///////////////////////////////////////////////////////////
  //
  // Trigger efficiency
  //

  /*if(Ntp->TriggerAccept("HLT_Mu17_Mu8_v")){
	  for(unsigned j=0;j<Ntp->NHLTTrigger_objs();j++){
		if(Ntp->HLTTrigger_objs_trigger(j).find("HLT_Mu17_Mu8_v") != string::npos){
			std::cout << "Trigger objects firing HLT_Mu17_Mu8_v:" << std::endl;
			for(unsigned k=0;k<Ntp->NHLTTrigger_objs(j);k++){
				std::cout << "id " << Ntp->HLTTrigger_objs_Id(j,k) << std::endl;
				std::cout << "pt " << Ntp->HLTTrigger_objs_Pt(j,k) << std::endl;
			}
		}
	  }
  }
  if(Ntp->TriggerAccept("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")){
	  for(unsigned j=0;j<Ntp->NHLTTrigger_objs();j++){
		if(Ntp->HLTTrigger_objs_trigger(j).find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") != string::npos){
			std::cout << "Trigger objects firing HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v:" << std::endl;
			for(unsigned k=0;k<Ntp->NHLTTrigger_objs(j);k++){
				std::cout << "id " << Ntp->HLTTrigger_objs_Id(j,k) << std::endl;
				std::cout << "pt " << Ntp->HLTTrigger_objs_Pt(j,k) << std::endl;
			}
		}
	  }
  }*/
  met.at(t).Fill(Ntp->MET_CorrT0pcT1_et(),w);

  if(Ntp->TriggerAccept("HLT_IsoMu24_eta2p1_v")
		  //&& Ntp->MET_CorrT0pcT1_et()<40.
		  ){
	  doubleTrigger("HLT_IsoMu24_eta2p1_v","HLT_Mu17_Mu8_v",vertex,"muon",t,w);
  }
  if(Ntp->TriggerAccept("HLT_Ele27_WP80_v")
		  //&& Ntp->MET_CorrT0pcT1_et()<40.
		  ){
	  doubleTrigger("HLT_Ele27_WP80_v","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v",vertex,"electron",t,w);
  }

  if(verbose)std::cout << "ZtoEMu_Fakerate::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

bool ZtoEMu_Fakerate::isGoodVtx(unsigned int i){
	if(fabs(Ntp->Vtx(i).z())>=24) return false;
	if(Ntp->Vtx(i).Perp()>=2) return false;
	if(Ntp->Vtx_ndof(i)<=4) return false;
	if(Ntp->Vtx_isFake(i)!=0) return false;
	return true;
}

bool ZtoEMu_Fakerate::jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx){
	for(unsigned i=0;i<vtx_track_idx.size();i++){
		if(vtx_track_idx.at(i)==leadingtrack_idx)return true;
	}
	return false;
}

double ZtoEMu_Fakerate::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double ZtoEMu_Fakerate::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

double ZtoEMu_Fakerate::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}

double ZtoEMu_Fakerate::vertexSignificance(TVector3 vec, unsigned int vertex){
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

bool ZtoEMu_Fakerate::passZVeto(unsigned int vertex, std::string object){
	if(object=="muon"){
		for(unsigned i=0;i<Ntp->NMuons();i++){
			if(Ntp->Muon_p4(i).Pt()<20) continue;
			if(fabs(Ntp->Muon_p4(i).Eta())>2.4) continue;
			if(!isFakeMuon(i,vertex)) continue;
			for(unsigned j=i+1;j<Ntp->NMuons();j++){
				if(Ntp->Muon_Charge(i)==Ntp->Muon_Charge(j)) continue;
				if(Ntp->Muon_p4(j).Pt()<20) continue;
				if(fabs(Ntp->Muon_p4(j).Eta())>2.4) continue;
				if(!isFakeMuon(j,vertex)) continue;
				return false;
			}
		}
		return true;
	}else if(object=="electron"){
		for(unsigned i=0;i<Ntp->NElectrons();i++){
			if(Ntp->Electron_p4(i).Et()<20) continue;
			if(fabs(Ntp->Electron_supercluster_eta(i))>2.5) continue;
			if(!isFakeElectron(i,vertex)) continue;
			for(unsigned j=i+1;j<Ntp->NElectrons();j++){
				if(Ntp->Electron_Charge(i)==Ntp->Electron_Charge(j)) continue;
				if(Ntp->Electron_p4(j).Et()<20) continue;
				if(fabs(Ntp->Electron_supercluster_eta(j))>2.5) continue;
				if(!isFakeElectron(j,vertex)) continue;
				return false;
			}
		}
		return true;
	}else{
		std::cout << "WARNING: no known object given" << std::endl;
		return false;
	}
}

bool ZtoEMu_Fakerate::matchTrigger(unsigned int i, double dr, std::string trigger, std::string object){
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
	for(unsigned j=0;j<Ntp->NHLTTrigger_objs();j++){
		if(Ntp->HLTTrigger_objs_trigger(j).find(trigger) != string::npos){
			for(unsigned k=0;k<Ntp->NHLTTrigger_objs(j);k++){
				if(Ntp->HLTTrigger_objs_Id(j,k)==id){
					triggerObj.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(j,k),
							Ntp->HLTTrigger_objs_Eta(j,k),
							Ntp->HLTTrigger_objs_Phi(j,k),
							Ntp->HLTTrigger_objs_E(j,k));
				}
				if(triggerObj.Pt()>0.
						&& particle.Pt()>0.
						&& particle.DeltaR(triggerObj)<dr) return true;
			}
		}
	}
	return false;
}

bool ZtoEMu_Fakerate::matchTruth(unsigned int i, double dr, std::string object){
	bool match = false;
	if(object=="muon"){
		for(unsigned j=0;j<Ntp->NMCParticles();j++){
			if(fabs(Ntp->MCParticle_pdgid(j))==13
					&& fabs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(j)))==23
					){
				if(Ntp->Muon_p4(i).DeltaR(Ntp->MCParticle_p4(j))<dr) match = true;
			}
		}
	}
	if(object=="electron"){
		for(unsigned j=0;j<Ntp->NMCParticles();j++){
			if(fabs(Ntp->MCParticle_pdgid(j))==11
					&& fabs(Ntp->MCParticle_pdgid(Ntp->MCParticle_midx(j)))==23
					){
				if(Ntp->Electron_p4(i).DeltaR(Ntp->MCParticle_p4(j))<dr) match = true;
			}
		}
	}
	return match;
}

bool ZtoEMu_Fakerate::matchFirstLeg(unsigned int i, double dr, std::string trigger, std::string object){
	TLorentzVector particle(0.,0.,0.,0.);
	TLorentzVector triggerObj(0.,0.,0.,0.);
	int id = 0;
	if(object=="muon"){
		particle = Ntp->Muon_p4(i);
		id = 83;
	}
	if(object=="electron"){
		particle = Ntp->Electron_p4(i);
		id = 82;
	}
	for(unsigned j=0;j<Ntp->NHLTTrigger_objs();j++){
		if(Ntp->HLTTrigger_objs_trigger(j).find(trigger) != string::npos){
			if(Ntp->NHLTTrigger_objs(j)==2){
				if(Ntp->HLTTrigger_objs_Id(j,0)==id && Ntp->HLTTrigger_objs_Id(j,1)==id){
					triggerObj.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(j,0),
							Ntp->HLTTrigger_objs_Eta(j,0),
							Ntp->HLTTrigger_objs_Phi(j,0),
							Ntp->HLTTrigger_objs_E(j,0));
				}
			}
			if(particle.Pt()>0.
					&& triggerObj.Pt()>0.
					&& particle.DeltaR(triggerObj)<dr
					){
				return true;
			}
		}
	}
	return false;
}

bool ZtoEMu_Fakerate::matchSecondLeg(unsigned int i, double dr, std::string trigger, std::string object){
	TLorentzVector particle(0.,0.,0.,0.);
	TLorentzVector triggerObj(0.,0.,0.,0.);
	int id = 0;
	bool flag = false;
	if(object=="muon"){
		particle = Ntp->Muon_p4(i);
		id = 83;
	}
	if(object=="electron"){
		particle = Ntp->Electron_p4(i);
		id = 82;
	}
	for(unsigned j=0;j<Ntp->NHLTTrigger_objs();j++){
		if(Ntp->HLTTrigger_objs_trigger(j).find(trigger) != string::npos){
			if(Ntp->NHLTTrigger_objs(j)==2){
				if(Ntp->HLTTrigger_objs_Id(j,0)==id && Ntp->HLTTrigger_objs_Id(j,1)==id){
					triggerObj.SetPtEtaPhiE(Ntp->HLTTrigger_objs_Pt(j,1),
							Ntp->HLTTrigger_objs_Eta(j,1),
							Ntp->HLTTrigger_objs_Phi(j,1),
							Ntp->HLTTrigger_objs_E(j,1));
					flag = true;
				}
			}
		}
		if(flag) break;
	}
	if(particle.Pt()>0.
			&& triggerObj.Pt()>0.
			&& particle.DeltaR(triggerObj)<dr
			){
		return true;
	}
	return false;
}

void ZtoEMu_Fakerate::doubleTrigger(std::string tagtrigger, std::string probetrigger, int vertex, std::string object, int t, double w){
	unsigned int ntag = 0;
	unsigned int nprobe = 0;
	std::vector<unsigned int> tags;
	std::vector<unsigned int> probes;
	bool pass = true;
	if(vertex>=0){
		if(object=="muon"){
			for(unsigned i=0;i<Ntp->NMuons();i++){
				if(Ntp->Muon_p4(i).Pt()<5) continue;
				if(fabs(Ntp->Muon_p4(i).Eta())>2.4) continue;
				if(!Ntp->isData() && Ntp->Muon_p4(i).Pt()<20 && !matchTruth(i,0.2,"muon")) continue;
				if(!isTightMuon(i,vertex)) continue;
				if(Muon_RelIso(i)>=0.12) continue;
				if(matchTrigger(i,0.2,tagtrigger,"muon")
						&& Ntp->Muon_p4(i).Pt()>26.
						&& matchFirstLeg(i,0.2,probetrigger,"muon")
						){
					ntag++;
					tags.push_back(i);
				}else{
					nprobe++;
					probes.push_back(i);
				}
			}
			ntags.at(t).Fill(ntag,w);
			nprobes.at(t).Fill(nprobe,w);
			for(unsigned i=0;i<tags.size();i++){
				tagmupt.at(t).Fill(Ntp->Muon_p4(tags.at(i)).Pt(),w);
				tagmueta.at(t).Fill(Ntp->Muon_p4(tags.at(i)).Eta(),w);
				mttagmu.at(t).Fill(sqrt(2*Ntp->Muon_p4(tags.at(i)).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(tags.at(i)).Px(),Ntp->Muon_p4(tags.at(i)).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
				for(unsigned j=0;j<probes.size();j++){
					if(Ntp->Muon_Charge(tags.at(i))!=Ntp->Muon_Charge(probes.at(j))
						&& (Ntp->Muon_p4(tags.at(i))+Ntp->Muon_p4(probes.at(j))).M()>=60.
						&& (Ntp->Muon_p4(tags.at(i))+Ntp->Muon_p4(probes.at(j))).M()<120.
						){
							probemupt.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
							probemueta.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Eta(),w);
							mtprobemu.at(t).Fill(sqrt(2*Ntp->Muon_p4(probes.at(j)).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Muon_p4(probes.at(j)).Px(),Ntp->Muon_p4(probes.at(j)).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
							drmumu.at(t).Fill(Ntp->Muon_p4(tags.at(i)).DeltaR(Ntp->Muon_p4(probes.at(j))),w);
							ptbalmumu.at(t).Fill((Ntp->Muon_p4(tags.at(i))+Ntp->Muon_p4(probes.at(j))).Pt(),w);
							mmumu.at(t).Fill((Ntp->Muon_p4(tags.at(i))+Ntp->Muon_p4(probes.at(j))).M(),w);
							muleg_denominator.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),fabs(Ntp->Muon_p4(probes.at(j)).Eta()),w);
							muleg_denominator_rebin.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),fabs(Ntp->Muon_p4(probes.at(j)).Eta()),w);
							if(fabs(Ntp->Muon_p4(probes.at(j)).Eta())<0.9) muleg_denom_09.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
							else if(fabs(Ntp->Muon_p4(probes.at(j)).Eta())<1.2) muleg_denom_12.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
							else if(fabs(Ntp->Muon_p4(probes.at(j)).Eta())<2.1) muleg_denom_21.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
							else if(fabs(Ntp->Muon_p4(probes.at(j)).Eta())<2.4) muleg_denom_24.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
							muleg_denom_all.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
							if(!matchSecondLeg(probes.at(j),0.2,probetrigger,"muon")) mmumufail.at(t).Fill((Ntp->Muon_p4(tags.at(i))+Ntp->Muon_p4(probes.at(j))).M(),w);
							if(matchSecondLeg(probes.at(j),0.2,probetrigger,"muon")){
								mmumupass.at(t).Fill((Ntp->Muon_p4(tags.at(i))+Ntp->Muon_p4(probes.at(j))).M(),w);
								muleg_numerator.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),fabs(Ntp->Muon_p4(probes.at(j)).Eta()),w);
								muleg_numerator_rebin.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),fabs(Ntp->Muon_p4(probes.at(j)).Eta()),w);
								if(fabs(Ntp->Muon_p4(probes.at(j)).Eta())<0.9) muleg_num_09.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
								else if(fabs(Ntp->Muon_p4(probes.at(j)).Eta())<1.2) muleg_num_12.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
								else if(fabs(Ntp->Muon_p4(probes.at(j)).Eta())<2.1) muleg_num_21.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
								else if(fabs(Ntp->Muon_p4(probes.at(j)).Eta())<2.4) muleg_num_24.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
								muleg_num_all.at(t).Fill(Ntp->Muon_p4(probes.at(j)).Pt(),w);
							}
						}
				}
			}
		}
		if(object=="electron"){
			for(unsigned i=0;i<Ntp->NElectrons();i++){
				if(Ntp->Electron_p4(i).Pt()<8) continue;
				if(fabs(Ntp->Electron_supercluster_eta(i))>2.5) continue;
				if(!Ntp->isData() && Ntp->Electron_p4(i).Pt()<20 && !matchTruth(i,0.2,"electron")) continue;
				if(!isMVATrigElectron(i)) continue;
				if(Electron_RelIso(i)>=0.15) continue;
				if(matchTrigger(i,0.2,tagtrigger,"electron")
						&& Ntp->Electron_p4(i).Pt()>29.
						&& matchFirstLeg(i,0.2,probetrigger,"electron")
						){
					ntag++;
					tags.push_back(i);
				}else{
					nprobe++;
					probes.push_back(i);
				}
			}
			ntags.at(t).Fill(ntag);
			nprobes.at(t).Fill(nprobe);
			for(unsigned i=0;i<tags.size();i++){
				tagept.at(t).Fill(Ntp->Electron_p4(tags.at(i)).Pt(),w);
				tageeta.at(t).Fill(Ntp->Electron_supercluster_eta(tags.at(i)),w);
				mttage.at(t).Fill(sqrt(2*Ntp->Electron_p4(tags.at(i)).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Electron_p4(tags.at(i)).Px(),Ntp->Electron_p4(tags.at(i)).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
				for(unsigned j=0;j<probes.size();j++){
					if(Ntp->Electron_Charge(tags.at(i))!=Ntp->Electron_Charge(probes.at(j))
							&& (Ntp->Electron_p4(tags.at(i))+Ntp->Electron_p4(probes.at(j))).M()>=60.
							&& (Ntp->Electron_p4(tags.at(i))+Ntp->Electron_p4(probes.at(j))).M()<120.
							){
						probeept.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
						probeeeta.at(t).Fill(Ntp->Electron_supercluster_eta(probes.at(j)),w);
						mtprobee.at(t).Fill(sqrt(2*Ntp->Electron_p4(probes.at(j)).Pt()*Ntp->MET_CorrT0pcT1_et()*(1-cosphi2d(Ntp->Electron_p4(probes.at(j)).Px(),Ntp->Electron_p4(probes.at(j)).Py(),Ntp->MET_CorrT0pcT1_ex(),Ntp->MET_CorrT0pcT1_ey()))),w);
						dree.at(t).Fill(Ntp->Electron_p4(tags.at(i)).DeltaR(Ntp->Electron_p4(probes.at(j))),w);
						ptbalee.at(t).Fill((Ntp->Electron_p4(tags.at(i))+Ntp->Electron_p4(probes.at(j))).Pt(),w);
						mee.at(t).Fill((Ntp->Electron_p4(tags.at(i))+Ntp->Electron_p4(probes.at(j))).M(),w);
						eleg_denominator.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),fabs(Ntp->Electron_supercluster_eta(probes.at(j))),w);
						eleg_denominator_rebin.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),fabs(Ntp->Electron_supercluster_eta(probes.at(j))),w);
						if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<0.8) eleg_denom_08.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
						else if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<1.4442) eleg_denom_14.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
						else if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<1.5) eleg_denom_15.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
						else if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<2.) eleg_denom_20.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
						else if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<2.5) eleg_denom_25.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
						eleg_denom_all.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
						if(!matchSecondLeg(probes.at(j),0.2,probetrigger,"electron")) meefail.at(t).Fill((Ntp->Electron_p4(tags.at(i))+Ntp->Electron_p4(probes.at(j))).M(),w);
						if(matchSecondLeg(probes.at(j),0.2,probetrigger,"electron")){
							meepass.at(t).Fill((Ntp->Electron_p4(tags.at(i))+Ntp->Electron_p4(probes.at(j))).M(),w);
							eleg_numerator.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),fabs(Ntp->Electron_supercluster_eta(probes.at(j))),w);
							eleg_numerator_rebin.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),fabs(Ntp->Electron_supercluster_eta(probes.at(j))),w);
							if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<0.8) eleg_num_08.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
							else if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<1.4442) eleg_num_14.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
							else if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<1.5) eleg_num_15.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
							else if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<2.) eleg_num_20.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
							else if(fabs(Ntp->Electron_supercluster_eta(probes.at(j)))<2.5) eleg_num_25.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
							eleg_num_all.at(t).Fill(Ntp->Electron_p4(probes.at(j)).Pt(),w);
						}
					}
				}
			}
		}
	}
}

//////////////////////////////
//
// Muon related functions
//

bool ZtoEMu_Fakerate::isTightMuon(unsigned int i){
	if(!Ntp->Muon_isGlobalMuon(i)) return false;
	if(!Ntp->Muon_isPFMuon(i)) return false;
	if(Ntp->Muon_normChi2(i)>=10.) return false;
	if(Ntp->Muon_hitPattern_numberOfValidMuonHits(i)<=0) return false;
	if(Ntp->Muon_numberOfMatchedStations(i)<=1) return false;
	if(Ntp->Muon_numberofValidPixelHits(i)<=0) return false;
	if(Ntp->Muon_trackerLayersWithMeasurement(i)<=5) return false;
	return true;
}

bool ZtoEMu_Fakerate::isTightMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.2) return false;
	if(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.5) return false;
	return true;
}

bool ZtoEMu_Fakerate::isHiggsMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	if(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
}

bool ZtoEMu_Fakerate::isLooseMuon(unsigned int i){
	if(!Ntp->Muon_isPFMuon(i)) return false;
	if(!(Ntp->Muon_isGlobalMuon(i) || Ntp->Muon_isTrackerMuon(i))) return false;
	return true;
}

bool ZtoEMu_Fakerate::isFakeMuon(unsigned int i){
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

bool ZtoEMu_Fakerate::isFakeMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isFakeMuon(i)) return false;
	if(dxy(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.2) return false;
	if(dz(Ntp->Muon_p4(i),Ntp->Muon_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
}

double ZtoEMu_Fakerate::Muon_RelIso(unsigned int i){
	return (Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i)))/Ntp->Muon_p4(i).Pt();
}

double ZtoEMu_Fakerate::Muon_AbsIso(unsigned int i){
	return Ntp->Muon_sumChargedHadronPt04(i)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(i)+Ntp->Muon_sumPhotonEt04(i)-0.5*Ntp->Muon_sumPUPt04(i));
}

//////////////////////////////
//
// Electron related functions
//

bool ZtoEMu_Fakerate::isTrigPreselElectron(unsigned int i){
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

bool ZtoEMu_Fakerate::isTrigNoIPPreselElectron(unsigned int i){
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

bool ZtoEMu_Fakerate::isMVATrigNoIPElectron(unsigned int i){
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

bool ZtoEMu_Fakerate::isMVANonTrigElectron(unsigned int i, unsigned int j){
	double mvapt = Ntp->Electron_p4(i).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(i));
	if(Ntp->Electron_numberOfMissedHits(i)>1) return false;
	if(vertexSignificance(Ntp->Electron_Poca(i),j)>=4) return false;
	if(mvapt>7. && mvapt<10.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.47) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.004) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.295) return false;
	}
	if(mvapt>=10.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=-0.34) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=-0.65) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(i)<=0.6) return false;
	}
	return true;
}

bool ZtoEMu_Fakerate::isHiggsElectron(unsigned int i, unsigned int j){
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

bool ZtoEMu_Fakerate::isMVATrigElectron(unsigned int i){
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

bool ZtoEMu_Fakerate::isTightElectron(unsigned int i){
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

bool ZtoEMu_Fakerate::isTightElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isTightElectron(i)) return false;
	if(dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	if(dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	return true;
}

bool ZtoEMu_Fakerate::isFakeElectron(unsigned int i){
	if(Ntp->Electron_p4(i).Pt()<=10) return false;
	if(fabs(Ntp->Electron_supercluster_eta(i))>=2.5) return false;
	if(Ntp->Electron_HasMatchedConversions(i)) return false;
	if(Ntp->Electron_numberOfMissedHits(i)>0) return false;
	if(Ntp->Electron_tkSumPt03(i)/Ntp->Electron_p4(i).Pt()>=0.2) return false;
	if(std::max(Ntp->Electron_ecalRecHitSumEt03(i)-1.,0.)/Ntp->Electron_p4(i).Pt()>=0.2) return false;
	if((Ntp->Electron_hcalDepth1TowerSumEt03(i)+Ntp->Electron_hcalDepth2TowerSumEt03(i))/Ntp->Electron_p4(i).Pt()>=0.2) return false;
	if(fabs(Ntp->Electron_supercluster_eta(i))<1.479){
		if(Ntp->Electron_sigmaIetaIeta(i)>=0.01) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.15) return false;
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.007) return false;
		if(Ntp->Electron_hadronicOverEm(i)>=0.12) return false;
	}
	if(fabs(Ntp->Electron_supercluster_eta(i))>=1.479 && fabs(Ntp->Electron_supercluster_eta(i))<2.5){
		if(Ntp->Electron_sigmaIetaIeta(i)>=0.03) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.10) return false;
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.009) return false;
		if(Ntp->Electron_hadronicOverEm(i)>=0.10) return false;
	}
	return true;
}

bool ZtoEMu_Fakerate::isFakeElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=Ntp->NVtx()) return false;
	if(!isFakeElectron(i)) return false;
	if(dz(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.1) return false;
	if(dxy(Ntp->Electron_p4(i),Ntp->Electron_Poca(i),Ntp->Vtx(j))>=0.02) return false;
	return true;
}

double ZtoEMu_Fakerate::Electron_RelIso(unsigned int i){
	return (Ntp->Electron_chargedHadronIso(i)+std::max((double)0.,Ntp->Electron_neutralHadronIso(i)+Ntp->Electron_photonIso(i)-Ntp->RhoIsolationAllInputTags()*Electron_Aeff_R04(Ntp->Electron_supercluster_eta(i))))/Ntp->Electron_p4(i).Pt();
}

double ZtoEMu_Fakerate::Electron_Aeff_R04(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.208;
	else if(eta>=1. && eta<1.479) return 0.209;
	else if(eta>=1.479 && eta<2.) return 0.115;
	else if(eta>=2. && eta<2.2) return 0.143;
	else if(eta>=2.2 && eta<2.3) return 0.183;
	else if(eta>=2.3 && eta<2.3) return 0.194;
	else if(eta>=2.4) return 0.261;
}

double ZtoEMu_Fakerate::Electron_Aeff_R03(double Eta){
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
// Finish function
//

void ZtoEMu_Fakerate::Finish(){
	// computing fake rates
	double muprob;
	double eprob;
	for(unsigned int i=1;i<=fakemu_rebin.at(0).GetNbinsX();i++){
		for(unsigned int j=1;j<=fakemu_rebin.at(0).GetNbinsY();j++){
			if(fakemu_rebin.at(0).GetBinContent(i,j)!=0){
				muprob = tightmu_rebin.at(0).GetBinContent(i,j)/fakemu_rebin.at(0).GetBinContent(i,j);
				mueff.at(0).SetBinContent(i,j,muprob);
			}
		}
	}
	for(unsigned int i=1;i<=fakee_rebin.at(0).GetNbinsX();i++){
		for(unsigned int j=1;j<=fakee_rebin.at(0).GetNbinsY();j++){
			if(fakee_rebin.at(0).GetBinContent(i,j)!=0){
				eprob = tighte_rebin.at(0).GetBinContent(i,j)/fakee_rebin.at(0).GetBinContent(i,j);
				eeff.at(0).SetBinContent(i,j,eprob);
			}
		}
	}

	// computing trigger efficiencies
	double mueff_data, mueff_mumu, mueff_tautau;
	double eeff_data, eeff_ee, eeff_tautau;
	for(unsigned i=1;i<=muleg_numerator_rebin.at(0).GetNbinsX();i++){
		for(unsigned j=1;j<=muleg_numerator_rebin.at(0).GetNbinsY();j++){
			if(muleg_denominator_rebin.at(0).GetBinContent(i,j)==0) mueff_data = 0;
			if(muleg_denominator_rebin.at(0).GetBinContent(i,j)>0) mueff_data = muleg_numerator_rebin.at(0).GetBinContent(i,j)/muleg_denominator_rebin.at(0).GetBinContent(i,j);
			if(muleg_denominator_rebin.at(12).GetBinContent(i,j)==0) mueff_tautau = 0;
			if(muleg_denominator_rebin.at(12).GetBinContent(i,j)>0) mueff_tautau = muleg_numerator_rebin.at(12).GetBinContent(i,j)/muleg_denominator_rebin.at(12).GetBinContent(i,j);
			if(muleg_denominator_rebin.at(13).GetBinContent(i,j)==0) mueff_mumu = 0;
			if(muleg_denominator_rebin.at(13).GetBinContent(i,j)>0) mueff_mumu = muleg_numerator_rebin.at(13).GetBinContent(i,j)/muleg_denominator_rebin.at(13).GetBinContent(i,j);

			muleg_eff.at(0).SetBinContent(i,j,mueff_data);
			muleg_eff.at(12).SetBinContent(i,j,mueff_tautau);
			muleg_eff.at(13).SetBinContent(i,j,mueff_mumu);
		}
	}
	for(unsigned i=1;i<=eleg_numerator_rebin.at(0).GetNbinsX();i++){
		for(unsigned j=1;j<=eleg_numerator_rebin.at(0).GetNbinsY();j++){
			if(eleg_denominator_rebin.at(0).GetBinContent(i,j)==0) eeff_data = 0;
			if(eleg_denominator_rebin.at(0).GetBinContent(i,j)>0) eeff_data = eleg_numerator_rebin.at(0).GetBinContent(i,j)/eleg_denominator_rebin.at(0).GetBinContent(i,j);
			if(eleg_denominator_rebin.at(12).GetBinContent(i,j)==0) eeff_tautau = 0;
			if(eleg_denominator_rebin.at(12).GetBinContent(i,j)>0) eeff_tautau = eleg_numerator_rebin.at(12).GetBinContent(i,j)/eleg_denominator_rebin.at(12).GetBinContent(i,j);
			if(eleg_denominator_rebin.at(14).GetBinContent(i,j)==0) eeff_ee = 0;
			if(eleg_denominator_rebin.at(14).GetBinContent(i,j)>0) eeff_ee = eleg_numerator_rebin.at(14).GetBinContent(i,j)/eleg_denominator_rebin.at(14).GetBinContent(i,j);

			eleg_eff.at(0).SetBinContent(i,j,eeff_data);
			eleg_eff.at(12).SetBinContent(i,j,eeff_tautau);
			eleg_eff.at(14).SetBinContent(i,j,eeff_ee);
		}
	}

	double muscale(0);
	double escale(0);
	for(unsigned i=1;i<=muleg_scale.at(0).GetNbinsX();i++){
		for(unsigned j=1;j<=muleg_scale.at(0).GetNbinsY();j++){
			if(muleg_eff.at(13).GetBinContent(i,j)>0) muscale = muleg_eff.at(0).GetBinContent(i,j)/muleg_eff.at(13).GetBinContent(i,j);
			muleg_scale.at(0).SetBinContent(i,j,muscale);
		}
	}
	for(unsigned i=1;i<=eleg_scale.at(0).GetNbinsX();i++){
		for(unsigned j=1;j<=eleg_scale.at(0).GetNbinsY();j++){
			if(eleg_eff.at(14).GetBinContent(i,j)>0) escale = eleg_eff.at(0).GetBinContent(i,j)/eleg_eff.at(14).GetBinContent(i,j);
			eleg_scale.at(0).SetBinContent(i,j,escale);
		}
	}

	//mueff.at(0).Write();
	//eeff.at(0).Write();
	//outfile->Write();
	//outfile->Close();
	Selection::Finish();
}
