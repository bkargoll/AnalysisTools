#include "ZtoEMu_Skim.h"
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

ZtoEMu_Skim::ZtoEMu_Skim(TString Name_, TString id_):
  Selection(Name_,id_)
  ,mu_ptlow(10)
  ,mu_pthigh(20)
  ,e_ptlow(10)
  ,e_pthigh(20)
  ,mu_eta(2.4)
  ,e_eta(2.5)
{
    //verbose=true;

	doHiggsObjects = false;
	doWWObjects = true;
}

ZtoEMu_Skim::~ZtoEMu_Skim(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "ZtoEMu_Skim::~ZtoEMu_Skim Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "ZtoEMu_Skim::~ZtoEMu_Skim()" << std::endl;
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

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");
  // Setup Extra Histograms
  NPV=HConfig.GetTH1D(Name+"_NPV","NPV",60,0.,60.,"Number of vertices");
  mupt=HConfig.GetTH1D(Name+"_mupt","mupt",40,0.,100.,"p_{T}^{#mu}");
  mueta=HConfig.GetTH1D(Name+"_mueta","mueta",20,-2.5,2.5,"#eta_{#mu}");
  ept=HConfig.GetTH1D(Name+"_ept","ept",40,0.,100.,"p_{T}^{e}");
  eeta=HConfig.GetTH1D(Name+"_eeta","eeta",20,-2.5,2.5,"#eta_{e}");
  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  ZtoEMu_Skim::Store_ExtraDist(){
	Extradist1d.push_back(&NPV);
	Extradist1d.push_back(&mupt);
	Extradist1d.push_back(&mueta);
	Extradist1d.push_back(&ept);
	Extradist1d.push_back(&eeta);
}

void  ZtoEMu_Skim::doEvent(){
  if(verbose)std::cout << "ZtoEMu_Skim::doEvent() START" << std::endl;
  unsigned int t;
  int id(Ntp->GetMCID());
  if(verbose)std::cout << "id: " << id << std::endl;
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" << std::endl; return;}

  // Apply Selection
    
  ///////////////////////////////////////////////
  //
  // Trigger passed?
  //
  if(verbose)std::cout << "Trigger passed?" << std::endl;
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
	  if(isGoodVtx(i)){
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
  
  // muon ID cuts
  for(unsigned i=0;i<Ntp->NMuons();i++){
	  if(Ntp->Muon_p4(i).Pt()>mu_ptlow
			  && fabs(Ntp->Muon_p4(i).Eta())<mu_eta
			  && vertex>=0
			  && (matchTrigger(i,0.2,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","muon") || matchTrigger(i,0.2,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","muon") || Ntp->GetMCID()==DataMCType::DY_emu_embedded)
			  ){
		  if(doHiggsObjects){
			  if(isHiggsMuon(i,vertex)
					  && ((fabs(Ntp->Muon_p4(i).Eta())<1.479 && Muon_RelIso(i)<0.15) || (fabs(Ntp->Muon_p4(i).Eta())>=1.479 && Muon_RelIso(i)<0.10))
							  ){
				  GoodMuons.push_back(i);
			  }else if(isFakeMuon(i,vertex)
					  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
					  ){
				  GoodMuons.push_back(i);
			  }
		  }else{
			  if(isTightMuon(i,vertex)
					  && Muon_RelIso(i)<0.12
					  ){
				  GoodMuons.push_back(i);
			  }else if(doWWObjects
					  && isFakeMuon(i,vertex)
					  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
					  ){
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
  bool matchRecoMuon = false;
  
  // electron ID cuts (eta-dependent MVA or simple cut-based)
  for(unsigned i=0;i<Ntp->NElectrons();i++){
	  hasMuonTrack = false;
	  matchRecoMuon = false;
	  if(Ntp->Electron_p4(i).Et()>e_ptlow
			  && fabs(Ntp->Electron_supercluster_eta(i))<e_eta
			  && vertex>=0
			  && (matchTrigger(i,0.2,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","electron") || matchTrigger(i,0.2,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v","electron") || Ntp->GetMCID()==DataMCType::DY_emu_embedded)
			  ){
		  // don't use electrons with tracks matching those of selected muons
		  for(unsigned j=0;j<GoodMuons.size();j++){
			  if(Ntp->Electron_Track_idx(i)==Ntp->Muon_Track_idx(GoodMuons.at(j))) hasMuonTrack = true;
		  }
		  if(hasMuonTrack) continue;
		  // no overlapping reco muons
		  for(unsigned j=0;j<Ntp->NMuons();j++){
			  if(j==muidx)continue;
			  if(Ntp->Electron_p4(i).DeltaR(Ntp->Muon_p4(j))<0.3) matchRecoMuon = true;
		  }
		  if(matchRecoMuon) continue;
		  if(doHiggsObjects){
			  if(isHiggsElectron(i,vertex)
					  && ((fabs(Ntp->Electron_supercluster_eta(i))<1.479 && Electron_RelIso(i)<0.15) || (fabs(Ntp->Electron_supercluster_eta(i))>=1.479 && Electron_RelIso(i)<0.10))
					  ){
				  GoodElectrons.push_back(i);
			  }else if(isFakeElectron(i,vertex)
					  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
					  ){
				  GoodElectrons.push_back(i);
			  }
		  }else{
			  if(isWWElectron(i,vertex)
					  && Electron_RelIso(i)<0.15
					  ){
				  GoodElectrons.push_back(i);
			  }else if(doWWObjects
					  && isFakeElectron(i,vertex)
					  && (Ntp->isData() || Ntp->GetMCID()==DataMCType::DY_ee || Ntp->GetMCID()==DataMCType::DY_mumu || Ntp->GetMCID()==DataMCType::DY_tautau || Ntp->GetMCID()==DataMCType::DY_ll)
					  ){
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
  bool passembed = false;
  bool leadingmu = false;
  value.at(ptthreshold)=0;
  if(muidx!=999 && eidx!=999){
	  value.at(ptthreshold)=1;
	  if(Ntp->Muon_p4(muidx).Pt()<=mu_ptlow || Ntp->Electron_p4(eidx).Et()<=e_ptlow) value.at(ptthreshold)=0;
	  if(Ntp->Muon_p4(muidx).Pt()<mu_pthigh && Ntp->Electron_p4(eidx).Et()<e_pthigh) value.at(ptthreshold)=0;
	  if(value.at(ptthreshold)==1 && Ntp->GetMCID()==DataMCType::DY_emu_embedded) passembed = true;
	  if(Ntp->Muon_p4(muidx).Pt()<mu_pthigh){
		  if(!Ntp->TriggerAccept("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) value.at(ptthreshold)=0;
	  }else if(Ntp->Electron_p4(eidx).Et()<e_pthigh){
		  if(!Ntp->TriggerAccept("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v")) value.at(ptthreshold)=0;
		  else leadingmu = true;
	  }
	  if(passembed) value.at(ptthreshold)=1;
  }
  pass.at(ptthreshold)=(value.at(ptthreshold)==cut.at(ptthreshold));
  
  //////////////////////////////////////////////////////////
  if(verbose) std::cout << "do weights" << std::endl;
  double wobs(1),w(1);
  if(!Ntp->isData() && Ntp->GetMCID()!=DataMCType::DY_emu_embedded){
    //w*=Ntp->PUWeight();
  }
  else{w=1;wobs=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs);
  if(verbose)std::cout << "status: " << status << std::endl;
  ///////////////////////////////////////////////////////////
  // Add plots
  if(verbose) std::cout << "add plots" << std::endl;
  NPV.at(t).Fill(Ntp->NVtx(),w);
  if(pass.at(TriggerOk)
		  && pass.at(PrimeVtx)
		  && pass.at(NMu)
		  && pass.at(NE)
		  && pass.at(ptthreshold)
		  ){
	  NPV.at(t).Fill(Ntp->NVtx(),w);
	  mupt.at(t).Fill(Ntp->Muon_p4(muidx).Pt(),w);
	  mueta.at(t).Fill(Ntp->Muon_p4(muidx).Eta(),w);
	  ept.at(t).Fill(Ntp->Electron_p4(eidx).Pt(),w);
	  eeta.at(t).Fill(Ntp->Electron_supercluster_eta(eidx),w);
  }
  if(verbose)std::cout << "ZtoEMu_Skim::doEvent() doEvent END" << std::endl;
}

//////////////////////
//
// Utility functions
//

bool ZtoEMu_Skim::isGoodVtx(unsigned int i){
	if(fabs(Ntp->Vtx(i).z())>=24) return false;
	if(Ntp->Vtx(i).Perp()>=2) return false;
	if(Ntp->Vtx_ndof(i)<=4) return false;
	if(Ntp->Vtx_isFake(i)!=0) return false;
	return true;
}

bool ZtoEMu_Skim::jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx){
	for(unsigned i=0;i<vtx_track_idx.size();i++){
		if(vtx_track_idx.at(i)==leadingtrack_idx)return true;
	}
	return false;
}

double ZtoEMu_Skim::cosphi2d(double px1, double py1, double px2, double py2){
	return (px1*px2+py1*py2)/(sqrt(pow(px1,2)+pow(py1,2))*sqrt(pow(px2,2)+pow(py2,2)));
}

double ZtoEMu_Skim::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

double ZtoEMu_Skim::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}

double ZtoEMu_Skim::vertexSignificance(TVector3 vec, unsigned int vertex){
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

bool ZtoEMu_Skim::matchTrigger(unsigned int i, double dr, std::string trigger, std::string object){
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

int ZtoEMu_Skim::matchTruth(TLorentzVector tvector){
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

bool ZtoEMu_Skim::matchTruth(TLorentzVector tvector, int pid, double dr){
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

bool ZtoEMu_Skim::isTightMuon(unsigned int idx){
	if(!Ntp->Muon_isGlobalMuon(idx)) return false;
	if(!Ntp->Muon_isPFMuon(idx)) return false;
	if(Ntp->Muon_normChi2(idx)>=10.) return false;
	if(Ntp->Muon_hitPattern_numberOfValidMuonHits(idx)<=0) return false;
	if(Ntp->Muon_numberOfMatchedStations(idx)<=1) return false;
	if(Ntp->Muon_numberofValidPixelHits(idx)<=0) return false;
	if(Ntp->Muon_trackerLayersWithMeasurement(idx)<=5) return false;
	return true;
}

bool ZtoEMu_Skim::isTightMuon(unsigned int idx, unsigned int vtx){
	if(vtx<0 || vtx>=Ntp->NVtx()) return false;
	if(!isTightMuon(idx)) return false;
	if(dxy(Ntp->Muon_p4(idx),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.2) return false;
	if(dz(Ntp->Muon_p4(idx),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.5) return false;
	return true;
}

bool ZtoEMu_Skim::isHiggsMuon(unsigned int idx, unsigned int vtx){
	if(vtx<0 || vtx>=Ntp->NVtx()) return false;
	if(!isTightMuon(idx)) return false;
	if(dxy(Ntp->Muon_p4(idx),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.02) return false;
	if(dz(Ntp->Muon_p4(idx),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.1) return false;
	return true;
}

bool ZtoEMu_Skim::isLooseMuon(unsigned int idx){
	if(!Ntp->Muon_isPFMuon(idx)) return false;
	if(!(Ntp->Muon_isGlobalMuon(idx) || Ntp->Muon_isTrackerMuon(idx))) return false;
	return true;
}

bool ZtoEMu_Skim::isFakeMuon(unsigned int idx){
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

bool ZtoEMu_Skim::isFakeMuon(unsigned int idx, unsigned int vtx){
	if(vtx<0 || vtx>=Ntp->NVtx()) return false;
	if(!isFakeMuon(idx)) return false;
	if(dxy(Ntp->Muon_p4(idx),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.2) return false;
	if(dz(Ntp->Muon_p4(idx),Ntp->Muon_Poca(idx),Ntp->Vtx(vtx))>=0.1) return false;
	return true;
}

double ZtoEMu_Skim::Muon_RelIso(unsigned int idx){
	return (Ntp->Muon_sumChargedHadronPt04(idx)+std::max(0.,Ntp->Muon_sumNeutralHadronEt04(idx)+Ntp->Muon_sumPhotonEt04(idx)-0.5*Ntp->Muon_sumPUPt04(idx)))/Ntp->Muon_p4(idx).Pt();
}

//////////////////////////////
//
// Electron related functions
//

bool ZtoEMu_Skim::isTrigPreselElectron(unsigned int idx){
	if(fabs(Ntp->Electron_supercluster_eta(idx))>2.5) return false;
	if(Ntp->Electron_numberOfMissedHits(idx)>0) return false;
	if(Ntp->Electron_Gsf_dr03TkSumPt(idx)/Ntp->Electron_p4(idx).Pt()>0.2) return false;
	if(Ntp->Electron_Gsf_dr03HcalTowerSumEt(idx)/Ntp->Electron_p4(idx).Pt()>0.2) return false;
	if(fabs(Ntp->Electron_supercluster_eta(idx))<1.479){
		if(Ntp->Electron_sigmaIetaIeta(idx)>0.014) return false;
		if(Ntp->Electron_hadronicOverEm(idx)>0.15) return false;
	}else{
		if(Ntp->Electron_sigmaIetaIeta(idx)>0.035) return false;
		if(Ntp->Electron_hadronicOverEm(idx)>0.1) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isTrigNoIPPreselElectron(unsigned int idx){
	if(fabs(Ntp->Electron_supercluster_eta(idx))>2.5) return false;
	if(Ntp->Electron_numberOfMissedHits(idx)>0) return false;
	if(Ntp->Electron_Gsf_dr03TkSumPt(idx)/Ntp->Electron_p4(idx).Pt()>0.2) return false;
	if(Ntp->Electron_Gsf_dr03HcalTowerSumEt(idx)/Ntp->Electron_p4(idx).Pt()>0.2) return false;
	if(fabs(Ntp->Electron_supercluster_eta(idx))<1.479){
		if(Ntp->Electron_sigmaIetaIeta(idx)>0.01) return false;
		if(Ntp->Electron_hadronicOverEm(idx)>0.12) return false;
		if(fabs(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(idx))>0.007) return false;
		if(fabs(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(idx))>0.15) return false;
	}else{
		if(Ntp->Electron_sigmaIetaIeta(idx)>0.03) return false;
		if(Ntp->Electron_hadronicOverEm(idx)>0.1) return false;
		if(fabs(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(idx))>0.009) return false;
		if(fabs(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(idx))>0.1) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isMVATrigNoIPElectron(unsigned int idx){
	double mvapt = Ntp->Electron_p4(idx).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(idx));
	if(mvaeta>2.5) return false;
	if(Ntp->Electron_HasMatchedConversions(idx)) return false;
	if(Ntp->Electron_numberOfMissedHits(idx)>0) return false;
	if(!isTrigNoIPPreselElectron(idx)) return false;
	if(mvaeta<1.479 && Electron_RelIso(idx)>0.15) return false;// TODO: correct isolation cuts?
	if(mvaeta>=1.479 && Electron_RelIso(idx)>0.10) return false;
	if(mvapt<20){
		if(mvaeta<0.8 && Ntp->Electron_MVA_TrigNoIP_discriminator(idx)<=-0.5375) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_TrigNoIP_discriminator(idx)<=-0.375) return false;
		if(mvaeta>=1.479 && Ntp->Electron_MVA_TrigNoIP_discriminator(idx)<=-0.025) return false;
	}
	if(mvapt>=20){
		if(mvaeta<0.8 && Ntp->Electron_MVA_TrigNoIP_discriminator(idx)<=0.325) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_TrigNoIP_discriminator(idx)<=0.775) return false;
		if(mvaeta>=1.479 && Ntp->Electron_MVA_TrigNoIP_discriminator(idx)<=0.775) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isMVANonTrigElectron(unsigned int idx, unsigned int vtx){
	double mvapt = Ntp->Electron_p4(idx).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(idx));
	if(mvapt<7.) return false;
	if(mvaeta>2.5) return false;
	if(Ntp->Electron_numberOfMissedHits(idx)>1) return false;
	if(vertexSignificance(Ntp->Electron_Poca(idx),vtx)>=4) return false;
	if(mvapt>7. && mvapt<10.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.47) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.004) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.295) return false;
	}
	if(mvapt>=10.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=-0.34) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=-0.65) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.6) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isHiggsElectron(unsigned int idx, unsigned int vtx){
	double mvapt = Ntp->Electron_p4(idx).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(idx));
	if(mvaeta>2.5) return false;
	if(Ntp->Electron_numberOfMissedHits(idx)>0) return false;
	if(Ntp->Electron_HasMatchedConversions(idx)) return false;
	if(dxy(Ntp->Electron_p4(idx),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>=0.02) return false;
	if(dz(Ntp->Electron_p4(idx),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>=0.1) return false;
	if(mvapt<20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.925) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.915) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.965) return false;
	}
	if(mvapt>=20.){
		if(mvaeta<0.8 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.905) return false;
		if(mvaeta>=0.8 && mvaeta<1.479 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.955) return false;
		if(mvaeta>=1.479 && mvaeta<2.5 && Ntp->Electron_MVA_NonTrig_discriminator(idx)<=0.975) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isWWElectron(unsigned int idx, unsigned int vtx){
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

bool ZtoEMu_Skim::isMVATrigElectron(unsigned int idx){
	double mvapt = Ntp->Electron_p4(idx).Pt();
	double mvaeta = fabs(Ntp->Electron_supercluster_eta(idx));
	if(mvapt<10.) return false;
	if(mvaeta>2.5) return false;
	if(Ntp->Electron_numberOfMissedHits(idx)>0) return false;
	if(Ntp->Electron_HasMatchedConversions(idx)) return false;
	if(!isTrigPreselElectron(idx)) return false;
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

bool ZtoEMu_Skim::isTightElectron(unsigned int idx){
	if(Ntp->Electron_HasMatchedConversions(idx)) return false;
	if(Ntp->Electron_numberOfMissedHits(idx)>0) return false;
	if(Electron_RelIso(idx)>=0.1) return false;
	if(fabs(1/Ntp->Electron_ecalEnergy(idx)-1/Ntp->Electron_trackMomentumAtVtx(idx))>=0.05) return false;
	if(fabs(Ntp->Electron_supercluster_eta(idx))<=1.479){
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(idx)>=0.004) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(idx)>=0.03) return false;
		if(Ntp->Electron_sigmaIetaIeta(idx)>=0.01) return false;
		if(Ntp->Electron_hadronicOverEm(idx)>=0.12) return false;
	}
	if(fabs(Ntp->Electron_supercluster_eta(idx))>1.479 && fabs(Ntp->Electron_supercluster_eta(idx))<2.5){
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(idx)>=0.005) return false;
		if(Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(idx)>=0.02) return false;
		if(Ntp->Electron_sigmaIetaIeta(idx)>=0.03) return false;
		if(Ntp->Electron_hadronicOverEm(idx)>=0.10) return false;
		if(Ntp->Electron_p4(idx).Pt()<20 && Electron_RelIso(idx)>=0.07) return false;
	}
	return true;
}

bool ZtoEMu_Skim::isTightElectron(unsigned int idx, unsigned int vtx){
	if(vtx<0 || vtx>=Ntp->NVtx()) return false;
	if(!isTightElectron(idx)) return false;
	if(dxy(Ntp->Electron_p4(idx),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>=0.02) return false;
	if(dz(Ntp->Electron_p4(idx),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>=0.1) return false;
	return true;
}

bool ZtoEMu_Skim::isFakeElectron(unsigned int idx){
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

bool ZtoEMu_Skim::isFakeElectron(unsigned int idx, unsigned int vtx){
	if(vtx<0 || vtx>=Ntp->NVtx()) return false;
	if(!isFakeElectron(idx)) return false;
	if(dz(Ntp->Electron_p4(idx),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>0.1) return false;
	if(dxy(Ntp->Electron_p4(idx),Ntp->Electron_Poca(idx),Ntp->Vtx(vtx))>0.02) return false;
	return true;
}

double ZtoEMu_Skim::Electron_RelIso(unsigned int idx){
	return (Ntp->Electron_chargedHadronIso(idx)+std::max((double)0.,Ntp->Electron_neutralHadronIso(idx)+Ntp->Electron_photonIso(idx)-Ntp->RhoIsolationAllInputTags()*Electron_Aeff_R04(Ntp->Electron_supercluster_eta(idx))))/Ntp->Electron_p4(idx).Pt();
}

double ZtoEMu_Skim::Electron_Aeff_R04(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.208;
	else if(eta>=1. && eta<1.479) return 0.209;
	else if(eta>=1.479 && eta<2.) return 0.115;
	else if(eta>=2. && eta<2.2) return 0.143;
	else if(eta>=2.2 && eta<2.3) return 0.183;
	else if(eta>=2.3 && eta<2.4) return 0.194;
	else if(eta>=2.4) return 0.261;
}

double ZtoEMu_Skim::Electron_Aeff_R03(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.13;
	else if(eta>=1. && eta<1.479) return 0.14;
	else if(eta>=1.479 && eta<2.) return 0.07;
	else if(eta>=2. && eta<2.2) return 0.09;
	else if(eta>=2.2 && eta<2.3) return 0.11;
	else if(eta>=2.3 && eta<2.4) return 0.11;
	else if(eta>=2.4) return 0.14;
}

bool ZtoEMu_Skim::isLooseElectron(unsigned int idx){
	if(fabs(Ntp->Electron_supercluster_eta(idx))<=1.479){ //barrel
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(idx)<0.007 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(idx)<0.15 &&
				Ntp->Electron_sigmaIetaIeta(idx)<0.01 &&
				Ntp->Electron_hadronicOverEm(idx)<0.12 &&
				fabs(1/Ntp->Electron_ecalEnergy(idx)-1/Ntp->Electron_trackMomentumAtVtx(idx))<0.05 &&
				Electron_RelIso(idx)<0.15 &&
				!Ntp->Electron_HasMatchedConversions(idx) &&
				Ntp->Electron_numberOfMissedHits(idx)<=1
				){
			return true;
		}
	}else if(fabs(Ntp->Electron_supercluster_eta(idx))>1.479 && fabs(Ntp->Electron_supercluster_eta(idx))<2.5){ //endcaps
		if(Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(idx)<0.009 &&
				Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(idx)<0.10 &&
				Ntp->Electron_sigmaIetaIeta(idx)<0.03 &&
				Ntp->Electron_hadronicOverEm(idx)<0.10 &&
				fabs(1/Ntp->Electron_ecalEnergy(idx)-1/Ntp->Electron_trackMomentumAtVtx(idx))<0.05 &&
				!Ntp->Electron_HasMatchedConversions(idx) &&
				Ntp->Electron_numberOfMissedHits(idx)<=1
				){
			if(Ntp->Electron_p4(idx).Pt()>=20.0 && Electron_RelIso(idx)<0.15){
				return true;
			}else if(Ntp->Electron_p4(idx).Pt()<20.0 && Electron_RelIso(idx)<0.10){
				return true;
			}
		}
	}
	return false;
}

//////////////////////////////
//
// Finish function
//

void ZtoEMu_Skim::Finish(){
	Selection::Finish();
}
