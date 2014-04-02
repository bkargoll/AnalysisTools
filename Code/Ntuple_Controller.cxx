//Ntuple_Controller.cxx IMPLEMENTATION FILE


#include "Ntuple_Controller.h"
#include "Tools.h"
#include "PDG_Var.h"

// External code
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

///////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::Ntuple_Controller(std::vector<TString> RootFiles):
  copyTree(false)
  ,ObjEvent(-1)
  ,verbose(false)
{
  // TChains the ROOTuple file
  TChain *chain = new TChain("t");
  std::cout << "Ntuple_Controller" << RootFiles.size() << std::endl;
  for(int i=0; i<RootFiles.size(); i++){
    chain->Add(RootFiles[i]);
  }
  TTree *tree = (TTree*)chain;  
  if(chain==0){
    std::cout << "failed" << std::endl;
  }
  std::cout << "Number of Events in Ntuple: " << chain->GetEntries() << std::endl;
  Ntp=new NtupleReader(tree);
  nbytes=0; 
  nb=0;
  std::cout << "Ntuple Configured" << std::endl;

  // Fit setup 

}

///////////////////////////////////////////////////////////////////////
//
// Function: Int_t Get_Entries()
//
// Purpose: To get the number of events in the Ntuple
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_Entries(){
  std::cout << Ntp->fChain->GetEntries() << std::endl;
  return Int_t(Ntp->fChain->GetEntries());
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the event _jentry
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Get_Event(int _jentry){
  jentry=_jentry;
  Int_t ientry = Ntp->LoadTree(jentry);
  nb = Ntp->fChain->GetEntry(jentry);   nbytes += nb;
}


///////////////////////////////////////////////////////////////////////
//
// Function: void Get_EventIndex()
//
// Purpose: To get the event index (jentry)
//
///////////////////////////////////////////////////////////////////////
Int_t Ntuple_Controller::Get_EventIndex(){
  return jentry;
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Get_Event(int _jentry)
//
// Purpose: To get the file name of the root file currently being 
//          accesses
//
///////////////////////////////////////////////////////////////////////
TString Ntuple_Controller::Get_File_Name(){
  return Ntp->fChain->GetCurrentFile()->GetName();
}

///////////////////////////////////////////////////////////////////////
//
// Function: void Branch_Setup(TString B_Name, int type)
//
// Purpose: To setup a branch
//
///////////////////////////////////////////////////////////////////////
void Ntuple_Controller::Branch_Setup(TString B_Name, int type){   
  Ntp->fChain->SetBranchStatus(B_Name,type);
}

///////////////////////////////////////////////////////////////////////
//
// destructor
//
///////////////////////////////////////////////////////////////////////
Ntuple_Controller::~Ntuple_Controller() {
  std::cout << "Ntuple_Controller::~Ntuple_Controller()" << std::endl;
  delete Ntp;
    std::cout << "Ntuple_Controller::~Ntuple_Controller() complete" << std::endl;
}


void Ntuple_Controller::CloneTree(TString n){
  if(!copyTree){
    std::cout << "Starting D3PD cloning " << std::endl;  
    newfile = new TFile(n+".root","recreate");
    SkimmedTree=Ntp->fChain->CloneTree(0);
    copyTree=true;
  }
}

void Ntuple_Controller::SaveCloneTree(){
  std::cout << "Ntuple_Controller::SaveCloneTree()"<< std::endl;
  if(copyTree){
    SkimmedTree->AutoSave();
    newfile->Close();
  }
  std::cout << "Ntuple_Controller::SaveCloneTree() Done"<< std::endl;
}

void Ntuple_Controller::ThinTree(){
  std::cout << "Ntuple_Controller::ThinTree" << std::endl;

  std::cout << "Ntuple_Controller::ThinTree complete" << std::endl;  
}

int Ntuple_Controller::SetupSystematics(TString sys){
  return Default;
}


void Ntuple_Controller::ConfigureObjects(){
  if(ObjEvent!=EventNumber()){
    ObjEvent=EventNumber();
    doElectrons();
    doPhotons();
    doJets();
    doMuons();
    doTaus();
    doMET();
  }
}

void Ntuple_Controller::doElectrons(){
  electrons.clear();
  electrons_default.clear();
}

void Ntuple_Controller::doPhotons(){
  photons.clear();
  photons_default.clear();

}

void Ntuple_Controller::doJets(){
  jets.clear();
  jets_default.clear();
}

void Ntuple_Controller::doMuons(){
  muons.clear();
  muons_default.clear();
}

void Ntuple_Controller::doTaus(){
  taus.clear();
  taus_default.clear();
}

void Ntuple_Controller::doMET(){
}


//Physics get Functions
int Ntuple_Controller::GetMCID(){
	// define signal dependent on analysis (i.e. analyst)
	#if defined(USE_cherepanov) || defined(USE_inugent)
	if((Ntp->DataMC_Type)==DataMCType::DY_ll_Signal && HistoC.hasID(DataMCType::DY_ll_Signal)){
	  for(int i=0;i<NMCSignalParticles();i++){
		  if(abs(MCSignalParticle_pdgid(i))==PDGInfo::Z0){
			  if(fabs(MCSignalParticle_p4(i).M()-PDG_Var::Z_mass())<3*PDG_Var::Z_width()){
				  return DataMCType::Signal;
			  }
		  }
	  }
	  return Ntp->DataMC_Type;
	}
	if(Ntp->DataMC_Type>100){
	  if(HistoC.hasID(Ntp->DataMC_Type%100)){
		  return Ntp->DataMC_Type%100;
	  }
	}
	#endif

    #if defined(USE_nehrkorn)
	if((Ntp->DataMC_Type)==DataMCType::DY_emu && HistoC.hasID(DataMCType::DY_emu))
		return DataMCType::Signal;
	#endif

	#if defined(USE_kargoll)
	if( (Ntp->DataMC_Type)==DataMCType::H_tautau && HistoC.hasID(DataMCType::H_tautau) ||
		(Ntp->DataMC_Type)==DataMCType::H_tautau_ggF && HistoC.hasID(DataMCType::H_tautau_ggF) ||
		(Ntp->DataMC_Type)==DataMCType::H_tautau_VBF && HistoC.hasID(DataMCType::H_tautau_VBF) ||
		(Ntp->DataMC_Type)==DataMCType::H_tautau_WHZHTTH && HistoC.hasID(DataMCType::H_tautau_WHZHTTH) )
			return DataMCType::Signal;
	#endif

	#if defined(USE_pistone)
	if( (Ntp->DataMC_Type)==DataMCType::Hpm_taunu && HistoC.hasID(DataMCType::Hpm_taunu) ||
		(Ntp->DataMC_Type)==DataMCType::HplusBWB && HistoC.hasID(DataMCType::HplusBWB) )
			return  DataMCType::Signal;
	#endif


  if(HConfig.hasID(Ntp->DataMC_Type))return Ntp->DataMC_Type;  
  return -999;
}

TMatrixF Ntuple_Controller::Vtx_Cov(unsigned int i){
  unsigned int dim=3;
  TMatrixF M(dim,dim);
  for(unsigned int j=0;j<dim;j++){
    for(unsigned int k=0;k<=j;k++){
      M[j][k]=Ntp->Vtx_Cov->at(i).at(j).at(k);
      M[k][j]=Ntp->Vtx_Cov->at(i).at(j).at(k);
    }
  }
  return M;
}

bool Ntuple_Controller::isVtxGood(unsigned int i){
  if(0<=i && i<NVtx()){
    if(Vtx_Track_idx(i).size()>4)return true;
  }
  return false;
}

bool Ntuple_Controller::isGoodVtx(unsigned int i){
	if(fabs(Vtx(i).z())>=24) return false;
	if(Vtx(i).Perp()>=2) return false;
	if(Vtx_ndof(i)<=4) return false;
	if(Vtx_isFake(i)!=0) return false;
	return true;
}

bool Ntuple_Controller::isGoodMuon(unsigned int i){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger 
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel       
  //  isGoodMuon_nooverlapremoval(i) with
  //  ΔR(μ,jet)>0.3 where jet is any jet passing the jet requirements not applied applied       
  if(isGoodMuon_nooverlapremoval(i)){
    unsigned int jet_idx=0;
    return !muonhasJetOverlap(i,jet_idx);
  }
  return false;
}

bool Ntuple_Controller::muonhasJetOverlap(unsigned int muon_idx,unsigned int &jet_idx){
  for(unsigned int j=0;j<NPFJets();j++){
    if(isGoodJet_nooverlapremoval(j)){
      if(Tools::dr(Muon_p4(muon_idx),PFJet_p4(j))>0.2 && Tools::dr(Muon_p4(muon_idx),PFJet_p4(j))<0.4){ jet_idx=j;return true;}
    }
  }
  return false;
}

bool Ntuple_Controller::muonhasJetMatch(unsigned int muon_idx,unsigned int &jet_idx){
  for(unsigned int j=0;j<NPFJets();j++){
    if(isGoodJet_nooverlapremoval(j)){
      if(Tools::dr(Muon_p4(muon_idx),PFJet_p4(j))<0.2){ jet_idx=j;return true;}
    }
  }
  return false;
}



bool Ntuple_Controller::isGoodMuon_nooverlapremoval(unsigned int i){
  //  Top Dilepton muon selection without Transverse IP cut and PT cut at 17GeV for our trigger and no overlpar removal applied
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel  
  //  GlobalMuon && TrackerMuon applied applied                                         
  //  pT > 20 GeV                                                                       
  //  fasb(eta) < 2.4                                                                   
  //  normChi2 < 10                                                                     
  //  n Tracker Hits > 10                                                               
  //  n Muon Hit > 0                                                                    
  //  Transverse IP of the muon wrt the beam spot (cm) < 0.02                           
  //  relIso < 0.20                                                                      
  //  fabs( muon.vertex().z() - PV.z()) not applied                           
  //  muon.innerTrack()->hitPattern().pixelLayersWithMeasurement() not applied
  //  numberOfMatchedStations() not applied                              
  if(Muon_isGlobalMuon(i) && Muon_isStandAloneMuon(i)){
    if(Muon_p4(i).Pt()>15.0){
      if(fabs(Muon_p4(i).Eta())<2.4){
	if(Muon_normChi2(i)<10.0){
	  if(Muon_innerTrack_numberofValidHits(i)>10){
	    if(Muon_hitPattern_numberOfValidMuonHits(i)>0){
	      return true;
	    }
	  }
	}
      }
    }
  }
  return false;
}

/////////////////////////////////////////////////////////////////////
//
// Official muon id code
//

bool Ntuple_Controller::isTightMuon(unsigned int i){
	if(!Muon_isGlobalMuon(i)) return false;
	if(!Muon_isPFMuon(i)) return false;
	if(Muon_normChi2(i)>=10.) return false;
	if(Muon_hitPattern_numberOfValidMuonHits(i)<=0) return false;
	if(Muon_numberOfMatchedStations(i)<=1) return false;
	if(Muon_numberofValidPixelHits(i)<=0) return false;
	if(Muon_trackerLayersWithMeasurement(i)<=5) return false;
	return true;
}

bool Ntuple_Controller::isTightMuon(unsigned int i, unsigned int j){
	if(j<0 || j>=NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Muon_p4(i),Muon_Poca(i),Vtx(j))>=0.2) return false;
	if(dz(Muon_p4(i),Muon_Poca(i),Vtx(j))>=0.5) return false;
	return true;
}

bool Ntuple_Controller::isLooseMuon(unsigned int i){
	if(!Muon_isPFMuon(i)) return false;
	if(!(Muon_isGlobalMuon(i) || Muon_isTrackerMuon(i))) return false;
	return true;
}

float Ntuple_Controller::Muon_RelIso(unsigned int i){
	return (Muon_sumChargedHadronPt04(i)+std::max(0.,Muon_sumNeutralHadronEt04(i)+Muon_sumPhotonEt04(i)-0.5*Muon_sumPUPt04(i)))/Muon_p4(i).Pt();
}

/////////////////////////////////////////////////////////////////////

bool Ntuple_Controller::isSelectedMuon(unsigned int i, unsigned int j, double impact_xy, double impact_z){
	if(j<0 || j>=NVtx()) return false;
	if(!isTightMuon(i)) return false;
	if(dxy(Muon_p4(i),Muon_Poca(i),Vtx(j))>=impact_xy) return false;
	if(dz(Muon_p4(i),Muon_Poca(i),Vtx(j))>=impact_z) return false;
	return true;
}

/////////////////////////////////////////////////////////////////////
//
// Official electron id code
//

bool Ntuple_Controller::isTrigPreselElectron(unsigned int i){
	if(fabs(Electron_supercluster_eta(i))>2.5) return false;
	if(fabs(Electron_supercluster_eta(i))<1.479){
		if(Electron_sigmaIetaIeta(i)>0.014) return false;
		if(Electron_hadronicOverEm(i)>0.15) return false;
		if(Electron_Gsf_dr03TkSumPt(i)/Electron_p4(i).Pt()>0.2) return false;
		if(Electron_Gsf_dr03HcalTowerSumEt(i)/Electron_p4(i).Pt()>0.2) return false;
		if(Electron_numberOfMissedHits(i)>0) return false;
	}else{
		if(Electron_sigmaIetaIeta(i)>0.035) return false;
		if(Electron_hadronicOverEm(i)>0.1) return false;
		if(Electron_Gsf_dr03TkSumPt(i)/Electron_p4(i).Pt()>0.2) return false;
		if(Electron_Gsf_dr03HcalTowerSumEt(i)/Electron_p4(i).Pt()>0.2) return false;
		if(Electron_numberOfMissedHits(i)>0) return false;
	}
	return true;
}

bool Ntuple_Controller::isMVATrigElectron(unsigned int i){
	// !!! make sure to also apply Electron_RelIso<0.15 in your analysis !!!
	double mvapt = Electron_p4(i).Pt();
	double mvaeta = fabs(Electron_supercluster_eta(i));
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_HasMatchedConversions(i)) return false;
	if(!isTrigPreselElectron(i)) return false;
	if(mvapt>10. && mvapt<20.){
		if(mvaeta<0.8 && Electron_MVA_Trig_discriminator(i)<=0.00) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_Trig_discriminator(i)<=0.10) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_Trig_discriminator(i)<=0.62) return false;
	}else if(mvapt>=20.){
		if(mvaeta<0.8 && Electron_MVA_Trig_discriminator(i)<=0.94) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_Trig_discriminator(i)<=0.85) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_Trig_discriminator(i)<=0.92) return false;
	}
	return true;
}

bool Ntuple_Controller::isTrigNoIPPreselElectron(unsigned int i){
	if(fabs(Electron_supercluster_eta(i))>2.5) return false;
	if(fabs(Electron_supercluster_eta(i))<1.479){
		if(Electron_sigmaIetaIeta(i)>0.01) return false;
		if(Electron_hadronicOverEm(i)>0.12) return false;
		if(fabs(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i))>0.007) return false;
		if(fabs(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i))>0.15) return false;
		if(Electron_Gsf_dr03TkSumPt(i)/Electron_p4(i).Pt()>0.2) return false;
		if(Electron_Gsf_dr03HcalTowerSumEt(i)/Electron_p4(i).Pt()>0.2) return false;
		if(Electron_numberOfMissedHits(i)>0) return false;
	}else{
		if(Electron_sigmaIetaIeta(i)>0.03) return false;
		if(Electron_hadronicOverEm(i)>0.1) return false;
		if(fabs(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i))>0.009) return false;
		if(fabs(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i))>0.1) return false;
		if(Electron_Gsf_dr03TkSumPt(i)/Electron_p4(i).Pt()>0.2) return false;
		if(Electron_Gsf_dr03HcalTowerSumEt(i)/Electron_p4(i).Pt()>0.2) return false;
		if(Electron_numberOfMissedHits(i)>0) return false;
	}
	return true;
}

bool Ntuple_Controller::isMVATrigNoIPElectron(unsigned int i){
	double mvapt = Electron_p4(i).Pt();
	double mvaeta = fabs(Electron_supercluster_eta(i));
	if(Electron_HasMatchedConversions(i)) return false;
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(!isTrigNoIPPreselElectron(i)) return false;
	if(mvapt<20){
		if(mvaeta<0.8 && Electron_MVA_TrigNoIP_discriminator(i)<=-0.5375) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_TrigNoIP_discriminator(i)<=-0.375) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_TrigNoIP_discriminator(i)<=-0.025) return false;
	}
	if(mvapt>=20){
		if(mvaeta<0.8 && Electron_MVA_TrigNoIP_discriminator(i)<=0.325) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_TrigNoIP_discriminator(i)<=0.775) return false;
	}
	return true;
}

bool Ntuple_Controller::isMVANonTrigElectron(unsigned int i, unsigned int j){
	// !!! make sure to also apply Electron_RelIso<0.4 in your analysis !!!
	double mvapt = Electron_p4(i).Pt();
	double mvaeta = fabs(Electron_supercluster_eta(i));
	if(Electron_numberOfMissedHits(i)>1) return false;
	if(vertexSignificance(Electron_Poca(i),j)>=4) return false;
	if(mvapt>7. && mvapt<10.){
		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=0.47) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=0.004) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.295) return false;
	}
	if(mvapt>=10.){
		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=-0.34) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=-0.65) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.6) return false;
	}
	return true;
}

bool Ntuple_Controller::isTightElectron(unsigned int i){
	if(Electron_HasMatchedConversions(i)) return false;
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_RelIso04(i)>=0.1) return false;
	if(fabs(1/Electron_ecalEnergy(i)-1/Electron_trackMomentumAtVtx(i))>=0.05) return false;
	if(fabs(Electron_supercluster_eta(i))<=1.479){
		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.004) return false;
		if(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.03) return false;
		if(Electron_sigmaIetaIeta(i)>=0.01) return false;
		if(Electron_hadronicOverEm(i)>=0.12) return false;
	}
	if(fabs(Electron_supercluster_eta(i))>1.479 && fabs(Electron_supercluster_eta(i))<2.5){
		if(Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(i)>=0.005) return false;
		if(Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(i)>=0.02) return false;
		if(Electron_sigmaIetaIeta(i)>=0.03) return false;
		if(Electron_hadronicOverEm(i)>=0.10) return false;
		if(Electron_p4(i).Pt()<20 && Electron_RelIso04(i)>=0.07) return false;
	}
	return true;
}

bool Ntuple_Controller::isTightElectron(unsigned int i, unsigned int j){
	if(j<0 || j>=NVtx()) return false;
	if(!isTightElectron(i)) return false;
	if(dxy(Electron_p4(i),Electron_Poca(i),Vtx(j))>=0.02) return false;
	if(dz(Electron_p4(i),Electron_Poca(i),Vtx(j))>=0.1) return false;
	return true;
}

float Ntuple_Controller::Electron_RelIso03(unsigned int i){
	return (Electron_chargedHadronIso(i)+std::max((float)0.,Electron_neutralHadronIso(i)+Electron_photonIso(i)-RhoIsolationAllInputTags()*Electron_Aeff_R03(Electron_supercluster_eta(i))))/Electron_p4(i).Pt();
}

float Ntuple_Controller::Electron_RelIso04(unsigned int i){
	return (Electron_chargedHadronIso(i)+std::max((float)0.,Electron_neutralHadronIso(i)+Electron_photonIso(i)-RhoIsolationAllInputTags()*Electron_Aeff_R04(Electron_supercluster_eta(i))))/Electron_p4(i).Pt();
}

float Ntuple_Controller::Electron_Aeff_R04(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.208;
	else if(eta>=1. && eta<1.479) return 0.209;
	else if(eta>=1.479 && eta<2.) return 0.115;
	else if(eta>=2. && eta<2.2) return 0.143;
	else if(eta>=2.2 && eta<2.3) return 0.183;
	else if(eta>=2.3 && eta<2.3) return 0.194;
	else if(eta>=2.4) return 0.261;
}

float Ntuple_Controller::Electron_Aeff_R03(double Eta){
	double eta=fabs(Eta);
	if(eta>=0. && eta<1.) return 0.13;
	else if(eta>=1. && eta<1.479) return 0.14;
	else if(eta>=1.479 && eta<2.) return 0.07;
	else if(eta>=2. && eta<2.2) return 0.09;
	else if(eta>=2.2 && eta<2.3) return 0.11;
	else if(eta>=2.3 && eta<2.3) return 0.11;
	else if(eta>=2.4) return 0.14;
}

/////////////////////////////////////////////////////////////////////

bool Ntuple_Controller::isSelectedElectron(unsigned int i, unsigned int j, double impact_xy, double impact_z){
	double mvapt = Electron_p4(i).Pt();
	double mvaeta = fabs(Electron_supercluster_eta(i));
	if(Electron_numberOfMissedHits(i)>0) return false;
	if(Electron_HasMatchedConversions(i)) return false;
	if(dxy(Electron_p4(i),Electron_Poca(i),Vtx(j))>=impact_xy) return false;
	if(dz(Electron_p4(i),Electron_Poca(i),Vtx(j))>=impact_z) return false;
	if(mvapt<20.){
		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=0.925) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=0.915) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.965) return false;
	}
	if(mvapt>=20.){
		if(mvaeta<0.8 && Electron_MVA_NonTrig_discriminator(i)<=0.905) return false;
		else if(mvaeta>=0.8 && mvaeta<1.479 && Electron_MVA_NonTrig_discriminator(i)<=0.955) return false;
		else if(mvaeta>=1.479 && mvaeta<2.5 && Electron_MVA_NonTrig_discriminator(i)<=0.975) return false;
	}
	return true;
}

bool Ntuple_Controller::isGoodJet(unsigned int i){
  //  Top Dilepton Jet selection with pt 15GeV
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel 
  //  isGoodJet_nooverlapremoval(i) with:
  //  deltaR jet-electron cleaning < 0.4 (2 selected lepton only) < 0.3 (e+jets only) 0.3
  //  deltaR jet-muon cleaning < 0.4(2 selected lepton only) < 0.3 (mu+jets only), 0.1 for PF and JET 0.3
  if(isGoodJet_nooverlapremoval(i)){
    unsigned int muon_idx;
    return !jethasMuonOverlap(i,muon_idx);
  }
  return false;
}

bool Ntuple_Controller::jethasMuonOverlap(unsigned int jet_idx,unsigned int &muon_idx){
  for(unsigned int j=0;j<NMuons();j++){
    if(isGoodMuon_nooverlapremoval(j) && Muon_RelIso(j)<0.2){
      if(Tools::dr(Muon_p4(j),PFJet_p4(jet_idx))<0.4){ muon_idx=j;return true;}
    }
  }
  return false;
}


bool Ntuple_Controller::isGoodJet_nooverlapremoval(unsigned int i){
  //  Top Dilepton Jet selection with pt 15GeV
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel
  //  Jet ID defined in isJetID(i)               
  //  cut         dilepton    l+jet
  //  corrected pT 30 GeV 30 GeV    
  //  residual correction (data) applied applied
  //  abs(eta) < 2.5 < 2.4                      
  //  jet ID applied applied     
  if(isJetID(i)){
    if(PFJet_p4(i).Pt()>15.0){
      if(fabs(PFJet_p4(i).Eta())<2.4){
	return true;
      }
    }
  }
  return false;
}

bool Ntuple_Controller::isJetID(unsigned int i){
  //  Top Dilepton Jet selection with pt and iso matching the muon and tau.
  //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel  
  //  Jet ID :
  //  number of constituents>1 (patJet->numberOfDaughters())
  //  CEF<0.99 (patJet->chargedEmEnergyFraction())
  //  NHF<0.99 (patJet->neutralHadronEnergyFraction())
  //  NEF<0.99 (patJet->neutralEmEnergyFraction())
  //  if |η|<2.4, CHF>0 (patJet->chargedHadronEnergyFraction())
  //  if |η|<2.4, NCH>0 (patJet->chargedMultiplicity()) 
  /////////////////////////////////////////////////////////////////////////
  // apply jet ID
  bool JetID_ok=false;
  if(PFJet_numberOfDaughters(i)>1){
    if(PFJet_chargedEmEnergyFraction(i)<0.99){
      if(PFJet_neutralHadronEnergyFraction(i)<0.99){
	if(PFJet_neutralEmEnergyFraction(i)<0.99){
	  if(fabs(PFJet_p4(i).Eta())<2.4){
	    if(PFJet_chargedHadronEnergyFraction(i)>0){
	      if(PFJet_chargedMultiplicity(i)>0){
		return true;
	      }
	    }
	  }
	  else{
	    return true;
	  }
	}
      }
    }
  }
  return false;
}



double Ntuple_Controller::TauSpinerGet(int SpinType){
#ifdef USE_TauSpinner
  if(!isData()){
    TauDecay taudecay;
    std::vector<SimpleParticle> tau_daughters, tau_daughters2;
    SimpleParticle tau, tau2;
    for(int i=0; i<NMCSignalParticles();i++){
      // check for signal Boson
      if(MCSignalParticle_pdgid(i)==25 || 
	 MCSignalParticle_pdgid(i)==36 || 
	 MCSignalParticle_pdgid(i)==22 || 
	 MCSignalParticle_pdgid(i)==23){
	if(MCSignalParticle_Tauidx(i).size()==2){
	  SimpleParticle X(MCSignalParticle_p4(i).Px(),MCSignalParticle_p4(i).Py(),MCSignalParticle_p4(i).Pz(),
			   MCSignalParticle_p4(i).E(),MCSignalParticle_pdgid(i));
	  bool tau1good(false),tau2good(false);
	  //first tau
	  unsigned int tauidx=MCSignalParticle_Tauidx(i).at(0);
	  if(verbose)std::cout  << "tau 1 indx " << tauidx  << " Number of Tau and Products " << NMCTauDecayProducts(tauidx) << std::endl;
	  if(tauidx<NMCTaus()){
	    for(int t=0;t<NMCTauDecayProducts(tauidx);t++){
	      int mypdgid=abs((int)MCTauandProd_pdgid(tauidx,t));
	      if(abs( mypdgid)==abs(PDGInfo::tau_plus) ){
		tau=SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),MCTauandProd_p4(tauidx,t).Py(),MCTauandProd_p4(tauidx,t).Pz(),
				   MCTauandProd_p4(tauidx,t).E(),MCTauandProd_pdgid(tauidx,t));
	      }
	      else if((taudecay.isTauFinalStateParticle(mypdgid) && mypdgid!=PDGInfo::gamma)){
		tau_daughters.push_back(SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),MCTauandProd_p4(tauidx,t).Py(),
						       MCTauandProd_p4(tauidx,t).Pz(),MCTauandProd_p4(tauidx,t).E(),
						       MCTauandProd_pdgid(tauidx,t)));
		tau1good=true;
	      }
	    }
	  }
	  // second tau
	  tauidx=MCSignalParticle_Tauidx(i).at(1);
	  if(verbose)std::cout  << "tau 2 indx " << tauidx  << " Number of Tau and Products " << NMCTauDecayProducts(tauidx) << std::endl;
	  if(tauidx<NMCTaus()){
	    for(int t=0;t<NMCTauDecayProducts(tauidx);t++){
	      int mypdgid=abs((int)MCTauandProd_pdgid(tauidx,t));
	      if(abs( mypdgid)==abs(PDGInfo::tau_plus)){
		    tau2=SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
					MCTauandProd_p4(tauidx,t).Py(),
					MCTauandProd_p4(tauidx,t).Pz(),
					MCTauandProd_p4(tauidx,t).E(),
					MCTauandProd_pdgid(tauidx,t));
		    tau2good=true;
	      }
	      if((taudecay.isTauFinalStateParticle(mypdgid) && mypdgid!=PDGInfo::gamma)){
		if(verbose) std::cout << "isDaughter" << std::endl;
		if(tau_daughters.size()>0)tau2good=true;
		tau_daughters2.push_back(SimpleParticle(MCTauandProd_p4(tauidx,t).Px(),
							MCTauandProd_p4(tauidx,t).Py(),
							MCTauandProd_p4(tauidx,t).Pz(),
							MCTauandProd_p4(tauidx,t).E(),
							MCTauandProd_pdgid(tauidx,t)));
	      }
	    }
	  }
	  if(tau1good && tau2good){
	    if(verbose)std::cout  << "Two Taus found: " << tau_daughters.size() << " " << tau_daughters2.size() << std::endl;
	    return TauSpinerInt.Get(SpinType,X,tau,tau_daughters,tau2,tau_daughters2);
	  }
	}
      }
    }
  }
#endif
  return 1.0;
}





bool Ntuple_Controller::hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &tau_idx){
  for(int i=0; i<NMCSignalParticles();i++){
    if(MCSignalParticle_pdgid(i)==parent_pdgid){
      for(int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	if(MCSignalParticle_Tauidx(i).at(j)>=NMCTaus()){
	  std::cout << "Warning INVALID Tau index... Skipping event! MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << std::endl;
	  return false;
	}
      }
      for(int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
	unsigned int tauidx=MCSignalParticle_Tauidx(i).at(j);
	if(verbose)std::cout << "MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << " " << Ntp->MCTau_JAK->size() << " " << Ntp->MCTauandProd_pdgid->size() << std::endl;
	if(MCTau_JAK(tauidx)==tau_jak){ tau_idx=tauidx;Boson_idx=i;return true;}
      }
    }
  }
  return false;
}

bool Ntuple_Controller::hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,unsigned int &tau1_idx, unsigned int &tau2_idx){
  for(int i=0; i<NMCSignalParticles();i++){
    if(MCSignalParticle_pdgid(i)==parent_pdgid){
      for(int j=0; j<MCSignalParticle_Tauidx(i).size();j++){
        if(MCSignalParticle_Tauidx(i).at(j)>=NMCTaus()){
	  std::cout << "Warning INVALID Tau index... Skipping event! MCSignalParticle_Tauidx: " << MCSignalParticle_Tauidx(i).at(j) << " Number of MC Taus: " << NMCTaus() << std::endl;
          return false;
        }
      }
      if(MCSignalParticle_Tauidx(i).size()==2){
	tau1_idx=MCSignalParticle_Tauidx(i).at(0);
        tau2_idx=MCSignalParticle_Tauidx(i).at(1);
	Boson_idx=i;
	return true;
      }
    }
  }
  return false;
}




bool Ntuple_Controller::TriggerAccept(TString n){
  unsigned int i=0;
  if(GetTriggerIndex(n,i))return TriggerAccept(i);
  return false;
}

unsigned int Ntuple_Controller::HLTPrescale(TString n){
  unsigned int i=0;
  if(GetTriggerIndex(n,i))return HLTPrescale(i);
  return 1;
}

unsigned int Ntuple_Controller::L1SEEDPrescale(TString n){
  unsigned int i=0;
  if(GetTriggerIndex(n,i))return L1SEEDPrescale(i);
  return 1;
}

bool Ntuple_Controller::GetTriggerIndex(TString n, unsigned int &i){
  for(i=0; i<Ntp->HTLTriggerName->size();i++){
    TString name=HTLTriggerName(i);
    if(name.Contains(n))return true;
  }
  return false;
}


TMatrixTSym<double> Ntuple_Controller::PFTau_TIP_primaryVertex_cov(unsigned int i){
  TMatrixTSym<double> V_cov(LorentzVectorParticle::NVertex);
  int l=0;
  for(unsigned int j=0;j<LorentzVectorParticle::NVertex;j++){
    for(unsigned int k=j;k<LorentzVectorParticle::NVertex;k++){
      //if(j==k) V_cov(i,j)=pow(0.0001,2.0);
      V_cov(j,k)=Ntp->PFTau_TIP_primaryVertex_cov->at(i).at(l);
      V_cov(k,j)=Ntp->PFTau_TIP_primaryVertex_cov->at(i).at(l);
      l++;
    }
  }
  return  V_cov;
}

TMatrixTSym<double> Ntuple_Controller::PFTau_TIP_secondaryVertex_cov(unsigned int i){
  TMatrixTSym<double> V_cov(LorentzVectorParticle::NVertex);
  int l=0;
  for(unsigned int j=0;j<LorentzVectorParticle::NVertex;j++){
    for(unsigned int k=j;k<LorentzVectorParticle::NVertex;k++){
      V_cov(j,k)=Ntp->PFTau_TIP_secondaryVertex_cov->at(i).at(l);
      V_cov(k,j)=Ntp->PFTau_TIP_secondaryVertex_cov->at(i).at(l);
      l++;
    }
  }
  return  V_cov;
}

LorentzVectorParticle Ntuple_Controller::PFTau_a1_lvp(unsigned int i){
  TMatrixT<double>    a1_par(LorentzVectorParticle::NLorentzandVertexPar,1);
  TMatrixTSym<double> a1_cov(LorentzVectorParticle::NLorentzandVertexPar);
  int l=0;
  if(Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar){
    for(int k=0; k<LorentzVectorParticle::NLorentzandVertexPar; k++){
      a1_par(k,0)=Ntp->PFTau_a1_lvp->at(i).at(k);
      for(int j=k; j<LorentzVectorParticle::NLorentzandVertexPar; j++){
	a1_cov(k,j)=Ntp->PFTau_a1_cov->at(i).at(l);
	a1_cov(j,k)=Ntp->PFTau_a1_cov->at(i).at(l);
	l++;
      } 
    }
  }
  return LorentzVectorParticle(a1_par,a1_cov,Ntp->PFTau_a1_pdgid->at(i).at(0),Ntp->PFTau_a1_charge->at(i).at(0),Ntp->PFTau_a1_B->at(i).at(0));
}

std::vector<TrackParticle> Ntuple_Controller::PFTau_daughterTracks(unsigned int i){
  std::vector<TrackParticle> daughter;
  for(unsigned int d=0;d<Ntp->PFTau_daughterTracks_poca->at(i).size();d++){
    TMatrixT<double>    a1_par(TrackParticle::NHelixPar,1);
    TMatrixTSym<double> a1_cov(TrackParticle::NHelixPar);
    int l=0;
    for(int k=0; k<TrackParticle::NHelixPar; k++){
      a1_par(k,0)=Ntp->PFTau_daughterTracks->at(i).at(d).at(k);
      for(int j=k; j<TrackParticle::NHelixPar; j++){
	a1_cov(k,j)=Ntp->PFTau_daughterTracks->at(i).at(d).at(l);
	a1_cov(j,k)=Ntp->PFTau_daughterTracks->at(i).at(d).at(l);
	l++;
      }
    }
    daughter.push_back(TrackParticle(a1_par,a1_cov,Ntp->PFTau_daughterTracks_pdgid->at(i).at(d),Ntp->PFTau_daughterTracks_M->at(i).at(d),Ntp->PFTau_daughterTracks_charge->at(i).at(d),Ntp->PFTau_daughterTracks_B->at(i).at(d)));
  }
  return daughter;
}

std::vector<TVector3> Ntuple_Controller::PFTau_daughterTracks_poca(unsigned int i){
  std::vector<TVector3> poca;
  for(unsigned int k=0;k<Ntp->PFTau_daughterTracks_poca->at(i).size();k++){
    poca.push_back(TVector3(Ntp->PFTau_daughterTracks_poca->at(i).at(k).at(0),Ntp->PFTau_daughterTracks_poca->at(i).at(k).at(1),Ntp->PFTau_daughterTracks_poca->at(i).at(k).at(2)));
  }
  return poca;
}
   

TMatrixTSym<double> Ntuple_Controller::PF_Tau_FlightLegth3d_TauFrame_cov(unsigned int i){
  TVector3 f=PFTau_FlightLength3d(i);
  TMatrixT<double> Res(5,1);
  Res(0,0)=f.X();
  Res(1,0)=f.Y();
  Res(2,0)=f.Z();
  Res(3,0)=f.Phi();
  Res(4,0)=f.Theta();
  TMatrixTSym<double> ResCov(5);
  TMatrixTSym<double> cov=PFTau_FlightLength3d_cov(i);
  for(unsigned int s=0;s<LorentzVectorParticle::NVertex;s++){
    for(unsigned int t=0;t<LorentzVectorParticle::NVertex;t++){
      ResCov(s,t)=cov(s,t);
    }
  }
  TMatrixT<double> Resp=MultiProngTauSolver::RotateToTauFrame(Res);
  TMatrixTSym<double> RespCov=ErrorMatrixPropagator::PropogateError(&MultiProngTauSolver::RotateToTauFrame,Res,ResCov);
  for(unsigned int s=0;s<LorentzVectorParticle::NVertex;s++){
    for(unsigned int t=0;t<LorentzVectorParticle::NVertex;t++){
      cov(s,t)=RespCov(s,t);
    }
  }
  return cov;
}

TVector3 Ntuple_Controller::PF_Tau_FlightLegth3d_TauFrame(unsigned int i){
  TVector3 f=PFTau_FlightLength3d(i);
  TMatrixT<double> Res(5,1);
  Res(0,0)=f.X();
  Res(1,0)=f.Y();
  Res(2,0)=f.Z();
  Res(3,0)=f.Phi();
  Res(4,0)=f.Theta();
  TMatrixT<double> Resp=MultiProngTauSolver::RotateToTauFrame(Res);
  return TVector3(Resp(0,0),Resp(1,0),Resp(2,0));
}

float Ntuple_Controller::dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs((-(poca.X()-vtx.X())*fourvector.Py()+(poca.Y()-vtx.Y())*fourvector.Px())/fourvector.Pt());
}

float Ntuple_Controller::dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx){
	return fabs(poca.Z()-vtx.Z()-((poca.X()-vtx.X())*fourvector.Px()+(poca.Y()-vtx.Y())*fourvector.Py())*fourvector.Pz()/pow(fourvector.Pt(),2));
}

float Ntuple_Controller::vertexSignificance(TVector3 vec, unsigned int vertex){
	if(vertex>=0 && vertex<NVtx()){
		const float elm[3] = {(vec.X()-Vtx(vertex).X()),(vec.Y()-Vtx(vertex).Y()),(vec.Z()-Vtx(vertex).Z())};
		TVectorF diff(3,elm);
		TMatrixF M(Vtx_Cov(vertex));
		if(M.IsValid()){
			double mag = diff.Norm2Sqr();
			double sim = M.Similarity(diff);
			return mag/sqrt(sim);
		}
	}
	return 999;
}
