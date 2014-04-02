//Ntuple_Controller.h HEADER FILE

#ifndef Ntuple_Controller_h
#define Ntuple_Controller_h


// Root include files
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TVectorT.h"

// Include files (C & C++ libraries)
#include<iostream>
#include <vector>
#include <string.h>

#include "NtupleReader.h"

#include "HistoConfig.h"
#ifdef USE_TauSpinner
#include "TauSpinerInterface.h"
#endif
#include "HistoConfig.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"
#include "TauDataFormat/TauNtuple/interface/TauDecay.h"


#include "SimpleFits/FitSoftware/interface/TrackParticle.h"
#include "SimpleFits/FitSoftware/interface/LorentzVectorParticle.h"
#include "SimpleFits/FitSoftware/interface/MultiProngTauSolver.h"
#include "SimpleFits/FitSoftware/interface/ErrorMatrixPropagator.h"
#include "SimpleFits/FitSoftware/interface/TauA1NuConstrainedFitter.h"
///////////////////////////////////////////////////////////////////////////////
//*****************************************************************************
//*
//*   Class: Ntuple_Controller
//*   
//*   Purpose: The purpose of this class is to provide a interface to the
//*            Ntuple
//*
//*   Designed by: Ian Nugent
//*
//*
//*****************************************************************************
///////////////////////////////////////////////////////////////////////////////


class Ntuple_Controller{
 private:
  NtupleReader *Ntp;
  TFile *newfile;
  TTree *SkimmedTree;
  int nbytes;
  int jentry;
  int nb;
  bool copyTree;

  bool verbose;

  int currentEvent;

  // Ntuple Access Functions
  virtual void Branch_Setup(TString B_Name, int type);
  virtual void Branch_Setup(){}

  // Functions to configure objects
  virtual void ConfigureObjects(); 
  void doElectrons();
  void doPhotons();
  void doJets();
  void doMuons();
  void doTaus();
  void doMET();
  unsigned int ObjEvent;

  // Object Variables
  std::vector<TLorentzVector> electrons_default;
  std::vector<TLorentzVector> photons_default;
  std::vector<TLorentzVector> jets_default;
  std::vector<TLorentzVector> muons_default;
  std::vector<TLorentzVector> taus_default;
  TLorentzVector              met_default;
  std::vector<TLorentzVector> electrons;
  std::vector<TLorentzVector> photons;
  std::vector<TLorentzVector> jets;
  std::vector<TLorentzVector> muons;
  std::vector<TLorentzVector> taus;
  TLorentzVector              met;

  // Systematic controls variables
  int theSys;
  HistoConfig HConfig;

  // Interfaces
#ifdef USE_TauSpinner  
  TauSpinerInterface TauSpinerInt;
#endif
  HistoConfig HistoC;

  // Fit Variables
  LorentzVectorParticle               theTau;
  std::vector<LorentzVectorParticle>  daughter;
  double                              LC_chi2;
  double                              ndof;
  bool                                fitStatus;

 public:
  // Constructor
  Ntuple_Controller(std::vector<TString> RootFiles);

  // Destructor
  ~Ntuple_Controller() ;

  //TauSpiner function
  double TauSpinerGet(int SpinType);
  void TauSpinerSetSignal(int signalcharge){
#ifdef USE_TauSpinner
TauSpinerInt.SetTauSignalCharge(signalcharge);
#endif
}
   enum beamspot{BS_x0,BS_y0,BS_z0,BS_sigmaZ,BS_dxdz,BS_dydz,BS_BeamWidthX,NBS_par};
   enum TrackQuality {
     undefQuality = -1, loose = 0, tight = 1, highPurity = 2,
     confirmed = 3, goodIterative = 4, looseSetWithPV = 5, highPuritySetWithPV = 6,
     qualitySize = 7
   };
  enum TrackPar{i_qoverp = 0, i_lambda, i_phi, i_dxy,i_dsz};

  // Ntuple Access Functions 
  virtual Int_t Get_Entries();
  virtual void Get_Event(int _jentry);
  virtual Int_t Get_EventIndex();
  virtual TString Get_File_Name();

  //Ntuple Cloning Functions
  virtual void CloneTree(TString n);
  virtual void SaveCloneTree();
  inline void AddEventToCloneTree(){if(copyTree)SkimmedTree->Fill();}

  // Systematic controls
  enum    Systematic {Default=0,NSystematics};

  int     SetupSystematics(TString sys_);
  void    SetSysID(int sysid){theSys=sysid;}


  // Data/MC switch and thin
  bool isData(){return Ntp->Event_isRealData;}
  void ThinTree();

  // Physics Varible Get Functions
  // Event Varibles
  int GetMCID();
  unsigned int RunNumber(){return Ntp->Event_RunNumber;}
  unsigned int EventNumber(){ return Ntp->Event_EventNumber;}
  int BunchCrossing(){ return Ntp->Event_bunchCrossing;}
  int OrbitNumber(){ return Ntp->Event_orbitNumber;}
  unsigned int LuminosityBlock(){return Ntp->Event_luminosityBlock;}
  int           PileupInfo_TrueNumInteractions_nm1(){return Ntp->PileupInfo_TrueNumInteractions_nm1;}
  int           PileupInfo_TrueNumInteractions_n0(){return Ntp->PileupInfo_TrueNumInteractions_n0;}
  int           PileupInfo_TrueNumInteractions_np1(){return Ntp->PileupInfo_TrueNumInteractions_np1;}
  double        PUWeight(){return Ntp->PUWeight;}
  double        PUWeight_p5(){return Ntp->PUWeight_p5;}
  double        PUWeight3D_m5(){return Ntp->PUWeight3D_m5;}
  double        PUWeight3D(){return Ntp->PUWeight3D;}
  double        PUWeight3D_p5(){return Ntp->PUWeight3D_p5;}
  double        PUWeight_m5(){return Ntp->PUWeight_m5;}

  TVectorT<double>      beamspot_par(){TVectorT<double> BS(NBS_par);for(unsigned int i=0;i<NBS_par;i++)BS(i)=Ntp->beamspot_par->at(i);return BS;}

  TMatrixTSym<double>   beamspot_cov(){
    TMatrixTSym<double> BS_cov(NBS_par);
    unsigned int l=0;
    for(unsigned int i=0;i<NBS_par;i++){
      for(unsigned int j=i;j<NBS_par;j++){
	BS_cov(i,j)=Ntp->beamspot_cov->at(l);
      }
    }
    return BS_cov;
  }
  
  float  beamspot_emittanceX(){return Ntp->beamspot_emittanceX;}
  float  beamspot_emittanceY(){return Ntp->beamspot_emittanceY;}
  float  beamspot_betaStar(){return Ntp->beamspot_betaStar;}

  // Vertex Information
  unsigned int NVtx(){return Ntp->Vtx_ndof->size();}
  TVector3     Vtx(unsigned int i){return TVector3(Ntp->Vtx_x->at(i),Ntp->Vtx_y->at(i),Ntp->Vtx_z->at(i));}
  float        Vtx_chi2(unsigned int i){return Ntp->Vtx_chi2->at(i);}
  float        Vtx_nTrk(unsigned int i){return Ntp->Vtx_nTrk->at(i);}
  float        Vtx_ndof(unsigned int i){return Ntp->Vtx_ndof->at(i);}
  TMatrixF     Vtx_Cov(unsigned int i);
  std::vector<int>  Vtx_Track_idx(unsigned int i){return Ntp->Vtx_Track_idx->at(i);}
  float Vtx_isFake(unsigned int i){return Ntp->Vtx_isFake->at(i);}
  TLorentzVector Vtx_TracksP4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->Vtx_TracksP4->at(i).at(j).at(1),Ntp->Vtx_TracksP4->at(i).at(j).at(2),Ntp->Vtx_TracksP4->at(i).at(j).at(3),Ntp->Vtx_TracksP4->at(i).at(j).at(0));}

  bool isVtxGood(unsigned int i);
  bool isGoodVtx(unsigned int i);

  // Muon information
  unsigned int   NMuons(){return Ntp->Muon_p4->size();}
  TLorentzVector Muon_p4(unsigned int i){return TLorentzVector(Ntp->Muon_p4->at(i).at(1),Ntp->Muon_p4->at(i).at(2),Ntp->Muon_p4->at(i).at(3),Ntp->Muon_p4->at(i).at(0));}
  TVector3       Muon_Poca(unsigned int i){return TVector3(Ntp->Muon_Poca->at(i).at(0),Ntp->Muon_Poca->at(i).at(1),Ntp->Muon_Poca->at(i).at(2));}
  bool           Muon_isGlobalMuon(unsigned int i){return Ntp->Muon_isGlobalMuon->at(i);}
  bool           Muon_isStandAloneMuon(unsigned int i){return Ntp->Muon_isStandAloneMuon->at(i);}
  bool           Muon_isTrackerMuon(unsigned int i){return Ntp->Muon_isTrackerMuon->at(i);}
  bool           Muon_isCaloMuon(unsigned int i){return Ntp->Muon_isCaloMuon->at(i);}
  bool           Muon_isIsolationValid(unsigned int i){return Ntp->Muon_isIsolationValid->at(i);}
  bool           Muon_isQualityValid(unsigned int i){return Ntp->Muon_isQualityValid->at(i);}
  bool           Muon_isTimeValid(unsigned int i){return Ntp->Muon_isTimeValid->at(i);}
  float          Muon_emEt03(unsigned int i){return Ntp->Muon_emEt03->at(i);}
  float          Muon_emVetoEt03(unsigned int i){return Ntp->Muon_emVetoEt03->at(i);}
  float          Muon_hadEt03(unsigned int i){return Ntp->Muon_hadEt03->at(i);}
  float          Muon_hadVetoEt03(unsigned int i){return Ntp->Muon_hadVetoEt03->at(i);}
  float          Muon_nJets03(unsigned int i){return Ntp->Muon_nJets03->at(i);}
  float          Muon_nTracks03(unsigned int i){return Ntp->Muon_nTracks03->at(i);}
  float          Muon_sumPt03(unsigned int i){return Ntp->Muon_sumPt03->at(i);}
  float          Muon_trackerVetoPt03(unsigned int i){return Ntp->Muon_trackerVetoPt03->at(i);}
  float          Muon_emEt05(unsigned int i){return Ntp->Muon_emEt05->at(i);}
  float          Muon_emVetoEt05(unsigned int i){return Ntp->Muon_emVetoEt05->at(i);}
  float          Muon_hadEt05(unsigned int i){return Ntp->Muon_hadEt05->at(i);}
  float          Muon_hadVetoEt05(unsigned int i){return Ntp->Muon_hadVetoEt05->at(i);}
  float          Muon_nJets05(unsigned int i){return Ntp->Muon_nJets05->at(i);}
  float          Muon_nTracks05(unsigned int i){return Ntp->Muon_nTracks05->at(i);}
  float          Muon_sumPt05(unsigned int i){return Ntp->Muon_sumPt05->at(i);}
  float          Muon_trackerVetoPt05(unsigned int i){return Ntp->Muon_trackerVetoPt05->at(i);}
  unsigned int   Muon_Track_idx(unsigned int i){return Ntp->Muon_Track_idx->at(i);}
  float          Muon_hitPattern_pixelLayerwithMeas(unsigned int i){return Ntp->Muon_hitPattern_pixelLayerwithMeas->at(i);}
  float          Muon_numberOfMatchedStations(unsigned int i){return Ntp->Muon_numberOfMatchedStations->at(i);}
  float          Muon_normChi2(unsigned int i){return Ntp->Muon_normChi2->at(i);}
  float          Muon_hitPattern_numberOfValidMuonHits(unsigned int i){return Ntp->Muon_hitPattern_numberOfValidMuonHits->at(i);}
  float          Muon_innerTrack_numberofValidHits(unsigned int i){return Ntp->Muon_innerTrack_numberofValidHits->at(i);}
  float          Muon_numberOfMatches(unsigned int i){return Ntp->Muon_numberOfMatches->at(i);}
  int            Muon_numberOfChambers(unsigned int i){return Ntp->Muon_numberOfChambers->at(i);}
  int            Muon_Charge(unsigned int i){return Ntp->Muon_charge->at(i);}

  bool           Muon_isPFMuon(unsigned int i){return Ntp->Muon_isPFMuon->at(i);}                                                     
  float          Muon_sumChargedHadronPt03(unsigned int i){return Ntp->Muon_sumChargedHadronPt03->at(i);}                             // sum-pt of charged Hadron					                               
  float          Muon_sumChargedParticlePt03(unsigned int i){return Ntp->Muon_sumChargedParticlePt03->at(i);}			      // sum-pt of charged Particles(inludes e/mu)			    
  float          Muon_sumNeutralHadronEt03(unsigned int i){return Ntp->Muon_sumNeutralHadronEt03->at(i);}			      // sum pt of neutral hadrons					    
  float          Muon_sumNeutralHadronEtHighThreshold03(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold03->at(i);}   // sum pt of neutral hadrons with a higher threshold		    
  float          Muon_sumPhotonEt03(unsigned int i){return Ntp->Muon_sumPhotonEt03->at(i);}					      // sum pt of PF photons					    
  float          Muon_sumPhotonEtHighThreshold03(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold03->at(i);}		      // sum pt of PF photons with a higher threshold		    
  float          Muon_sumPUPt03(unsigned int i){return Ntp->Muon_sumPUPt03->at(i);}						      // sum pt of charged Particles not from PV (for Pu corrections)  
																      								    
  float          Muon_sumChargedHadronPt04(unsigned int i){return Ntp->Muon_sumChargedHadronPt04->at(i);}			      // sum-pt of charged Hadron					    
  float          Muon_sumChargedParticlePt04(unsigned int i){return Ntp->Muon_sumChargedParticlePt04->at(i);}			      // sum-pt of charged Particles(inludes e/mu)			    
  float          Muon_sumNeutralHadronEt04(unsigned int i){return Ntp->Muon_sumNeutralHadronEt04->at(i);}			      // sum pt of neutral hadrons					    
  float          Muon_sumNeutralHadronEtHighThreshold04(unsigned int i){return Ntp->Muon_sumNeutralHadronEtHighThreshold04->at(i);}   // sum pt of neutral hadrons with a higher threshold		    
  float          Muon_sumPhotonEt04(unsigned int i){return Ntp->Muon_sumPhotonEt04->at(i);}					      // sum pt of PF photons					    
  float          Muon_sumPhotonEtHighThreshold04(unsigned int i){return Ntp->Muon_sumPhotonEtHighThreshold04->at(i);}		      // sum pt of PF photons with a higher threshold		    
  float          Muon_sumPUPt04(unsigned int i){return Ntp->Muon_sumPUPt04->at(i);}						      // sum pt of charged Particles not from PV (for Pu corrections)  

  int            Muon_numberofValidPixelHits(unsigned int i){return Ntp->Muon_numberofValidPixelHits->at(i);}
  int            Muon_trackerLayersWithMeasurement(unsigned int i){return Ntp->Muon_trackerLayersWithMeasurement->at(i);}


  TrackParticle Muon_TrackParticle(unsigned int i){
    TMatrixT<double>    mu_par(TrackParticle::NHelixPar,1);
    TMatrixTSym<double> mu_cov(TrackParticle::NHelixPar);
    unsigned int l=0;
    for(int k=0; k<TrackParticle::NHelixPar; k++){
      mu_par(k,0)=Ntp->Muon_par->at(i).at(k);
      for(int j=k; j<TrackParticle::NHelixPar; j++){
        mu_cov(k,j)=Ntp->Muon_cov->at(i).at(l);
        l++;
      }
    }
    return TrackParticle(mu_par,mu_cov,Ntp->Muon_pdgid->at(i),Ntp->Muon_M->at(i),Ntp->Muon_charge->at(i),Ntp->Muon_B->at(i));
  }

  bool           isGoodMuon(unsigned int i);
  bool           isGoodMuon_nooverlapremoval(unsigned int i);

  bool			 isTightMuon(unsigned int i);
  bool			 isTightMuon(unsigned int i, unsigned int j);
  bool           isSelectedMuon(unsigned int i, unsigned int j, double impact_xy, double impact_z);
  bool			 isLooseMuon(unsigned int i);
  float          Muon_RelIso(unsigned int i);


  //Base Tau Information (PF)
   unsigned int      NPFTaus(){return Ntp->PFTau_p4->size();}
   TLorentzVector    PFTau_p4(unsigned int i){return TLorentzVector(Ntp->PFTau_p4->at(i).at(1),Ntp->PFTau_p4->at(i).at(2),Ntp->PFTau_p4->at(i).at(3),Ntp->PFTau_p4->at(i).at(0));}
   TVector3          PFTau_Poca(unsigned int i){return TVector3(Ntp->PFTau_Poca->at(i).at(0),Ntp->PFTau_Poca->at(i).at(1),Ntp->PFTau_Poca->at(i).at(2));}
   bool PFTau_isTightIsolation(unsigned int i){return Ntp->PFTau_isTightIsolation->at(i);}
   bool PFTau_isMediumIsolation(unsigned int i){return  Ntp->PFTau_isMediumIsolation->at(i);}
   bool PFTau_isLooseIsolation(unsigned int i){return  Ntp->PFTau_isLooseIsolation->at(i);}
   bool PFTau_isTightIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isTightIsolationDBSumPtCorr->at(i);}
   bool PFTau_isMediumIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isMediumIsolationDBSumPtCorr->at(i);}
   bool PFTau_isLooseIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isLooseIsolationDBSumPtCorr->at(i);}
   bool PFTau_isVLooseIsolationDBSumPtCorr(unsigned int i){return  Ntp->PFTau_isVLooseIsolationDBSumPtCorr->at(i);}
   bool PFTau_isHPSAgainstElectronsLoose(unsigned int i){return  Ntp->PFTau_isHPSAgainstElectronsLoose->at(i);}
   bool PFTau_isHPSAgainstElectronsMedium(unsigned int i){return  Ntp->PFTau_isHPSAgainstElectronsMedium->at(i);}
   bool PFTau_isHPSAgainstElectronsTight(unsigned int i){return  Ntp->PFTau_isHPSAgainstElectronsTight->at(i);}
   bool PFTau_isHPSAgainstMuonLoose(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonLoose->at(i);}
   bool PFTau_isHPSAgainstMuonMedium(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonMedium->at(i);}
   bool PFTau_isHPSAgainstMuonTight(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonTight->at(i);}
   bool PFTau_isHPSAgainstMuonLoose2(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonLoose2->at(i);}
   bool PFTau_isHPSAgainstMuonMedium2(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonMedium2->at(i);}
   bool PFTau_isHPSAgainstMuonTight2(unsigned int i){return  Ntp->PFTau_isHPSAgainstMuonTight2->at(i);}
   bool PFTau_isHPSByDecayModeFinding(unsigned int i){return  Ntp->PFTau_isHPSByDecayModeFinding->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3LooseElectronRejection->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3MediumElectronRejection->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3TightElectronRejection->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMVA3VTightElectronRejection->at(i);}
   bool PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits->at(i);}
   bool PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits->at(i);}
   float PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(unsigned int i){return Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits->at(i);}
   bool PFTau_HPSPFTauDiscriminationByLooseIsolationMVA(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMediumIsolationMVA(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA->at(i);}
   bool PFTau_HPSPFTauDiscriminationByTightIsolationMVA(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA->at(i);}
   bool PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByLooseIsolationMVA2->at(i);}
   bool PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByMediumIsolationMVA2->at(i);}
   bool PFTau_HPSPFTauDiscriminationByTightIsolationMVA2(unsigned int i){return  Ntp->PFTau_HPSPFTauDiscriminationByTightIsolationMVA2->at(i);}
   int PFTau_hpsDecayMode(unsigned int i){return  Ntp->PFTau_hpsDecayMode->at(i);}
   int PFTau_Charge(unsigned int i){return  Ntp->PFTau_Charge->at(i);}
   std::vector<int> PFTau_Track_idx(unsigned int i){return  Ntp->PFTau_Track_idx->at(i);}
   TVector3 PFTau_TIP_primaryVertex_pos(unsigned int i){return  TVector3(Ntp->PFTau_TIP_primaryVertex_pos->at(i).at(0),Ntp->PFTau_TIP_primaryVertex_pos->at(i).at(1),Ntp->PFTau_TIP_primaryVertex_pos->at(i).at(2));}
   TMatrixTSym<double> PFTau_TIP_primaryVertex_cov(unsigned int i);
   bool PFTau_TIP_hassecondaryVertex(unsigned int i){if(Ntp->PFTau_TIP_secondaryVertex_pos->at(i).size()==3)return true; return false;}
   TVector3 PFTau_TIP_secondaryVertex_pos(unsigned int i){return  TVector3(Ntp->PFTau_TIP_secondaryVertex_pos->at(i).at(0),Ntp->PFTau_TIP_secondaryVertex_pos->at(i).at(1),Ntp->PFTau_TIP_secondaryVertex_pos->at(i).at(2));}
   TMatrixTSym<double> PFTau_TIP_secondaryVertex_cov(unsigned int i);
   float PFTau_TIP_secondaryVertex_vtxchi2(unsigned int i){if(Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).size()==1) return  Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).at(0); return 0;}
   float PFTau_TIP_secondaryVertex_vtxndof(unsigned int i){if(Ntp->PFTau_TIP_secondaryVertex_vtxndof->at(i).size()==1) return  Ntp->PFTau_TIP_secondaryVertex_vtxndof->at(i).at(0);  return 0;}
   bool PFTau_TIP_hasA1Momentum(unsigned int i){if(Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar)return true; return false;}
   LorentzVectorParticle PFTau_a1_lvp(unsigned int i);
   TLorentzVector PFTau_3PS_A1_LV(unsigned int i){return PFTau_a1_lvp(i).LV();}
   std::vector<TrackParticle> PFTau_daughterTracks(unsigned int i);
   std::vector<TVector3> PFTau_daughterTracks_poca(unsigned int i);   
   TMatrixTSym<double> PFTau_FlightLength3d_cov(unsigned int i){return  PFTau_TIP_secondaryVertex_cov(i)+PFTau_TIP_primaryVertex_cov(i);}
   TVector3 PFTau_FlightLength3d(unsigned int i){return PFTau_TIP_secondaryVertex_pos(i)-PFTau_TIP_primaryVertex_pos(i);}
   TMatrixTSym<double> PF_Tau_FlightLegth3d_TauFrame_cov(unsigned int i);
   TVector3 PF_Tau_FlightLegth3d_TauFrame(unsigned int i);
   double   PFTau_FlightLength_significance(unsigned int i){float e=PFTau_FlightLength_error(0); if(e>0) return PFTau_FlightLength(i)/e; return 0;}
   double   PFTau_FlightLength_error(unsigned int i){return PF_Tau_FlightLegth3d_TauFrame_cov(i)(LorentzVectorParticle::vz,LorentzVectorParticle::vz);}
   double   PFTau_FlightLength(unsigned int i){return PFTau_FlightLength3d(i).Mag();}
   
   bool ThreeProngTauFit(unsigned int i, unsigned int j,LorentzVectorParticle &theTau,std::vector<LorentzVectorParticle> &daughter,double &LC_chi2){
     if(Ntp->PFTau_TIP_secondaryVertex_vtxchi2->at(i).size()==1 && Ntp->PFTau_a1_lvp->at(i).size()==LorentzVectorParticle::NLorentzandVertexPar){
       // Tau Solver
       TVector3 pv=PFTau_TIP_primaryVertex_pos(i);
       TMatrixTSym<double> pvcov=PFTau_TIP_primaryVertex_cov(i);
       LorentzVectorParticle a1=PFTau_a1_lvp(i);
       TauA1NuConstrainedFitter TauA1NU(j,a1,pv,pvcov);
       bool fitStatus= TauA1NU.Fit();
       if(fitStatus){
	 theTau=TauA1NU.GetMother();
	 daughter=TauA1NU.GetReFitDaughters();
	 return  true;
       }
     }
     return false;
   }
   
   ////////////////////////////////////////////////
   // wrapper for backwards compatibility to KFit do not use in new code!!!
   int NKFTau(){return NPFTaus();}
   int KFTau_discriminatorByQC(unsigned int i,unsigned int j=0){
     TLorentzVector a1=PFTau_a1_lvp(i).LV();
     if(0.8<a1.M() && a1.M()<1.5)return true;
     return false;
   }

   bool KFTau_discriminatorByKFit(unsigned int i, unsigned int j=0){
     LorentzVectorParticle theTau;
     std::vector<LorentzVectorParticle> daughter;
     double LC_chi2;
     return ThreeProngTauFit(i,j,theTau,daughter,LC_chi2);
   }

   TLorentzVector KFTau_TauFit_p4(unsigned int i, unsigned int j=0){
     LorentzVectorParticle theTau;
     std::vector<LorentzVectorParticle> daughter;
     double LC_chi2;
     if(ThreeProngTauFit(i,j,theTau,daughter,LC_chi2))return theTau.LV();return TLorentzVector(0,0,0,0);
   }
   int KFTau_MatchedHPS_idx(unsigned int i){return i;}
   bool isGoodKFTau(unsigned int i,unsigned int j=0){
     LorentzVectorParticle theTau;
     std::vector<LorentzVectorParticle> daughter;
     double LC_chi2;
     if(ThreeProngTauFit(i,j,theTau,daughter,LC_chi2)){
       if(KFTau_discriminatorByQC(i,j)){
	 if(PFTau_hpsDecayMode(KFTau_MatchedHPS_idx(i)) == 10 && PFTau_isHPSByDecayModeFinding(i)){
	   if(PFTau_isVLooseIsolationDBSumPtCorr(KFTau_MatchedHPS_idx(i))){
	     return true;
	   }
	 }
       }
     }
     return false;
   }
   
   TLorentzVector   KFTau_TauVis_p4(unsigned int i,unsigned int j=0){return PFTau_a1_lvp(i).LV();}
   TLorentzVector   KFTau_Neutrino_p4(unsigned int i,unsigned int j=0){return (KFTau_TauFit_p4(i,j)-PFTau_a1_lvp(i).LV());}
   int      KFTau_nKinTaus(){return NKFTau();}
   int      KFTau_indexOfFitInfo(unsigned int i){return i;}
   TVector3 KFTau_Fit_TauPrimVtx(unsigned int i){return PFTau_TIP_primaryVertex_pos(i);}
   float    KFTau_Fit_ndf(unsigned int i,unsigned int j=0){return PFTau_TIP_secondaryVertex_vtxndof(i);}
   float    KFTau_Fit_chi2(unsigned int i,unsigned int j=0){
     LorentzVectorParticle theTau;
     std::vector<LorentzVectorParticle> daughter;
     double LC_chi2;
     if(ThreeProngTauFit(i,j,theTau,daughter,LC_chi2)){
       return PFTau_TIP_secondaryVertex_vtxchi2(i)+LC_chi2;
     }
     return 999;
   }
   float    KFTau_Fit_Chi2Prob(unsigned int i,unsigned int j=0){return TMath::Prob(KFTau_Fit_chi2(i,j),(int)KFTau_Fit_ndf(i,j));}
   int      KFTau_Fit_charge(unsigned int i){if(Ntp->PFTau_a1_charge->at(i).size()>0) return Ntp->PFTau_a1_charge->at(i).at(0); return 0;}
   int      KFTau_Fit_csum(unsigned int i,unsigned int j=0){return 0;}
   int      KFTau_Fit_iterations(unsigned int i, unsigned int j=0){return 0;}
   double   KFTau_RotatedVtx(unsigned int i, unsigned int j=0){return 0;}
   
   double   KFTau_Fit_RefitVisibleMass(unsigned int i, unsigned int j=0){return PFTau_a1_lvp(i).LV().M();}
   double   KFTau_Fit_SV_PV_significance(unsigned int i, unsigned int j=0){return PFTau_FlightLength_significance(i);}
   double   KFTau_Fit_SV_PV_error(unsigned int i, unsigned int j=0){return PFTau_FlightLength_error(i);}
   double   KFTau_Fit_SV_PV_distance(unsigned int i, unsigned int j=0){return PFTau_FlightLength(i);}

   /////////////////////////////////////////////////////////////////////////////////


      // Jet Information
   unsigned int       NPFJets(){return Ntp->PFJet_p4->size();}
   TLorentzVector     PFJet_p4(unsigned int i){return TLorentzVector(Ntp->PFJet_p4->at(i).at(1),Ntp->PFJet_p4->at(i).at(2),Ntp->PFJet_p4->at(i).at(3),Ntp->PFJet_p4->at(i).at(0));}
   float              PFJet_chargedEmEnergy(unsigned int i){return Ntp->PFJet_chargedEmEnergy->at(i);}
   float              PFJet_chargedHadronEnergy(unsigned int i){return Ntp->PFJet_chargedHadronEnergy->at(i);}
   float              PFJet_chargedHadronMultiplicity(unsigned int i){return Ntp->PFJet_chargedHadronMultiplicity->at(i);}
   float              PFJet_chargedMuEnergy(unsigned int i){return Ntp->PFJet_chargedMuEnergy->at(i);}
   float              PFJet_chargedMultiplicity(unsigned int i){return Ntp->PFJet_chargedMultiplicity->at(i);}
   float              PFJet_electronEnergy(unsigned int i){return Ntp->PFJet_electronEnergy->at(i);}
   float              PFJet_electronMultiplicity(unsigned int i){return Ntp->PFJet_electronMultiplicity->at(i);}
   float              PFJet_HFEMEnergy(unsigned int i){return Ntp->PFJet_HFEMEnergy->at(i);}
   float              PFJet_HFEMMultiplicity(unsigned int i){return Ntp->PFJet_HFEMMultiplicity->at(i);}
   float              PFJet_HFHadronEnergy(unsigned int i){return Ntp->PFJet_HFHadronEnergy->at(i);}
   float              PFJet_HFHadronMultiplicity(unsigned int i){return Ntp->PFJet_HFHadronMultiplicity->at(i);}
   float              PFJet_muonEnergy(unsigned int i){return Ntp->PFJet_muonEnergy->at(i);}
   float              PFJet_muonMultiplicity(unsigned int i){return Ntp->PFJet_muonMultiplicity->at(i);}
   float              PFJet_neutralEmEnergy(unsigned int i){return Ntp->PFJet_neutralEmEnergy->at(i);}
   float              PFJet_neutralHadronEnergy(unsigned int i){return Ntp->PFJet_neutralHadronEnergy->at(i);}
   float              PFJet_neutralHadronMultiplicity(unsigned int i){return Ntp->PFJet_neutralHadronMultiplicity->at(i);}
   float              PFJet_photonEnergy(unsigned int i){return Ntp->PFJet_photonEnergy->at(i);}
   float              PFJet_photonMultiplicity(unsigned int i){return Ntp->PFJet_photonMultiplicity->at(i);}
   float              PFJet_jetArea(unsigned int i){return Ntp->PFJet_jetArea->at(i);}
   float              PFJet_maxDistance(unsigned int i){return Ntp->PFJet_maxDistance->at(i);}
   int                PFJet_nConstituents(unsigned int i){return Ntp->PFJet_nConstituents->at(i);}
   float              PFJet_pileup(unsigned int i){return Ntp->PFJet_pileup->at(i);}
   float              PFJet_etaetaMoment(unsigned int i){return Ntp->PFJet_etaetaMoment->at(i);}
   float              PFJet_etaphiMoment(unsigned int i){return Ntp->PFJet_etaphiMoment->at(i);}
   std::vector<int>   PFJet_Track_idx(unsigned int i){return Ntp->PFJet_Track_idx->at(i);}
   int                PFJet_MatchedHPS_idx(unsigned int i){return Ntp->PFJet_MatchedHPS_idx->at(i);}
   int                PFJet_numberOfDaughters(unsigned int i){return Ntp->PFJet_numberOfDaughters->at(i);}
   float              PFJet_chargedEmEnergyFraction(unsigned int i){return Ntp->PFJet_chargedEmEnergyFraction->at(i);}
   float              PFJet_chargedHadronEnergyFraction(unsigned int i){return Ntp->PFJet_chargedHadronEnergyFraction->at(i);}
   float              PFJet_neutralHadronEnergyFraction(unsigned int i){return Ntp->PFJet_neutralHadronEnergyFraction->at(i);}
   float              PFJet_neutralEmEnergyFraction(unsigned int i){return Ntp->PFJet_neutralEmEnergyFraction->at(i);}
   bool               isGoodJet(unsigned int i);
   bool               isGoodJet_nooverlapremoval(unsigned int i);
   bool               isJetID(unsigned int i);
   float              PFJet_nTrk(unsigned int i){return Ntp->PFJet_nTrk->at(i);}
   TLorentzVector     PFJet_TracksP4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->PFJet_TracksP4->at(i).at(j).at(1),Ntp->PFJet_TracksP4->at(i).at(j).at(2),Ntp->PFJet_TracksP4->at(i).at(j).at(3),Ntp->PFJet_TracksP4->at(i).at(j).at(0));}
   int                PFJet_nTracks(unsigned int i){return Ntp->PFJet_TracksP4->at(i).size();}

   float              PFJet_PUJetID_discr(unsigned int i){return Ntp->PFJet_PUJetID_discr->at(i);}
   float              PFJet_PUJetID_looseWP(unsigned int i){return Ntp->PFJet_PUJetID_looseWP->at(i);}
   float              PFJet_PUJetID_mediumWP(unsigned int i){return Ntp->PFJet_PUJetID_mediumWP->at(i);}
   float              PFJet_PUJetID_tightWP(unsigned int i){return Ntp->PFJet_PUJetID_tightWP->at(i);}

   float              PFJet_partonFlavour(unsigned int i){return Ntp->PFJet_partonFlavour->at(i);}
   float              PFJet_bDiscriminator(unsigned int i){return Ntp->PFJet_bDiscriminator->at(i);}
   //float              PFJet_BTagWeight(unsigned int i){return Ntp->PFJet_BTagWeight->at(i);} // not implemented at the moment

   //MET information
   float              MET_Uncorr_et(){return Ntp->MET_Uncorr_et;}
   float              MET_Uncorr_phi(){return Ntp->MET_Uncorr_phi;}
   float              MET_Uncorr_ex(){return Ntp->MET_Uncorr_et*cos(Ntp->MET_Uncorr_phi);}
   float              MET_Uncorr_ey(){return Ntp->MET_Uncorr_et*sin(Ntp->MET_Uncorr_phi);}
   float              MET_Uncorr_significance(){return Ntp->MET_Uncorr_significance;}
   float              MET_Uncorr_significance_xx(){return Ntp->MET_Uncorr_significance_xx;}
   float              MET_Uncorr_significance_xy(){return Ntp->MET_Uncorr_significance_xy;}
   float              MET_Uncorr_significance_yy(){return Ntp->MET_Uncorr_significance_yy;}

   float              MET_CorrT0rt_et(){return Ntp->MET_CorrT0rt_et;}
   float              MET_CorrT0rt_phi(){return Ntp->MET_CorrT0rt_phi;}
   float              MET_CorrT0rt_ex(){return Ntp->MET_CorrT0rt_et*cos(Ntp->MET_CorrT0rt_phi);}
   float              MET_CorrT0rt_ey(){return Ntp->MET_CorrT0rt_et*sin(Ntp->MET_CorrT0rt_phi);}
   float              MET_CorrT0rtT1_et(){return Ntp->MET_CorrT0rtT1_et;}
   float              MET_CorrT0rtT1_phi(){return Ntp->MET_CorrT0rtT1_phi;}
   float              MET_CorrT0rtT1_ex(){return Ntp->MET_CorrT0rtT1_et*cos(Ntp->MET_CorrT0rtT1_phi);}
   float              MET_CorrT0rtT1_ey(){return Ntp->MET_CorrT0rtT1_et*sin(Ntp->MET_CorrT0rtT1_phi);}
   float              MET_CorrT0pc_et(){return Ntp->MET_CorrT0pc_et;}
   float              MET_CorrT0pc_phi(){return Ntp->MET_CorrT0pc_phi;}
   float              MET_CorrT0pc_ex(){return Ntp->MET_CorrT0pc_et*cos(Ntp->MET_CorrT0pc_phi);}
   float              MET_CorrT0pc_ey(){return Ntp->MET_CorrT0pc_et*sin(Ntp->MET_CorrT0pc_phi);}
   float              MET_CorrT0pcT1_et(){return Ntp->MET_CorrT0pcT1_et;}
   float              MET_CorrT0pcT1_phi(){return Ntp->MET_CorrT0pcT1_phi;}
   float              MET_CorrT0pcT1_ex(){return Ntp->MET_CorrT0pcT1_et*cos(Ntp->MET_CorrT0pcT1_phi);}
   float              MET_CorrT0pcT1_ey(){return Ntp->MET_CorrT0pcT1_et*sin(Ntp->MET_CorrT0pcT1_phi);}
   float              MET_CorrT0rtTxy_et(){return Ntp->MET_CorrT0rtTxy_et;}
   float              MET_CorrT0rtTxy_phi(){return Ntp->MET_CorrT0rtTxy_phi;}
   float              MET_CorrT0rtTxy_ex(){return Ntp->MET_CorrT0rtTxy_et*cos(Ntp->MET_CorrT0rtTxy_phi);}
   float              MET_CorrT0rtTxy_ey(){return Ntp->MET_CorrT0rtTxy_et*sin(Ntp->MET_CorrT0rtTxy_phi);}
   float              MET_CorrT0rtT1Txy_et(){return Ntp->MET_CorrT0rtT1Txy_et;}
   float              MET_CorrT0rtT1Txy_phi(){return Ntp->MET_CorrT0rtT1Txy_phi;}
   float              MET_CorrT0rtT1Txy_ex(){return Ntp->MET_CorrT0rtT1Txy_et*cos(Ntp->MET_CorrT0rtT1Txy_phi);}
   float              MET_CorrT0rtT1Txy_ey(){return Ntp->MET_CorrT0rtT1Txy_et*sin(Ntp->MET_CorrT0rtT1Txy_phi);}
   float              MET_CorrT0pcTxy_et(){return Ntp->MET_CorrT0pcTxy_et;}
   float              MET_CorrT0pcTxy_phi(){return Ntp->MET_CorrT0pcTxy_phi;}
   float              MET_CorrT0pcTxy_ex(){return Ntp->MET_CorrT0pcTxy_et*cos(Ntp->MET_CorrT0pcTxy_phi);}
   float              MET_CorrT0pcTxy_ey(){return Ntp->MET_CorrT0pcTxy_et*sin(Ntp->MET_CorrT0pcTxy_phi);}
   float              MET_CorrT0pcT1Txy_et(){return Ntp->MET_CorrT0pcT1Txy_et;}
   float              MET_CorrT0pcT1Txy_phi(){return Ntp->MET_CorrT0pcT1Txy_phi;}
   float              MET_CorrT0pcT1Txy_ex(){return Ntp->MET_CorrT0pcT1Txy_et*cos(Ntp->MET_CorrT0pcT1Txy_phi);}
   float              MET_CorrT0pcT1Txy_ey(){return Ntp->MET_CorrT0pcT1Txy_et*sin(Ntp->MET_CorrT0pcT1Txy_phi);}
   float              MET_CorrT1_et(){return Ntp->MET_CorrT1_et;}
   float              MET_CorrT1_phi(){return Ntp->MET_CorrT1_phi;}
   float              MET_CorrT1_ex(){return Ntp->MET_CorrT1_et*cos(Ntp->MET_CorrT1_phi);}
   float              MET_CorrT1_ey(){return Ntp->MET_CorrT1_et*sin(Ntp->MET_CorrT1_phi);}
   float              MET_CorrT1Txy_et(){return Ntp->MET_CorrT1Txy_et;}
   float              MET_CorrT1Txy_phi(){return Ntp->MET_CorrT1Txy_phi;}
   float              MET_CorrT1Txy_ex(){return Ntp->MET_CorrT1Txy_et*cos(Ntp->MET_CorrT1Txy_phi);}
   float              MET_CorrT1Txy_ey(){return Ntp->MET_CorrT1Txy_et*sin(Ntp->MET_CorrT1Txy_phi);}
   float              MET_CorrCaloT1_et(){return Ntp->MET_CorrCaloT1_et;}
   float              MET_CorrCaloT1_phi(){return Ntp->MET_CorrCaloT1_phi;}
   float              MET_CorrCaloT1_ex(){return Ntp->MET_CorrCaloT1_et*cos(Ntp->MET_CorrCaloT1_phi);}
   float              MET_CorrCaloT1_ey(){return Ntp->MET_CorrCaloT1_et*sin(Ntp->MET_CorrCaloT1_phi);}
   float              MET_CorrCaloT1T2_et(){return Ntp->MET_CorrCaloT1T2_et;}
   float              MET_CorrCaloT1T2_phi(){return Ntp->MET_CorrCaloT1T2_phi;}
   float              MET_CorrCaloT1T2_ex(){return Ntp->MET_CorrCaloT1T2_et*cos(Ntp->MET_CorrCaloT1T2_phi);}
   float              MET_CorrCaloT1T2_ey(){return Ntp->MET_CorrCaloT1T2_et*sin(Ntp->MET_CorrCaloT1T2_phi);}

   float              MET_CorrMVA_et(){return Ntp->MET_CorrMVA_et;}
   float              MET_CorrMVA_phi(){return Ntp->MET_CorrMVA_phi;}
   float              MET_CorrMVA_ex(){return Ntp->MET_CorrMVA_et*cos(Ntp->MET_CorrMVA_phi);}
   float              MET_CorrMVA_ey(){return Ntp->MET_CorrMVA_et*sin(Ntp->MET_CorrMVA_phi);}
   float              MET_CorrMVA_significance(){return Ntp->MET_CorrMVA_significance;}
   float              MET_CorrMVA_significance_xx(){return Ntp->MET_CorrMVA_significance_xx;}
   float              MET_CorrMVA_significance_xy(){return Ntp->MET_CorrMVA_significance_xy;}
   float              MET_CorrMVA_significance_yy(){return Ntp->MET_CorrMVA_significance_yy;}

   float              MET_CorrMVAMuTau_et(){return Ntp->MET_CorrMVAMuTau_et;}
   float              MET_CorrMVAMuTau_phi(){return Ntp->MET_CorrMVAMuTau_phi;}
   float              MET_CorrMVAMuTau_ex(){return Ntp->MET_CorrMVAMuTau_et*cos(Ntp->MET_CorrMVAMuTau_phi);}
   float              MET_CorrMVAMuTau_ey(){return Ntp->MET_CorrMVAMuTau_et*sin(Ntp->MET_CorrMVAMuTau_phi);}
   float              MET_CorrMVAMuTau_significance(){return Ntp->MET_CorrMVAMuTau_significance;}
   float              MET_CorrMVAMuTau_significance_xx(){return Ntp->MET_CorrMVAMuTau_significance_xx;}
   float              MET_CorrMVAMuTau_significance_xy(){return Ntp->MET_CorrMVAMuTau_significance_xy;}
   float              MET_CorrMVAMuTau_significance_yy(){return Ntp->MET_CorrMVAMuTau_significance_yy;}

   //Track Information
   unsigned int      NTracks(){return Ntp->Track_p4->size();}
   TLorentzVector    Track_p4(unsigned int i){return TLorentzVector(Ntp->Track_p4->at(i).at(1),Ntp->Track_p4->at(i).at(2),Ntp->Track_p4->at(i).at(3),Ntp->Track_p4->at(i).at(0));}
   TVector3          Track_Poca(unsigned int i){return TVector3(Ntp->Track_Poca->at(i).at(0),Ntp->Track_Poca->at(i).at(1),Ntp->Track_Poca->at(i).at(2));}
   int               Track_charge(unsigned int i){return Ntp->Track_charge->at(i);}
   float             Track_chi2(unsigned int i){return Ntp->Track_chi2->at(i);}
   float             Track_ndof(unsigned int i){return Ntp->Track_ndof->at(i);}
   unsigned short    Track_numberOfLostHits(unsigned int i){return Ntp->Track_numberOfLostHits->at(i);}
   unsigned short    Track_numberOfValidHits(unsigned int i){return Ntp->Track_numberOfValidHits->at(i);}
   unsigned int      Track_qualityMask(unsigned int i){return Ntp->Track_qualityMask->at(i);}

   TrackParticle Track_TrackParticle(unsigned int i){
     TMatrixT<double>    track_par(TrackParticle::NHelixPar,1);
     TMatrixTSym<double> track_cov(TrackParticle::NHelixPar);
     unsigned int l=0;
     for(int k=0; k<TrackParticle::NHelixPar; k++){
       track_par(k,0)=Ntp->Track_par->at(i).at(k);
       for(int j=k; j<TrackParticle::NHelixPar; j++){
	 track_cov(k,j)=Ntp->Track_cov->at(i).at(l);
	 l++;
       }
     }
     return TrackParticle(track_par,track_cov,Ntp->Track_pdgid->at(i),Ntp->Track_M->at(i),Ntp->Track_charge->at(i),Ntp->Track_B->at(i));
   }

   // MC Informaation
   // Singal particles (Z0,W+/-,H0,H+/-)
   unsigned int               NMCSignalParticles(){return Ntp->MCSignalParticle_p4->size();}
   TLorentzVector             MCSignalParticle_p4(unsigned int i){return TLorentzVector(Ntp->MCSignalParticle_p4->at(i).at(1),Ntp->MCSignalParticle_p4->at(i).at(2),Ntp->MCSignalParticle_p4->at(i).at(3),Ntp->MCSignalParticle_p4->at(i).at(0));}
   int                        MCSignalParticle_pdgid(unsigned int i){return Ntp->MCSignalParticle_pdgid->at(i);}
   int                        MCSignalParticle_charge(unsigned int i){return Ntp->MCSignalParticle_charge->at(i);}
   TVector3                   MCSignalParticle_Poca(unsigned int i){return TVector3(Ntp->MCSignalParticle_Poca->at(i).at(0),Ntp->MCSignalParticle_Poca->at(i).at(1),Ntp->MCSignalParticle_Poca->at(i).at(2));}
   std::vector<unsigned int>  MCSignalParticle_Tauidx(unsigned int i){return Ntp->MCSignalParticle_Tauidx->at(i);}
   unsigned int               NMCParticles(){return Ntp->MC_p4->size();}
   TLorentzVector             MCParticle_p4(unsigned int i){return TLorentzVector(Ntp->MC_p4->at(i).at(1),Ntp->MC_p4->at(i).at(2),Ntp->MC_p4->at(i).at(3),Ntp->MC_p4->at(i).at(0));}
   int                        MCParticle_pdgid(unsigned int i){return Ntp->MC_pdgid->at(i);}
   int                        MCParticle_charge(unsigned int i){return Ntp->MC_charge->at(i);}
   int                        MCParticle_midx(unsigned int i){return Ntp->MC_midx->at(i);}
   std::vector<int>           MCParticle_childpdgid(unsigned int i){return Ntp->MC_childpdgid->at(i);}

   // Tau decays (Tau is first element of vector)
   int NMCTaus(){return Ntp->MCTauandProd_p4->size();}
   TLorentzVector MCTau_p4(unsigned int i){return MCTauandProd_p4(i,0);}
   int MCTau_pdgid(unsigned int i){return MCTauandProd_pdgid(i,0);}
   int MCTau_charge(unsigned int i){return MCTauandProd_charge(i,0);}
   unsigned int MCTau_JAK(unsigned int i){return Ntp->MCTau_JAK->at(i);}
   unsigned int MCTau_DecayBitMask(unsigned int i){return Ntp->MCTau_DecayBitMask->at(i);}
   //Tau and decay products
   int NMCTauDecayProducts(unsigned int i){if(0<=i && i<NMCTaus()) return Ntp->MCTauandProd_p4->at(i).size(); return 0;}
   TLorentzVector MCTauandProd_p4(unsigned int i, unsigned int j){return TLorentzVector(Ntp->MCTauandProd_p4->at(i).at(j).at(1),Ntp->MCTauandProd_p4->at(i).at(j).at(2),Ntp->MCTauandProd_p4->at(i).at(j).at(3),Ntp->MCTauandProd_p4->at(i).at(j).at(0));}
   int MCTauandProd_pdgid(unsigned int i, unsigned int j){return Ntp->MCTauandProd_pdgid->at(i).at(j);}
   unsigned int MCTauandProd_midx(unsigned int i, unsigned int j){return Ntp->MCTauandProd_midx->at(i).at(j);}
   int MCTauandProd_charge(unsigned int i, unsigned int j){return Ntp->MCTauandProd_charge->at(i).at(j);}
   TVector3 MCTauandProd_Vertex(unsigned int i, unsigned int j){
     return TVector3(Ntp->MCTauandProd_Vertex->at(i).at(j).at(0),Ntp->MCTauandProd_Vertex->at(i).at(j).at(1),Ntp->MCTauandProd_Vertex->at(i).at(j).at(2));
   }
   bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,TauDecay::JAK tau_jak, unsigned int &idx);
   bool hasSignalTauDecay(PDGInfo::PDGMCNumbering parent_pdgid,unsigned int &Boson_idx,unsigned int &tau1_idx, unsigned int &tau2_idx);

   bool jethasMuonOverlap(unsigned int jet_idx,unsigned int &muon_idx);
   bool muonhasJetOverlap(unsigned int muon_idx,unsigned int &jet_idx);
   bool muonhasJetMatch(unsigned int muon_idx,unsigned int &jet_idx);


   // Electrons
   unsigned int       NElectrons(){return Ntp->Electron_p4->size();}
   TLorentzVector     Electron_p4(unsigned int i){return TLorentzVector(Ntp->Electron_p4->at(i).at(1),Ntp->Electron_p4->at(i).at(2),Ntp->Electron_p4->at(i).at(3),Ntp->Electron_p4->at(i).at(0));}
   TVector3           Electron_Poca(unsigned int i){return TVector3(Ntp->Electron_Poca->at(i).at(0),Ntp->Electron_Poca->at(i).at(1),Ntp->Electron_Poca->at(i).at(2));}
   int   Electron_Charge(unsigned int i){return Ntp->Electron_charge->at(i);}
   float   Electron_Gsf_deltaEtaEleClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaEtaEleClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaEtaSeedClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaEtaSeedClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaEtaSuperClusterTrackAtVtx(unsigned int i){return Ntp->Electron_Gsf_deltaEtaSuperClusterTrackAtVtx->at(i);}
   float   Electron_Gsf_deltaPhiEleClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaPhiEleClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaPhiSeedClusterTrackAtCalo(unsigned int i){return Ntp->Electron_Gsf_deltaPhiSeedClusterTrackAtCalo->at(i);}
   float   Electron_Gsf_deltaPhiSuperClusterTrackAtVtx(unsigned int i){return Ntp->Electron_Gsf_deltaPhiSuperClusterTrackAtVtx->at(i);}
   float   Electron_Gsf_dr03EcalRecHitSumE(unsigned int i){return Ntp->Electron_Gsf_dr03EcalRecHitSumE->at(i);}
   float   Electron_Gsf_dr03HcalDepth1TowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth1TowerSumEt->at(i);}
   float   Electron_Gsf_dr03HcalDepth1TowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth1TowerSumEtBc->at(i);}
   float   Electron_Gsf_dr03HcalDepth2TowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth2TowerSumEt->at(i);}
   float   Electron_Gsf_dr03HcalDepth2TowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalDepth2TowerSumEtBc->at(i);}
   float   Electron_Gsf_dr03HcalTowerSumEt(unsigned int i){return Ntp->Electron_Gsf_dr03HcalTowerSumEt->at(i);}
   float   Electron_Gsf_dr03HcalTowerSumEtBc(unsigned int i){return Ntp->Electron_Gsf_dr03HcalTowerSumEtBc->at(i);}
   float   Electron_Gsf_dr03TkSumPt(unsigned int i){return Ntp->Electron_Gsf_dr03TkSumPt->at(i);}
   bool    Electron_Gsf_passingCutBasedPreselection(unsigned int i){return Ntp->Electron_Gsf_passingCutBasedPreselection->at(i);}
   bool    Electron_Gsf_passingMvaPreselection(unsigned int i){return Ntp->Electron_Gsf_passingMvaPreselection->at(i);}
   int     Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits(unsigned int i){return Ntp->Electron_gsftrack_trackerExpectedHitsInner_numberOfLostHits->at(i);}
   float   Electron_supercluster_e(unsigned int i){return Ntp->Electron_supercluster_e->at(i);}
   float   Electron_supercluster_phi(unsigned int i){return Ntp->Electron_supercluster_phi->at(i);}
   float   Electron_supercluster_eta(unsigned int i){return Ntp->Electron_supercluster_eta->at(i);}
   float   Electron_supercluster_centroid_x(unsigned int i){return Ntp->Electron_supercluster_centroid_x->at(i);}
   float   Electron_supercluster_centroid_y(unsigned int i){return Ntp->Electron_supercluster_centroid_y->at(i);}
   float   Electron_supercluster_centroid_z(unsigned int i){return Ntp->Electron_supercluster_centroid_z->at(i);}
   unsigned int Electron_Track_idx(unsigned int i){return Ntp->Electron_Track_idx->at(i);}

   float    Electron_ecalRecHitSumEt03(unsigned int i){return Ntp->Electron_ecalRecHitSumEt03->at(i);}
   float    Electron_hcalDepth1TowerSumEt03(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEt03->at(i);}
   float    Electron_hcalDepth1TowerSumEtBc03(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEtBc03->at(i);}
   float    Electron_hcalDepth2TowerSumEt03(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEt03->at(i);}
   float    Electron_hcalDepth2TowerSumEtBc03(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEtBc03->at(i);}
   float    Electron_tkSumPt03(unsigned int i){return Ntp->Electron_tkSumPt03->at(i);}
   float    Electron_ecalRecHitSumEt04(unsigned int i){return Ntp->Electron_ecalRecHitSumEt04->at(i);}
   float    Electron_hcalDepth1TowerSumEt04(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEt04->at(i);}
   float    Electron_hcalDepth1TowerSumEtBc04(unsigned int i){return Ntp->Electron_hcalDepth1TowerSumEtBc04->at(i);}
   float    Electron_hcalDepth2TowerSumEt04(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEt04->at(i);}
   float    Electron_hcalDepth2TowerSumEtBc04(unsigned int i){return Ntp->Electron_hcalDepth2TowerSumEtBc04->at(i);}
   float    Electron_tkSumPt04(unsigned int i){return Ntp->Electron_tkSumPt04->at(i);}
   
   float    Electron_chargedHadronIso(unsigned int i){return Ntp->Electron_chargedHadronIso->at(i);}
   float    Electron_neutralHadronIso(unsigned int i){return Ntp->Electron_neutralHadronIso->at(i);}
   float    Electron_photonIso(unsigned int i){return Ntp->Electron_photonIso->at(i);}
   
   float    Electron_sigmaIetaIeta(unsigned int i){return Ntp->Electron_sigmaIetaIeta->at(i);}
   float    Electron_hadronicOverEm(unsigned int i){return Ntp->Electron_hadronicOverEm->at(i);}
   float    Electron_fbrem(unsigned int i){return Ntp->Electron_fbrem->at(i);}
   float    Electron_eSuperClusterOverP(unsigned int i){return Ntp->Electron_eSuperClusterOverP->at(i);}
   float    Electron_ecalEnergy(unsigned int i){return Ntp->Electron_ecalEnergy->at(i);}
   float    Electron_trackMomentumAtVtx(unsigned int i){return Ntp->Electron_trackMomentumAtVtx->at(i);}
   float    Electron_numberOfMissedHits(unsigned int i){return Ntp->Electron_numberOfMissedHits->at(i);}
   bool     Electron_HasMatchedConversions(unsigned int i){return Ntp->Electron_HasMatchedConversions->at(i);}

   float    Electron_MVA_Trig_discriminator(unsigned int i){return Ntp->Electron_MVA_Trig_discriminator->at(i);}
   float    Electron_MVA_TrigNoIP_discriminator(unsigned int i){return Ntp->Electron_MVA_TrigNoIP_discriminator->at(i);}
   float    Electron_MVA_NonTrig_discriminator(unsigned int i){return Ntp->Electron_MVA_NonTrig_discriminator->at(i);}
   float    RhoIsolationAllInputTags(){return Ntp->RhoIsolationAllInputTags;}

   float	Electron_RegEnergy(unsigned int i){return Ntp->Electron_RegEnergy->at(i);}
   float	Electron_RegEnergyError(unsigned int i){return Ntp->Electron_RegEnergyError->at(i);}

   TrackParticle Electron_TrackParticle(unsigned int i){
     TMatrixT<double>    e_par(TrackParticle::NHelixPar,1);
     TMatrixTSym<double> e_cov(TrackParticle::NHelixPar);
     unsigned int l=0;
     for(int k=0; k<TrackParticle::NHelixPar; k++){
       e_par(k,0)=Ntp->Electron_par->at(i).at(k);
       for(int j=k; j<TrackParticle::NHelixPar; j++){
	 e_cov(k,j)=Ntp->Electron_cov->at(i).at(l);
	 l++;
       }
     }
     return TrackParticle(e_par,e_cov,Ntp->Electron_pdgid->at(i),Ntp->Electron_M->at(i),Ntp->Electron_charge->at(i),Ntp->Electron_B->at(i));
   }

   bool isTrigPreselElectron(unsigned int i);
   bool isTrigNoIPPreselElectron(unsigned int i);
   bool isMVATrigElectron(unsigned int i);
   bool isMVATrigNoIPElectron(unsigned int i);
   bool isMVANonTrigElectron(unsigned int i, unsigned int j);
   bool isTightElectron(unsigned int i);
   bool isTightElectron(unsigned int i, unsigned int j);
   float Electron_RelIso03(unsigned int i);
   float Electron_RelIso04(unsigned int i);
   float Electron_Aeff_R04(double Eta);
   float Electron_Aeff_R03(double Eta);
   bool isSelectedElectron(unsigned int i, unsigned int j, double impact_xy, double impact_z);

   // Trigger
   bool         TriggerAccept(TString n);
   unsigned int HLTPrescale(TString n);
   unsigned int L1SEEDPrescale(TString n);
   bool         GetTriggerIndex(TString n, unsigned int &i);
   unsigned int NHLTTriggers(){return Ntp->HTLTriggerName->size();}
   std::string  HTLTriggerName(unsigned int i){return Ntp->HTLTriggerName->at(i);}
   bool         TriggerAccept(unsigned int i){return Ntp->TriggerAccept->at(i);}
   bool         TriggerError(unsigned int i){return Ntp->TriggerError->at(i);}
   bool         TriggerWasRun(unsigned int i){return Ntp->TriggerWasRun->at(i);}
   unsigned int HLTPrescale(unsigned int i){return Ntp->HLTPrescale->at(i);}
   unsigned int NHLTL1GTSeeds(unsigned int i){return Ntp->NHLTL1GTSeeds->at(i);}
   unsigned int L1SEEDPrescale(unsigned int i){return Ntp->L1SEEDPrescale->at(i);}
   bool         L1SEEDInvalidPrescale(unsigned int i){return Ntp->L1SEEDInvalidPrescale->at(i);}
   float        MuonTriggerMatch(unsigned int i, unsigned int j){if(j<Ntp->TriggerMatchMuon->at(i).size()) return Ntp->TriggerMatchMuon->at(i).at(j); return 999;}
   //float        ElectronTriggerMatch(unsigned int i, unsigned int j){if(j<Ntp->ElectronTriggerMatch->at(i).size()) return Ntp->ElectronTriggerMatch->at(i).at(j);return 999;}
   float        JetTriggerMatch(unsigned int i, unsigned int j){if(j<Ntp->TriggerMatchJet->at(i).size()) return Ntp->TriggerMatchJet->at(i).at(j);return 999;}
   float        TauTriggerMatch(unsigned int i, unsigned int j){if(j<Ntp->TriggerMatchTau->at(i).size()) return Ntp->TriggerMatchTau->at(i).at(j);return 999;}
   unsigned int NHLTTriggerObject(unsigned int i){return Ntp->HLTTrigger_objs_Eta->at(i).size();}
   TLorentzVector HLTTriggerObject_p4(unsigned int i, unsigned int j){
     TLorentzVector L(0,0,0,0); 
     if(j<Ntp->HLTTrigger_objs_Eta->at(i).size())L.SetPtEtaPhiM(Ntp->HLTTrigger_objs_Pt->at(i).at(j),Ntp->HLTTrigger_objs_Eta->at(i).at(j), Ntp->HLTTrigger_objs_Phi->at(i).at(j),0.0);
     return L;
   }
   int          NHLTTrigger_objs(){return Ntp->HLTTrigger_objs_Pt->size();}
   int          NHLTTrigger_objs(unsigned int i){return Ntp->HLTTrigger_objs_Pt->at(i).size();}
   float        HLTTrigger_objs_Pt(unsigned int i, unsigned int j){return Ntp->HLTTrigger_objs_Pt->at(i).at(j);}
   float        HLTTrigger_objs_Eta(unsigned int i, unsigned int j){return Ntp->HLTTrigger_objs_Eta->at(i).at(j);}
   float        HLTTrigger_objs_Phi(unsigned int i, unsigned int j){return Ntp->HLTTrigger_objs_Phi->at(i).at(j);}
   float        HLTTrigger_objs_E(unsigned int i,unsigned int j){return Ntp->HLTTrigger_objs_E->at(i).at(j);}
   int          HLTTrigger_objs_Id(unsigned int i,unsigned int j){return Ntp->HLTTrigger_objs_Id->at(i).at(j);}
   std::string  HLTTrigger_objs_trigger(unsigned int i){return Ntp->HLTTrigger_objs_trigger->at(i);}

   // helper functions
   float        dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   float        dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
   float        vertexSignificance(TVector3 vec, unsigned int vertex);

};

#endif

