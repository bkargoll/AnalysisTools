#ifndef HToTaumuTauh_h
#define HToTaumuTauh_h

#include <TH1.h>
#include <TString.h>
#include <cmath>
#include <vector>

#include "../Selection.h"

class TLorentzVector;
class TVector3;

class HToTaumuTauh : public Selection {

 public:
  HToTaumuTauh(TString Name_, TString id_);
  virtual ~HToTaumuTauh();

  virtual void Configure();

  enum cuts {
	  TriggerOk=0,
	  PrimeVtx,
	  NMuId,
	  NMuKin,
	  DiMuonVeto,
	  NTauId,
	  NTauIso,
	  NTauKin,
	  TriLeptonVeto,
	  OppCharge,
	  MT,
	  BJetVeto,
	  CatCut1,
	  CatCut2,
	  CatCut3,
	  CatCut4,
	  CatCut5,
	  NCuts};

  	// cuts in categories
	enum cuts_VBFTight {
		VbfTight_NJet	= CatCut1,
		VbfTight_DeltaEta,
		VbfTight_NJetRapGap,
		VbfTight_JetInvM,
		VbfTight_HiggsPt,
		VbfTight_NCuts
	};
	enum cuts_VBFLoose {
		VbfLoose_NJet	= CatCut1,
		VbfLoose_DeltaEta,
		VbfLoose_NJetRapGap,
		VbfLoose_JetInvM,
		VbfLoose_NotVbfTight,
		VbfLoose_NCuts
	};
	enum cuts_OneJetLow {
		OneJetLow_NJet = CatCut1,
		OneJetLow_NotVbf,
		OneJetLow_TauPt,
		OneJetLow_NCuts
	};
	enum cuts_OneJetHigh {
		OneJetHigh_NJet = CatCut1,
		OneJetHigh_NotVbf,
		OneJetHigh_TauPt,
		OneJetHigh_HiggsPt,
		OneJetHigh_NCuts
	};
	enum cuts_OneJetBoost {
		OneJetBoost_NJet = CatCut1,
		OneJetBoost_NotVbf,
		OneJetBoost_TauPt,
		OneJetBoost_HiggsPt,
		OneJetBoost_NCuts
	};
	enum cuts_ZeroJetHigh {
		ZeroJetHigh_NJet = CatCut1,
		ZeroJetHigh_TauPt,
		ZeroJetHigh_NCuts
	};
	enum cuts_ZeroJetLow {
		ZeroJetLow_NJet = CatCut1,
		ZeroJetLow_TauPt,
		ZeroJetLow_NCuts
	};
	enum cuts_NoCategory {
		NoCategory_NCuts = CatCut1
	};

 protected:
  virtual void Setup();
  virtual void doEvent();
  virtual void Store_ExtraDist();
  virtual void  Finish();

  // Selection Variables
  std::vector<TH1D> NCatFired;
  std::vector<TH1D> CatFired;

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> VtxZ;
  std::vector<TH1D> VtxRho;
  std::vector<TH1D> VtxPhi;
  std::vector<TH1D> VtxNdof;
  std::vector<TH1D> VtxIsfake;

  std::vector<TH1D> MuDxy;
  std::vector<TH1D> MuDz;
  std::vector<TH1D> MuRelIso;
  std::vector<TH1D> MuPt;
  std::vector<TH1D> MuEta;
  std::vector<TH1D> MuPhi;

  std::vector<TH1D> MuSelPt;
  std::vector<TH1D> MuSelEta;
  std::vector<TH1D> MuSelPhi;
  std::vector<TH1D> MuSelDxy;
  std::vector<TH1D> MuSelDz;
  std::vector<TH1D> MuSelRelIso;
  std::vector<TH1D> MuSelFakesTauID;
  std::vector<TH1D> MuSelDrHlt;

  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPhi;
  std::vector<TH1D> TauDecayMode;
  std::vector<TH1D> TauIso;

  std::vector<TH1D> TauSelPt;
  std::vector<TH1D> TauSelEta;
  std::vector<TH1D> TauSelPhi;
  std::vector<TH1D> TauSelDrHlt; // todo: not filled at the moment
  std::vector<TH1D> TauSelDecayMode;
  std::vector<TH1D> TauSelIso;

  std::vector<TH1D> MuVetoDPtSelMuon;
  std::vector<TH1D> MuVetoInvM;
  std::vector<TH1D> MuVetoPtPositive;
  std::vector<TH1D> MuVetoPtNegative;
  std::vector<TH1D> MuVetoDRTau;
  std::vector<TH1D> MuVetoDeltaR;

  std::vector<TH1D> NMuonTriLepVeto;
  std::vector<TH1D> NElecTriLepVeto;

  std::vector<TH1D> MuCharge;
  std::vector<TH1D> TauCharge;

  std::vector<TH1D> MuTauDR;
  std::vector<TH1D> MuTauDPhi;
  std::vector<TH1D> MuTauDEta;
  std::vector<TH1D> MuTauDPt;
  std::vector<TH1D> MuTauRelDPt;
  std::vector<TH2D> MuPtVsTauPt;

  std::vector<TH1D> MetPt;
  std::vector<TH1D> MetPhi;

  std::vector<TH1D> MetLepMuDr;
  std::vector<TH1D> MetLepTauDr;
  std::vector<TH1D> MetLepNMu;
  std::vector<TH1D> MetLepNTau;
  std::vector<TH1D> MetLepNMuMinusNMu;
  std::vector<TH1D> MetLepNTauMinusNTau;

  std::vector<TH1D> NJetsKin;
  std::vector<TH1D> JetKin1Pt;
  std::vector<TH1D> JetKin1Eta;
  std::vector<TH1D> JetKin1Phi;
  std::vector<TH1D> JetKin1IsLooseId;
  std::vector<TH1D> JetKin2IsLooseId;
  std::vector<TH1D> JetKin2Pt;
  std::vector<TH1D> JetKin2Eta;
  std::vector<TH1D> JetKin2Phi;
  std::vector<TH1D> NJetsId;
  std::vector<TH1D> Jet1Pt;
  std::vector<TH1D> Jet1Eta;
  std::vector<TH1D> Jet1Phi;
  std::vector<TH1D> Jet1IsB;
  std::vector<TH1D> Jet2Pt;
  std::vector<TH1D> Jet2Eta;
  std::vector<TH1D> Jet2Phi;
  std::vector<TH1D> Jet2IsB;

  std::vector<TH1D> NBJets;
  std::vector<TH1D> BJet1Pt;
  std::vector<TH1D> BJet1Eta;
  std::vector<TH1D> BJet1Phi;

  std::vector<TH1D> HiggsPt;
  std::vector<TH1D> HiggsPhi;
  std::vector<TH1D> JetsDEta;
  std::vector<TH1D> JetsInEtaGap;
  std::vector<TH1D> JetsInvM;

  std::vector<TH1D> MetPhiMet10GeV;
  std::vector<TH1D> MtMet10GeV;
  std::vector<TH1D> HiggsPtMet10GeV;
  std::vector<TH1D> HiggsPhiMet10GeV;
  std::vector<TH1D> MetPhiMet20GeV;
  std::vector<TH1D> MtMet20GeV;
  std::vector<TH1D> HiggsPtMet20GeV;
  std::vector<TH1D> HiggsPhiMet20GeV;

  std::vector<TH1D> MtAfterMuon;
  std::vector<TH1D> MtAfterDiMuonVeto;
  std::vector<TH1D> MtAfterTau;
  std::vector<TH1D> MtAfterTriLepVeto;
  std::vector<TH1D> MtAfterOppCharge;
  std::vector<TH1D> MtAfterBJetVeto;
  std::vector<TH1D> MtOnlyTau;
  std::vector<TH1D> MtOnlyTriLepVeto;
  std::vector<TH1D> MtOnlyOppCharge;
  std::vector<TH1D> MtOnlyBJet;
  std::vector<TH1D> MtMuPlusOnly;
  std::vector<TH1D> MtMuMinusOnly;
  std::vector<TH1D> MtMuPlusOnlyBGSubt;
  std::vector<TH1D> MtMuMinusOnlyBGSubt;
  std::vector<TH1D> Mt1ProngOnly;
  std::vector<TH1D> Mt3ProngOnly;
  std::vector<TH1D> Mt3ProngSV;

  std::vector<TH1D> Cat0JetLowQcdShapeRegion;
  std::vector<TH1D> Cat0JetHighLowQcdShapeRegion;
  std::vector<TH1D> Cat1JetLowQcdShapeRegion;
  std::vector<TH1D> Cat1JetHighQcdShapeRegion;
  std::vector<TH1D> Cat1JetBoostQcdShapeRegion;
  std::vector<TH1D> CatVBFLooseQcdShapeRegion;
  std::vector<TH1D> CatVBFTightQcdShapeRegion;
  std::vector<TH1D> CatInclusiveQcdShapeRegion;

  unsigned verbose;

  // cut values
  double cMu_dxy, cMu_dz, cMu_relIso, cMu_pt, cMu_eta, cMu_dRHltMatch;
  double cTau_pt, cTau_eta, cTau_rawIso, cMuTau_dR, cTau_dRHltMatch;
  double cMuTriLep_pt, cMuTriLep_eta, cEleTriLep_pt, cEleTriLep_eta;
  std::vector<TString> cTriggerNames;
  double cCat_jetPt, cCat_jetEta, cCat_bjetPt, cCat_bjetEta, cCat_btagDisc, cCat_splitTauPt, cJetClean_dR;

  // flag for category to run
  TString categoryFlag;
  // flag for WJets background source
  TString wJetsBGSource;
  // flag for data-driven QCD shape (set to false for yield estimation!)
  bool qcdShapeFromData;

  // map to hold WJets yields for each category
  std::map<TString, double> wJetsYieldMap;
  //std::map<TString, double> wJetsYieldScaleMap;

  // map to hold QCD yields for each category
  std::map<TString, double> qcdYieldMap;

  // variables to hold selected objects (to be used e.g. for sync Ntuple)
  int selVertex;
  int selMuon;
  int selTau;
  std::vector<int> selJets, selBJets;
  double selMjj, selJetdeta;
  int selNjetingap;
  double w; // event weight
  unsigned int t; // index of histogram
  bool isWJetMC; // for Wjets background method

  // booleans for different analysis stages
  void setStatusBooleans(bool resetAll = false);
  bool passedVertex;
  bool passedMuId;
  bool passedMu;
  bool passedTauIdIso;
  bool passedTau;
  bool passedObjects;
  bool passedDiMuonVeto;
  bool passedFullInclusiveSelNoBVeto;
  bool passedFullInclusiveSel;
  bool passedFullInclusiveSelNoMt;
  bool passedFullInclusiveSelNoMtNoOS;
  bool passedFullInclusiveNoTauNoMuNoCharge;
  bool passedObjectsFailDiMuonVeto;
  bool passed_VBFTight;
  bool passed_VBFLoose;
  bool passed_VBF;
  bool passed_OneJetHigh;
  bool passed_OneJetLow;
  bool passed_OneJetBoost;
  bool passed_ZeroJetHigh;
  bool passed_ZeroJetLow;
  bool passed_NoCategory;
  bool passed_VBFTightRelaxed;
  bool passed_VBFLooseRelaxed;

  bool hasRelaxedIsoTau, hasAntiIsoMuon;

  // function definitions
  bool selectMuon_Id(unsigned i, unsigned vertex);
  bool selectMuon_Kinematics(unsigned i);

  bool selectMuon_antiIso(unsigned i, unsigned vertex);

  bool selectMuon_diMuonVeto(unsigned i, unsigned i_vtx);
  bool selectMuon_triLeptonVeto(unsigned i, int selectedMuon, unsigned i_vtx);

  bool selectElectron_triLeptonVeto(unsigned i, unsigned i_vtx, std::vector<int>);

  bool selectPFTau_Id(unsigned i);
  bool selectPFTau_Id(unsigned i, std::vector<int> muonCollection);
  bool selectPFTau_Iso(unsigned i);
  bool selectPFTau_Kinematics(unsigned i);

  bool selectPFTau_relaxedIso(unsigned i, std::vector<int> muonCollection);

  bool selectPFJet_Cleaning(unsigned i, int selectedMuon, int selectedTau);
  bool selectPFJet_Kinematics(unsigned i);
  bool selectPFJet_Id(unsigned i);
  bool selectBJet(unsigned i, int selectedMuon, int selectedTau);

  // categories
  std::vector<float> cut_VBFTight, cut_VBFLoose;
  std::vector<float> cut_OneJetHigh, cut_OneJetLow, cut_OneJetBoost;
  std::vector<float> cut_ZeroJetHigh, cut_ZeroJetLow;
  std::vector<float> cut_NoCategory;
  // relaxed categories for background methods
  std::vector<float> cut_VBFTightRelaxed, cut_VBFLooseRelaxed;

  bool migrateCategoryIntoMain(TString thisCategory, std::vector<float> categoryValueVector, std::vector<float> categoryPassVector, int categoryNCuts);

  void configure_VBFTight();
  bool category_VBFTight(unsigned NJets, double DEta, int NJetsInGap, double Mjj, double higgsPt);

  void configure_VBFLoose();
  bool category_VBFLoose(unsigned NJets, double DEta, int NJetsInGap, double Mjj, bool passedVBFTight);

  void configure_OneJetHigh();
  bool category_OneJetHigh(unsigned NJets, double TauPt, double higgsPt, bool passedVBF);

  void configure_OneJetLow();
  bool category_OneJetLow(unsigned NJets, double TauPt, bool passedVBF);

  void configure_OneJetBoost();
  bool category_OneJetBoost(unsigned NJets, double TauPt, double higgsPt, bool passedVBF);

  void configure_ZeroJetHigh();
  bool category_ZeroJetHigh(unsigned NJets, double TauPt);

  void configure_ZeroJetLow();
  bool category_ZeroJetLow(unsigned NJets, double TauPt);

  void configure_NoCategory();
  bool category_NoCategory();

  bool helperCategory_VBFLooseRelaxed_WYieldShape(bool useRelaxedForPlots, unsigned NJets, double DEta, int NJetsInGap, double Mjj);
  bool helperCategory_VBFTightRelaxed_WYield(bool useRelaxedForPlots, unsigned NJets, double DEta, int NJetsInGap, double Mjj, double higgsPt);

 private:
  // everything is in protected to be accessible by derived classes

};

#endif
