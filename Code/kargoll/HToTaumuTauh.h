#ifndef HToTaumuTauh_h
#define HToTaumuTauh_h

#include <TBenchmark.h>
#include <TH1.h>
#include <TString.h>
#include <cmath>
#include <vector>

#include "Selection.h"
#include "ReferenceScaleFactors.h"
#include "../DataFormats/SVFitObject.h"
#include "SVFitStorage.h"
#include "UncertaintyValue.h"

class TLorentzVector;
class TVector3;
class CategoryStandalone;

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

 protected:
  virtual void Setup();
  virtual void Store_ExtraDist();
  virtual void Finish();

  virtual void doEvent();

  virtual void doEventSetup();
  virtual void doSelection(bool runAnalysisCuts);
  virtual void doPlotting();

  // Histograms
  std::vector<TH1D> h_MuSelPt;
  std::vector<TH1D> h_MuSelEta;
  std::vector<TH1D> h_MuSelPhi;
  std::vector<TH1D> h_MuSelDxy;
  std::vector<TH1D> h_MuSelDz;
  std::vector<TH1D> h_MuSelRelIso;
  std::vector<TH1D> h_MuSelFakesTauID;
  std::vector<TH1D> h_MuSelDrHlt;

  std::vector<TH1D> h_TauSelPt;
  std::vector<TH1D> h_TauSelEta;
  std::vector<TH1D> h_TauSelPhi;
  std::vector<TH1D> h_TauSelDecayMode;
  std::vector<TH1D> h_TauSelIso;
  std::vector<TH1D> h_TauSelMass;

  std::vector<TH1D> h_MuCharge;
  std::vector<TH1D> h_TauCharge;

  std::vector<TH1D> h_MetPt;
  std::vector<TH1D> h_MetPhi;

  std::vector<TH1D> h_NJetsKin;
  std::vector<TH1D> h_NJetsId;
  std::vector<TH1D> h_Jet1Pt;
  std::vector<TH1D> h_Jet1Eta;
  std::vector<TH1D> h_Jet1Phi;
  std::vector<TH1D> h_Jet1IsB;
  std::vector<TH1D> h_Jet2Pt;
  std::vector<TH1D> h_Jet2Eta;
  std::vector<TH1D> h_Jet2Phi;
  std::vector<TH1D> h_Jet2IsB;

  std::vector<TH1D> h_NBJets;
  std::vector<TH1D> h_BJet1Pt;
  std::vector<TH1D> h_BJet1Eta;
  std::vector<TH1D> h_BJet1Phi;

  std::vector<TH1D> h_HiggsPt;
  std::vector<TH1D> h_HiggsPt_1JetSel;
  std::vector<TH1D> h_HiggsPhi;
  std::vector<TH1D> h_JetsDEta;
  std::vector<TH1D> h_JetsInEtaGap;
  std::vector<TH1D> h_JetsInvM;

  std::vector<TH1D> h_QcdShapeRegion;

  std::vector<TH1D> h_SVFitMass;
  std::vector<TH1D> h_visibleMass;

  std::vector<TH1D> h_shape_VisM;
  std::vector<TH1D> h_shape_SVfitM;

  std::vector<TH1D> h_shape_VisM_ZLScaleUp;
  std::vector<TH1D> h_shape_VisM_ZLScaleDown;
  std::vector<TH1D> h_shape_SVfitM_ZLScaleUp;
  std::vector<TH1D> h_shape_SVfitM_ZLScaleDown;

  std::vector<TH1D> h_shape_VisM_TauESUp;
  std::vector<TH1D> h_shape_VisM_TauESDown;
  std::vector<TH1D> h_shape_SVfitM_TauESUp;
  std::vector<TH1D> h_shape_SVfitM_TauESDown;

  std::vector<TH1D> h_SVFitTimeReal;
  std::vector<TH1D> h_SVFitTimeCPU;

  std::vector<TH1D> h_SVFitStatus;

  std::vector<TH1D> h_SVFitMassResol;
  std::vector<TH1D> h_visibleMassResol;

  std::vector<TH1D> h_BGM_Mt;
  std::vector<TH1D> h_BGM_MtSideband;
  std::vector<TH1D> h_BGM_MtExtrapolation;
  std::vector<TH1D> h_BGM_MtSS;
  std::vector<TH1D> h_BGM_MtSidebandSS;
  std::vector<TH1D> h_BGM_MtExtrapolationSS;
  std::vector<TH1D> h_BGM_MtSidebandInclusive;
  std::vector<TH1D> h_BGM_MtExtrapolationInclusive;
  std::vector<TH1D> h_BGM_MtSidebandSSInclusive;
  std::vector<TH1D> h_BGM_MtExtrapolationSSInclusive;
  std::vector<TH1D> h_BGM_MtAntiIso;
  std::vector<TH1D> h_BGM_MtAntiIsoSS;
  std::vector<TH1D> h_BGM_QcdAbcd;
  std::vector<TH1D> h_BGM_QcdAbcdInclusive;
  std::vector<TH1D> h_BGM_QcdOSMuIso;
  std::vector<TH1D> h_BGM_QcdOSTauIso;
  std::vector<TH1D> h_BGM_QcdSSMuIso;
  std::vector<TH1D> h_BGM_QcdSSTauIso;
  std::vector<TH1D> h_BGM_QcdEff;

  std::vector<TH1D> h_recMinusGenTauMass_recDecayModeEq0;
  std::vector<TH1D> h_recMinusGenTauMass_recDecayModeEq1;
  std::vector<TH1D> h_recMinusGenTauMass_recDecayModeEq10;
  std::vector<TH1D> h_recTauPtDivGenTauPt_recDecayModeEq0;
  std::vector<TH1D> h_recTauPtDivGenTauPt_recDecayModeEq1;
  std::vector<TH1D> h_recTauPtDivGenTauPt_recDecayModeEq10;


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
  // flag for data-driven QCD yield via efficiency method (only for categories where appropriate, for others ABCD method is used)
  bool qcdUseEfficiencyMethod;
  // flag to use embedding
  bool useEmbedding;


  // object corrections to use
  TString correctTaus;
  TString correctMuons;
  TString correctElecs;
  TString correctJets;

  // store if category has passed or not
  bool catPassed;

  // event variables
  double w; // event weight
  double wobs; // second event weight
  unsigned int t; // index of histogram
  int64_t id;
  int idStripped;

  // variables to hold selected objects (to be used for plotting, etc.)
  int selVertex;
  int selMuon;
  int selTau;
  double higgs_GenPtWeight;
  double higgs_GenPt;
  double higgsPhi;
  int diMuonNeg, diMuonPos;

  // booleans about event status
  bool status; // true only if full selection passed
  bool isWJetMC; // for Wjets background method
  bool isQCDShapeEvent, isWJetShapeEvent;
  bool isSignal;
  std::vector<bool> originalPass; // not affected by QCD/WJet shape methods

  // variables to hold information for categorization
  unsigned nJets_;
  double tauPt_;
  double higgsPt_;
  double jetdEta_;
  int nJetsInGap_;
  double mjj_;

  // selected objects at various stages
  std::vector<int> selectedMuonsId;
  std::vector<int> selectedMuons;
  std::vector<int> antiIsoMuons;
  std::vector<int> selectedTausId;
  std::vector<int> selectedTausIso;
  std::vector<int> selectedTaus;
  std::vector<int> diMuonVetoMuonsPositive;	// muons selected for the dimuon veto
  std::vector<int> diMuonVetoMuonsNegative;	// muons selected for the dimuon veto
  std::vector<int> relaxedIsoTaus;
  std::vector<int> triLepVetoMuons;
  std::vector<int> triLepVetoElecs;
  std::vector<int> selectedJetsClean;
  std::vector<int> selectedJetsKin;
  std::vector<int> selectedJets;
  std::vector<int> selectedLooseJets;
  std::vector<int> selectedBJets;

  // instance of reference scale factor class
  ReferenceScaleFactors* RSF;

  // instances of SVFitStorage class
  SVFitStorage svfitstorage;
  SVFitStorage svfitstorTauESUp;
  SVFitStorage svfitstorTauESDown;

  // timing information
  TBenchmark* clock;

  // booleans for different analysis stages
  enum passedStages{
	Vtx,
	VtxMuId,
	VtxMu,
	VtxTauIdIso,
	VtxTau,
	Objects,
	ObjectsDiMuonVeto,
	ObjectsFailDiMuonVeto,
	FullInclusiveSelNoMt,
	FullInclusiveSelNoMtNoOS,
	FullInclusiveNoTauNoMuNoCharge,
	FullInclusiveSelNoBVeto,
	FullInclusiveSel
  };
  bool getStatusBoolean(passedStages stage, const std::vector<bool>& passVec);
  bool getStatusBoolean(passedStages stage){return getStatusBoolean(stage, pass);};
  void resetPassBooleans();

  bool passed_VBFTightRelaxed;
  bool passed_VBFLooseRelaxed;

  bool hasRelaxedIsoTau, hasAntiIsoMuon;

  // relaxed categories for background methods
  std::vector<float> cut_VBFTightRelaxed, cut_VBFLooseRelaxed;

  // function definitions
  bool selectMuon_Id(unsigned i, unsigned vertex);
  bool selectMuon_Kinematics(unsigned i);

  bool selectMuon_antiIso(unsigned i, unsigned vertex);

  bool selectMuon_diMuonVeto(unsigned i, unsigned i_vtx);
  bool selectMuon_triLeptonVeto(int i, int selectedMuon, unsigned i_vtx);

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

  bool selectPFJet_Relaxed(unsigned i, int selectedMuon, int selectedTau);

  void calculateJetVariables(const std::vector<int>& jetCollection);

  bool migrateCategoryIntoMain(TString thisCategory, std::vector<float> categoryValueVector, std::vector<float> categoryPassVector, unsigned categoryNCuts);

  double getWJetsMCPrediction();
  UncDouble yield_DdBkg_WJets(int flag = Standard);
  void applyDdBkg_WJets();
  UncDouble yield_DdBkg_QCDAbcd(int flag = Standard);
  UncDouble yield_DdBkg_QCDEff();
  void applyDdBkg_QCD();

  enum bg_flag{
	  Standard = 1,		// Opposite sign, including category selection
	  SameSign,			// Same sign, including category selection
	  Inclusive,		// Opposite sign, without category selection
	  SameSignInclusive	// Same sign, without category selection
  };

 private:
  // everything is in protected to be accessible by derived classes

};

#endif
