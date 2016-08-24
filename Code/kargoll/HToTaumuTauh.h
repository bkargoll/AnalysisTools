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

  std::vector<TH1D> h_TauSelPt;

  std::vector<TH1D> h_MetPt;

  std::vector<TH1D> h_NJetsId;

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

  std::vector<TH1D> h_Tau3p_FlightLengthSignificance;

  std::vector<TH1D> h_Tau3p_Plus_E_Resol;
  std::vector<TH1D> h_Tau3p_Plus_E_AbsResol;
  std::vector<TH1D> h_Tau3p_Minus_E_Resol;
  std::vector<TH1D> h_Tau3p_Minus_E_AbsResol;
  std::vector<TH1D> h_Tau3p_Zero_E_Resol;
  std::vector<TH1D> h_Tau3p_Zero_E_AbsResol;
  std::vector<TH1D> h_Tau3p_True_E_Resol;
  std::vector<TH1D> h_Tau3p_True_E_AbsResol;
  std::vector<TH1D> h_Tau3p_False_E_Resol;
  std::vector<TH1D> h_Tau3p_False_E_AbsResol;

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

  // number of tau prongs: 1, 3, or 0 (i.e. all)
  int cTau_Prongs;
  // flight length significance cut for 3prong taus
  int cTau_flightLength;


  // object corrections to use
  TString correctTaus;
  TString correctMuons;
  TString correctElecs;
  TString correctJets;

  // systematic shifts
  TString tauESShift;

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
  double flightLengthSig_;

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
  TString svFitSuffix;

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
