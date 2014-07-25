#ifndef HToTaumuTauh_h
#define HToTaumuTauh_h

#include <TH1.h>
#include <TString.h>
#include <cmath>
#include <vector>

#include "../Selection.h"

class TLorentzVector;
class TVector3;

// small struct needed to allow sorting indices by some value
struct sortIdxByValue {
    bool operator()(const std::pair<int,double> &left, const std::pair<int,double> &right) {
        return left.second > right.second;
    }
};

class HToTaumuTauh : public Selection {

 public:
  HToTaumuTauh(TString Name_, TString id_);
  virtual ~HToTaumuTauh();

  virtual void  Configure();

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
  virtual void doEvent();
  virtual void Store_ExtraDist();
  virtual void  Finish();

  // Selection Variables
  std::vector<TH1D> NCatFired;
  std::vector<TH1D> CatFired;

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NVtxFullSelection;
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
  std::vector<TH1D> MuSelFakesTauID;
  std::vector<TH1D> MuSelDrHlt;

  std::vector<TH1D> TauPt;
  std::vector<TH1D> TauEta;
  std::vector<TH1D> TauPhi;

  std::vector<TH1D> TauSelPt;
  std::vector<TH1D> TauSelEta;
  std::vector<TH1D> TauSelPhi;
  std::vector<TH1D> TauSelDrHlt; // todo: not filled at the moment
  std::vector<TH1D> TauSelDecayMode;

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

  std::vector<TH1D> MetPt;
  std::vector<TH1D> MetPhi;

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

  std::vector<TH1D> Cat0JetLowMt;
  std::vector<TH1D> Cat0JetLowMtSideband;
  std::vector<TH1D> Cat0JetLowMtExtrapolation;
  std::vector<TH1D> Cat0JetHighMt;
  std::vector<TH1D> Cat0JetHighMtSideband;
  std::vector<TH1D> Cat0JetHighMtExtrapolation;
  std::vector<TH1D> Cat1JetLowMt;
  std::vector<TH1D> Cat1JetLowMtSideband;
  std::vector<TH1D> Cat1JetLowMtExtrapolation;
  std::vector<TH1D> Cat1JetHighMt;
  std::vector<TH1D> Cat1JetHighMtSideband;
  std::vector<TH1D> Cat1JetHighMtExtrapolation;
  std::vector<TH1D> Cat1JetBoostMt;
  std::vector<TH1D> Cat1JetBoostMtSideband;
  std::vector<TH1D> Cat1JetBoostMtExtrapolation;
  std::vector<TH1D> CatVBFLooseMt;
  std::vector<TH1D> CatVBFLooseMtSideband;
  std::vector<TH1D> CatVBFLooseRelaxMt;
  std::vector<TH1D> CatVBFLooseRelaxMtExtrapolation;
  std::vector<TH1D> CatVBFTightMt;
  std::vector<TH1D> CatVBFTightMtSideband;
  std::vector<TH1D> CatVBFTightRelaxMt;
  std::vector<TH1D> CatVBFTightRelaxMtExtrapolation;
  std::vector<TH1D> CatInclusiveMt;
  std::vector<TH1D> CatInclusiveMtSideband;
  std::vector<TH1D> CatInclusiveMtExtrapolation;

  std::vector<TH2D> Cat0JetLowQcdAbcd;
  std::vector<TH2D> Cat0JetHighQcdAbcd;
  std::vector<TH2D> Cat1JetLowQcdAbcd;
  std::vector<TH2D> Cat1JetHighQcdAbcd;
  std::vector<TH2D> Cat1JetBoostQcdAbcd;
  std::vector<TH2D> CatVBFLooseQcdAbcd;
  std::vector<TH2D> CatVBFTightQcdAbcd;
  std::vector<TH2D> CatInclusiveQcdAbcd;

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

  // map to hold WJets yields for each category
  std::map<TString, double> wJetsYieldMap;
  //std::map<TString, double> wJetsYieldScaleMap;

  // variables to hold selected objects (to be used e.g. for sync Ntuple)
  int selVertex;
  int selMuon;
  int selTau;
  std::vector<int> selJets, selBJets;
  double selMjj, selJetdeta;
  int selNjetingap;


  // function definitions
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);

  double matchTrigger(unsigned int i_obj, std::vector<TString> trigger, std::string objectType);

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

  std::vector<int> sortPFjets();
  bool selectPFJet_Cleaning(unsigned i, int selectedMuon, int selectedTau);
  bool selectPFJet_Kinematics(unsigned i);
  bool selectPFJet_Id(unsigned i);
  bool selectBJet(unsigned i, int selectedMuon, int selectedTau);

  inline double transverseMass(double pt1, double phi1, double pt2, double phi2){
	  return sqrt(2 * pt1 * pt2 * (1 - cos(phi1 - phi2)));
  }

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

  bool helperCategory_VBFLooseRelaxed(bool useRelaxedForPlots, unsigned NJets, double DEta, int NJetsInGap, double Mjj);
  bool helperCategory_VBFTightRelaxed(bool useRelaxedForPlots, unsigned NJets, double DEta, int NJetsInGap, double Mjj, double higgsPt);

 private:
  // everything is in protected to be accessible by derived classes

};

#endif
