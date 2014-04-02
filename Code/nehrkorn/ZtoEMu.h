#ifndef ZtoEMu_h
#define ZtoEMu_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"

class ZtoEMu : public Selection {

 public:
  ZtoEMu(TString Name_, TString id_);
  virtual ~ZtoEMu();

  virtual void Configure();
  virtual void Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
		 NMu,
		 NE,
		 ptthreshold,
		 drEMu,
		 diMuonVeto,
		 triLeptonVeto,
		 charge,
		 jetVeto,
		 MtMu,
	     ptBalance,
	     ZMassmax,
	     ZMassmin,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables

  std::vector<TH1D> RelIsoE;
  std::vector<TH1D> RelIsoMu;
  std::vector<TH1D> EPt;
  std::vector<TH1D> MuPt;
  std::vector<TH1D> mtMu;
  std::vector<TH1D> mtE;
  std::vector<TH1D> etaMu;
  std::vector<TH1D> etaE;
  std::vector<TH1D> jetsum;
  std::vector<TH1D> NJets;
  std::vector<TH1D> NJetsLoose;
  std::vector<TH1D> NJetsMedium;
  std::vector<TH1D> NJetsTight;
  std::vector<TH1D> NJetsOwn;
  std::vector<TH1D> PUJetId;
  std::vector<TH1D> chargesum;
  std::vector<TH1D> drmue;
  std::vector<TH1D> deltaphi;
  std::vector<TH1D> ptbal;
  std::vector<TH1D> chargesumsigned;
  std::vector<TH1D> FirstJetPt;
  std::vector<TH1D> SecondJetPt;
  std::vector<TH1D> jeteta;
  
  std::vector<TH1D> invmass_zmass;
  std::vector<TH1D> invmass_ptbalance;
  std::vector<TH1D> invmass_mtmu;
  std::vector<TH1D> invmass_jetveto;
  std::vector<TH1D> invmass_vetos;
  std::vector<TH1D> invmass_only_object_id;
  
  std::vector<TH1D> nm0_met;
  std::vector<TH1D> nm0_jetsum;
  std::vector<TH1D> nm0_onejet;
  std::vector<TH1D> nm0_mtmu;
  std::vector<TH1D> nm0_ptbalance;
  
  std::vector<TH1D> NPV;
  std::vector<TH1D> num_interactions;
  std::vector<TH1D> evtweight;
  
  std::vector<TH1D> met;
  std::vector<TH1D> met_xycorr;
  std::vector<TH1D> met_uncorr;
  std::vector<TH1D> onejet;
  std::vector<TH1D> mte_mtmu;
  std::vector<TH1D> leadingjet;
  std::vector<TH1D> subleadingjet;
  std::vector<TH1D> sumjets;
  std::vector<TH1D> NbJets;
  std::vector<TH1D> NbJetsVtxL;
  std::vector<TH1D> NbJetsVtxM;
  std::vector<TH1D> NbJetsVtxT;

  // comparison of generators

  std::vector<TH1D> zpt;
  std::vector<TH1D> zeta;
  std::vector<TH1D> zmass;
  std::vector<TH1D> leadingjet_pt;
  std::vector<TH1D> subleadingjet_pt;
  std::vector<TH1D> leadingjet_eta;
  std::vector<TH1D> subleadingjet_eta;
  std::vector<TH1D> jetsumcustom;

  double mu_ptlow,mu_pthigh,mu_eta,e_ptlow,e_pthigh,e_eta,jet_pt,jet_eta,jet_sum,zmin,zmax,mtmu,ptbalance;
  int n_mu,n_e;
  bool doHiggsObjects;
  bool doWWObjects;
  bool doOldJetVeto;
  
  double csvl,csvm,csvt;

  double calculatePzeta(int muiterator, int eiterator);
  double calculatePzetaDQM(int muiterator, int eiterator);
  double cosphi2d(double px1, double py1, double px2, double py2);
  double cosphi3d(TVector3 vec1, TVector3 vec2);
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  bool jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx);
  bool isGoodVtx(unsigned int i);
  double vertexSignificance(TVector3 vec, unsigned int vertex);
  bool matchTrigger(unsigned int i, double dr, std::string trigger, std::string object);
  int matchTruth(TLorentzVector tvector);
  bool matchTruth(TLorentzVector tvector, int pid, double dr);
  int findBin(TGraphAsymmErrors* graph, double xval);
  
  bool isTightMuon(unsigned int idx);
  bool isTightMuon(unsigned int idx, unsigned int vtx);
  bool isHiggsMuon(unsigned int idx, unsigned int vtx);
  bool isLooseMuon(unsigned int idx);
  bool isFakeMuon(unsigned int idx);
  bool isFakeMuon(unsigned int idx, unsigned int vtx);
  double Muon_RelIso(unsigned int idx);
  
  bool isTrigPreselElectron(unsigned int idx);
  bool isTrigNoIPPreselElectron(unsigned int idx);
  bool isMVATrigElectron(unsigned int idx);
  bool isMVATrigNoIPElectron(unsigned int idx);
  bool isMVANonTrigElectron(unsigned int idx, unsigned int vtx);
  bool isHiggsElectron(unsigned int idx, unsigned int vtx);
  bool isWWElectron(unsigned int idx, unsigned int vtx);
  bool isTightElectron(unsigned int idx);
  bool isTightElectron(unsigned int idx, unsigned int vtx);
  bool isLooseElectron(unsigned int idx);
  bool isFakeElectron(unsigned int idx);
  bool isFakeElectron(unsigned int idx, unsigned int vtx);
  double Electron_RelIso(unsigned int idx);
  double Electron_Aeff_R04(double Eta);
  double Electron_Aeff_R03(double Eta);
  
  double MuonIDeff(unsigned int idx);
  double MuonIDerrUp(unsigned int idx);
  double MuonIDerrDown(unsigned int idx);
  double MuonHiggsIDeff(unsigned int idx);
  double MuonTriggerEff(unsigned int idx);
  double MuonTriggerErr(unsigned int idx);
  double ElectronIDeff(unsigned int idx, std::string id);
  double ElectronIDerr(unsigned int idx, std::string id);
  double ElectronTrigIDeff(unsigned int idx);
  double ElectronTrigIDerr(unsigned int idx);
  double ElectronNonTrigIDeff(unsigned int idx);
  double ElectronNonTrigIDerr(unsigned int idx);
  double ElectronHiggsIDeff(unsigned int idx);
  double ElectronTriggerEff(unsigned int idx);
  double ElectronTriggerErr(unsigned int idx);
  
  double TriggerEff(unsigned int muid, unsigned int eid, TString path);
  double SingleEle(unsigned int idx);
  double DoubleEleLeading(unsigned int idx);
  double DoubleEleTrailing(unsigned int idx);
  double SingleMu(unsigned int idx);
  double DoubleMuLeading(unsigned int idx);
  double DoubleMuTrailing(unsigned int idx);

  double ElectronMassScale(unsigned int idx);
  double ZPtReweight(double zpt);
  double rundependentJetPtCorrection(double jeteta, int runnumber);

  //double JECuncertainty(unsigned int i, TString datamc);

  double Fakerate(TLorentzVector vec, TH2D *fakeRateHist, std::string type);
  double FakerateWW(unsigned int idx, std::string type);
  
  TFile* FRFile;
  TFile* EmbEffFile;
  TH2D* ElectronFakeRate;
  TH2D* MuonFakeRate;
  TH2D* EmbEff;
  double fakeRate;
  double fakeRateMu;
  double fakeRateE;
  
  TFile* MuIdEffFile;
  TFile* MuIsoEffFile;
  TFile* ETrigIdEffFile;
  TFile* ENonTrigIdEffFile;
  TFile* TriggerEfficiencies;
  TFile* FakeRates;

  TH2D* ElectronTrigEff;
  TH2D* ElectronNonTrigEff;
  TGraphAsymmErrors* MuIdEff09;
  TGraphAsymmErrors* MuIdEff12;
  TGraphAsymmErrors* MuIdEff21;
  TGraphAsymmErrors* MuIdEff24;
  TGraphAsymmErrors* MuIsoEff09;
  TGraphAsymmErrors* MuIsoEff12;
  TGraphAsymmErrors* MuIsoEff21;
  TGraphAsymmErrors* MuIsoEff24;

  TGraphAsymmErrors* SingleEle15;
  TGraphAsymmErrors* SingleEle25;
  TGraphAsymmErrors* DoubleEleLead15;
  TGraphAsymmErrors* DoubleEleLead25;
  TGraphAsymmErrors* DoubleEleTrail15;
  TGraphAsymmErrors* DoubleEleTrail25;
  TGraphAsymmErrors* SingleMu08;
  TGraphAsymmErrors* SingleMu12;
  TGraphAsymmErrors* SingleMu21;
  TGraphAsymmErrors* SingleMu25;
  TGraphAsymmErrors* DoubleMuLead12;
  TGraphAsymmErrors* DoubleMuLead21;
  TGraphAsymmErrors* DoubleMuLead25;
  TGraphAsymmErrors* DoubleMuTrail12;
  TGraphAsymmErrors* DoubleMuTrail21;
  TGraphAsymmErrors* DoubleMuTrail25;

  TGraphAsymmErrors* EleFake1;
  TGraphAsymmErrors* EleFake15;
  TGraphAsymmErrors* EleFake2;
  TGraphAsymmErrors* EleFake25;
  TGraphAsymmErrors* MuFake1;
  TGraphAsymmErrors* MuFake15;
  TGraphAsymmErrors* MuFake2;
  TGraphAsymmErrors* MuFake25;

  TF1* gause;
  TF1* gausmu;
  TH1D* eres;
  TH1D* mures;
  double eleres;
  double muonres;

};
#endif

