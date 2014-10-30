#ifndef ZtoEMu_Skim_h
#define ZtoEMu_Skim_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "ReferenceScaleFactors.h"

class ZtoEMu_Skim : public Selection {

 public:
  ZtoEMu_Skim(TString Name_, TString id_);
  virtual ~ZtoEMu_Skim();

  virtual void Configure();
  virtual void Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
		 NMu,
		 NE,
		 ptthreshold,
		 mll,
		 charge,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  ReferenceScaleFactors *RSF;

  double mu_ptlow,mu_pthigh,mu_eta,e_ptlow,e_pthigh,e_eta,jet_pt,jet_eta,singlejet,zmin,zmax,mtmu,ptbalance,mmin;
  int n_mu,n_e;
  bool doHiggsObjects;
  bool doWWObjects;
  bool useMadgraphZ;
  TString mucorr, ecorr, jetcorr;

  double csvl,csvm,csvt;

  double calculatePzeta(int muiterator, int eiterator);
  double calculatePzetaDQM(int muiterator, int eiterator);
  double cosphi2d(double px1, double py1, double px2, double py2);
  double cosphi3d(TVector3 vec1, TVector3 vec2);
  int findBin(TGraphAsymmErrors* graph, double xval);
  
  bool isFakeMuon(unsigned int idx, TString corr="");
  bool isFakeMuon(unsigned int idx, unsigned int vtx, TString corr="");
  bool isWWElectron(unsigned int idx, unsigned int vtx, TString corr="");
  bool isFakeElectron(unsigned int idx, TString corr="");
  bool isFakeElectron(unsigned int idx, unsigned int vtx, TString corr="");

  double ZPtReweight(double zpt);
  double PowhegReweight(double zpt);
  double ZPtRelUnc(double zpt);
  double ZPtMadgraphRelUnc(double zpt);

  double Fakerate(double pt, double eta, TH2D *fakeRateHist);
  double FakerateError(double pt, double eta, TH2D *fakeRateHist);
  
  TFile* FRFile;
  TFile* ZptCorrFile;
  TH1D* ZptCorrection;
  TH2D* ElectronFakeRate35;
  TH2D* ElectronFakeRate20;
  TH2D* ElectronFakeRate50;
  TH2D* MuonFakeRate15;
  TH2D* MuonFakeRate5;
  TH2D* MuonFakeRate30;
  double fakeRate;
  double fakeRateMu;
  double fakeRateE;

};
#endif

