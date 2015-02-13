#ifndef Tvariable_Base_h
#define Tvariable_Base_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "ReferenceScaleFactors.h"
#include "PDFweights.h"

class Tvariable_Base : public Selection {

 public:
	Tvariable_Base(TString Name_, TString id_);
  virtual ~Tvariable_Base();

  virtual void Configure();
  virtual void Finish();

  enum cuts {TriggerOk=0,
	  	 PrimeVtx,
		 NMu,
		 NE,
		 ptthreshold,
		 mll,
		 triLeptonVeto,
		 charge,
		 oneJet,
		 MET,
	     ptBalance,
	     ZMassmax,
	     ZMassmin,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();
  TString categoryFlag;

  // Selection Variables
  ReferenceScaleFactors *RSF;

  std::vector<TH1D> LeadPt;
  std::vector<TH1D> TrailPt;
  std::vector<TH1D> etaLead;
  std::vector<TH1D> etaTrail;
  std::vector<TH1D> NJets;
  std::vector<TH1D> NJetsLoose;
  std::vector<TH1D> NJetsMedium;
  std::vector<TH1D> NJetsTight;
  std::vector<TH1D> PUJetId;
  std::vector<TH1D> drleadtrail;
  std::vector<TH1D> deltaphi;
  std::vector<TH1D> ptbal;

  std::vector<TH1D> invmass_zmass;
  std::vector<TH1D> invmass_ptbalance;
  std::vector<TH1D> invmass_met;
  std::vector<TH1D> invmass_jetveto;
  std::vector<TH1D> invmass_vetos;
  std::vector<TH1D> invmass_only_object_id;

  std::vector<TH1D> invmass_zmass_m;
  std::vector<TH1D> invmass_ptbalance_m;
  std::vector<TH1D> invmass_met_m;
  std::vector<TH1D> invmass_jetveto_m;
  std::vector<TH1D> invmass_vetos_m;
  std::vector<TH1D> invmass_only_object_id_m;

  std::vector<TH1D> invmass_trilepton_only;
  std::vector<TH1D> invmass_charge_only;
  std::vector<TH1D> invmass_jetveto_only;
  std::vector<TH1D> invmass_met_only;
  std::vector<TH1D> invmass_ptbal_only;

  std::vector<TH1D> nm0_met;
  std::vector<TH1D> nm0_mvamet;
  std::vector<TH1D> nm0_onejet;
  std::vector<TH1D> nm0_ptbalance;
  std::vector<TH1D> nmm_met;
  std::vector<TH1D> nmm_mvamet;
  std::vector<TH1D> nmm_onejet;
  std::vector<TH1D> nmm_ptbalance;

  std::vector<TH1D> NPV;
  std::vector<TH1D> NPV3d;
  std::vector<TH1D> NPVfine;

  std::vector<TH1D> met_corr;
  std::vector<TH1D> met_uncorr;
  std::vector<TH1D> onejet;
  std::vector<TH1D> onejet_eta;
  std::vector<TH1D> NbJets;
  std::vector<TH1D> NbJetsVtxL;
  std::vector<TH1D> NbJetsVtxM;
  std::vector<TH1D> NbJetsVtxT;

  std::vector<TH1D> ptbal_zoom;
  std::vector<TH1D> ht_pseudo;
  std::vector<TH1D> zmass_zoom;
  std::vector<TH1D> invmass_high;
  std::vector<TH1D> ptsum;
  std::vector<TH1D> ptsum_nm0;
  std::vector<TH1D> mvamet;
  std::vector<TH1D> invmass_ptbalance_widerange;

  std::vector<TH1D> pdf_w0;
  std::vector<TH1D> pdf_w1;

  double mu_ptlow,e_ptlow,mu_eta,e_eta,lead_pt,trail_pt,lead_eta,trail_eta,mmin,jet_pt,jet_eta,singlejet,met,ptbalance,zmin,zmax;
  double csvl,csvm,csvt;
  double normunc_dy,normunc_tt,normunc_tw,normunc_diboson,normunc_qcd;
  int n_mu,n_e;
  bool doWWObjects;
  bool useMadgraphZ;
  bool doPDFuncertainty;
  bool doTriggerUncertainty,doPileupUncertainty;
  bool doElectronIdUncertainty,doElectronScaleUncertainty,doElectronResUncertainty;
  bool doMuonIdUncertainty,doMuonScaleUncertainty,doMuonResUncertainty;
  bool doJECUncertainty, doJERUncertainty;
  bool doFakeRateUncertainty;
  bool doMetUncertainty;
  bool upwardUncertainty,systValid;
  TString mucorr, ecorr, jetcorr;

  TString pdfname1;
  TString pdfname2;
  PDFweights* pdf;
  int nPDFmembers;

  double calculatePzeta(int muiterator, int eiterator);
  double calculatePzetaDQM(int muiterator, int eiterator);
  double cosphi2d(double px1, double py1, double px2, double py2);
  double cosphi3d(TVector3 vec1, TVector3 vec2);
  int findBin(TGraphAsymmErrors* graph, double xval);

  bool isFakeMuon(unsigned int idx);
  bool isFakeMuon(unsigned int idx, unsigned int vtx);
  bool isWWElectron(unsigned int idx, unsigned int vtx);
  bool isFakeElectron(unsigned int idx);
  bool isFakeElectron(unsigned int idx, unsigned int vtx);

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
  double fakeRateMu1, fakeRateMu2;
  double fakeRateE1, fakeRateE2;

 private:
  // everything is protected for inheritance reasons
};
#endif

