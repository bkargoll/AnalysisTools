#ifndef ZToTaumuTauh_h
#define ZToTaumuTauh_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "ReferenceScaleFactors.h"

class ZToTaumuTauh : public Selection {

 public:
  ZToTaumuTauh(TString Name_, TString id_);
  virtual ~ZToTaumuTauh();

  virtual void  Configure();
  virtual void  Finish();

  enum cuts {TriggerOk=0,
	  PrimeVtx,
	  NMuId,
	  NMuKin,
	  NMuIso,
	  NTauId,
	  NTauKin,
	  NTauIso,
	  ChargeSum,
	  MT_MuMET,
	  TauDecayMode,
	  TauFLSigma,
	  NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

  ReferenceScaleFactors *RSF;

  // cut values
  double cMu_dxy, cMu_dz, cMu_relIso, cMu_pt, cMu_eta, cMu_dRHltMatch;
  double cTau_pt, cTau_eta, cMuTau_dR, cTau_IsoRaw, cTau_dRHltMatch;
  std::vector<TString> cTriggerNames;

  double OneProngNoPiWeight;
  //dummy values
  int TriggerOkDummy, selVertexDummy, selMuon_IsoDummy, selMuon_AntiIsoDummy, selTauDummy, ChargeSumDummy;
  double MTDummy, MvisDummy, TauFLSigmaDummy;

  int Charge;
  double SB_lowerLimit, SB_upperLimit;
  bool Scaleby_Counting;
  TString tau_corr;
  bool verbose, Use_Embedded;

  // Histograms

  std::vector<std::vector<TH1D>* > Extradist1d_OS, Extradist1d_SS;

  std::vector<TH1D> NVtx;
  std::vector<TH1D> NGoodVtx;
  std::vector<TH1D> NTrackperVtx;
  std::vector<TH1D> NMtTauMET;
  std::vector<TH1D> NMvis, Mvis_SignalOnly, Mvis_SignalOnly_genMu, Mvis_SignalOnly_genA1, Mvis_SignalOnly_genTaumu, Mvis_SignalOnly_genTauh;
  std::vector<TH1D> NSignal_SB_WJets;
  std::vector<TH1D> NSB_Data;
  std::vector<TH1D> NSB;
  std::vector<TH1D> Mu_pt, Mu_eta, Mu_phi, Tau_pt, Tau_eta, Tau_phi, MET_phi;
  std::vector<TH1D> TauFL, TauFLSigned, TauFLSigmaSigned, TauFLSigmaUnsigned;
  std::vector<TH1D> A1mass, A1mass10GeV;

  std::vector<TH1D> Mvis3Prong, Mvis1Prong, MvisIncl;
  std::vector<TH1D> MTMuMET3Prong, MTMuMET1Prong, MTMuMETIncl;

  std::vector<TH1D> dR_selTauh_genTauh, dR_selMu_genMu;
  std::vector<TH1D> POCAPV_Mag;
  std::vector<TH1D> Phi_SVPV, Phi_genTauh, Theta_SVPV, Theta_genTauh, dPhi_SVPV_genTauh, dTheta_SVPV_genTauh, Angle_SVPV_genTauh;
  std::vector<TH1D> Phi_POCAPV, Phi_genTaumu, Theta_POCAPV, Theta_genTaumu, dPhi_POCAPV_genTaumu, dTheta_POCAPV_genTaumu;
  std::vector<TH1D> dPhi_MinusSVPV_genTaumu, dTheta_MinusSVPV_genTaumu, Angle_MinusSVPV_genTaumu;
  std::vector<TH1D> GJ_Tauh, GJ_Taumu;
  std::vector<TH1D> dPhi_DiTauGen, Pt_DiTauGen, Pt_ZGen, M_ZGen, M_DiTauPtBalance, dM_DiTau, dPt_GenTaumuPtBalance, dP_GenTaumuPtBalance, dP_GenTauh;
  std::vector<TH2D> dP_GenTauMuPtBalance_vs_dPTauh, Pt_vs_dPhi_DiTauGen;
  std::vector<TH2D> TauFLSigmaCut_vs_Res, TauFLSigma_vs_Res;

  std::vector<TH1D> NQCD;
  std::vector<TH1D> QCD_MT_MuMET_A, QCD_MT_MuMET_B, QCD_MT_MuMET_C, QCD_MT_MuMET_D;
  std::vector<TH1D> QCD_MET_A, QCD_MET_B, QCD_MET_C, QCD_MET_D;

  bool selectMuon_Id(unsigned i, unsigned vertex);
  bool selectMuon_Kinematics(unsigned i);
  bool selectMuon_Isolation(unsigned i);
  bool selectMuon_AntiIsolation(unsigned i);

  bool selectPFTau_Id(unsigned i);
  bool selectPFTau_Id(unsigned i, std::vector<int> muonCollection);
  bool selectPFTau_Isolation(unsigned i);
  bool selectPFTau_Kinematics(unsigned i);
  double Reconstruct_hadronicTauEnergy(unsigned i);

 private:

};
#endif
