#ifndef ZtoEMu_Skim_h
#define ZtoEMu_Skim_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "TGraphAsymmErrors.h"

class ZtoEMu_Skim : public Selection {

 public:
  ZtoEMu_Skim(TString Name_, TString id_);
  virtual ~ZtoEMu_Skim();

  virtual void  Configure();
  virtual void Finish();

  enum cuts {TriggerOk=0,
	     PrimeVtx,
		 NMu,
		 NE,
		 ptthreshold,
	     NCuts};

 protected:
  virtual void doEvent();
  virtual void Store_ExtraDist();

 private:
  // Selection Variables
  
  std::vector<TH1D> NPV;
  std::vector<TH1D> mupt;
  std::vector<TH1D> mueta;
  std::vector<TH1D> ept;
  std::vector<TH1D> eeta;

  double mu_ptlow,mu_pthigh,mu_eta,e_ptlow,e_pthigh,e_eta;
  int n_mu,n_e;

  bool doHiggsObjects;
  bool doWWObjects;

  double cosphi2d(double px1, double py1, double px2, double py2);
  double dxy(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  double dz(TLorentzVector fourvector, TVector3 poca, TVector3 vtx);
  bool jetFromVtx(std::vector<int> vtx_track_idx, int leadingtrack_idx);
  bool isGoodVtx(unsigned int i);
  double vertexSignificance(TVector3 vec, unsigned int vertex);
  bool matchTrigger(unsigned int i, double dr, std::string trigger, std::string object);
  int matchTruth(TLorentzVector tvector);
  bool matchTruth(TLorentzVector tvector, int pid, double dr);
  
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

};
#endif

