#ifndef MuTauSync_h
#define MuTauSync_h

#include "Selection.h"
#include <vector>
#include "TString.h"
#include "HToTaumuTauh.h"

class MuTauSync: public HToTaumuTauh {

public:
	MuTauSync(TString Name_, TString id_);
	virtual ~MuTauSync();

	virtual void Configure();

protected:
	virtual void doEvent();
	virtual void Finish();

private:
	// tree to save synch nTuple
	TFile* syncFile;
	TTree* syncTree_VBFHiggs;
	TTree* syncTree_TauPlusXA;

	void defineBranches(TTree* tree);

	// member variables to hold synch nTuple values
	// see here for definition: https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/HiggsTauTau/src/HTTSync.cc#L48-L201
	unsigned int run;
	unsigned int lumi;
	unsigned int evt;
	 // Event Variables
	int npv;
	float npu;
	float rho;
	 // Event Weights
	float mcweight;
	float puweight;
	float  trigweight_1;
	float trigweight_2;
	float idweight_1;
	float idweight_2;
	float isoweight_1;
	float isoweight_2;
	float fakeweight;
	float effweight;
	float weight;
	float embeddedWeight;
	float signalWeight;
	 // SV Fit variables
	float mvis;
	float m_sv;
	float pt_sv;
	float eta_sv;
	float phi_sv;
	float m_sv_Up;
	float m_sv_Down;
	 // First lepton : muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
	float pt_1;
	float phi_1;
	float eta_1;
	float m_1;
	int q_1;
	float iso_1;
	float mva_1;
	float d0_1;
	float dZ_1;
	bool passid_1;
	bool passiso_1;
	float mt_1;
	 // Second lepton : hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT) Tau for Tau Tau
	float pt_2;
	float phi_2;
	float eta_2;
	float m_2;
	int q_2;
	float iso_2;
	float d0_2;
	float dZ_2;
	float pt_tt;
	float byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
	float againstElectronMVA3raw_2;
	float  byIsolationMVA2raw_2;
	float againstMuonLoose2_2;
	float againstMuonMedium2_2;
	float againstMuonTight2_2;
	float mva_2;
	bool passid_2;
	bool passiso_2;
	float mt_2;
	 // Met related variables
	float met;
	float metphi;
	float l1met;
	float l1metphi;
	float l1metcorr;
	float calomet;
	float calometphi;
	float calometcorr;
	float calometphicorr;
	float mvamet;
	float mvametphi;
	float pzetavis;
	float pzetamiss;
	// met covariance matrices
	float metcov00;
	float metcov01;
	float metcov10;
	float metcov11;
	// mva met covariance matrices
	float mvacov00;
	float mvacov01;
	float mvacov10;
	float mvacov11;
	// First jet: leading jet after applying Jet energy corrections (excluding hadronic Tau)
	float jpt_1;
	float jeta_1;
	float jphi_1;
	float jptraw_1;
	float jptunc_1;
	float jmva_1;
	float jlrm_1;
	int jctm_1;
	int jpass_1;
	 // Second Jet : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
	float jpt_2;
	float jeta_2;
	float jphi_2;
	float jptraw_2;
	float jptunc_2;
	float jmva_2;
	float jlrm_2;
	int jctm_2;
	int jpass_2;
	 // B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
	float bpt;
	float beta;
	float bphi;
	 // Di Jet kinematic variables for VBF selection ==> Two leading pT Jets
	float mjj;
	float jdeta;
	int njetingap;
	float  mva;
	 // variables that go into the VBF MVA
	float  jdphi;
	float dijetpt;
	float dijetphi;
	float  hdijetphi;
	float  visjeteta;
	float  ptvis;
	 // number of btags passing btag id (pt > 20)
	int nbtag;
	 // number of jets passing jet id ( pt > 30 )
	int njets;
	int njetspt20;
	 // mva output for e+mu channel
	float mva_gf;
	float mva_vbf;


};
#endif
