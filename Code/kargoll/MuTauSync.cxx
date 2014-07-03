#include "MuTauSync.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

MuTauSync::MuTauSync(TString Name_, TString id_):
  HToTaumuTauh(Name_,id_)
{
	std::cout << "Setting up the class MuTauSync" << std::endl;
	// always run without category for sync exercise
	categoryFlag = "NoCategory";
}

MuTauSync::~MuTauSync(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "MuTauSync::~MuTauSync Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "MuTauSync::~MuTauSync()" << std::endl;
}

void  MuTauSync::Configure(){
	HToTaumuTauh::Configure();

	// implement changes between main analysis and sync exercise
	cut.at(OppCharge) = 999; // set to 999 to disable opp. charge cut
	title.at(OppCharge) += " DISABLED";
	cut.at(MT) = 999; // set to 999 to disable mt cut
	title.at(MT) += " DISABLED";

	// set up synch ntuple
	syncFile = new TFile("MuTauSyncTree.root", "RECREATE");
	syncTree = new TTree("syncTree", "syncTree");
	defineBranches();
}

void MuTauSync::doEvent(){
	HToTaumuTauh::doEvent();

	//std::cout << "Fill the sync tree!" << std::endl;
	//////// fill synch tree
	bool passSynchSelection = pass.at(TriggerOk) && pass.at(PrimeVtx) && pass.at(NMuId) && pass.at(NMuKin);
	passSynchSelection = passSynchSelection && pass.at(DiMuonVeto) && pass.at(NTauId) && pass.at(NTauIso) && pass.at(NTauKin) && pass.at(TriLeptonVeto);

	//std::cout << "Did this event pass? " << passSynchSelection << std::endl;
	//std::cout << "Selected Tau is " << selTau << std::endl;

	if (passSynchSelection){
		run = Ntp->RunNumber();
		lumi = Ntp->LuminosityBlock();
		evt = Ntp->EventNumber();
		 // Event Variables
		npv = Ntp->NVtx();
		npu = Ntp->PileupInfo_TrueNumInteractions_n0();
		rho = Ntp->RhoIsolationAllInputTags();
		 // Event Weights
		mcweight = -10;
		puweight = Ntp->PUWeightFineBins();
		trigweight_1 = -10;
		trigweight_2 = -10;
		idweight_1 = -10;
		idweight_2 = -10;
		isoweight_1 = -10;
		isoweight_2 = -10;
		fakeweight = -10;
		effweight = -10;
		weight = -10;
		embeddedWeight = -10;
		signalWeight = -10;
		 // SV Fit variables
		mvis = (Ntp->Muon_p4(selMuon) + Ntp->PFTau_p4(selTau)).M();
		m_sv = -10;
		pt_sv = -10;
		eta_sv = -10;
		phi_sv = -10;
		m_sv_Up = -10;
		m_sv_Down = -10;
		 // First lepton : muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
		pt_1 = Ntp->Muon_p4(selMuon).Pt();
		phi_1 = Ntp->Muon_p4(selMuon).Phi();
		eta_1 = Ntp->Muon_p4(selMuon).Eta();
		m_1 = Ntp->Muon_p4(selMuon).M();
		q_1 = Ntp->Muon_Charge(selMuon);
		iso_1 = Ntp->Muon_RelIso(selMuon);
		mva_1 = 0; // to be set to 0 for muons
		d0_1 = Ntp->dxy(Ntp->Muon_p4(selMuon),Ntp->Muon_Poca(selMuon),Ntp->Vtx(selVertex));
		dZ_1 = Ntp->dz(Ntp->Muon_p4(selMuon),Ntp->Muon_Poca(selMuon),Ntp->Vtx(selVertex));
		passid_1 = true; // passes obviously, as it was selected
		passiso_1 = true;
		mt_1 = transverseMass(pt_1,phi_1,Ntp->MET_CorrMVAMuTau_et(),Ntp->MET_CorrMVAMuTau_phi());
		 // Second lepton : hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT) Tau for Tau Tau
		pt_2 = Ntp->PFTau_p4(selTau).Pt();
		phi_2 = Ntp->PFTau_p4(selTau).Phi();
		eta_2 = Ntp->PFTau_p4(selTau).Eta();
		m_2 = Ntp->PFTau_p4(selTau).M();
		q_2 = Ntp->PFTau_Charge(selTau);
		iso_2 = -10; //MVA iso for hadronic Tau
		d0_2 = Ntp->dxy(Ntp->PFTau_p4(selTau),Ntp->PFTau_Poca(selTau),Ntp->Vtx(selVertex));
		dZ_2 = Ntp->dz(Ntp->PFTau_p4(selTau),Ntp->PFTau_Poca(selTau),Ntp->Vtx(selVertex));
		pt_tt = -10; //(Ntp->Muon_p4(selMuon) + Ntp->PFTau_p4(selTau) + met_p4).Pt()
		byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = Ntp->PFTau_HPSPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits(selTau);
		againstElectronMVA3raw_2 = -10;
		 byIsolationMVA2raw_2 = -10;
		againstMuonLoose2_2 = Ntp->PFTau_isHPSAgainstMuonLoose2(selTau);
		againstMuonMedium2_2 = Ntp->PFTau_isHPSAgainstMuonMedium2(selTau);
		againstMuonTight2_2 = Ntp->PFTau_isHPSAgainstMuonTight2(selTau);
		mva_2 = -10;
		passid_2 = true;
		passiso_2 = true;
		mt_2 = transverseMass(pt_2,phi_2,Ntp->MET_CorrMVAMuTau_et(),Ntp->MET_CorrMVAMuTau_phi());
		 // Met related variables
		met = Ntp->MET_Uncorr_et();
		metphi = Ntp->MET_Uncorr_phi();
		l1met = -10;
		l1metphi = -10;
		l1metcorr = -10;
		calomet = -10; // raw calo met
		calometphi = -10;
		calometcorr = Ntp->MET_CorrCaloT1T2_et();
		calometphicorr = Ntp->MET_CorrCaloT1T2_phi();
		mvamet = Ntp->MET_CorrMVAMuTau_et();
		mvametphi = Ntp->MET_CorrMVAMuTau_phi();
		pzetavis = -10;
		pzetamiss = -10;
		// met covariance matrices
		metcov00 = Ntp->MET_Uncorr_significance_xx();
		metcov01 = Ntp->MET_Uncorr_significance_xy();
		metcov10 = Ntp->MET_Uncorr_significance_xy();
		metcov11 = Ntp->MET_Uncorr_significance_yy();
		// mva met covariance matrices
		mvacov00 = Ntp->MET_CorrMVAMuTau_significance_xx();
		mvacov01 = Ntp->MET_CorrMVAMuTau_significance_xy();
		mvacov10 = Ntp->MET_CorrMVAMuTau_significance_xy();
		mvacov11 = Ntp->MET_CorrMVAMuTau_significance_yy();
		// First jet: leading jet after applying Jet energy corrections (excluding hadronic Tau)
		if (selJets.size() > 0){
			jpt_1 = Ntp->PFJet_p4(selJets.at(0)).Pt();
			jeta_1 = Ntp->PFJet_p4(selJets.at(0)).Eta();
			jphi_1 = Ntp->PFJet_p4(selJets.at(0)).Phi();
			jptraw_1 = -10;
			jptunc_1 = -10;
			jmva_1 = Ntp->PFJet_PUJetID_discr(selJets.at(0));
			jlrm_1 = -10;
			jctm_1 = -10;
			jpass_1 = int(Ntp->PFJet_PUJetID_looseWP(selJets.at(0)));
		}
		else {
			jpt_1 = -20;
			jeta_1 = -20;
			jphi_1 = -20;
			jptraw_1 = -20;
			jptunc_1 = -20;
			jmva_1 = -20;
			jlrm_1 = -20;
			jctm_1 = -20;
			jpass_1 = -20;
		}
		 // Second Jet : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
		if (selJets.size() > 1){
			jpt_2 = Ntp->PFJet_p4(selJets.at(1)).Pt();
			jeta_2 = Ntp->PFJet_p4(selJets.at(1)).Eta();
			jphi_2 = Ntp->PFJet_p4(selJets.at(1)).Phi();
			jptraw_2 = -10;
			jptunc_2 = -10;
			jmva_2 = Ntp->PFJet_PUJetID_discr(selJets.at(1));
			jlrm_2 = -10;
			jctm_2 = -10;
			jpass_2 = int(Ntp->PFJet_PUJetID_looseWP(selJets.at(1)));
		}
		else {
			jpt_2 = -20;
			jeta_2 = -20;
			jphi_2 = -20;
			jptraw_2 = -20;
			jptunc_2 = -20;
			jmva_2 = -20;
			jlrm_2 = -20;
			jctm_2 = -20;
			jpass_2 = -20;
		}
		 // B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
		if (selBJets.size() > 0){
			bpt = Ntp->PFJet_p4(selBJets.at(0)).Pt();
			beta = Ntp->PFJet_p4(selBJets.at(0)).Eta();
			bphi = Ntp->PFJet_p4(selBJets.at(0)).Phi();
		}
		else{
			bpt = -20;
			beta = -20;
			bphi = -20;
		}
		 // Di Jet kinematic variables for VBF selection ==> Two leading pT Jets
		mjj = selMjj;
		jdeta = fabs(selJetdeta);
		njetingap = selNjetingap;
		 mva = -10;
		 // variables that go into the VBF MVA
		 jdphi = -10;
		dijetpt = -10;
		dijetphi = -10;
		 hdijetphi = -10;
		 visjeteta = -10;
		 ptvis = -10;
		 // number of btags passing btag id (pt > 20)
		nbtag = selBJets.size();
		 // number of jets passing jet id ( pt > 30 )
		njets = selJets.size();
		njetspt20 = -10;
		 // mva output for e+mu channel
		mva_gf = -10;
		mva_vbf = 10;


		syncTree->Fill();
	}
}

void MuTauSync::Finish(){
	HToTaumuTauh::Finish();

	syncTree->SetDirectory(syncFile);
	syncTree->GetCurrentFile()->Write();
	syncFile->Close();
}

void MuTauSync::defineBranches(){
	// copied from here: https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/HiggsTauTau/src/HTTSync.cc#L48-L201
	syncTree->Branch("run" ,			&run,"run/I" );//Run
	syncTree->Branch("lumi" ,			&lumi,"lumi/I" );//Lumi
	syncTree->Branch("evt" ,			&evt ,"evt/I" );//Evt

	    //Event Variables
	syncTree->Branch("npv" ,			&npv ,"npv/I" );//NPV
	syncTree->Branch("npu" ,			&npu ,"npu/I" );//NPU
	syncTree->Branch("rho" ,			&rho,"rho/F" );//Rho

	    //Event Weights
	syncTree->Branch("mcweight" ,		&mcweight ,"mcweight/F");//MC Weight (xs/nevents * additional wieght (ie pt weight for gghiggs))
	syncTree->Branch("puweight" ,		&puweight ,"puweight/F");//Pielup Weight

	syncTree->Branch("trigweight_1" ,	&trigweight_1,"trigweight_1/F");//Effieiency Scale factor (all components multiplied in)
	syncTree->Branch("trigweight_2" ,	&trigweight_2 ,"trigweight_2/F");//Effieiency Scale factor (all components multiplied in)
	syncTree->Branch("idweight_1" ,		&idweight_1 ,"idweight_1/F");//Effieiency Scale factor (all components multiplied in)
	syncTree->Branch("idweight_2" ,		&idweight_2 ,"idweight_2/F");//Effieiency Scale factor (all components multiplied in)
	syncTree->Branch("isoweight_1" ,	&isoweight_1,"isoweight_1/F");//Effieiency Scale factor (all components multiplied in)
	syncTree->Branch("isoweight_2" ,	&isoweight_2 ,"isoweight_2/F");//Effieiency Scale factor (all components multiplied in)
	syncTree->Branch("fakeweight" ,		&fakeweight ,"fakeweight/F");//Effieiency Scale factor (all components multiplied in)

	syncTree->Branch("effweight" ,		&effweight ,"effweight/F");//Effieiency Scale factor (all components multiplied in)
	syncTree->Branch("weight" ,			&weight,"weight/F" );//mcweight*puweight*effweight
	syncTree->Branch("embeddedWeight" ,	&embeddedWeight ,"embeddedWeight/F" );
	syncTree->Branch("signalWeight" ,	&signalWeight ,"signalWeight/F" );

	    //SV Fit variables
	syncTree->Branch("mvis" ,			&mvis ,"mvis/F" );//SV Fit using integration method
	syncTree->Branch("m_sv" ,			&m_sv,"m_sv/F" );//SV Fit using integration method
	syncTree->Branch("pt_sv" ,			&pt_sv,"pt_sv/F" );//SV Fit using integration method
	syncTree->Branch("eta_sv" ,			&eta_sv ,"eta_sv/F" );//SV Fit using integration method
	syncTree->Branch("phi_sv" ,			&phi_sv,"phi_sv/F" );//SV Fit using integration method
	syncTree->Branch("m_sv_Up" ,		&m_sv_Up ,"m_sv_Up/F" );//High Energy scale shape
	syncTree->Branch("m_sv_Down" ,		&m_sv_Down ,"m_sv_Down/F" );//Low Energy Scale Shape

	    ///First lepton : muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
	syncTree->Branch("pt_1" ,			&pt_1 ,"pt_1/F" ); //pT
	syncTree->Branch("phi_1" ,			&phi_1 ,"phi_1/F" ); //Phi
	syncTree->Branch("eta_1" ,			&eta_1 ,"eta_1/F" ); //Eta
	syncTree->Branch("m_1" ,			&m_1 ,"m_1/F" ); //Mass
	syncTree->Branch("q_1" ,			&q_1 ,"q_1/I" ); //Mass
	syncTree->Branch("iso_1" ,			&iso_1 ,"iso_1/F" ); //Delta Beta iso value
	syncTree->Branch("mva_1" ,			&mva_1 ,"mva_1/F" );//MVA id (when using electron) 0 otherwise
	syncTree->Branch("d0_1" ,			&d0_1,"d0_1/F" );//d0 with respect to primary vertex
	syncTree->Branch("dZ_1" ,			&dZ_1 ,"dZ_1/F" );//dZ with respect to primary vertex
	syncTree->Branch("passid_1" ,		&passid_1 ,"passid_1/B" );//Whether it passes id (not necessarily iso)
	syncTree->Branch("passiso_1" ,		&passiso_1,"passiso_1/B");//Whether it passes iso (not necessarily id)
	syncTree->Branch("mt_1" ,			&mt_1 ,"mt_1/F" );//mT of first lepton wrt to MVA met

	    ///Second lepton : hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT) Tau for Tau Tau
	syncTree->Branch("pt_2" ,&pt_2 ,"pt_2/F" );//pT
	syncTree->Branch("phi_2" ,&phi_2 ,"phi_2/F" );//Phi
	syncTree->Branch("eta_2" ,&eta_2 ,"eta_2/F" );//Eta
	syncTree->Branch("m_2" ,&m_2 ,"m_2/F" );//Mass (visible mass for hadronic Tau)
	syncTree->Branch("q_2" ,&q_2 ,"q_2/I" ); //Mass
	syncTree->Branch("iso_2" ,&iso_2,"iso_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	syncTree->Branch("d0_2" ,&d0_2 ,"d0_2/F" );//d0 with respect to primary vertex
	syncTree->Branch("dZ_2" ,&dZ_2 ,"dZ_2/F" );//dZ with respect to primary vertex

	syncTree->Branch("pt_tt" ,&pt_tt ,"pt_tt/F" );//pT

	syncTree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2" ,&byCombinedIsolationDeltaBetaCorrRaw3Hits_2 ,"byCombinedIsolationDeltaBetaCorrRaw3Hits_2/F" );
	syncTree->Branch("againstElectronMVA3raw_2" ,&againstElectronMVA3raw_2 ,"againstElectronMVA3raw_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	syncTree->Branch("byIsolationMVA2raw_2" ,&byIsolationMVA2raw_2 ,"byIsolationMVA2raw_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	syncTree->Branch("againstMuonLoose2_2" ,&againstMuonLoose2_2 ,"againstMuonLoose2_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	syncTree->Branch("againstMuonMedium2_2" ,&againstMuonMedium2_2 ,"againstMuonMedium2_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	syncTree->Branch("againstMuonTight2_2" ,&againstMuonTight2_2 ,"againstMuonTight2_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon


	syncTree->Branch("mva_2" ,&mva_2 ,"mva_2/F" );//MVA id (for anti electron id)
	syncTree->Branch("passid_2" ,&passid_2 ,"passid_2/B" );//Whether it passes id (not necessarily iso)
	syncTree->Branch("passiso_2" ,&passiso_2 ,"passiso_2/B");//Whether it passes iso (not necessarily id)
	syncTree->Branch("mt_2" ,&mt_2,"mt_2/F" );//mT of 2nd lepton wrt to MVA met

	    //Met related variables
	syncTree->Branch("met" ,&met,"met/F" ); //pfmet
	syncTree->Branch("metphi" ,&metphi ,"metphi/F" ); //pfmet Phi

	syncTree->Branch("l1met" ,&l1met ,"l1met/F" ); //l1met
	syncTree->Branch("l1metphi" ,&l1metphi ,"l1metphi/F" ); //pfmet Phi
	syncTree->Branch("l1metcorr" ,&l1metcorr ,"l1metcorr/F" ); //l1met

	syncTree->Branch("calomet" ,&calomet ,"calomet/F" );
	syncTree->Branch("calometphi" ,&calometphi ,"calometphi/F" );
	syncTree->Branch("calometcorr" ,&calometcorr ,"calometcorr/F" );
	syncTree->Branch("calometphicorr" ,&calometphicorr ,"calometphicorr/F" );

	syncTree->Branch("mvamet" ,&mvamet ,"mvamet/F" ); //mvamet
	syncTree->Branch("mvametphi" ,&mvametphi ,"mvametphi/F" ); //mvamet Phi
	syncTree->Branch("pzetavis" ,&pzetavis ,"pzetavis/F" ); //pZeta Visible
	syncTree->Branch("pzetamiss" ,&pzetamiss,"pzetamiss/F"); //pZeta Missing
	    //MET covariance matrices
	syncTree->Branch("metcov00" ,&metcov00 ,"metcov00/F"); //pf met covariance matrix 00
	syncTree->Branch("metcov01" ,&metcov01 ,"metcov01/F"); //pf met covariance matrix 01
	syncTree->Branch("metcov10" ,&metcov10 ,"metcov10/F"); //pf met covariance matrix 10
	syncTree->Branch("metcov11" ,&metcov11 ,"metcov11/F"); //pf met covariance matrix 11
	    //MVAMet covariance matrices
	syncTree->Branch("mvacov00" ,&mvacov00 ,"mvacov00/F"); //mva met covariance matrix 00
	syncTree->Branch("mvacov01" ,&mvacov01 ,"mvacov01/F"); //mva met covariance matrix 01
	syncTree->Branch("mvacov10" ,&mvacov10 ,"mvacov10/F"); //mva met covariance matrix 10
	syncTree->Branch("mvacov11" ,&mvacov11,"mvacov11/F"); //mva met covariance matrix 11

	    //First Jet : leading jet after applying Jet energy corrections (excluding hadronic Tau)
	syncTree->Branch("jpt_1" ,&jpt_1 ,"jpt_1/F" );//Jet Pt after corrections
	syncTree->Branch("jeta_1" ,&jeta_1 ,"jeta_1/F" );//Jet Eta
	syncTree->Branch("jphi_1" ,&jphi_1 ,"jphi_1/F" );//Jet Phi
	syncTree->Branch("jptraw_1" ,&jptraw_1 ,"jptraw_1/F" );//Jet Raw Pt (before corrections)
	syncTree->Branch("jptunc_1" ,&jptunc_1 ,"jptunc_1/F" );//Jet Unc (relative to Jet corrected pT)
	syncTree->Branch("jmva_1" ,&jmva_1 ,"jmva_1/F" );//Jet MVA id value
	syncTree->Branch("jlrm_1" ,&jlrm_1 ,"jlrm_1/F" );//Jet MVA id value
	syncTree->Branch("jctm_1" ,&jctm_1 ,"jctm_1/I" );//Jet MVA id value
	syncTree->Branch("jpass_1" ,&jpass_1 ,"jpass_1/B" );//Whether Jet pass PU Id Loose WP

	    //Second Jet : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
	syncTree->Branch("jpt_2" ,&jpt_2 ,"jpt_2/F" );//Jet Pt after corrections
	syncTree->Branch("jeta_2" ,&jeta_2 ,"jeta_2/F" );//Jet Eta
	syncTree->Branch("jphi_2" ,&jphi_2 ,"jphi_2/F" );//Jet Phi
	syncTree->Branch("jptraw_2" ,&jptraw_2 ,"jptraw_2/F" );//Jet Raw Pt (before corrections)
	syncTree->Branch("jptunc_2" ,&jptunc_2,"jptunc_2/F" );//Jet Unc (relative to Jet corrected pT)
	syncTree->Branch("jmva_2" ,&jmva_2 ,"jmva_2/F" );//Jet MVA id value
	syncTree->Branch("jlrm_2" ,&jlrm_2 ,"jlrm_2/F" );//Jet MVA id value
	syncTree->Branch("jctm_2" ,&jctm_2 ,"jctm_2/I" );//Jet MVA id value
	syncTree->Branch("jpass_2" ,&jpass_2 ,"jpass_2/B" );//Whether jet passes PU Id Loose WP

	    //B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
	syncTree->Branch("bpt" ,&bpt ,"bpt/F" );//Corrected BTag Pt
	syncTree->Branch("beta" ,&beta,"beta/F" );//Btag Eta
	syncTree->Branch("bphi" ,&bphi ,"bphi/F" );//Btag Phi

	    //Di Jet kinematic variables for VBF selection ==> Two leading pT Jets
	syncTree->Branch("mjj" ,&mjj ,"mjj/F" );//Mass Di Jet system
	syncTree->Branch("jdeta" ,&jdeta ,"jdeta/F" );//|jeta_1-jeta_2|
	syncTree->Branch("njetingap" ,&njetingap ,"njetingap/I");//# of Jets between two jets
	syncTree->Branch("mva" ,&mva ,"mva/F" );//VBF MVA value

	    //Variables that go into the VBF MVA
	syncTree->Branch("jdphi" ,&jdphi ,"jdphi/F" );//Delta Phi between two leading jets
	syncTree->Branch("dijetpt" ,&dijetpt ,"dijetpt/F" );//Pt of the di jet system
	syncTree->Branch("dijetphi" ,&dijetphi ,"dijetphi/F" );//Phi of the di jet system
	syncTree->Branch("hdijetphi" ,&hdijetphi ,"hdijetphi/F" );//Phi of the di jet system - Higgs system phi
	syncTree->Branch("visjeteta" ,&visjeteta ,"visjeteta/F");//TMath::Min(eta_vis - jeta,eta_vis,jeta2);
	syncTree->Branch("ptvis" ,&ptvis,"ptvis/F" );//Pt Vis

	    //number of btags passing btag id ( pt > 20 )
	syncTree->Branch("nbtag" ,&nbtag ,"nbtag/I");

	    //number of jets passing jet id ( pt > 30 )
	syncTree->Branch("njets" ,&njets ,"njets/I");
	syncTree->Branch("njetspt20" ,&njetspt20 ,"njetspt20/I");

	syncTree->Branch("mva_gf" ,&mva_gf ,"mva_gf/F");
	syncTree->Branch("mva_vbf" ,&mva_vbf ,"mva_vbf/F");
}
