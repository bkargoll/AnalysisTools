#include "MuTauSync.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>
#include "SimpleFits/FitSoftware/interface/Logger.h"

MuTauSync::MuTauSync(TString Name_, TString id_):
  HToTaumuTauh(Name_,id_)
{
	Logger(Logger::Info) << "Setting up the class MuTauSync" << std::endl;
	// always run without category for sync exercise
	categoryFlag = "NoCategory";
}

MuTauSync::~MuTauSync(){
  for(unsigned int j=0; j<Npassed.size(); j++){
    Logger(Logger::Info) << "Selection Summary before: "
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  Logger(Logger::Info) << "Done." << std::endl;
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
	syncTree_TauPlusXA = new TTree("syncTree_TauPlusXA", "syncTree_TauPlusXA");
	syncTree_VBFHiggs = new TTree("syncTree_VBFHiggs", "syncTree_VBFHiggs");
	defineBranches(syncTree_TauPlusXA);
	defineBranches(syncTree_VBFHiggs);
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
		trigweight_1 = (!Ntp->isData()) ? RSF->HiggsTauTau_MuTau_Trigger_Mu_ScaleMCtoData(Ntp->Muon_p4(selMuon)) : -10;
		trigweight_2 = (!Ntp->isData()) ? RSF->HiggsTauTau_MuTau_Trigger_Tau_ScaleMCtoData(Ntp->PFTau_p4(selTau, "")) : -10;
		idweight_1 = (!Ntp->isData()) ? RSF->HiggsTauTau_MuTau_Id_Mu(Ntp->Muon_p4(selMuon)) : -10;
		idweight_2 = -10;
		isoweight_1 = (!Ntp->isData()) ? RSF->HiggsTauTau_MuTau_Iso_Mu(Ntp->Muon_p4(selMuon)) : -10;
		isoweight_2 = -10;
		fakeweight = -10;
		effweight = trigweight_1 * trigweight_2 * idweight_1 * isoweight_1;
		weight = w;
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
		d0_1 = Ntp->dxySigned(Ntp->Muon_p4(selMuon),Ntp->Muon_Poca(selMuon),Ntp->Vtx(selVertex));
		dZ_1 = Ntp->dzSigned(Ntp->Muon_p4(selMuon),Ntp->Muon_Poca(selMuon),Ntp->Vtx(selVertex));
		passid_1 = true; // passes obviously, as it was selected
		passiso_1 = true;
		mt_1 = Ntp->transverseMass(pt_1,phi_1,Ntp->MET_CorrMVAMuTau_et(),Ntp->MET_CorrMVAMuTau_phi());
		 // Second lepton : hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT) Tau for Tau Tau
		pt_2 = Ntp->PFTau_p4(selTau, "").Pt();
		phi_2 = Ntp->PFTau_p4(selTau, "").Phi();
		eta_2 = Ntp->PFTau_p4(selTau, "").Eta();
		m_2 = Ntp->PFTau_p4(selTau, "").M();
		q_2 = Ntp->PFTau_Charge(selTau);
		iso_2 = -10; //MVA iso for hadronic Tau
		d0_2 = Ntp->dxy(Ntp->PFTau_p4(selTau, ""),Ntp->PFTau_Poca(selTau),Ntp->Vtx(selVertex));
		dZ_2 = Ntp->dz(Ntp->PFTau_p4(selTau, ""),Ntp->PFTau_Poca(selTau),Ntp->Vtx(selVertex));
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
		mt_2 = Ntp->transverseMass(pt_2,phi_2,Ntp->MET_CorrMVAMuTau_et(),Ntp->MET_CorrMVAMuTau_phi());
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
		if (selectedJetsClean.size() > 0){
			jpt_1 = Ntp->PFJet_p4(selectedJetsClean.at(0)).Pt();
			jeta_1 = Ntp->PFJet_p4(selectedJetsClean.at(0)).Eta();
			jphi_1 = Ntp->PFJet_p4(selectedJetsClean.at(0)).Phi();
			jptraw_1 = -10;
			jptunc_1 = -10;
			jmva_1 = Ntp->PFJet_PUJetID_discr(selectedJetsClean.at(0));
			jlrm_1 = -10;
			jctm_1 = -10;
			jpass_1 = int(Ntp->PFJet_PUJetID_looseWP(selectedJetsClean.at(0)));
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
		if (selectedJetsClean.size() > 1){
			jpt_2 = Ntp->PFJet_p4(selectedJetsClean.at(1)).Pt();
			jeta_2 = Ntp->PFJet_p4(selectedJetsClean.at(1)).Eta();
			jphi_2 = Ntp->PFJet_p4(selectedJetsClean.at(1)).Phi();
			jptraw_2 = -10;
			jptunc_2 = -10;
			jmva_2 = Ntp->PFJet_PUJetID_discr(selectedJetsClean.at(1));
			jlrm_2 = -10;
			jctm_2 = -10;
			jpass_2 = int(Ntp->PFJet_PUJetID_looseWP(selectedJetsClean.at(1)));
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
		if (selectedBJets.size() > 0){
			bpt = Ntp->PFJet_p4(selectedBJets.at(0)).Pt();
			beta = Ntp->PFJet_p4(selectedBJets.at(0)).Eta();
			bphi = Ntp->PFJet_p4(selectedBJets.at(0)).Phi();
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
		nbtag = selectedBJets.size();
		 // number of jets passing jet id ( pt > 30 )
		njets = selectedJetsClean.size();
		njetspt20 = -10;
		 // mva output for e+mu channel
		mva_gf = -10;
		mva_vbf = 10;

		if(Ntp->isData()) 	syncTree_TauPlusXA->Fill();
		else				syncTree_VBFHiggs->Fill();
	}
}

void MuTauSync::Finish(){
	HToTaumuTauh::Finish();

	syncTree_TauPlusXA->SetDirectory(syncFile);
	syncTree_TauPlusXA->GetCurrentFile()->Write();
	syncTree_VBFHiggs->SetDirectory(syncFile);
	syncTree_VBFHiggs->GetCurrentFile()->Write();
	syncFile->Close();
}

void MuTauSync::defineBranches(TTree* tree){
	// copied from here: https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/HiggsTauTau/src/HTTSync.cc#L48-L201
	tree->Branch("run" ,			&run,"run/I" );//Run
	tree->Branch("lumi" ,			&lumi,"lumi/I" );//Lumi
	tree->Branch("evt" ,			&evt ,"evt/I" );//Evt

	    //Event Variables
	tree->Branch("npv" ,			&npv ,"npv/I" );//NPV
	tree->Branch("npu" ,			&npu ,"npu/I" );//NPU
	tree->Branch("rho" ,			&rho,"rho/F" );//Rho

	    //Event Weights
	tree->Branch("mcweight" ,		&mcweight ,"mcweight/F");//MC Weight (xs/nevents * additional wieght (ie pt weight for gghiggs))
	tree->Branch("puweight" ,		&puweight ,"puweight/F");//Pielup Weight

	tree->Branch("trigweight_1" ,	&trigweight_1,"trigweight_1/F");//Effieiency Scale factor (all components multiplied in)
	tree->Branch("trigweight_2" ,	&trigweight_2 ,"trigweight_2/F");//Effieiency Scale factor (all components multiplied in)
	tree->Branch("idweight_1" ,		&idweight_1 ,"idweight_1/F");//Effieiency Scale factor (all components multiplied in)
	tree->Branch("idweight_2" ,		&idweight_2 ,"idweight_2/F");//Effieiency Scale factor (all components multiplied in)
	tree->Branch("isoweight_1" ,	&isoweight_1,"isoweight_1/F");//Effieiency Scale factor (all components multiplied in)
	tree->Branch("isoweight_2" ,	&isoweight_2 ,"isoweight_2/F");//Effieiency Scale factor (all components multiplied in)
	tree->Branch("fakeweight" ,		&fakeweight ,"fakeweight/F");//Effieiency Scale factor (all components multiplied in)

	tree->Branch("effweight" ,		&effweight ,"effweight/F");//Effieiency Scale factor (all components multiplied in)
	tree->Branch("weight" ,			&weight,"weight/F" );//mcweight*puweight*effweight
	tree->Branch("embeddedWeight" ,	&embeddedWeight ,"embeddedWeight/F" );
	tree->Branch("signalWeight" ,	&signalWeight ,"signalWeight/F" );

	    //SV Fit variables
	tree->Branch("mvis" ,			&mvis ,"mvis/F" );//SV Fit using integration method
	tree->Branch("m_sv" ,			&m_sv,"m_sv/F" );//SV Fit using integration method
	tree->Branch("pt_sv" ,			&pt_sv,"pt_sv/F" );//SV Fit using integration method
	tree->Branch("eta_sv" ,			&eta_sv ,"eta_sv/F" );//SV Fit using integration method
	tree->Branch("phi_sv" ,			&phi_sv,"phi_sv/F" );//SV Fit using integration method
	tree->Branch("m_sv_Up" ,		&m_sv_Up ,"m_sv_Up/F" );//High Energy scale shape
	tree->Branch("m_sv_Down" ,		&m_sv_Down ,"m_sv_Down/F" );//Low Energy Scale Shape

	    ///First lepton : muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
	tree->Branch("pt_1" ,			&pt_1 ,"pt_1/F" ); //pT
	tree->Branch("phi_1" ,			&phi_1 ,"phi_1/F" ); //Phi
	tree->Branch("eta_1" ,			&eta_1 ,"eta_1/F" ); //Eta
	tree->Branch("m_1" ,			&m_1 ,"m_1/F" ); //Mass
	tree->Branch("q_1" ,			&q_1 ,"q_1/I" ); //Mass
	tree->Branch("iso_1" ,			&iso_1 ,"iso_1/F" ); //Delta Beta iso value
	tree->Branch("mva_1" ,			&mva_1 ,"mva_1/F" );//MVA id (when using electron) 0 otherwise
	tree->Branch("d0_1" ,			&d0_1,"d0_1/F" );//d0 with respect to primary vertex
	tree->Branch("dZ_1" ,			&dZ_1 ,"dZ_1/F" );//dZ with respect to primary vertex
	tree->Branch("passid_1" ,		&passid_1 ,"passid_1/B" );//Whether it passes id (not necessarily iso)
	tree->Branch("passiso_1" ,		&passiso_1,"passiso_1/B");//Whether it passes iso (not necessarily id)
	tree->Branch("mt_1" ,			&mt_1 ,"mt_1/F" );//mT of first lepton wrt to MVA met

	    ///Second lepton : hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT) Tau for Tau Tau
	tree->Branch("pt_2" ,&pt_2 ,"pt_2/F" );//pT
	tree->Branch("phi_2" ,&phi_2 ,"phi_2/F" );//Phi
	tree->Branch("eta_2" ,&eta_2 ,"eta_2/F" );//Eta
	tree->Branch("m_2" ,&m_2 ,"m_2/F" );//Mass (visible mass for hadronic Tau)
	tree->Branch("q_2" ,&q_2 ,"q_2/I" ); //Mass
	tree->Branch("iso_2" ,&iso_2,"iso_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	tree->Branch("d0_2" ,&d0_2 ,"d0_2/F" );//d0 with respect to primary vertex
	tree->Branch("dZ_2" ,&dZ_2 ,"dZ_2/F" );//dZ with respect to primary vertex

	tree->Branch("pt_tt" ,&pt_tt ,"pt_tt/F" );//pT

	tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2" ,&byCombinedIsolationDeltaBetaCorrRaw3Hits_2 ,"byCombinedIsolationDeltaBetaCorrRaw3Hits_2/F" );
	tree->Branch("againstElectronMVA3raw_2" ,&againstElectronMVA3raw_2 ,"againstElectronMVA3raw_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	tree->Branch("byIsolationMVA2raw_2" ,&byIsolationMVA2raw_2 ,"byIsolationMVA2raw_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	tree->Branch("againstMuonLoose2_2" ,&againstMuonLoose2_2 ,"againstMuonLoose2_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	tree->Branch("againstMuonMedium2_2" ,&againstMuonMedium2_2 ,"againstMuonMedium2_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon
	tree->Branch("againstMuonTight2_2" ,&againstMuonTight2_2 ,"againstMuonTight2_2/F" );//MVA iso for hadronic Tau, Delta Beta for muon


	tree->Branch("mva_2" ,&mva_2 ,"mva_2/F" );//MVA id (for anti electron id)
	tree->Branch("passid_2" ,&passid_2 ,"passid_2/B" );//Whether it passes id (not necessarily iso)
	tree->Branch("passiso_2" ,&passiso_2 ,"passiso_2/B");//Whether it passes iso (not necessarily id)
	tree->Branch("mt_2" ,&mt_2,"mt_2/F" );//mT of 2nd lepton wrt to MVA met

	    //Met related variables
	tree->Branch("met" ,&met,"met/F" ); //pfmet
	tree->Branch("metphi" ,&metphi ,"metphi/F" ); //pfmet Phi

	tree->Branch("l1met" ,&l1met ,"l1met/F" ); //l1met
	tree->Branch("l1metphi" ,&l1metphi ,"l1metphi/F" ); //pfmet Phi
	tree->Branch("l1metcorr" ,&l1metcorr ,"l1metcorr/F" ); //l1met

	tree->Branch("calomet" ,&calomet ,"calomet/F" );
	tree->Branch("calometphi" ,&calometphi ,"calometphi/F" );
	tree->Branch("calometcorr" ,&calometcorr ,"calometcorr/F" );
	tree->Branch("calometphicorr" ,&calometphicorr ,"calometphicorr/F" );

	tree->Branch("mvamet" ,&mvamet ,"mvamet/F" ); //mvamet
	tree->Branch("mvametphi" ,&mvametphi ,"mvametphi/F" ); //mvamet Phi
	tree->Branch("pzetavis" ,&pzetavis ,"pzetavis/F" ); //pZeta Visible
	tree->Branch("pzetamiss" ,&pzetamiss,"pzetamiss/F"); //pZeta Missing
	    //MET covariance matrices
	tree->Branch("metcov00" ,&metcov00 ,"metcov00/F"); //pf met covariance matrix 00
	tree->Branch("metcov01" ,&metcov01 ,"metcov01/F"); //pf met covariance matrix 01
	tree->Branch("metcov10" ,&metcov10 ,"metcov10/F"); //pf met covariance matrix 10
	tree->Branch("metcov11" ,&metcov11 ,"metcov11/F"); //pf met covariance matrix 11
	    //MVAMet covariance matrices
	tree->Branch("mvacov00" ,&mvacov00 ,"mvacov00/F"); //mva met covariance matrix 00
	tree->Branch("mvacov01" ,&mvacov01 ,"mvacov01/F"); //mva met covariance matrix 01
	tree->Branch("mvacov10" ,&mvacov10 ,"mvacov10/F"); //mva met covariance matrix 10
	tree->Branch("mvacov11" ,&mvacov11,"mvacov11/F"); //mva met covariance matrix 11

	    //First Jet : leading jet after applying Jet energy corrections (excluding hadronic Tau)
	tree->Branch("jpt_1" ,&jpt_1 ,"jpt_1/F" );//Jet Pt after corrections
	tree->Branch("jeta_1" ,&jeta_1 ,"jeta_1/F" );//Jet Eta
	tree->Branch("jphi_1" ,&jphi_1 ,"jphi_1/F" );//Jet Phi
	tree->Branch("jptraw_1" ,&jptraw_1 ,"jptraw_1/F" );//Jet Raw Pt (before corrections)
	tree->Branch("jptunc_1" ,&jptunc_1 ,"jptunc_1/F" );//Jet Unc (relative to Jet corrected pT)
	tree->Branch("jmva_1" ,&jmva_1 ,"jmva_1/F" );//Jet MVA id value
	tree->Branch("jlrm_1" ,&jlrm_1 ,"jlrm_1/F" );//Jet MVA id value
	tree->Branch("jctm_1" ,&jctm_1 ,"jctm_1/I" );//Jet MVA id value
	tree->Branch("jpass_1" ,&jpass_1 ,"jpass_1/B" );//Whether Jet pass PU Id Loose WP

	    //Second Jet : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
	tree->Branch("jpt_2" ,&jpt_2 ,"jpt_2/F" );//Jet Pt after corrections
	tree->Branch("jeta_2" ,&jeta_2 ,"jeta_2/F" );//Jet Eta
	tree->Branch("jphi_2" ,&jphi_2 ,"jphi_2/F" );//Jet Phi
	tree->Branch("jptraw_2" ,&jptraw_2 ,"jptraw_2/F" );//Jet Raw Pt (before corrections)
	tree->Branch("jptunc_2" ,&jptunc_2,"jptunc_2/F" );//Jet Unc (relative to Jet corrected pT)
	tree->Branch("jmva_2" ,&jmva_2 ,"jmva_2/F" );//Jet MVA id value
	tree->Branch("jlrm_2" ,&jlrm_2 ,"jlrm_2/F" );//Jet MVA id value
	tree->Branch("jctm_2" ,&jctm_2 ,"jctm_2/I" );//Jet MVA id value
	tree->Branch("jpass_2" ,&jpass_2 ,"jpass_2/B" );//Whether jet passes PU Id Loose WP

	    //B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
	tree->Branch("bpt" ,&bpt ,"bpt/F" );//Corrected BTag Pt
	tree->Branch("beta" ,&beta,"beta/F" );//Btag Eta
	tree->Branch("bphi" ,&bphi ,"bphi/F" );//Btag Phi

	    //Di Jet kinematic variables for VBF selection ==> Two leading pT Jets
	tree->Branch("mjj" ,&mjj ,"mjj/F" );//Mass Di Jet system
	tree->Branch("jdeta" ,&jdeta ,"jdeta/F" );//|jeta_1-jeta_2|
	tree->Branch("njetingap" ,&njetingap ,"njetingap/I");//# of Jets between two jets
	tree->Branch("mva" ,&mva ,"mva/F" );//VBF MVA value

	    //Variables that go into the VBF MVA
	tree->Branch("jdphi" ,&jdphi ,"jdphi/F" );//Delta Phi between two leading jets
	tree->Branch("dijetpt" ,&dijetpt ,"dijetpt/F" );//Pt of the di jet system
	tree->Branch("dijetphi" ,&dijetphi ,"dijetphi/F" );//Phi of the di jet system
	tree->Branch("hdijetphi" ,&hdijetphi ,"hdijetphi/F" );//Phi of the di jet system - Higgs system phi
	tree->Branch("visjeteta" ,&visjeteta ,"visjeteta/F");//TMath::Min(eta_vis - jeta,eta_vis,jeta2);
	tree->Branch("ptvis" ,&ptvis,"ptvis/F" );//Pt Vis

	    //number of btags passing btag id ( pt > 20 )
	tree->Branch("nbtag" ,&nbtag ,"nbtag/I");

	    //number of jets passing jet id ( pt > 30 )
	tree->Branch("njets" ,&njets ,"njets/I");
	tree->Branch("njetspt20" ,&njetspt20 ,"njetspt20/I");

	tree->Branch("mva_gf" ,&mva_gf ,"mva_gf/F");
	tree->Branch("mva_vbf" ,&mva_vbf ,"mva_vbf/F");
}
