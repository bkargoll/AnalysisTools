// stolen from crovelli (https://github.com/crovelli/Utili/blob/master/macro/duplicatesRemoval.C)
// and then adapted

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include <TEntryList.h>
#include <iostream>
#include <set>

using namespace std;

void removeDuplicatesFromTTree() {

	// minimal set of variables to be read from the input tree
	UInt_t RunNumber;
	UInt_t LumiNumber;
	UInt_t EventNumber;
	TBranch *b_RunNumber = 0;
	TBranch *b_LumiNumber = 0;
	TBranch *b_EventNumber = 0;

	// opening the input tree root file
	printf("open file...\n");
	TFile *inFile = TFile::Open("SVFitCache_unclean.root");
	if (inFile == 0) {
		printf("Error: cannot open the input root file");
		return;
	}

	// preparing the output file
	TString outputFileName = "outputFile.root";
	TFile fileOut(outputFileName.Data(), "RECREATE");
	fileOut.cd();

	std::vector<TString> treeNames;
	treeNames.push_back("Tree_TauPlusX");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-110_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-115_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-125_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-130_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-135_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-140_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-145_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-150_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-155_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-160_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-110_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-115_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-125_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-130_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-135_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-140_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-145_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-150_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-155_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-160_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-110_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-115_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-120_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-125_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-130_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-135_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-140_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-145_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-150_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-155_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_WH_ZH_TTH_HToTauTau_M-160_8TeV-pythia6-tauola");
	treeNames.push_back("Tree_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball");
	treeNames.push_back("Tree_WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola");
	treeNames.push_back("Tree_WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola");
	treeNames.push_back("Tree_WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola");
	treeNames.push_back("Tree_ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola");
	treeNames.push_back("Tree_ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola");
	treeNames.push_back("Tree_ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola");
	treeNames.push_back("Tree_TT_CT10_TuneZ2star_8TeV-powheg-tauola");
	treeNames.push_back("Tree_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola");
	treeNames.push_back("Tree_Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola");
	treeNames.push_back("Tree_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball");
	treeNames.push_back("Tree_GluGluToHToTauTau_M-120_8TeV-powheg-pythia6");
	treeNames.push_back("Tree_VBF_HToTauTau_M-120_8TeV-powheg-pythia6");

	int nTrees = treeNames.size();
	printf("loop over %i trees...\n", nTrees);
	for (int itree = 0; itree < nTrees; itree++) {
		// counting variables
		unsigned n_events(0), n_stored(0), n_removed(0);

		printf("load tree...\n");
		TTree *fChain = (TTree *) inFile->Get(treeNames.at(itree));
		if (fChain == 0) {
			printf("%s could not be loaded. Skip.\n", treeNames.at(itree).Data());
			continue;
		}
		n_events = fChain->GetEntries();
		//printf("tree has %i entries\n", n_events);

		// reading my input tree
		//printf("set branches...\n");
		fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
		fChain->SetBranchAddress("LumiNumber", &LumiNumber, &b_LumiNumber);
		fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);

		// preparing the clean tree after duplicated events removal
		//printf("create output tree...\n");
		TTree *tCleaned = fChain->CloneTree(0);
		tCleaned->SetName(treeNames.at(itree));

		// needed to remove duplicates
		//printf("create set and list...\n");
		set<Long64_t> eventIds;
		Long64_t myProduct;
		TEntryList *tlist = new TEntryList(fChain);
		tlist->SetName("listCleaned");

		// checking the original tree based on run and event numbers
		Long64_t nb;
		//printf("start loop...\n");
		for (Long64_t jentry = 0; jentry < n_events; jentry++) {
			nb = fChain->GetEntry(jentry);

			//if (jentry % 5000 == 0)
			//	cout << jentry << endl;

			Long64_t major = (Long64_t) RunNumber << 13 + LumiNumber;
			Long64_t minor = EventNumber;
			myProduct = (major << 31) + EventNumber;

			// if we have not seen this event yet, add it to the set and to the entry list
			if (eventIds.count(myProduct) == 0) {
				eventIds.insert(myProduct);
				tlist->Enter(jentry, fChain);
				tCleaned->Fill();
				n_stored++;
			} else {
				//cout << "this event is duplicated: " << myProduct << " " << RunNumber << " " << LumiNumber<< " " << EventNumber << endl;
				n_removed++;
			}
		}

		// saving the trees
		tCleaned->Write();
		//tlist->Write();
		printf("%s: %i - %i = %i \n", treeNames.at(itree).Data(), n_events, n_removed, n_stored);
	}

	fileOut.Close();
	printf("Cleaned trees have been saved to %s.\nDone.\n", outputFileName.Data());

	// delete [] index;
}
