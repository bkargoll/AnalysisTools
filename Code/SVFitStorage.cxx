#include "SVFitStorage.h"
#include "Parameters.h"
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TTreeIndex.h"
#include "TFile.h"
#include "SVFitObject.h"
#include <iostream>

SVFitStorage::SVFitStorage():
	treeName_("invalid"),
	isConfigured_(false),
	intreeLoaded_(false){

	TString thelib= getenv ("DATAFORMATS_LIB");
	gSystem->Load(thelib.Data());

	inputFileName = "SVFitInput_temp_";

	// create output tree
	// make sure to set name and title to proper values in Configure(...)
	outtree_= new TTree("temp", "temp");
}

SVFitStorage::~SVFitStorage(){
	if (isConfigured_)
		SaveTree();

	delete outtree_;
	// files opened by TChain are closed by ~TChain automatically
}

void SVFitStorage::Configure(TString datasetName){
	if (isConfigured_){
		std::cout << "WARNING: SVFitStorage object has been configured before. Make sure to configure only once. Abort config..." << std::endl;
		return;
	}

	treeName_ = "Tree_" +  datasetName;

	// setup output tree
	outtree_->SetName(treeName_);
	outtree_->SetTitle(treeName_);
	outtree_->Branch("RunNumber", &RunNumber_);
	outtree_->Branch("LumiNumber", &LumiNumber_);
	outtree_->Branch("EventNumber", &EventNumber_);
	outtree_->Branch("svfit", &svfit_);

	// setup input tree
	LoadTree();

	isConfigured_ = true;
}

void SVFitStorage::LoadTree(){
	TString key = "InputAuxiliaryFile:";
	int nfiles = GetFile(key);
	if (nfiles == 0) {
		std::cout << "Key not found: " << key << std::endl;
		std::cout << "WARNING: No input SVFit storage specified. Will calculate all SVFit values while running." << std::endl;
	} else {
		TDirectory *gdirectory_save = gDirectory;
		TChain *chain = new TChain(treeName_);
		for (int i = 0; i < nfiles; i++) {
			TString name = inputFileName;
			name += i;
			name += ".root";
			chain->Add(name);
		}
		intree_ = (TTree*) chain;

		//Set branches
		intree_->SetBranchAddress("RunNumber", &RunNumber_, &b_RunNumber_);
		intree_->SetBranchAddress("LumiNumber", &LumiNumber_, &b_LumiNumber_);
		intree_->SetBranchAddress("EventNumber", &EventNumber_, &b_EventNumber_);
		intree_->SetBranchAddress("svfit", &svfit_, &b_svfit_);

		//Build index tree
		if (RunNumber_ > 262143 || LumiNumber_ > 4095) {
			std::cout << "ERROR: RunNumber must be smaller than 262144 and LumiNumber smaller than 4096" << std::endl;
			std::cout << "       Trying to access run " << RunNumber_ << " and lumi " << LumiNumber_ << std::endl;
			std::cout << "       This is going to crash..." << std::endl;
			return;
		}
		intree_->BuildIndex("(RunNumber<<13) + LumiNumber", "EventNumber"); // works only for run<262144 and lumi<4096 (true in 2012)
		index_ = (TTreeIndex*) intree_->GetTreeIndex();
		gDirectory = gdirectory_save;
		gDirectory->cd();

		intreeLoaded_ = true;
	}
}


void SVFitStorage::SaveTree(){
	if ( !isConfigured_ ){
		std::cout << "ERROR: SVFitStorage must be configured before SaveTree can be called." << std::endl;
		return;
	}
	//Load File name
	Parameters Par; // assumes configured in Analysis.cxx
	TString outputFileDCache;
	Par.GetString("OutputAuxiliaryFile:", outputFileDCache);
	TString outputFileLocal = "MySVFIT.root";
	//Save output
	TDirectory *gdirectory_save = gDirectory;
	TFile *outfile_ = TFile::Open(outputFileLocal, "RECREATE");
	if (!outfile_) {
		std::cout << "ERROR:  " << outputFileDCache << " not saved" << std::endl;
		return;
	}
	outfile_->cd();
	outtree_->Write(treeName_);
	outfile_->Close();
	gDirectory = gdirectory_save;
	gDirectory->cd();
	std::cout << "SVFit_Tree saved to " << outputFileLocal << std::endl;
	//Store file on the grid
	StoreFile(outputFileLocal, outputFileDCache);
	std::cout << outputFileLocal << " saved to the grid " << outputFileDCache << std::endl;
}

void SVFitStorage::SaveEvent(Int_t RunNumber, Int_t LumiNumber, Int_t EventNumber,SVFitObject svfit){
	if (!isConfigured_) {
		std::cout << "ERROR: SVFitStorage must be configured before SaveTree can be called." << std::endl;
		return;
	}

	//Fill event
	RunNumber_ = RunNumber;
	LumiNumber_ = LumiNumber;
	EventNumber_ = EventNumber;
	svfit_ = svfit;
	outtree_->Fill();
}

const SVFitObject& SVFitStorage::GetEvent(Int_t RunNumber, Int_t LumiNumber, Int_t EventNumber){
	if (!isConfigured_) {
		std::cout << "ERROR: SVFitStorage must be configured before GetEvent can be called." << std::endl;
		svfit_ = SVFitObject(); // invalid object
		return svfit_;
	}
	if (!intreeLoaded_) {
		svfit_ = SVFitObject(); // invalid object
		return svfit_;
	}

	//Get tree entry using index and then get svfit
	std::cout << (RunNumber << 13) + LumiNumber << "  " << EventNumber << " gives " << ((RunNumber << 13) + LumiNumber)<<31 + EventNumber << std::endl;
	Long64_t local = index_->GetEntryNumberWithIndex((RunNumber << 13) + LumiNumber, EventNumber);
	// make sure that entry exists, otherwise return invalid result
	if (local >= 0) {
		b_svfit_->GetEntry(local);
	} else {
		svfit_ = SVFitObject(); // invalid object
	}
	return svfit_;
}


