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
#include "TFile.h"
#include "SVFitObject.h"
#include <iostream>
#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "TChainIndex.h"
#include "TTreeIndex.h"

SVFitStorage::SVFitStorage():
	treeName_("invalid"),
	isConfigured_(false),
	intreeLoaded_(false){

	TString thelib= getenv ("DATAFORMATS_LIB");
	gSystem->Load(thelib.Data());

	inputFileName = "SVFitInput_temp_";

	// create ipnut and output tree
	// make sure to set name and title to proper values in Configure(...)
	outtree_= new TTree("temp", "temp");
	intree_ = new TChain("temp", "temp");

	// allocate memory for objects to be stored in tree
	svfit_ = new SVFitObject();
}

SVFitStorage::~SVFitStorage(){
	if (isConfigured_)
		SaveTree();

	delete outtree_;
	delete intree_;

	delete svfit_;
}

void SVFitStorage::Configure(TString datasetName){
	if (isConfigured_){
		Logger(Logger::Warning) << "SVFitStorage object has been configured before. Make sure to configure only once. Abort config..." << std::endl;
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
		Logger(Logger::Warning) << "Key not found: " << key <<
				"\n\tNo input SVFit storage specified. Will calculate all SVFit values while running." << std::endl;
	} else {
		intree_->SetName(treeName_);
		intree_->SetTitle(treeName_);
		TDirectory *gdirectory_save = gDirectory;
		for (int i = 0; i < nfiles; i++) {
			TString name = inputFileName;
			name += i;
			name += ".root";
			intree_->Add(name);
		}
		//intree_->GetEntries();
		//intree_ = (TTree*) chain;
		if (intree_->LoadTree(0) < 0)
			Logger(Logger::Error) << "Input TChain was not loaded correctly." << std::endl;

		//Set branches
		intree_->SetBranchAddress("RunNumber", &RunNumber_, &b_RunNumber_);
		intree_->SetBranchAddress("LumiNumber", &LumiNumber_, &b_LumiNumber_);
		intree_->SetBranchAddress("EventNumber", &EventNumber_, &b_EventNumber_);
		intree_->SetBranchAddress("svfit", &svfit_, &b_svfit_);

		Logger(Logger::Debug) << "Building the index ..." << std::endl;
		int a = intree_->BuildIndex("(RunNumber<<13) + LumiNumber", "EventNumber"); // works only for run<262144 and lumi<4096 (true in 2012)
		Logger(Logger::Debug) << "BuildIndex done with return value " << a << std::endl;
		index_ = (TChainIndex*) intree_->GetTreeIndex();
		gDirectory = gdirectory_save;
		gDirectory->cd();

		intreeLoaded_ = true;
	}
}


void SVFitStorage::SaveTree(){
	if ( !isConfigured_ ){
		Logger(Logger::Error) << "SVFitStorage must be configured before SaveTree can be called." << std::endl;
		return;
	}
	if ( outtree_->GetEntries() < 1){
		Logger(Logger::Info) << "Output tree contains no events, and thus no output file is created." << std::endl;
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
		Logger(Logger::Error) << outputFileDCache << " not saved" << std::endl;
		return;
	}
	outfile_->cd();
	outtree_->Write(treeName_);
	outfile_->Close();
	gDirectory = gdirectory_save;
	gDirectory->cd();
	Logger(Logger::Info) << "SVFit_Tree saved to " << outputFileLocal << std::endl;
	//Store file on the grid
	StoreFile(outputFileLocal, outputFileDCache);
	Logger(Logger::Info) << outputFileLocal << " saved to the grid " << outputFileDCache << std::endl;
}

void SVFitStorage::SaveEvent(Int_t RunNumber, Int_t LumiNumber, Int_t EventNumber, SVFitObject* svfit){
	if (!isConfigured_) {
		Logger(Logger::Error) << "SVFitStorage must be configured before SaveTree can be called." << std::endl;
		return;
	}

	//Fill event
	RunNumber_ = RunNumber;
	LumiNumber_ = LumiNumber;
	EventNumber_ = EventNumber;
	svfit_ = svfit;
	outtree_->Fill();
}

SVFitObject* SVFitStorage::GetEvent(UInt_t RunNumber, UInt_t LumiNumber, UInt_t EventNumber){
	if (!isConfigured_) {
		Logger(Logger::Error) << "SVFitStorage must be configured before GetEvent can be called." << std::endl;
		*svfit_ = SVFitObject(); // invalid object
		return svfit_;
	}
	if (!intreeLoaded_) {
		Logger(Logger::Debug) << "SVFitStorage not configured, thus GetEvent does not work." << std::endl;
		*svfit_ = SVFitObject(); // invalid object
		return svfit_;
	}
	if (RunNumber > 262143 || LumiNumber > 4095) {
		Logger(Logger::Error) << "RunNumber must be smaller than 262144 and LumiNumber smaller than 4096" <<
				"\n\tTrying to access run " << RunNumber << " and lumi " << LumiNumber <<
				"\n\tExpect undefined behavior!" << std::endl;
		*svfit_ = SVFitObject(); // invalid object
		return svfit_;
	}

	//Get tree entry using index and then get svfit
	Logger(Logger::Debug) << "Try to access run " << RunNumber << ", lumi " << LumiNumber << ", Event " << EventNumber <<
			"\n\ti.e. major = " << ((RunNumber << 13) + LumiNumber) << ", minor = " << EventNumber << std::endl;
	Long64_t local = index_->GetEntryNumberWithIndex((RunNumber << 13) + LumiNumber, EventNumber);

	// make sure that entry exists, otherwise return invalid result
	if (local < 0) {
		Logger(Logger::Debug) << "Index " << local << " not available in input tree." << std::endl;
		*svfit_ = SVFitObject(); // invalid object
		return svfit_;
	}
	else {
		Logger(Logger::Debug) << "Accessing index " << local << std::endl;
		intree_->GetEntry(local);
	}

	// check if correct event was loaded
	if (RunNumber_ != RunNumber) {
		Logger(Logger::Error) << "Event " << EventNumber << " was not loaded correctly."<< std::endl;
		*svfit_ = SVFitObject(); // rather return invalid object than wrong one
	}
	return svfit_;
}


