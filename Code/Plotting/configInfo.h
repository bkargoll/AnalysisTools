#ifndef configInfo_h
#define configInfo_h

#include "userConfig.h"

/////////////////////////////////////////////////
//
// configInfo: struct to hold user specific information
// to be used in plotting macro
//
/////////////////////////////////////////////////

struct configInfo{
	TString fileName;
	TFile* file;
	bool isFileLumiScaled;
	TString prefix;
	double lumi;

	
	// constructors
	// default constructor
	configInfo() {};
	// main constructor
	configInfo(TString inputFile, bool IsFileLumiScaled, TString Prefix, double Lumi){
		fileName = inputFile;
		file = new TFile(fileName, "READ");
		isFileLumiScaled = IsFileLumiScaled;
		prefix = Prefix;
		lumi = Lumi;
	}
};

// check validity of configuration
bool testConfigInfo(configInfo config, bool verbose){
	bool isZombie = config.file->IsZombie();
	bool isOpen = config.file->IsOpen();

	if(verbose){
		std::cout << "### Info on used configInfo object:" << std::endl;
		std::cout << "    Input file name: " << config.fileName << std::endl;
		if(!isZombie && isOpen)	std::cout << "    Input file is open and valid." << std::endl;
		std::cout << "    Histograms will be accessed by using prefix " << config.prefix << std::endl;
		std::cout << "    Input will be scaled by lumi: " << !config.isFileLumiScaled << std::endl;
		if(!config.isFileLumiScaled) std::cout << "    For scaling will use lumi = " << config.lumi <<  std::endl;
	}
	if(!isOpen) {
		std::cout << "ERROR: Input file has not been opened!" << std::endl;
		return false;
	}
	else if(isZombie) {
		std::cout << "ERROR: Input file is in zombie state." << std::endl;
		return false;
	}
	if(!config.isFileLumiScaled && config.lumi < 0.1){
		std::cout << "ERROR: Input histos need to be scaled, but lumi is missing." << std::endl;
		return false;
	}

	return true;
}

#endif
