#ifndef sample_h
#define sample_h

/////////////////////////////////////////////////
//
// sample: struct to hold information on each MC sample to be plotted
// to be used in plotting macro
//
/////////////////////////////////////////////////

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "Riostream.h"
#include <fstream>
#include <istream>

#include "configInfo.h"
#include "userConfig.h"

// container to hold all information for one sample
struct sample{
	int id;
	double xsec;
	int nEvt;
	int color;
	TString legName;
	std::vector<double> mcScale;
	double syst;
	std::vector<TString> identifier;

	//// constructors
	// default constructor
	sample(){};
	// create sample from raw numbers
	sample(TString LegName, double Xsec, int NEvt, double lumi, int Color, TString Identifier){
		legName = LegName;
		xsec = Xsec;
		nEvt = NEvt;
		color = Color;
		identifier.push_back(Identifier);
		mcScale.push_back(lumi * Xsec / NEvt);
		//todo: define syst
	}
	// create sample using dataMCType and SkimSummary info
	sample(int dataMcType, TString LegName, double Xsec, double lumi, int Color, TString Identifier, std::map<int,int> nEventsMap){
		id = dataMcType;
		legName = LegName;
		xsec = Xsec;
		nEvt = nEventsMap[id];
		color = Color;
		identifier.push_back(Identifier);
		mcScale.push_back(lumi * xsec / nEvt);
		//todo: define syst
	}
	// create slim sample without scaling information
	sample(TString LegName, int Color, TString Identifier){
		legName = LegName;
		color = Color;
		identifier.push_back(Identifier);
		mcScale.push_back(1);
		//todo: define syst
	}
	// copy constructor
	sample(const sample &original){
		id = original.id;
		xsec = original.xsec;
		nEvt = original.nEvt;
		color = original.color;
		legName = original.legName;
		mcScale = original.mcScale;
		syst = original.syst;
		identifier = original.identifier;
	}
	sample(const sample &original, TString LegName, int Color){
		id = original.id;
		xsec = original.xsec;
		nEvt = original.nEvt;
		color = Color;
		legName = LegName;
		mcScale = original.mcScale;
		syst = original.syst;
		identifier = original.identifier;
	}
	// add sample to existing sample
	sample operator+=(const sample& rhs){
		if(rhs.identifier.size() != 1){
			std::cout << "ERROR: sample to add must contain exactly one element, size = " << rhs.identifier.size() << std::endl;
			return sample();
		}
		identifier.push_back(rhs.identifier.at(0));
		if(rhs.mcScale.size() == 1) mcScale.push_back(rhs.mcScale.at(0));
		//todo: define syst
		return *this;
	}
};

// print sample, for debugging purposes
void printSample(sample& s, configInfo conf){
	printf("**** Printing sample %s ****\n", s.legName.Data());
	if(s.identifier.size() == 0){
		printf("  Empty sample. Why are you even looking at this?\n");
		return;
	}

	if( !conf.isFileLumiScaled && (s.identifier.size() != s.mcScale.size()) ){
			printf("  WARNING: sizes of vectors in sample differ! mcScale info is needed!\n");
			printf("  size(identifier) = %i, size(mcScale) = %i\n", s.identifier.size(), s.mcScale.size());
			return;
	}

	printf("  will be drawn with color = %i\n", s.color);
	printf("  sample consists of %i subsamples:\n", s.identifier.size());
	for(unsigned i = 0; i<s.identifier.size(); i++){
		printf("    subsample %i has identifier %s \n", i+1, s.identifier.at(i).Data());
		if(!conf.isFileLumiScaled) printf("      and mcScale = %.6f\n", s.mcScale.at(i));
	}
	//todo: add syst
	return;
}

// function to read number of events from SkimSummary
std::map<int, int> readSkimSummary(TString skimsummaryFile){
	if(verbose) std::cout << "--> readSkimSummary(...)" << std::endl;
	ifstream input_file;
	char *file_=(char*)skimsummaryFile.Data();
	input_file.open(file_, std::ios::in);
	if (!(input_file)){
		std::cout << "\nERROR: Opening SkimSummary file "<< skimsummaryFile <<" fhas failed.\n" << std::endl;
	}
	std::cout << "\nOpened SkimSummary SkimEff file: "<< skimsummaryFile <<".\n" << std::endl;
	std::string s;
	int a=0;
	std::map<int, int> nEvtMap;
	while(getline(input_file, s)){
		a++;
		if(a>1000) break;
		std::stringstream line(s);
		TString tmp;
		int id;
		int nevents;
		float neventserr;
		float nevents_sel;
		float neventserr_sel;
		float noweight;
		float noweight_sel;
		line >> tmp >> id >> tmp >> nevents >> tmp >> neventserr >> tmp >> nevents_sel >> tmp >> neventserr_sel >> tmp >> noweight >> tmp >> noweight_sel;
		// only consider main data-mc-types, thus removing all but 2 last digits
		nEvtMap[id%100] += nevents;
	}
	return nEvtMap;
}

// test validity of sample vector
bool testSamples(std::vector<sample> samples, configInfo conf){
	if(samples.size() < 1){
		std::cout << "ERROR: No sample specified." << std::endl;
		return false;
	}
	for(unsigned s = 0; s < samples.size(); s++){
		sample sam = samples.at(s);
		if(sam.legName == ""){
			std::cout << "ERROR: Sample has no LegName." << std::endl;
			return false;
		}
		if(sam.identifier.at(0) == ""){
			std::cout << "ERROR: Sample has no identifier." << std::endl;
			return false;
		}
		if(!conf.isFileLumiScaled){
			if(sam.mcScale.at(0) == 0){
				std::cout << "ERROR: Sample has no mcScale, but needs to be scaled." << std::endl;
				return false;
			}
		}else{
			for(unsigned ss = 0; ss < sam.mcScale.size(); ss++){
				if(sam.mcScale.at(ss) != 1){
					sam.mcScale.at(ss) = 1;
					std::cout << "WARNING: Sample " << sam.identifier.at(ss) << " is already scaled but additional scale was given. Scale set to 1." << std::endl;
				}
			}
		}

	}
	return true;
}

#endif

