#include "Tvariable_MuMu.h"

Tvariable_MuMu::Tvariable_MuMu(TString Name_, TString id_):
Tvariable_Base(Name_,id_)
{
	std::cout << "Setting up the class Tvariable_MuMu" << std::endl;
	// run ZtoMuMu category
	categoryFlag = "ZtoMuMu";
}

Tvariable_MuMu::~Tvariable_MuMu()
{
	for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "Tvariable_MuMu::~Tvariable_MuMu Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts+1) << std::endl;
	}
	std::cout << "Tvariable_MuMu::~Tvariable_MuMu()" << std::endl;
}
