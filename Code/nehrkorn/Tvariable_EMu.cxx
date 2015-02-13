#include "Tvariable_EMu.h"

Tvariable_EMu::Tvariable_EMu(TString Name_, TString id_):
Tvariable_Base(Name_,id_)
{
	std::cout << "Setting up the class Tvariable_EMu" << std::endl;
	// run ZtoEMu category
	categoryFlag = "ZtoEMu";
}

Tvariable_EMu::~Tvariable_EMu()
{
	for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "Tvariable_EMu::~Tvariable_EMu Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts+1) << std::endl;
	}
	std::cout << "Tvariable_EMu::~Tvariable_EMu()" << std::endl;
}
