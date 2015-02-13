#include "Tvariable_EE.h"

Tvariable_EE::Tvariable_EE(TString Name_, TString id_):
Tvariable_Base(Name_,id_)
{
	std::cout << "Setting up the class Tvariable_EE" << std::endl;
	// run ZtoEE category
	categoryFlag = "ZtoEE";
}

Tvariable_EE::~Tvariable_EE()
{
	for(unsigned int j=0; j<Npassed.size(); j++){
	    std::cout << "Tvariable_EE::~Tvariable_EE Selection Summary before: "
		 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
		 << Npassed.at(j).GetBinContent(NCuts+1) << " +/- " << Npassed.at(j).GetBinError(NCuts+1) << std::endl;
	}
	std::cout << "Tvariable_EE::~Tvariable_EE()" << std::endl;
}
