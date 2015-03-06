#include "SkimConfig.h"

#include <cstdlib>
#include <algorithm>                                     
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include "Riostream.h"
#include <fstream>
#include <istream>
#include <strstream>
#include <cstdlib>
#include "TH1D.h"
#include <math.h>

#include "HistoConfig.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"

std::vector<int64_t> SkimConfig::SkimIDs;
std::vector<float>   SkimConfig::NEvents;
std::vector<float>   SkimConfig::NEventsErr;
std::vector<float>   SkimConfig::NEvents_sel;
std::vector<float>   SkimConfig::NEventsErr_sel;
std::vector<float>   SkimConfig::NEvents_noweight;
std::vector<float>   SkimConfig::NEvents_noweight_sel;
bool SkimConfig::loaded=false;
bool SkimConfig::converted=false;

SkimConfig::SkimConfig()
{
}

bool SkimConfig::Load()
{
  return SkimConfig::Load("InputData/SkimSummary.log");
}

bool SkimConfig::Load(TString Name_)
{
  std::cout << "SkimConfig::LoadSkimEff("<< Name_ <<")" << std::endl;
  if(loaded) return false;
  SkimIDs.clear();
  NEvents.clear();
  NEventsErr.clear();
  NEvents_sel.clear();
  NEventsErr_sel.clear();
  NEvents_noweight.clear();
  NEvents_noweight_sel.clear();
  ifstream input_file;
  char *file_=(char*)Name_.Data();
  input_file.open(file_, std::ios::in);
  if (!(input_file)){
    std::cout << "\nERROR: Opening SkimEff file "<< Name_ <<" for SkimConfig has failed.\n" << std::endl;
    return false;
  }
  std::cout << "\nOpened SkimConfig SkimEff file: "<< Name_ <<".\n" << std::endl;
  loaded =true;
  std::string s;
  int a=0;
  while(getline(input_file, s)){
    a++;
    if(a>10000){
    	std::cout << "ERROR: SkimConfig has more than 10000 lines. This seems weird. Breaking..." << std::endl;
    	break;
    }

    std::stringstream line(s);
    TString tmp;
    int64_t id;
    float nevents;
    float neventserr;
    float nevents_sel;
    float neventserr_sel;
    float noweight;
    float noweight_sel;
    line >> tmp >> id >> tmp >> nevents >> tmp >> neventserr >> tmp >> nevents_sel >> tmp >> neventserr_sel >> tmp >> noweight >> tmp >> noweight_sel;

    // adapt for masses in DataMCType:
    // in case a DataMCType of length 8 (i.e. without mass digits) is given,
    // 3 zeroes are added for the mass digits
    if( id >= 100000 && id < 1000000000){
		// shift DataMCType 3 digits to the left
		int64_t jakid = id - (id%100);
		jakid *= 1000;
		id = jakid + (id%100);
    }

    SkimIDs.push_back(id);
    NEvents.push_back(nevents);
    NEventsErr.push_back(neventserr);
    NEvents_sel.push_back(nevents_sel);
    NEventsErr_sel.push_back(neventserr_sel);
    NEvents_noweight.push_back(noweight);
    NEvents_noweight_sel.push_back(noweight_sel);
  }
  input_file.close();

  CovertToHistoFormat();

  return true;
}

SkimConfig::~SkimConfig(){

}


double SkimConfig::GetNEvents(int64_t id){
  for(unsigned int i=0; i<SkimIDs.size();i++){
    if(id==SkimIDs.at(i))return NEvents.at(i);
  }
  return 0;
}

void SkimConfig::SaveEfficiency(TString Name,std::vector<int64_t> ids,std::vector<TH1D> NPassed, std::vector<TH1D> NPassed_noweight){
  if(!loaded){
    std::cout << "Error SkimConfig::SaveEfficiency input not loaded -> no skim summary " << std::endl;
    return;
  }
  if(!converted){
    if(!CovertToHistoFormat()) {
      std::cout << "Error SkimConfig::SaveEfficiency Input not mapped to Histograms " << std::endl;
      return;
    }
  }


  ofstream output;
  output.open(Name+"SkimEff.dat", std::ios::out);
  int nbins=NPassed.at(0).GetNbinsX();
  for(unsigned int i=0; i<ids.size();i++){
    double Eff(0),Eff_w(0);
    if(NPassed.at(i).GetBinContent(1)>0) Eff=NPassed.at(i).GetBinContent(nbins)/NPassed.at(i).GetBinContent(1);
    if(NPassed_noweight.at(i).GetBinContent(1)>0) Eff_w=NPassed_noweight.at(i).GetBinContent(nbins)/NPassed_noweight.at(i).GetBinContent(1);
    (output) << "ID= "             << ids.at(i)  << setprecision (15)
	     << "AllEvt= "         << NPassed.at(i).GetBinContent(1)           
	     << "AllEvtErr= "      << NPassed.at(i).GetBinError(1)
	     << "SelEvt= "         << NPassed.at(i).GetBinContent(nbins)       
	     << "SelEvtErr= "      << NPassed.at(i).GetBinError(nbins) 
	     << "AllEvtnoweight= " << NPassed_noweight.at(i).GetBinContent(1) 
	     << "SelEvtnoweight= " << NPassed_noweight.at(i).GetBinContent(nbins)
	     << "Eff(weight)= "    << Eff
	     << "Eff(noweight)= "  << Eff_w
	     << std::endl;
  }
}

void SkimConfig::ApplySkimEfficiency(std::vector<int64_t> ids,std::vector<TH1D> &NPassed, std::vector<TH1D> &NPassed_noweight){
  if(!loaded){
    std::cout << "Error SkimConfig::ApplySkimEfficiency input not loaded -> no skim summary " << std::endl;
    return;
  }
  if(!converted){
    if(!CovertToHistoFormat()) {
      std::cout << "Error SkimConfig::ApplySkimEfficiency Input not mapped to Histograms " << std::endl;
      return;
    }
  }
  for(unsigned int i=0; i<ids.size();i++){
    NPassed.at(i).SetBinContent(0,0);
    NPassed.at(i).SetBinError(0,0);
    NPassed_noweight.at(i).SetBinContent(0,0);
    NPassed_noweight.at(i).SetBinError(0,0);
  }
  for(unsigned int i=0; i<SkimIDs.size();i++){
    NPassed.at(i).SetBinContent(0,NEvents.at(i));
    NPassed.at(i).SetBinError(0,NEventsErr.at(i));
    NPassed_noweight.at(i).SetBinContent(0,NEvents_noweight.at(i));
    NPassed_noweight.at(i).SetBinError(0,NEvents_noweight.at(i)*NEvents_noweight.at(i));
  }
}




void SkimConfig::CheckNEvents(std::vector<int64_t> ids, std::vector<float> nevts){
  if(!loaded){
    std::cout << "Error SkimConfig::CheckNEvents input not loaded -> no skim summary " << std::endl;  
    return;
  }
  if(!converted){
    if(!CovertToHistoFormat()) {
      std::cout << "Error SkimConfig::CheckNEvents Input not mapped to Histograms " << std::endl;
      return;
    }
  }

  ofstream output;
  output.open("DQM.tex", ios::out);
  
  output << "\\section{Data Quality Table}" << std::endl;
  output << "\\begin{landscape}" << std::endl;
  output << "\\begin{table}[t]" << std::endl;
  output << "\\tiny " << std::endl;
  output << "\\begin{center}" << std::endl;
  output << "\\begin{tabular}{|p{3.0cm}|p{3.0cm}|p{3.0cm}|p{3.0cm}|p{3.0cm}|} \\hline" << endl; 
  output << "ID & DQM Status & Events Read & Events Expected & Ratio \\\\  \\hline " << std::endl;
  for(unsigned int i=0; i<ids.size();i++){
    if(fabs(nevts.at(i)-NEvents_noweight_sel.at(i))>0.0001){
      output << ids.at(i) << " & \\textcolor{red}{FAILED} & ";
    }
    else{
      output << ids.at(i) << " & \\textcolor{green}{PASSED} & ";
    }
    output <<  nevts.at(i) << " & " << NEvents_noweight_sel.at(i) << " & " <<   nevts.at(i)/NEvents_noweight_sel.at(i) 
	   << " \\\\  " << std::endl;
  }
  output << " \\hline"                                      << std::endl;
  output << "\\end{tabular}"                                << std::endl;
  output << "\\caption[Data Quality]{Data Quality}"         << std::endl;
  output << "\\end{center}"                                 << std::endl;
  output << "\\end{table}"                                  << std::endl;
  output << "\\normalsize"                                  << std::endl;
  output << "\\end{landscape}"                              << std::endl;
  output.close();

  return;
}



bool SkimConfig::CovertToHistoFormat(){
  if(!loaded){ std::cout << "SkimConfig::CovertToHistoFormat() ERROR SkimConfig has not been loaded" << std::endl; return false;}
  if(converted) return true;

  converted=true;
  HistoConfig H;
  std::vector<bool>    IDFlag(SkimIDs.size(),false);
  std::vector<int64_t>     SkimIDs_new;
  std::vector<float>   NEvents_new;
  std::vector<float>   NEventsErr_new;
  std::vector<float>   NEvents_sel_new;
  std::vector<float>   NEventsErr_sel_new;
  std::vector<float>   NEvents_noweight_new;
  std::vector<float>   NEvents_noweight_sel_new;

  for(unsigned int i=0; i<H.GetNHisto(); i++){
    SkimIDs_new.push_back(H.GetID(i));
    NEvents_new.push_back(0);
    NEventsErr_new.push_back(0);
    NEvents_sel_new.push_back(0);
    NEventsErr_sel_new.push_back(0);
    NEvents_noweight_new.push_back(0);
    NEvents_noweight_sel_new.push_back(0);
  }

  for(unsigned int i=0;i<SkimIDs_new.size();i++){
    for(unsigned int j=0;j<SkimIDs.size();j++){
      if(SkimIDs.at(j)==SkimIDs_new.at(i)){
	std::cout << "SkimConfig::CovertToHistoFormat() Found Full Mask: " << SkimIDs.at(j) << std::endl;
	IDFlag.at(j)=true;
	NEvents_new.at(i)=NEvents.at(j);
	NEventsErr_new.at(i)=NEventsErr.at(j);
	NEvents_sel_new.at(i)=NEvents_sel.at(j);
	NEventsErr_sel_new.at(i)=NEventsErr_sel.at(j);
	NEvents_noweight_new.at(i)=NEvents_noweight.at(j);
	NEvents_noweight_sel_new.at(i)=NEvents_noweight_sel.at(j);
      }
    }
  }
 
  for(unsigned int i=0;i<SkimIDs_new.size();i++){
    for(unsigned int j=0;j<SkimIDs.size();j++){
      if(!IDFlag.at(j)){
	if(SkimIDs.at(j)%100000==SkimIDs_new.at(i)){
	  std::cout << "SkimConfig::CovertToHistoFormat() Found Master Decay: " << SkimIDs.at(j) << std::endl;
	  IDFlag.at(j)=true;
	  NEvents_new.at(i)+=NEvents.at(j);
	  NEventsErr_new.at(i)+=sqrt(NEventsErr.at(j)*NEventsErr.at(j)+NEventsErr_new.at(i)*NEventsErr_new.at(i));
	  NEvents_sel_new.at(i)+=NEvents_sel.at(j);
	  NEventsErr_sel_new.at(i)+=sqrt(NEventsErr_sel.at(j)*NEventsErr_sel.at(j)+NEventsErr_sel_new.at(i)*NEventsErr_sel_new.at(i));
	  NEvents_noweight_new.at(i)+=NEvents_noweight.at(j);
	  NEvents_noweight_sel_new.at(i)+=NEvents_noweight_sel.at(j);
	}
      }
    }
  }

  for(unsigned int i=0;i<SkimIDs_new.size();i++){
    if(SkimIDs_new.at(i)==DataMCType::Signal){
      for(unsigned int j=0;j<SkimIDs.size();j++){
	std::cout << "SkimConfig::CovertToHistoFormat() Found Signal Match " << SkimIDs.at(j) << std::endl;
	if(SkimIDs.at(j)==DataMCType::DY_ll_Signal){
	  IDFlag.at(j)=true;
	  NEvents_new.at(i)=NEvents.at(j);
	  NEventsErr_new.at(i)=NEventsErr.at(j);
	  NEvents_sel_new.at(i)=NEvents_sel.at(j);
	  NEventsErr_sel_new.at(i)=NEventsErr_sel.at(j);
	  NEvents_noweight_new.at(i)=NEvents_noweight.at(j);
	  NEvents_noweight_sel_new.at(i)=NEvents_noweight_sel.at(j);
	}
      }
    }
  }

  for(unsigned int j=0;j<SkimIDs.size();j++){
    if(!IDFlag.at(j)){
      std::cout << "SkimConfig::CovertToHistoFormat() WARNING unmatched DataMCType: " << SkimIDs.at(j) << std::endl;
    }
  }

  // now set as default
  SkimIDs=SkimIDs_new;
  NEvents=NEvents_new;
  NEventsErr=NEventsErr_new;
  NEvents_sel=NEvents_sel_new;
  NEventsErr_sel=NEventsErr_sel_new;
  NEvents_noweight=NEvents_noweight_new;
  NEvents_noweight_sel=NEvents_noweight_sel_new;
  return true;
}
