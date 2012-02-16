#include "Parameters.h"

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
#include <math.h>

// Static var
TString Parameters::file="Tools/Par.dat";

Parameters::Parameters(){
}

Parameters::Parameters(TString f){
  SetFile(f);
}

Parameters::~Parameters(){
}


void Parameters::SetFile(TString f){
  file=f;
}

TString Parameters::GetFile(){
  return file;
}


void Parameters::GetString(TString p, TString &v, TString dv){
  return GetParameter(p,v,dv);
}

void Parameters::GetBool(TString p, bool &v, bool dv){
  TString vb;
  TString dvb="false";
  if(dv)dvb="true";
  GetParameter(p,vb,dvb);
  vb.ToLower();
  v=false;
  if(vb==true) v=true;
}

void Parameters::GetInt(TString p, int &v, int dv){
  return GetParameter(p,v,dv);
}

void Parameters::GetDouble(TString p, double &v, double dv){
  return GetParameter(p,v,dv);
}


void Parameters::GetVectorString(TString p, std::vector<TString> &v, TString dv){
  v.clear();
  // Open File
  ifstream input_file;
  char *file_=(char*)file.Data();
  input_file.open(file_, std::ios::in);
  if (!(input_file)){
    std::cout << "\nERROR: Opening xml file "<< file <<" for Parameters has failed.\n" << std::endl;
    return;
  }
  //std::cout << "\nOpened Parameters xml file: "<< file <<".\n" << std::endl;

  std::string s;
  unsigned int a=0;
  while(getline(input_file, s)){
    a++;
    if(a>250) break;
    std::stringstream line(s);
    TString par;
    TString val;
    line >> par >> val;
    par.ToLower();
    p.ToLower();
    if(p.Contains(par) && par.Contains(p)){
      v.push_back(val);
    }
  }
  input_file.close();
  for(int i=0; i<v.size(); i++){
    for(int j=i+1; j<v.size(); j++){
      if(v.at(i)==v.at(j)){
        v.erase(v.begin()+j);
        j--;
      }
    }
  }
  if(dv!="" && v.size()==0) v.push_back(dv);
  for(int i=0; i<v.size();i++){
    std::cout << "Parameters::GetVectorString File=" << file  << " Found: " <<  p << "=" << v.at(i) << std::endl;
  }
  return;
}


template<typename T>
void Parameters::GetParameter(TString p, T &v,T dv){

  // Open File
  ifstream input_file;
  char *file_=(char*)file.Data();
  input_file.open(file, std::ios::in);
  if (!(input_file)){
    std::cout << "\nERROR: Opening xml file "<< file <<" for Parameters has failed.\n" << std::endl;
    return;
  }
  //std::cout << "\nOpened Parameters xml file: "<< file <<".\n" << std::endl;

  std::string s;
  unsigned int a=0;
  while(getline(input_file, s)){
    a++;
    if(a>250) break;
    std::stringstream line(s); 
    TString par;
    T val;
    line >> par >> val; 
    par.ToLower();
    p.ToLower();
    if(p.Contains(par) && par.Contains(p)){
      v=val;
      std::cout << "Parameters::GetParameter File=" << file << " Found: " <<  p << "=" << v << std::endl;
      return;
    }
  }
  v=dv;
  std::cout << "Parameters::GetParameter File=" << file << " Not Found: " <<  p << "=" << v << std::endl;
  input_file.close();
  return;
}


void Parameters::GetVectorStringDouble(TString p, std::vector<TString> &v1, std::vector<double> &v2){
  v1.clear();
  v2.clear();
  // Open File
  ifstream input_file;
  char *file_=(char*)file.Data();
  input_file.open(file_, std::ios::in);
  if (!(input_file)){
    std::cout << "\nERROR: Opening xml file "<< file <<" for Parameters has failed.\n" << std::endl;
    return;
  }
  //std::cout << "\nOpened Parameters xml file: "<< file <<".\n" << std::endl;

  std::string s;
  unsigned int a=0;
  while(getline(input_file, s)){
    a++;
    if(a>250) break;
    std::stringstream line(s);
    TString par;
    TString val1;
    double val2;
    line >> par >> val1 >> val2;
    par.ToLower();
    p.ToLower();
    if(p.Contains(par) && par.Contains(p)){
      v1.push_back(val1);
      v2.push_back(val2);
    }
  }
  input_file.close();
  for(int i=0; i<v1.size(); i++){
    for(int j=i+1; j<v1.size(); j++){
      if(v1.at(i)==v1.at(j)){
        v1.erase(v1.begin()+j);
	v2.erase(v2.begin()+j);
        j--;
      }
    }
  }
  for(int i=0; i<v1.size();i++){
    std::cout << "Parameters::GetVectorStringDouble File=" << file  << " Found: " <<  p << "=" << v1.at(i) << " " << v2.at(i) << std::endl;
  }
  return;
}
