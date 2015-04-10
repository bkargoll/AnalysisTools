#include "SimpleFits/FitSoftware/interface/Logger.h"
#include "DataStorage.h"
#include "Parameters.h"
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <fstream>
#include "TString.h"
#include "TSystem.h"
#include <iostream>
#include <time.h>

DataStorage::DataStorage(){
  mydir= getenv ("PWD");
  inputFileName = "dataStorageInput_temp_";
}

DataStorage::~DataStorage(){
	for (unsigned int i_file = 0; i_file < filesToDelete.size(); i_file++){
		TString cmd = "rm " + filesToDelete.at(i_file);
		Logger(Logger::Info) << "calling shell command: " << cmd << std::endl;
		system(cmd.Data());
	}
}

int DataStorage::GetFile(TString key){
  Parameters Par; // assumes configured in Analysis.cxx
  TString gridsite="none";
  std::vector<TString> Files;
  Par.GetVectorString(key,Files);
  Par.GetString("GRIDSite:",gridsite);

  if(gridsite=="none" || Files.size()==0) return 0;
  
  for(unsigned int i=0;i<Files.size();i++){
    TString inFile=inputFileName; inFile+=i; inFile+=".root";
    TString cmd1= "srmcp srm://" + gridsite + ":8443/" + Files.at(i) + " file:////" + mydir + "/" + inFile;
    Logger(Logger::Info) << "calling shell command: "<< cmd1 << std::endl;
    system(cmd1.Data());
    filesToDelete.push_back(mydir + "/" + inFile);
    ifstream f(inFile);
    if(!f) return 0;
  }
  return Files.size();
}

void DataStorage::StoreFile(TString File, TString savedFile){
  // put timestamp in output file
  savedFile.ReplaceAll("{DATE}", timeStamp(true, false));
  savedFile.ReplaceAll("{TIME}", timeStamp(false, true));

  // compile path of output file
  TString gridsite;
  Parameters Par; // assumes configured in Analysis.cxx
  Par.GetString("GRIDSite:",gridsite);
  TString outfile = " srm://" + gridsite + ":8443" + savedFile;
  TString outpath = gSystem->DirName(outfile);

  // create output folder, if neccessary
  TString command;
  command = "srmls -count=1 " + outpath;
  Logger(Logger::Info) << "calling shell command: " << command << std::endl;
  int isDir =  system( (command).Data() );
  if (isDir != 0)	{command = "srmmkdir " + outpath; Logger(Logger::Info) << "calling shell command: " << command << std::endl; system( (command).Data() ); }

  TString cmd1 = "srmcp file:////" + mydir + "/" + File + outfile;
  Logger(Logger::Info) << "calling shell command: " << cmd1 << std::endl;
  system(cmd1.Data());
}

// get a TString from the current time
// stolen from http://www.cplusplus.com/reference/ctime/strftime/
TString DataStorage::timeStamp(bool Date, bool Time){
	if (!Date && !Time){
		return "";
	}
	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];

	time (&rawtime);
	timeinfo = localtime (&rawtime);

	// output pattern
	TString pattern = "";
	if (Date)			pattern += "%F";
	if (Date && Time)	pattern += "_";
	if (Time) 			pattern += "%H-%M";

	strftime (buffer,80,pattern.Data(),timeinfo);
	return TString(buffer);
}
