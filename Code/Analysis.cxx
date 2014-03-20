#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>
#include <cstdlib>
#include <stdlib.h> 

#include <string.h>
#include <istream>
#include <strstream>

#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

#include "HistoConfig.h"
#include "SkimConfig.h"
#include "Ntuple_Controller.h"
#include "Selection_Base.h"
#include "Selection_Factory.h"
#include "DoubleEventRemoval.h"
#include "Parameters.h"
#include "Plots.h"

int main() {
  std::cout << "\n\n\n================" << std::endl;
  cout << "Program is Starting" << endl;
  system("date");
  float neventsproc(0);
  time_t startTime,endTime,beforeLoop,afterLoop;
  time (&startTime);
 
  //TApplication App(argv[0],&argc,argv);
  gROOT->ProcessLine("#include <vector>");
  gSystem->Load("$ROOTSYS/lib/libPhysics.so") ; 
  gROOT->SetBatch(kTRUE);

  ////////////////////////////////////////////////
  // Get Input
  TString base = std::getenv("PWD");
  base+="/";
  std::cout << "Working Dir: " << base << std::endl;
  Parameters Par(base+"Input.txt");
  std::vector<TString> Files, UncertType, UncertList, Analysis;
  std::vector<double> UncertW;
  bool thin, skim;
  int mode, runtype;
  TString mode_str, runType_str, histofile, skimfile,PlotStyle,PlotLabel;
  double Lumi;
  Par.GetVectorString("File:",Files);
  Par.GetBool("Thin:",thin,"False");
  Par.GetBool("Skim:",skim,"False");
  Par.GetString("Mode:",mode_str,"RECONSTRUCT");// RECONSTRUCT/ANALYSIS
  Par.GetString("RunType:",runType_str,"LOCAL");// GRID/LOCAL
  Par.GetString("SkimInfo:",skimfile,"");
  Par.GetString("HistoFile:",histofile,"InputData/Histo.txt");
  Par.GetVectorString("Analysis:",Analysis);
  Par.GetVectorString("UncertType:",UncertType,"<default>");
  Par.GetVectorStringDouble("UncertList:",UncertList,UncertW);
  Par.GetDouble("Lumi:",Lumi,1);
  Par.GetString("PlotStyle:",PlotStyle,"style1");
  Par.GetString("PlotLabel:",PlotLabel,"none");
  /////////////////////////////////////////////////
  // Check Input
  HistoConfig H; 
  SkimConfig SC;
  if(!H.Load(base+histofile)){ 
    cout << "FATAL: No Histogram File or Histograms!!!" << endl;
    exit(6);
  }
  if(!SC.Load(base+skimfile)){
    cout << "FATAL: Invalid Skim File!!!" << endl;
    exit(6);
  }
  if(Analysis.size()==0){
    cout << "FATAL: No Analysis Selected!!!!" << endl;
    exit(6);
  }
  if(UncertType.size()==0){
    cout << "FATAL: Uncerainty Configuration!!!!" << endl;
    exit(6);
  }
  if(Files.size()==0){
    cout << "FATAL: No Input Files!!!!" << endl;
    exit(6);
  }
  mode_str.ToUpper();
  if(mode_str=="ANALYSIS"){
    std::cout << "Using Mode: ANALYSIS" << std::endl;
    mode=Selection_Base::ANALYSIS;
  }
  else if(mode_str=="RECONSTRUCT"){
    std::cout << "Using Mode: RECONSTRUCT" << std::endl;
    mode=Selection_Base::RECONSTRUCT;
  }
  else{
    cout << "FATAL: No Mode Type!!!!" << endl;
    exit(6);
  }
  runType_str.ToUpper();
  if(runType_str=="GRID"){
    std::cout << "Using RunType: GRID" << std::endl;
    runtype=Selection_Base::GRID;
  }
  else if(runType_str=="LOCAL"){
    std::cout << "Using RunType: Local" << std::endl;
    runtype=Selection_Base::Local;
  }
  else{
    cout << "FATAL: No Mode Type!!!!" << endl;
    exit(6);
  }
  Plots P;
  P.Set_Plot_Type(PlotStyle,PlotLabel);
  //////////////////////////////////////////////////
  // Configure Analysis
  Selection_Factory SF;
  vector<Selection_Base*> selections;
  for(int i=0;i<UncertType.size();i++){
    cout << "Configuring systematic " << UncertType.at(i) << endl;
    for(int j=0; j<Analysis.size();j++){
      cout << "Configuring Selection " << Analysis.at(j) << endl;
      selections.push_back(SF.Factory(Analysis.at(j),UncertType.at(i),mode,runtype,Lumi));
    }
  }
  cout << "Selection modules Setup NSelection= " << selections.size() << endl;

  //////////////////////////////////////////////////
  // Run Analysis on Ntuples
  if(mode==Selection_Base::ANALYSIS){
    cout << "Setuping up Ntuple control" << endl;
    Ntuple_Controller Ntp(Files);
    cout << "Ntuple control Setup" << endl;
    if(thin)Ntp.ThinTree();
    if(skim)Ntp.CloneTree("SKIMMED_NTUP");
    for(int j=0; j<selections.size();j++){
      selections[j]->Set_Ntuple(&Ntp);
    }

    //////////////////////////////////////
    // Event Loop
    cout << "Starting Event Loop" << endl;
    Int_t nentries = Ntp.Get_Entries();
    std::vector<TString> ListOfFilesRead;
    std::vector<int> EventsReadFromFile; 
    Int_t RunNumber(-999);
    int num=10000;
    cout << " Number of events in Ntuples " << nentries << endl;
    Int_t i=0; 
    int k=0;
    int p=0;
    DoubleEventRemoval DER;
    bool isok(false); 
    time (&beforeLoop);
    neventsproc=(float)nentries;
    for(i=0;i<nentries ;i++){
      TString fileName=Ntp.Get_File_Name();
      if(i==0 || fileName!=ListOfFilesRead[k]){
	ListOfFilesRead.push_back(fileName);
	EventsReadFromFile.push_back(0);
	k=ListOfFilesRead.size()-1;
      }
      EventsReadFromFile[k]++;
      num++;
      if(num>=10000){
	cout << "Starting event:" << i << " out of " << nentries << endl;
	num=0;
      }
      Ntp.Get_Event(i);
      if(Ntp.isData()){
	if(!DER.CheckDoubleEvents(Ntp.RunNumber(), Ntp.EventNumber())) continue;
	if(RunNumber!=Ntp.RunNumber()){
	  RunNumber=Ntp.RunNumber();
	  //cout << "RunNumber: " << RunNumber << endl;
	}
      }
      bool passed=false;
      for(int j=0; j<selections.size();j++){
	selections.at(j)->Event();
	if(selections.at(j)->Passed()) passed=true;
      }
      if(skim && passed){
	Ntp.AddEventToCloneTree();
	p++;
      }
    }
    time (&afterLoop);
    if(skim) Ntp.SaveCloneTree();
    if(ListOfFilesRead.at(0)==ListOfFilesRead.at(ListOfFilesRead.size()-1)) ListOfFilesRead.pop_back();
    if(i==nentries && Files.size()==ListOfFilesRead.size()) 
      cout << " Number of events Read: " << i << " out of " << nentries << " SUCCESSFULL" <<  endl;
    else 
      cout << " Number of events Read: " << i << " out of " << nentries <<  " Number of files Read: " 
	   << ListOfFilesRead.size() << " out of " << Files.size() << " FAILED" <<  endl;

    cout << " Number of skimmed events in SkimmedD3PD: " << p << endl;
    cout << " Number of files Read: " << ListOfFilesRead.size() << endl;
    for(k=0; k<ListOfFilesRead.size();k++){
      cout << " Number of events: " << EventsReadFromFile[k] << " from file: " << ListOfFilesRead[k] << endl;
    }
    cout << "Event Loop done" << endl;
  }
  ///////////////////////////////////////////
  // Reconstruct Analysis from Stored Histograms
  else if(mode==Selection_Base::RECONSTRUCT){
    cout << "Reconstructing histograms" << endl;
    for(int j=0; j<selections.size();j++){
      selections[j]->LoadResults(Files);
    }
    std::cout << "Loading Files Complete" << std::endl;
    if(runtype==Selection_Base::Local){
      Selection_Factory SF;
      for(int j=0; j<selections.size();j++){
	if(selections.at(j)->Get_SysType()=="default"){
	  Selection_Base *selectionsys;
	  for(int i=0;i<UncertList.size();i++){
	    TString n=selections.at(j)->Get_Name();
	    cout << "Adding Systematic Uncertainty " << UncertType.at(i) << endl;
	    Selection_Base *s=SF.Factory(selections.at(j)->Get_Analysis(),UncertList.at(i),mode,runtype,Lumi);
	    s->LoadResults(Files);
	    selections[j]->EvaluateSystematics(s,1.0);
	    delete s;
	  }
	}
      }
    }
  }
  cout << "Finishing" << endl;
  for(int j=0; j<selections.size();j++){
    cout << "Finishing: " << selections[j]->Get_Name() << endl; 
    selections[j]->Finish();
  }
  cout << "Finished" << endl;

  cout << "Cleaning up" << endl;
  for(int j=0; j<selections.size();j++){
    delete selections[j];
    selections[j]=NULL;
  }
  selections.clear();

  if(mode!=Selection_Base::RECONSTRUCT){
    double dt=difftime (afterLoop,beforeLoop);
    int hour=(int)(dt/(60*60));
    int min=(int)(dt/(60))-60*hour;
    int sec=(int)(dt-min*60-hour*60*60);
    cout << "Time in Ntuple Loop: " << hour << "hr " << min << "min " << sec << "s" << endl;
    float msec=1000*dt/neventsproc;
    cout << "Time in per event in Ntuple Loop: " << msec << "ms" << endl;
  }

  time (&endTime);
  double dt=difftime (endTime,startTime);
  int hour=(int)(dt/(60*60));
  int min=(int)(dt/(60))-60*hour;
  int sec=(int)(dt-min*60-hour*60*60);

  cout << "Run Time: " << hour << "hr " << min << "min " << sec << "s" << endl;
  cout << "Program is Finished" << endl;
  system("date");
}   

