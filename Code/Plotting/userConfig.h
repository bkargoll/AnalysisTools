#ifndef userConfig_h
#define userConfig_h

/////////////////////////////////////////////////
//
// userConfig: main place for user input
// define all needed stuff here
//
/////////////////////////////////////////////////

#include "configInfo.h"

// flags for testing
bool verbose = false;
bool testPlotting = false;

// draw signal stacked on top (or overlaying)
bool signaltop = false;

// create configInfo struct
double lumi = 19712.;
//TString infile = "/net/scratch_cms/institut_3b/kargoll/Taus30GeV/Categories/workdirAnalysis_Nov_04_2014/LOCAL_COMBINED_inclusive_default_LumiScaled.root";
//TString ident = "inclusive_default_";
TString infile = "/net/scratch_cms/institut_3b/kargoll/Taus30GeV/Backgrounds/workdirAnalysis_Oct_31_2014/LOCAL_COMBINED_htotaumutauhbackgrounds_default_LumiScaled.root";
TString ident = "htotaumutauhbackgrounds_default_";
bool isLumiScaled = true;
// define user-specific config
configInfo conf(infile, isLumiScaled, ident, lumi);


#endif
