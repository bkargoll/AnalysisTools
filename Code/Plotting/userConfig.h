#ifndef userConfig_h
#define userConfig_h

/////////////////////////////////////////////////
//
// userConfig: main place for user input
// define all needed stuff here
//
/////////////////////////////////////////////////

#include "configInfo.h" // todo: remove configInfo struct completely?

// flags for testing
bool verbose = true;
bool testPlotting = false;

// draw signal stacked on top (or overlaying)
bool signaltop = false;

// create configInfo struct
double lumi = 19712.;
TString infile = "/net/scratch_cms/institut_3b/kargoll/FirstRunWithQCD/Backgrounds/workdirAnalysis_Oct_23_2014/LOCAL_COMBINED_htotaumutauhbackgrounds_default_LumiScaled.root";
bool isLumiScaled = true;
// define user-specific config
configInfo conf(infile, isLumiScaled, "htotaumutauhbackgrounds_default_", lumi);



#endif
