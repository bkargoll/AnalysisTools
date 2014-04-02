#! /bin/bash
echo 'Run all the categories' 
cd OneJetBoost; source OneJetBoost.sh; cd ..
cd OneJetHigh; source OneJetHigh.sh; cd ..
cd OneJetLow; source OneJetLow.sh; cd ..
cd VBFLoose; source VBFLoose.sh; cd ..
cd VBFTight; source VBFTight.sh; cd ..
cd ZeroJetHigh; source ZeroJetHigh.sh; cd ..
cd ZeroJetLow; source ZeroJetLow.sh; cd ..
echo 'All categories have been executed. Done.'