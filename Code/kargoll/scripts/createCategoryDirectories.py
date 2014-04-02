#!/usr/bin/env python

import os

channels = ["VBFTight", "VBFLoose", "OneJetHigh", "OneJetLow", "OneJetBoost", "ZeroJetHigh", "ZeroJetLow"]

for chan in channels:
    os.mkdir(chan)
    os.system("cp Input.txt "+chan)
    os.system("cp Set_1/Set_1.sh "+chan+"/"+chan+".sh")
    os.chdir(chan)
    os.mkdir("EPS")
    os.system("ln -s "+os.path.abspath("../Code/InputData")+" InputData")
 
    ############ modify Input.txt
    file = open("Input.txt","r")
    input = file.readlines()
    file.close()
        
    for i in range(len(input)):
        if ("Skim" in input[i]) and ("True" in input[i]):
            input[i] = "Skim: False\n"
        if "RECONSTRUCT" in input[i]:
            input[i] = "Mode: ANALYSIS\n"
        if "Analysis: " in input[i]:
            input[i] = "Analysis: "+chan+"\n"
        if ("File:" in input[i]) and not ("HistoFile" in input[i]):
            input[i] = "File: /user/kargoll/Analysis/workdirAnalysis_Mar_08_2014/SKIMMED_NTUP_combined.root\n"
    
    printFile = True
    file = open("Input.txt","w")
    for line in input:
        if ("File:" in line) and not ("HistoFile" in line):
            if printFile:
                file.write(line)
                printFile = False
        else:
            file.write(line)
    file.close()
    
    ############ Modify shell script
    file = open(chan+".sh","r")
    data = file.readlines()
    file.close()
    
    file = open(chan+".sh","w")
    for line in data:
        if not (("Set_1-get.sh" in line) or ("Set_1-clean.sh" in line)):
            file.write(line.replace('Set_1', chan))
    
    os.chdir("..")
             
    
print "Directory structure is set up."