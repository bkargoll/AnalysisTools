#! /usr/bin/perl
use Cwd;
use POSIX;
use POSIX qw(strftime);

#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];

$UserID= POSIX::cuserid();
$UserIDCern=$UserID;
$UserDir="";
if($UserID eq "nugent"){
    $UserIDCern="inugent";
    $UserDir="--inugent";
}
if($UserID eq "cherepanov"){
    $UserIDCern="cherepan";
    $UserDir="--cherepanov";
}
if($UserID eq "nehrkorn"){
    $UserIDCern="anehrkor";
    $UserDir="--nehrkorn";
}
if($UserID eq "kargoll"){
    $UserIDCern="bkargoll";
    $UserDir="--kargoll";
}
if($UserID eq "pistone"){
    $UserDir="--pistone";
}



#Default vaules
$InputDir="/net/scratch_cms/institut_3b/$UserID/Test";
$OutputDir="/net/scratch_cms/institut_3b/$UserID";
#$OutputDir="~/";
$CodeDir="../Code";
$set="Analysis_";
$CMSSWRel="5_3_14_patch2";
$PileupVersion="V08-03-17";
$tag="03-00-12";
$TauReco="5_2_3_patch3_Dec_08_2012";
$BTag="NO";
$Cleaning ="NO";
$maxdata=100;
$maxmc=100;

if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\nThis code requires one input option. The systax is:./todo_Grid.pl [OPTION]");
    printf("\nPlease choose from the following options:\n");
    printf("\n./todo.pl --help                                   Prints this message\n");
    printf("\n./todo.pl --TauNtuple <TauNtupleDir>               Setups up TauNtuple and gives instructions for submission");
    printf("\n                                                   <TauNtupleDir> location of CMSSW for TauNtuple.");
    printf("\n                                                   Optional Commmads: ");
    printf("\n                                                     --CMSSWRel <CMSSW release #> The CMSSW release you want to use Default: $CMSSWRel");
    printf("\n                                                     --PileupVersion <Version # ie V08-03-17>");
    printf("\n                                                     --tag <tag for TauNuptle/SkimProduction/TiggerFiliter/...>");
    printf("\n                                                     --TauReco <tag for Recommended TauReco> current option: na, 4_4_Y_08_03_2012 or 5_2_3_patch3_Dec_08_2012");
    printf("\n                                                     --BTag <YES/NO> Default: NO");
    printf("\n./todo.pl --Local <Input.txt>                      INTENTED FOR SMALL SCALE TESTS ONLY");  
    printf("\n                                                   Configure a directory to run locally. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --InputDir <InputDir>   Default value: $InputDir"); 
    printf("\n                                                       Note: the root files must be inside a $InputDir/<dir>/ ");
    printf("\n                                                       where dir contains user, data or mc as a substring");
    printf("\n                                                     --OutputDir <OutputDir> Default value: $OutputDir");
    printf("\n                                                     --CodeDir  <CodeDir>    Default value: $CodeDir");
    printf("\n                                                     --SetName <SetName>     Default value: $set ");
    printf("\n                                                     --NMaxData <Max Number of data files per job >     Default value: $maxdata ");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job >     Default value: $maxmc ");
    printf("\n./todo.pl --DCache <Input.txt> <ListofDS.txt>      INTENTED FOR REGULAR USE (DEFAULT)");
    printf("\n                                                   Configure a directory to run from. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   <ListoDS.txt> list of DCache Dataset directories you want to run on.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --OutputDir <OutputDir> ");
    printf("\n                                                     --CodeDir <CodeDir>");
    printf("\n                                                     --SetName <SetName> "); 
    printf("\n                                                     --NMaxData <Max Number of data files per job >     Default value: $maxdata ");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job >     Default value: $maxmc ");
    printf("\n  ");
    printf("\n./todo.pl --GRID <Input.txt> <ListofDS.txt>        ALTERNATIVE FOR REGULAR USE");
    printf("\n                                                   Configure a directory to run from. <InputPar.txt> name of file that");
    printf("\n                                                   contains input command template.");
    printf("\n                                                   <ListoDS.txt> list of DCache Dataset directories you want to run on.");
    printf("\n                                                   Optional commands:  ");
    printf("\n                                                     --OutputDir <OutputDir> ");
    printf("\n                                                     --CodeDir <CodeDir>");
    printf("\n                                                     --SetName <SetName> ");
    printf("\n                                                     --NMaxData <Max Number of data files per job > Default value: $maxdata");
    printf("\n                                                     --NMaxMC <Max Number of MC files per job > Default value: $maxmc ");
    printf("\n                                                     --BuildRoot <ROOT Version> builds custom version for root instead of copying lib+include");
    printf("\n                                                     --ROOTSYS <ROOTSYS> the current ROOTSYS variable if --BuildRoot is not defined");
    printf("\n                                                     --GRIDSite <site> the grid site you wish to run on. Default=grid-srm.physik.rwth-aachen.de \n ");
    exit(0);  
} 

######################################
$InputFile=$ARGV[1];
$buildRoot=0;
$hasroot=0;
$gridsite="grid-srm.physik.rwth-aachen.de";
for($l=2;$l<$numArgs; $l++){
    if($ARGV[$l] eq "--InputDir"){
	$l++;
	$InputDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--OutputDir"){
	$l++;
	$OutputDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--CodeDir"){
	$l++;
	$CodeDir=$ARGV[$l];
    }
    if($ARGV[$l] eq "--SetName"){
	$l++;
	$set=$ARGV[$l];
    }
    if($ARGV[$l] eq "--CMSSWRel"){
        $l++;
        $CMSSWRel=$ARGV[$l];
    }
    if($ARGV[$l] eq "--PileupVersion"){
       $l++;
        $PileupVersion=$ARGV[$l];
    }
    if($ARGV[$l] eq "--tag"){
        $l++;
        $tag=$ARGV[$l];
    }
    if($ARGV[$l] eq "--TauReco"){
        $l++;
        $TauReco=$ARGV[$l];
    }
    if($ARGV[$l] eq "--BTag"){
	$l++;
	$BTag=$ARGV[$l];
    }
    if($ARGV[$l] eq "--Cleaning"){
	$l++;
	$Cleaning =$ARGV[$l];
    }
    if($ARGV[$l] eq  "--NMaxData" ){
	$l++;
	$maxdata=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--NMaxMC" ){
	$l++;
	$maxmc=$ARGV[$l];
    }
    if($ARGV[$l] eq  "--BuildRoot"){
	$l++;
	$buildrootversion=$ARGV[$l];
	$buildRoot=1;
    } 
    if($ARGV[$l] eq  "--ROOTSYS"){
        $l++;
        $MYROOTSYS=$ARGV[$l];
	$hasroot=1;
    }
    if($ARGV[$l] eq  "--GRIDSite"){
	$l++;
	$gridsite=$ARGV[$l];;
    }
}

$time= strftime("%h_%d_%Y",localtime);
$temp= $set . $time;
$set=$temp;

if( $ARGV[0] eq "--TauNtuple"){
    $currentdir=getcwd;
    if($ARGV[1] ne ""){
	$basedir=$ARGV[1];
    }
    else{
	printf("\nWorkingDir for CMSSW is required. Please follow the syntax:./todo.pl --TauNtuple <TauNtupleDir> ");
	printf("\nFor more details use: ./todo --help\n"); 
	exit(0);
    }

    printf("\nWorkingDir for CMSSW: $basedir");
    printf("\nCurrentDir is: $currentdir");
    printf("\nUsing CMSSW Release: $CMSSWRel");

    # setup CMSSW
    system(sprintf("rm Install_TauNtuple_$CMSSWRel-$time; rm Setup_$CMSSWRel-$time "));
    system(sprintf("echo \"rm $basedir/TauNtuple_$CMSSWRel-$time -rf; mkdir $basedir/TauNtuple_$CMSSWRel-$time\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cd $basedir/TauNtuple_$CMSSWRel-$time\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"export SCRAM_ARCH=\\\"slc5_amd64_gcc462\\\"\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"export SCRAM_ARCH=\\\"slc5_amd64_gcc462\\\"\" >> Setup_$CMSSWRel-$time"));
    
    system(sprintf("echo \"cd $basedir/TauNtuple_$CMSSWRel-$time\" >> Setup_$CMSSWRel-$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Setup_$CMSSWRel-$time"));
    
    system(sprintf("echo \"cmsrel CMSSW_$CMSSWRel\" >> Install_TauNtuple_$CMSSWRel-$time")); 
    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Install_TauNtuple_$CMSSWRel-$time")); 
    system(sprintf("echo \"cmsenv\" >> Install_TauNtuple_$CMSSWRel-$time")); 

    system(sprintf("echo \"cd CMSSW_$CMSSWRel/src\" >> Setup_$CMSSWRel-$time"));
    system(sprintf("echo \"cmsenv\" >> Setup_$CMSSWRel-$time"));

    system(sprintf("echo \"git init \" >> Setup_$CMSSWRel-$time"));

	# for SimpleFits
	system(sprintf("echo \"git cms-addpkg Validation/EventGenerator \" >> Install_TauNtuple_$CMSSWRel-$time"));
	
    # MVA-MET recipe for CMSSW_5_3_14_patchX (https://twiki.cern.ch/twiki/bin/viewauth/CMS/MVAMet#Installation)
    system(sprintf("echo \"git cms-addpkg PhysicsTools/PatAlgos \" >> Install_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"git cms-merge-topic cms-analysis-tools:5_3_14-updateSelectorUtils\" >> Install_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"git cms-merge-topic cms-analysis-tools:5_3_13_patch2-testNewTau\" >> Install_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"git cms-merge-topic -u TaiSakuma:53X-met-131120-01\" >> Install_TauNtuple_$CMSSWRel-$time"));
	system(sprintf("echo \"git-cms-merge-topic -u cms-met:53X-MVaNoPuMET-20131217-01\" >> Install_TauNtuple_$CMSSWRel-$time"));

	# ATTENTION: CMSSW_X_Y_Z/src must be completely EMPTY in order to run "git cms-addpkg"
    system(sprintf("echo \"mkdir data \" >> Install_TauNtuple_$CMSSWRel-$time"));

   #Tau Package recomendation for 5_3_X (X>=12): nothing to do
    
    # embedding: seems that there is nothing to do

    # Electron Tools
    system(sprintf("echo \"git cms-addpkg EgammaAnalysis/ElectronTools\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cd EgammaAnalysis/ElectronTools/data/\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cat download.url | xargs wget\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cd ../../../\" >> Install_TauNtuple_$CMSSWRel-$time"));
    
    # PUJetID  (documentation: https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID
    #			recipe: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolJetMetAnalysis#Jet_Recipe
    #			hypernews discussion: https://hypernews.cern.ch/HyperNews/CMS/get/met/333/1/1.html)
    system(sprintf("echo \"git clone https://github.com/violatingcp/Jets_Short.git\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cp -r Jets_Short/* .\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"rm -rf Jets_Short\" >> Install_TauNtuple_$CMSSWRel-$time"));

	# Ntuple code
    system(sprintf("echo \"git clone https://github.com/inugent/TauDataFormat TauDataFormat; cd TauDataFormat; git checkout; cd ../; \" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"git clone https://github.com/inugent/SkimProduction SkimProduction; cd SkimProduction; git checkout; cd ../; \" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"cp SkimProduction/CRAB/*.root data/ \" >> Install_TauNtuple_$CMSSWRel-$time"));
    # SimpleFits
    system(sprintf("echo \"git clone https://github.com/inugent/SimpleFits SimpleFits; cd SimpleFits; git checkout; cd ../; \" >> Install_TauNtuple_$CMSSWRel-$time"));

    # Setup CRAB
    system(sprintf("echo \"export VO_CMS_SW_DIR=\\\"/net/software_cms\\\"\" >> Install_TauNtuple_$CMSSWRel-$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Install_TauNtuple_$CMSSWRel-$time "));
    system(sprintf("echo \"source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh\" >> Install_TauNtuple_$CMSSWRel-$time"));

    system(sprintf("echo \"export VO_CMS_SW_DIR=\\\"/net/software_cms\\\"\" >> Setup_$CMSSWRel-$time"));
    system(sprintf("echo \"source /cvmfs/cms.cern.ch/cmsset_default.sh\" >> Setup_$CMSSWRel-$time "));
    system(sprintf("echo \"source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh\" >> Setup_$CMSSWRel-$time"));

    #build
    system(sprintf("echo \"scram b -j 4 \" >> Install_TauNtuple_$CMSSWRel-$time"));

    # print Instructions
    printf("\n\nInstructions");
    printf("\ngit config --global credential.helper 'cache --timeout=3600'");
    printf("\ngit config --global credential.helper cache");
    printf("\ngit config --global user.github $UserIDCern ");
    #printf("\nkinit $UserIDCern ");
    printf("\nsource Install_TauNtuple_$CMSSWRel-$time"); 
    printf("\nYou can now go to the Production Directory and submit jobs.");
    printf("\ncd  $basedir/SkimProduction/CRAB/");
    printf("\n1) Read the instruction for todo.pl: ./todo.pl --help");
    printf("\n2) Setup the jobs: ./todo.pl --Submit HTL_Tau_Ntuple_cfg.py Input.dat ");
    printf("\n3) Setup the jobs and submit 1 Test job: ./todo.pl --Submit HTL_Tau_Ntuple_cfg.py Input.dat --n 1");
    printf("\n4) Setup the jobs and start production:./todo.pl --Submit HTL_Tau_Ntuple_cfg.py Input.dat --n -1\n");
    printf("\nIf you want to use this release in a new window: source Setup_$CMSSWRel-$time \n");
}


if( $ARGV[0] eq "--Local" ){
    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Input directory will be:  $InputDir \n");
    printf("Code Repository is:       $CodeDir \n");
    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    system(sprintf("cd $OutputDir"));
    system(sprintf("rm -rf $OutputDir/workdir$set \n"));
    system(sprintf("mkdir $OutputDir/workdir$set "));
    printf("Cleaning complete \n");
    
    # create directory stucture
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/i386_linux "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));
    
    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source config \\\$\@ \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Combine")); 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Setup Condor Combine scripts
    system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 			
    system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"rm Set*/*.o; rm Set*/*.e; rm Set*/*.log; \" >> $OutputDir/workdir$set/Submit")) ;

    opendir(DIR,"$InputDir/");
    @dirs = grep {( /user/ || /data/ || /mc/)} readdir(DIR);
    closedir DIR;
    
    $B=0;
    $max=1;
    for($l=0;$l<2; $l++){
	printf("\n\nStarting Loop $l  maxdata = $maxdata maxmc = $maxmc \n");
	$A=$maxdata+$maxmc+10;
	foreach $subdir (@dirs){
	    printf("$subdir");
	    if(($l==0 && ($subdir =~ m/data/)) || ($l==1 && !($subdir =~ m/data/))){
		if($l==0){
		    $max=$maxdata;
		}
		else{
		    $max=$maxmc;
		}
		printf("Accessing Directory  $subdir \n");
		opendir(SUBDIR,"$InputDir/$subdir/");
		@files = grep { /root/ } readdir(SUBDIR);

		$nfiles = @files;
		$idx=0;
		foreach $file (@files){
		    $idx++;
		    printf("$file Set= $B  Index=  $A   Max.= $max N Files = $nfiles Current File = $idx \n");
		    if($A > $max){
			$A=1;
			$B++;
			# Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/  \" >>  $OutputDir/workdir$set/Input.txt ")) ;
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B \" >> $OutputDir/workdir$set/Submit")) ;
			system(sprintf("echo \"condor_submit  Condor_Set_$B  \" >> $OutputDir/workdir$set/Submit")) ;

			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;

            # Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe  | tee Set_$B.output \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));

			# Setup Input.txt
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Input.txt ")); 
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));

			# Setup Condor scripts
			system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"executable   = Set_$B.sh      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"output       = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"error        = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"log          = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));		
			system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
		    }
		    system(sprintf("echo \"File: $InputDir/$subdir/$file  \" >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    $A++;
		}
	    }
	}
    }

    # Finish Submit script
    system(sprintf("echo \"cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"chmod +x  Auto_Combine_submit  \" >> $OutputDir/workdir$set/Submit"));
    #system(sprintf("echo \"nohup  ./Auto_Combine_submit & \" >> $OutputDir/workdir$set/Submit"));
    GenerateAutoCombineScript("Auto_Combine_submit");
    system(sprintf("mv Auto_Combine_submit $OutputDir/workdir$set/"));

    # print Instructions  
    printf("\n\nInstructions");
    printf("\nPlease make sure you have run:");
    printf("\ngit config --global credential.helper 'cache --timeout=3600'");
    printf("\ngit config --global credential.helper cache");
    printf("\nNow you can run the analysis using dcache.");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
    printf("\nTo compile the code in the workdir: source compile ");
    printf("\nTo submit jobs to the batch queue: source Submit ");
    printf("\nTo combine jobs submitted to the batch queue: source Combine \n");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile; cd $OutputDir/workdir$set/Set_1; source Set_1 | tee log; cd ..\n");

} 



if( $ARGV[0] eq "--DCache" ){

    $TempDataSetFile=$ARGV[2];
    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Code Repository is:       $CodeDir \n");
    printf("List of dcache dir:       $TempDataSetFile \n");

    # Open ListofFile.txt
    @DataSets;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile!");
    while ($item = <DAT>) {
	chomp($item);
	print("File: $item \n");
	push(@DataSets,$item);
    }
    close(DAT);

    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    system(sprintf("cd $OutputDir"));
    system(sprintf("rm -rf $OutputDir/workdir$set \n"));
    system(sprintf("mkdir $OutputDir/workdir$set "));
    printf("Cleaning complete \n");
    
    # creat directory stucture
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/i386_linux "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));
    
    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source config \\\$\@ \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Setup Condor Combine scripts
    system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 			
    system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"rm Set*/*.o; rm Set*/*.e; rm Set*/*.log; \" >> $OutputDir/workdir$set/Submit")) ;
    $B=0;
    for($l=0;$l<2; $l++){
	system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Condor_Combine"));
	$max=1;
	foreach $DS (@DataSets){
	    if(($l==0 && ($DS =~ m/data/)) || ($l==1 && !($DS =~ m/data/))){
		if($l==0){
		    $max=$maxdata;
		}
		else{
		    $max=$maxmc;
		}
		printf("\n\nStarting Loop $l \n");
		$A=$maxdata+$maxmc+10;

		# find the root files for the current DataSet (DS)
		printf("Accessing Directory  $DS \n");
		system(sprintf("touch junk"));
		system(sprintf("uberftp grid-ftp.physik.rwth-aachen.de \"cd $DS; ls */  \" | grep root >& junk "));
		system(sprintf("cat junk | awk '{print \$8}' >& junk1")); 

		# Get list of files in dcache dir
		@files=();
		open(DAT, "junk1");
		while ($item = <DAT>) {
		    chomp($item);
		    push(@files,$item);
		}
		close(DAT);
		system(sprintf("rm junk"));
		system(sprintf("rm junk1"));
		$nfiles = @files;
		$idx=0;
		foreach $file (@files){
		    $idx++;
		    printf("$DS/$file Set= $B  Index=  $A   Max.= $max N Files = $nfiles Current File = $idx \n");
		    if($A > $max ){
			$A=1;
			$B++;

			# Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/ \" >>  $OutputDir/workdir$set/Input.txt ")) ;
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B \" >> $OutputDir/workdir$set/Submit")) ;
			system(sprintf("echo \"condor_submit  Condor_Set_$B  \" >> $OutputDir/workdir$set/Submit")) ;


			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;
			
			# Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-get.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ; 
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"mkdir  /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *   /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe 2>&1 | tee Set_$B.output \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *  $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"rm -r  /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));			
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));

            # Setup Set_$B_get.sh and Set_$B_clean.sh
			system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"mkdir /user/scratch/$UserID \" >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"cd /user/scratch/$UserID \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));

			system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
			system(sprintf("echo \"cd /user/scratch/$UserID \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh")); 

			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			# Setup Input.txt
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Input.txt ")); 
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));
			
			# Setup Condor scripts
			system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"executable   = Set_$B.sh      \"  >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"output       = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"error        = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"log          = Set_$B-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B")); 
			system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
			system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
		    }
		    ($a,$b,$c)=split('/',$file);
		    $myfile=$file;
		    if($a =~ m/root/){
			$myfile=$a;
		    }
		    if($b =~ m/root/){
                        $myfile=$b;
			system(sprintf("echo \"notification = Error        \" >> $OutputDir/workdir$set/Set_$B/Condor_Set_$B"));
                    }
		    if($c =~ m/root/){
                        $myfile=$c;
                    }
		    system(sprintf("echo \"dccp  dcap://grid-dcap-extern.physik.rwth-aachen.de/$DS/$file . \"  >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
		    system(sprintf("echo \"File:  /user/scratch/$UserID/$myfile \"     >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    system(sprintf("echo \"rm -rf /user/scratch/$UserID/$myfile  \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
		    $A++;
		}
	    }
	}
    }

    # Finish Submit script
    system(sprintf("echo \"cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"chmod +x  Auto_Combine_submit  \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"nohup  ./Auto_Combine_submit & \" >> $OutputDir/workdir$set/Submit"));
    GenerateAutoCombineScript("Auto_Combine_submit");
    system(sprintf("mv Auto_Combine_submit $OutputDir/workdir$set/"));

    # print Instructions
    printf("\n\nInstructions");
    printf("\nPlease make sure you have run:");
    printf("\nvoms-proxy-init");
    printf("\ngrid-proxy-init");
    printf("\ngit config --global credential.helper 'cache --timeout=3600'");
    printf("\ngit config --global credential.helper cache");
    printf("\nNow you can run the analysis using dcache.");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
    printf("\nTo compile the code in the workdir: source compile ");
    printf("\nTo submit jobs to the batch queue: source Submit ");
    printf("\nTo combine jobs submitted to the batch queue: source Combine \n");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile; cd $OutputDir/workdir$set/Set_1; source Set_1 | tee log; cd ..\n");
} 

if( $ARGV[0] eq "--GRID" ){

    if( $hasroot==0 ){
	printf("ERROR no root version defined ABORTING!!! Please read: ./todo.pl --help \n");
	exit(0);
    }

    $TempDataSetFile=$ARGV[2];
    # Print out input parameters
    printf("Active directory will be: $OutputDir/workdir$set \n");
    printf("Code Repository is:       $CodeDir \n");

    # Open ListofFile.txt
    @DataSets;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile!");
    while ($item = <DAT>) {
	chomp($item);
	print("File: $item \n");
	push(@DataSets,$item);
    }
    close(DAT);

    # Clean Directory in case workdir$set exists
    printf("Cleaning Directories \n");
    my $dir = getcwd;
    system(sprintf("cd $OutputDir"));
    system(sprintf("if [ -d  $OutputDir/workdir$set/ ]; then \n rm -rf $OutputDir/workdir$set \n fi "));
    system(sprintf("mkdir $OutputDir/workdir$set ")); 
    printf("Cleaning complete \n");
    
    #create directory stucture
    system(sprintf("cp $dir/subs $OutputDir/workdir$set "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code "));
    system(sprintf("mkdir  $OutputDir/workdir$set/Code/i386_linux "));
    system(sprintf("cp -r $CodeDir/* $OutputDir/workdir$set/Code/ "));
    system(sprintf("mkdir $OutputDir/workdir$set/EPS "));
    system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/InputData "));

    # Setup file to retrieve jobs
    system(sprintf("cp $dir/CheckandGet.sh $OutputDir/workdir$set/"));
    
    # generate compile script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd  $OutputDir/workdir$set/Code/\" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"source config \\\$\@ \"   >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"gmake all \" >> $OutputDir/workdir$set/compile "));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/compile")) ;
 
    # Generate Combine script 
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Combine")) ;
    system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Combine"));
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Combine")) ; 
    system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe \" >> $OutputDir/workdir$set/Combine")) ;

    # Generate Combine Input
    system(sprintf("cp $InputFile $OutputDir/workdir$set/Input.txt "));
    system(sprintf("echo \"Mode: RECONSTRUCT\" >> $OutputDir/workdir$set/Input.txt"));
    system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Input.txt"));

    # Setup Condor Combine scripts
    system(sprintf("echo \"universe     = vanilla      \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"rank         = memory       \"  >> $OutputDir/workdir$set/Condor_Combine"));
    system(sprintf("echo \"executable   = Combine      \"  >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"output       = Combine-Condor_\\\$(cluster)_\\\$(proccess).o  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"error        = Combine-Condor_\\\$(cluster)_\\\$(proccess).e  \" >> $OutputDir/workdir$set/Condor_Combine")); 
    system(sprintf("echo \"log          = Combine-Condor_\\\$(cluster)_\\\$(proccess).log  \" >> $OutputDir/workdir$set/Condor_Combine")); 			
    system(sprintf("echo \"queue = 1 \" >> $OutputDir/workdir$set/Condor_Combine"));

    # Start Submit script
    system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Submit")) ; 
    system(sprintf("echo \"cd $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit")) ;
    system(sprintf("echo \"if [ -f $OutputDir/workdir$set/Set*/*.log ]; then \n rm $OutputDir/workdir$set/Set*/*.o; rm $OutputDir/workdir$set/Set*/*.e; rm $OutputDir/workdir$set/Set*/*.log; rm $OutputDir/workdir$set/Set*/*.tar; \n fi  \" >> $OutputDir/workdir$set/Submit"));


    system(sprintf("echo \"echo 'Creating tarballs and installing on the GRID... '\" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"if [ -f workdir$set.tar  ]; then \n rm workdir$set.tar \n fi \n tar -cf workdir$set.tar root Code;\" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"srmrm  srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set.tar \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"srmcp  file:////$OutputDir/workdir$set/workdir$set.tar srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set.tar \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"echo 'Complete '\" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"rm   $OutputDir/workdir$set/jobs_submitted ; touch $OutputDir/workdir$set/jobs_submitted\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"rm   $OutputDir/workdir$set/jobs_complete ; touch $OutputDir/workdir$set/jobs_submitted\" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"echo 'Setting jobs to use mode GRID instead of LOCAL' \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"./subs 'RunType: LOCAL' 'RunType: GRID' Set*/Input.txt  \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"hasdir=\\\$(srmls -recursion_depth=1 srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/ | grep /pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/  | wc -l) \n if [ \\\$hasdir == 0 ]; then \" >> $OutputDir/workdir$set/Submit "));
    system(sprintf("echo \"srmmkdir srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set \n fi \" >> $OutputDir/workdir$set/Submit"));

    $B=0;
    for($l=0;$l<2; $l++){
	system(sprintf("echo \"notification = Error \" >> $OutputDir/workdir$set/Condor_Combine"));
	$max=1;
	foreach $DS (@DataSets){
	    if(($l==0 && ($DS =~ m/data/)) || ($l==1 && !($DS =~ m/data/))){
		if($l==0){
		    $max=$maxdata;
		}
		else{
		    $max=$maxmc;
		}
		printf("\n\nStarting Loop $l \n");
		$A=$maxdata+$maxmc+10;

		# find the root files for the current DataSet (DS)
		printf("Accessing Directory  $DS \n");
		system(sprintf("touch junk"));
		#system(sprintf("uberftp grid-ftp.physik.rwth-aachen.de \"cd $DS; ls */  \" | grep root >& junk "));
		#system(sprintf("cat junk | awk '{print \$8}' >& junk1")); 
		system(sprintf("srmls -recursion_depth=3  srm://$gridsite:8443$DS >& junk"));
		system(sprintf("cat junk | grep root | awk '{print \$2}' >& junk1"));
		
		# Get list of files in dcache dir
		@files=(); 
		open(DAT, "junk1");
		while ($item = <DAT>) {
		    chomp($item);
		    push(@files,$item);
		}
		close(DAT);

                @filessrm=();
                open(DAT, "junk");
                while ($item = <DAT>) {
                    chomp($item);
                    push(@filessrm,$item);
                }
                close(DAT);

		system(sprintf("rm junk"));
		system(sprintf("rm junk1"));
		$nfiles = @files;
		$idx=0;
		$nfilessrm = @filessrm;
		
		
		if($nfilessrm>999){
		    printf("Error more then 999 files in Directory! SRM does not handel this correctly -> ABORTING!");
		    exit(0);
		}

		foreach $file (@files){
		    $idx++;
		    printf("$DS/$file Set= $B  Index=  $A   Max.= $max N Files = $nfiles Current File = $idx \n");
		    if($A > $max ){
			$A=1;
			$B++;


			# Add Set information to Combining scripts and Input.txt
			system(sprintf("echo \"File: $OutputDir/workdir$set/Set_$B/ \" >>  $OutputDir/workdir$set/Input.txt ")) ;
			system(sprintf("echo \"cd $OutputDir/workdir$set/  \" >> $OutputDir/workdir$set/Submit")) ;
			system(sprintf("echo \"echo 'Submitting Set_$B/GRIDJob.jdl' \" >> $OutputDir/workdir$set/Submit "));
			system(sprintf("echo \"isSkim=\\\$(srmls -recursion_depth=2 srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/ | grep SKIMMED_NTUP_$B.root  | wc -l) \n if [ \\\$isSkim == 1 ]; then \" >> $OutputDir/workdir$set/Submit "));
			system(sprintf("echo \"srmrm  srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/SKIMMED_NTUP_$B.root \n fi \" >> $OutputDir/workdir$set/Submit "));
			system(sprintf("echo \"glite-ce-job-submit -a -r grid-ce.physik.rwth-aachen.de:8443/cream-pbs-short $OutputDir/workdir$set/Set_$B/GRIDJob.jdl  |  grep -A 5 \\\"glite-ce-job-submit Success\\\" | grep https | awk '{print \\\$1 \\\" $OutputDir/workdir$set/Set_$B \\\" }' >>  $OutputDir/workdir$set/jobs_submitted \" >> $OutputDir/workdir$set/Submit")) ;


			# Create and configure Set_$B dir
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B ")) ;
			system(sprintf("ln -s $OutputDir/workdir$set/Code/InputData $OutputDir/workdir$set/Set_$B/InputData "));
			system(sprintf("mkdir $OutputDir/workdir$set/Set_$B/EPS ")) ;
			
			# Setup Set_$B.sh
			system(sprintf("echo \"#! /bin/bash\" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;
			system(sprintf("echo \"echo 'Starting Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Code/; source config \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-get.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ; 
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"mkdir  /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *   /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"$OutputDir/workdir$set/Code/Analysis.exe | tee Set_$B.output \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cp -r *  $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"source $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"rm -r  /user/scratch/$UserID/workdir$set-Set_$B  \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));			
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh"));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B/ \" >> $OutputDir/workdir$set/Set_$B/Set_$B.sh")) ;

            # Setup Set_$B-GRID.sh 
            system(sprintf("echo \"echo 'Starting Job'; ls; echo \\\$PWD \" > $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo 'srmcp  srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set.tar file:////\$PWD/workdir$set.tar ; tar -xf workdir$set.tar; ' >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh  "));
			system(sprintf("echo \"cp -r \\\$PWD/Code/InputData/ \\\$PWD/InputData \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("chmod +x $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo \"source \\\$PWD/Set_$B-getGRID.sh \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
            system(sprintf("echo \"cd Code/; chmod +x  \\\$PWD/config  \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
            system(sprintf("echo \"echo \\\$PWD | awk '{split(\\\$1,a,\\\"Code\\\"); print \\\"#! /bin/bash  \\n source \\\" \\\$1 \\\"/config --useRoot \\\" a[1] \\\"root/ \\\"}'  > \\\$PWD/junk; source \\\$PWD/junk; \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo \" cd ../; echo \\\$PWD; ls  \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo \"chmod +x \\\$PWD/subs; \\\$PWD/subs '/user/scratch/$UserID' \\\$PWD Input.txt \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo \"export LD_LIBRARY_PATH=\\\$LD_LIBRARY_PATH:\\\$ROOTSYS/lib/\" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			system(sprintf("echo \"echo 'System Configured.'; printenv; ls ; echo \\\$PWD ; \\\$PWD/Code/Analysis.exe \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
			#
			system(sprintf("echo \"if [ -f SKIMMED_NTUP.root ]; then \n \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh "));
			system(sprintf("echo \"srmcp  file:////\\\$PWD/SKIMMED_NTUP.root srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/SKIMMED_NTUP_$B.root \n fi \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh "));
			#
			system(sprintf("echo \"echo 'Completed Job' \" >> $OutputDir/workdir$set/Set_$B/Set_$B-GRID.sh"));
            #
			system(sprintf("echo \"#! /bin/bash \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
			system(sprintf("echo \"isSkim=\\\$(srmls -recursion_depth=2 srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/ | grep SKIMMED_NTUP_$B.root  | wc -l) \n if [ \\\$isSkim == 1 ]; then \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));
			system(sprintf("echo \"cd $OutputDir/workdir$set/Set_$B; srmcp  srm://$gridsite:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$UserIDCern/workdir$set/SKIMMED_NTUP_$B.root  file:////\\\$PWD/SKIMMED_NTUP.root; cd $OutputDir/workdir$set/ \n fi \" >> $OutputDir/workdir$set/Set_$B/GRIDRetrieve.sh "));

			
            # Setup Set_$B_get.sh Set_$B-getGRID and Set_$B_clean.sh
			system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"mkdir /user/scratch/$UserID \" >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
			system(sprintf("echo \"cd /user/scratch/$UserID \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh")); 
            system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
            system(sprintf("echo \"cd /user/scratch/$UserID \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));
            system(sprintf("echo \"#! /bin/bash\"         >> $OutputDir/workdir$set/Set_$B/Set_$B-getGRID.sh"));

			# Setup Input.txt
			system(sprintf("cp   $InputFile $OutputDir/workdir$set/Set_$B/Input.txt ")); 
			system(sprintf("echo \"Mode: ANALYSIS\" >> $OutputDir/workdir$set/Set_$B/Input.txt")); 
			system(sprintf("echo \"RunType: LOCAL\" >> $OutputDir/workdir$set/Set_$B/Input.txt"));
			
			# Setup .jdl file
			system(sprintf("echo '['  >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
            system(sprintf("echo ' Executable    = \"Set_$B-GRID.sh\"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ' Arguments     = \"42\"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ' StdOutput     = \"out\"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ' StdError      = \"err\"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ' InputSandbox  = {\"Set_$B/Set_$B-GRID.sh\",\"Set_$B/Set_$B-getGRID.sh\",\"subs\",\"Set_$B/Input.txt\"}; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ' OutputSandbox = {\"out\",\"err\",\"*_ANALYSIS_*.root\"}; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ' Requirements  = \"\"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ' VirtualOrganisation = \"dcms\"; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ' Rank          = -other.GlueCEStateEstimatedResponseTime; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ' Retry         = 1; ' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			system(sprintf("echo ']' >> $OutputDir/workdir$set/Set_$B/GRIDJob.jdl"));
			
		    }
		    ($x,$y)=split($DS,$file);
		    ($a,$b,$c,$d,$e)=split('/',$y);
		    $file=$y;
		    if($a =~ m/root/){
			$myfile=$a;
		    }
		    if($b =~ m/root/){
                        $myfile=$b;
                    }
		    if($c =~ m/root/){
                        $myfile=$c;
                    }
                    if($d=~ m/root/){
                        $myfile=$d;
                    }

                    if($e =~ m/root/){
                        $myfile=$e;
                    }

		    system(sprintf("echo \"dccp  dcap://grid-dcap-extern.physik.rwth-aachen.de/$DS$file . \"  >> $OutputDir/workdir$set/Set_$B/Set_$B-get.sh"));
		    system(sprintf("echo 'srmcp  srm://$gridsite:8443/$DS$file file:////\$PWD/$myfile '  >> $OutputDir/workdir$set/Set_$B/Set_$B-getGRID.sh"));
		    system(sprintf("echo \"File:  /user/scratch/$UserID/$myfile \"     >> $OutputDir/workdir$set/Set_$B/Input.txt")) ;
		    $A++;
		}
	    }
	}
    }


    printf("Setting up root....\n");
    if($buildRoot==1){
        printf("Building custom root $buildrootversion... Enjoy your coffee!!! ");
	system(sprintf("wget ftp://root.cern.ch/root/root_v$buildrootversion.source.tar.gz"));
	system(sprintf("gzip -dc root_v$buildrootversion.source.tar.gz | tar -xf -"));
	system(sprintf("mkdir $OutputDir/workdir$set/root"));
	system(sprintf("cd root_v$buildrootversion; ./configure --enable-python --enable-roofit --enable-minuit2 --disable-xrootd --disable-sqlite --disable-python  --disable-mysql --prefix=$OutputDir/workdir$set/root; make & make install "));
    }
    else{
        printf("Copying local root $MYROOTSYS ");
        system(sprintf("mkdir $OutputDir/workdir$set/root/"));
        system(sprintf("cp -r $MYROOTSYS/* $OutputDir/workdir$set/root/"));
    }
    
    # Finish Submit script
    system(sprintf("echo \"cd  $OutputDir/workdir$set/ \" >> $OutputDir/workdir$set/Submit"));
    system(sprintf("echo \"rm -rf /user/scratch/$UserID/$myfile  \"    >> $OutputDir/workdir$set/Set_$B/Set_$B-clean.sh"));

    # print Instructions
    printf("\n\nInstructions");
    printf("\nPlease make sure you have run:");
    printf("\nvoms-proxy-init");
    printf("\ngrid-proxy-init");
    printf("\ngit config --global credential.helper 'cache --timeout=3600'");
    printf("\ngit config --global credential.helper cache");
    printf("\nNow you can run the analysis using dcache.");
    printf("\nTo go to the Test workdir: cd  $OutputDir/workdir$set ");
    printf("\nTo compile the code in the workdir: source compile --useRoot $OutputDir/workdir$set/root/ $UserDir ");
    printf("\nTo submit jobs to the GRID: source Submit ");
    printf("\nTo check the status of the GRID jobs and download finished jobs: source CheckandGet.sh");
    printf("\nTo additionally print the details of all the jobs: source CheckandGet.sh  --detailed");
    printf("\nTo test a single job: cd  $OutputDir/workdir$set; source compile  --useRoot $OutputDir/workdir$set/root/ $UserDir; cd $OutputDir/workdir$set/Set_1; source Set_1.sh; cd ..\n");
    
} 

sub GenerateAutoCombineScript{
    $file=$_[0];
    system(sprintf("echo \"#! /bin/bash\" >> $file"));
    system(sprintf("echo \"for i in {1..2881} \" >>  $file"));
    system(sprintf("echo \"do \" >>  $file"));
    system(sprintf("echo \"  nDir=\\\$(ls | grep Set_ | wc -l) \" >>  $file"));
    system(sprintf("echo \"  nFin=\\\$(grep \\\"Completed Job\\\" Set_*/Set_*.o | wc -l) \" >>  $file"));
    system(sprintf("echo \"  if [ \\\$nFin != \\\$nDir ]; then\" >>  $file"));
    system(sprintf("echo \"    sleep 300\" >>  $file"));
    system(sprintf("echo \"  else\" >>  $file"));
    system(sprintf("echo \"    sleep 180\" >>  $file"));
    system(sprintf("echo \"    nSuc=\\\$(grep \\\"Program is Finished\\\" Set_*/Set_*.output  | wc -l )  \" >>  $file"));
    system(sprintf("echo \"    if [ \\\$nFin == \\\$nSuc ]; then\" >>  $file"));
    system(sprintf("echo \"      echo \\\"All Jobs Succeeded: Submitting Condor_Combine \\\"\" >>  $file"));
    system(sprintf("echo \"      condor_submit Condor_Combine\" >>  $file"));
    system(sprintf("echo \"    else \" >>  $file"));
    system(sprintf("echo \"      echo \\\"Failed Jobs: Aborting combine submission \\\"\" >>  $file"));
    system(sprintf("echo \"    fi\" >>  $file"));
    system(sprintf("echo \"    break\" >>  $file"));
    system(sprintf("echo \"  fi\" >>  $file"));
    system(sprintf("echo \"  if [ \\\$i == 2880 ]; then\" >>  $file"));
    system(sprintf("echo \"    echo \\\"Out of time (48hr): Submit Condor_Combine manually when jobs complete \\\" \" >>  $file"));
    system(sprintf("echo \"    break\" >>  $file")); 
    system(sprintf("echo \"  fi\" >>  $file")); 
    system(sprintf("echo \"done\" >>  $file")); 
}

