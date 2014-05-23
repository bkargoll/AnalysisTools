include make.inc
#################################################
### Defining Compiling variables              ###
#################################################

CXX           = g++ -Wno-deprecated
#-m32  -Wno-deprecated
LD            = g++ -Wno-deprecated
#-m32  -Wno-deprecated 
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

SHAREDLIBFLAGS =-L./CommonUtils/lib -lTauNtuple 
SHAREDCXXFLAGS =-I./CommonUtils 
SUBDIRS = CommonUtils/ 

CXXFLAGS      = -I$(ROOTSYS)/include -c -I/include $(SHAREDCXXFLAGS) -c -I./ -c
LDFLAGS       =  
ROOTLIBS      =  -L$(ROOTSYS)/lib -L/usr/lib/ -L/lib/i686/  -lCore -lCint -lHist -lGraf  -lGraf3d -lGpad -lTree -lRint -lReflexDict -lReflex -lPostscript -lMatrix -lPhysics -lGui -LObj -lMinuit2 -lThread -rdynamic -Wl,--rpath $(ROOTSYS)/lib

# -ldl -lpthread

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS) $(SHAREDLIBFLAGS)

TARGETS       =	Ntuple_Controller \
		NtupleReader \
		Parameters \
		Selection_Base \
		Selection_Factory \
		Selection \
		Tables \
		Plots \
		DoubleEventRemoval \
		PDG_Var \
		Tools \
		SkimConfig \
		HistoConfig \
		Example \
		TauSolver \
		TauSpinExample \
		TauSpinerInterface

# now include conditional options set by flags in make.inc
ifdef USE_TauSpinner
	SUBDIRS += TauSpiner/
	SHAREDLIBFLAGS +=	-L./TauSpiner/lhapdf/lib -lLHAPDF \
				-L./TauSpiner/tauola++/1.1.3/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaTauSpinner \
				-L./TauSpiner/hepmc/install/lib -lHepMC
        SHAREDCXXFLAGS +=	-I./TauSpiner/hepmc/install/include \
				-I./TauSpiner/tauola++/1.1.3/include/Tauola/ \
				-I./TauSpiner/tauola++/1.1.3/include/TauSpinner/ \
				-I./TauSpiner/tauola++/1.1.3/include/ \
				-I./TauSpiner/lhapdf/include/ \
				-I./TauSpiner/pythia8/176/include/
	DEFS += -DUSE_TauSpinner=1
endif

ifdef USE_cherepanov
	TARGETS += Ztotautau_hadmu_ControlSample Tau_momentum_calculation Validation
	DEFS += -DUSE_cherepanov=1
endif

ifdef USE_inugent
	TARGETS += TauLifeTime TriggerStudy TriggerStudyMC ZDouble3prong Ztotautau_ControlSample Ztomumu_ControlSample ChargedHiggs_dilepontic ChargedHiggs_tauplusjet 
	DEFS += -DUSE_inugent=1
endif

ifdef USE_kargoll
	TARGETS += HToTaumuTauh HToTaumuTauhSkim MuTauSync OneJetBoost OneJetHigh OneJetLow VBFLoose VBFTight ZeroJetHigh ZeroJetLow
	DEFS += -DUSE_kargoll=1
endif 

ifdef USE_nehrkorn
	TARGETS += ZtoEMu ZtoEMu_Skim ZtoEMu_Fakerate ZtoEMu_Wjets ZtoEMu_QCD
	DEFS += -DUSE_nehrkorn=1
endif

ifdef USE_pistone
        DEFS += -DUSE_pistone=1
endif



HDRS          = $(TARGETS:=.h) 		
SRCS          = $(TARGETS:=.cxx)  Analysis.cxx
OBJS          = $(SRCS:.cxx=.o)

PROGRAM       = Analysis.exe

$(PROGRAM):     $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@echo $(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@echo "done"

VPATH = utilities:i386_linux
vpath %.cxx inugent
vpath %.h inugent
vpath %.cxx cherepanov
vpath %.h cherepanov
vpath %.cxx kargoll
vpath %.h kargoll
vpath %.cxx nehrkorn
vpath %.h nehrkorn

#	@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o ../$(PROGRAM)

##########################################################
###                                                    ###
##########################################################

$(OBJS): %.o : %.cxx
	$(CXX) $(CXXFLAGS) $(DEFS) $< -o i386_linux/$@ 

.PHONY: clean cleanall all install sharedlib 

install: Analysis.exe


clean:
	@rm i386_linux/*.o
	@rm Analysis.exe

cleanall:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone clean && cd ../../; \
	done
	@rm i386_linux/*.o
	@rm Analysis.exe

all: sharedlib install

sharedlib:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone && cd ../ && pwd; \
	done
	$(install)

