/*
 * VBFTightStandalone.cxx
 *
 *  Created on: May 08, 2015
 *      Author: kargoll
 */

#include "VBFTightStandalone.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

// incomplete constructor, allows to retrieve cut values
VBFTightStandalone::VBFTightStandalone() {
	nCuts_ = NCuts;
	cutValues_.resize(nCuts_);
}

VBFTightStandalone::VBFTightStandalone(unsigned nJets, double deltaEta, int nJetRapGap, double jetInvM, double higgsPt)
{
	nCuts_ = NCuts;
	cutValues_.resize(nCuts_);
	eventValues_.resize(nCuts_);
	passCut_.resize(nCuts_);

	eventValues_.at(NJet)		= nJets;
	eventValues_.at(DeltaEta)	= deltaEta;
	eventValues_.at(NJetRapGap)= nJetRapGap;
	eventValues_.at(JetInvM)	= jetInvM;
	eventValues_.at(HiggsPt)	= higgsPt;
}

VBFTightStandalone::~VBFTightStandalone() {
}

void VBFTightStandalone::setCutValues() {
	cutValues_.at(NJet)		= 2;
	cutValues_.at(DeltaEta)	= 4.0;
	cutValues_.at(NJetRapGap)	= 0;
	cutValues_.at(JetInvM)		= 700.0;
	cutValues_.at(HiggsPt)		= 100.0;
}

void VBFTightStandalone::execute(){
	passCut_.at(NJet) = (eventValues_.at(NJet) >= cutValues_.at(NJet));

	if(passCut_.at(NJet)){
		passCut_.at(DeltaEta) = (fabs(eventValues_.at(DeltaEta)) > cutValues_.at(DeltaEta));

		passCut_.at(NJetRapGap) = (eventValues_.at(NJetRapGap) <= cutValues_.at(NJetRapGap));

		passCut_.at(JetInvM) = (eventValues_.at(JetInvM) > cutValues_.at(JetInvM));
	}
	else{
		passCut_.at(DeltaEta) = true;
		passCut_.at(NJetRapGap) = true;
		passCut_.at(JetInvM) = true;
	}

	passCut_.at(HiggsPt) = (eventValues_.at(HiggsPt) > cutValues_.at(HiggsPt));
}
