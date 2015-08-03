/*
 * VBFLooseStandalone.cxx
 *
 *  Created on: May 8, 2015
 *      Author: kargoll
 */

#include "VBFLooseStandalone.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

// incomplete constructor, allows to retrieve cut values
VBFLooseStandalone::VBFLooseStandalone() {
	nCuts_ = NCuts;
	cutValues_.resize(nCuts_);
}

VBFLooseStandalone::VBFLooseStandalone(unsigned nJets, double deltaEta, int nJetRapGap, double jetInvM, bool notVBFTight)
{
	nCuts_ = NCuts;
	cutValues_.resize(nCuts_);
	eventValues_.resize(nCuts_);
	passCut_.resize(nCuts_);

	eventValues_.at(NJet)		= nJets;
	eventValues_.at(DeltaEta)	= deltaEta;
	eventValues_.at(NJetRapGap) = nJetRapGap;
	eventValues_.at(JetInvM)	= jetInvM;
	eventValues_.at(NotVbfTight)= notVBFTight;
}

VBFLooseStandalone::~VBFLooseStandalone() {
}

void VBFLooseStandalone::setCutValues() {
	cutValues_.at(NJet)		= 2;
	cutValues_.at(DeltaEta)	= 3.5;
	cutValues_.at(NJetRapGap)	= 0;
	cutValues_.at(JetInvM)		= 500.0;
	cutValues_.at(NotVbfTight)	= true;
}

void VBFLooseStandalone::execute(){
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

	passCut_.at(NotVbfTight) = ( eventValues_.at(NotVbfTight) == cutValues_.at(NotVbfTight) );
}
