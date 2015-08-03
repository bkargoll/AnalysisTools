/*
 * RelaxedVBFLooseStandalone.cxx
 *
 *  Created on: May 8, 2015
 *      Author: kargoll
 */

#include "RelaxedVBFLooseStandalone.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

// incomplete constructor, allows to retrieve cut values
RelaxedVBFLooseStandalone::RelaxedVBFLooseStandalone():
	VBFLooseStandalone() {}

RelaxedVBFLooseStandalone::RelaxedVBFLooseStandalone(unsigned nJets, float deltaEta, unsigned nJetRapGap, float jetInvM):
	VBFLooseStandalone(nJets, deltaEta, nJetRapGap, jetInvM, true) { }

RelaxedVBFLooseStandalone::~RelaxedVBFLooseStandalone() {
}

void RelaxedVBFLooseStandalone::setCutValues() {
	cutValues_.at(NJet)			= 2;
	cutValues_.at(DeltaEta)		= 2.0;
	cutValues_.at(NJetRapGap)	= 0;
	cutValues_.at(JetInvM)		= 200.0;
	cutValues_.at(NotVbfTight)	= true;
}
