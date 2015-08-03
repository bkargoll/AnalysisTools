/*
 * RelaxedVBFTightStandalone.cxx
 *
 *  Created on: May 8, 2015
 *      Author: kargoll
 */

#include "RelaxedVBFTightStandalone.h"
#include "SimpleFits/FitSoftware/interface/Logger.h"

// incomplete constructor, allows to retrieve cut values
RelaxedVBFTightStandalone::RelaxedVBFTightStandalone():
	VBFTightStandalone() {}

RelaxedVBFTightStandalone::RelaxedVBFTightStandalone(unsigned nJets, float deltaEta, unsigned nJetRapGap, float jetInvM, float higgsPt):
	VBFTightStandalone(nJets, deltaEta, nJetRapGap, jetInvM, higgsPt)
		{ }

RelaxedVBFTightStandalone::~RelaxedVBFTightStandalone() {
}

void RelaxedVBFTightStandalone::setCutValues() {
	cutValues_.at(NJet)			= 2;
	cutValues_.at(DeltaEta)		= 2.0;
	cutValues_.at(NJetRapGap)	= 0;
	cutValues_.at(JetInvM)		= 200.0;
	cutValues_.at(HiggsPt)		= 100.0;
}
