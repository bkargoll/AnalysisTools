/*
 * RelaxedVBFTightStandalone.h
 *
 *  Created on: May 8, 2015
 *      Author: kargoll
 */

#ifndef RelaxedVBFTightStandalone_H_
#define RelaxedVBFTightStandalone_H_

#include "VBFTightStandalone.h"

class RelaxedVBFTightStandalone: public VBFTightStandalone {
public:
	RelaxedVBFTightStandalone();
	RelaxedVBFTightStandalone(unsigned nJets, float deltaEta, unsigned nJetRapGap, float jetInvM, float higgsPt);
	virtual ~RelaxedVBFTightStandalone();

private:
	virtual void setCutValues();

};

#endif /* RelaxedVBFTightStandalone_H_ */
