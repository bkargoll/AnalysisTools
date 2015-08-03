/*
 * RelaxedVBFLooseStandalone.h
 *
 *  Created on: May 8, 2015
 *      Author: kargoll
 */

#ifndef RelaxedVBFLooseStandalone_H_
#define RelaxedVBFLooseStandalone_H_

#include "VBFLooseStandalone.h"

class RelaxedVBFLooseStandalone: public VBFLooseStandalone {
public:
	RelaxedVBFLooseStandalone();
	RelaxedVBFLooseStandalone(unsigned nJets, float deltaEta, unsigned nJetRapGap, float jetInvM);
	virtual ~RelaxedVBFLooseStandalone();

private:
	virtual void setCutValues();
};

#endif /* RelaxedVBFLooseStandalone_H_ */
