/*
 * VBFLooseStandalone.h
 *
 *  Created on: May 8, 2015
 *      Author: kargoll
 */

#ifndef VBFLOOSESTANDALONE_H_
#define VBFLOOSESTANDALONE_H_

#include "CategoryStandalone.h"

class VBFLooseStandalone: public CategoryStandalone {
public:
	VBFLooseStandalone();
	VBFLooseStandalone(unsigned nJets, double deltaEta, int nJetRapGap, double jetInvM, bool notVBFTight);
	virtual ~VBFLooseStandalone();

	// enumerator for the cuts in this category
	enum cuts_VBFLoose {
		NJet	= 0,
		DeltaEta,
		NJetRapGap,
		JetInvM,
		NotVbfTight,
		NCuts
	};

private:
	virtual void execute();
	virtual void setCutValues();
};

#endif /* VBFLOOSESTANDALONE_H_ */
