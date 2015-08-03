/*
 * VBFTightStandalone.h
 *
 *  Created on: May 08, 2015
 *      Author: kargoll
 */

#ifndef VBFTIGHTStandalone_H_
#define VBFTIGHTStandalone_H_

#include "CategoryStandalone.h"

class VBFTightStandalone : public CategoryStandalone {
public:
	VBFTightStandalone();
	VBFTightStandalone(unsigned nJets, double deltaEta, int nJetRapGap, double jetInvM, double higgsPt);
	virtual ~VBFTightStandalone();

	// enumerator for the cuts in this category
	enum cuts_VBFTightStandalone {
		NJet = 0,
		DeltaEta,
		NJetRapGap,
		JetInvM,
		HiggsPt,
		NCuts
	};

private:
	virtual void execute();
	virtual void setCutValues();
};

#endif /* VBFTIGHTStandalone_H_ */
