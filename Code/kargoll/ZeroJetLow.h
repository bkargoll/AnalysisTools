/*
 * ZeroJetLow.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ZEROJETLOW_H_
#define ZEROJETLOW_H_

#include "Category.h"

class ZeroJetLow: public Category {
public:
	ZeroJetLow(TString Name_, TString id_);
	virtual ~ZeroJetLow();

	// enumerator for the cuts in this category
	enum cuts_ZeroJetLow {
		NJet = HToTaumuTauh::CatCut1,
		TauPt,
		NCuts
	};

protected:
	void categoryConfiguration();
	bool categorySelection();
};

#endif /* ZEROJETLOW_H_ */
