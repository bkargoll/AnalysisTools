/*
 * OneJetLow.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ONEJETLOW_H_
#define ONEJETLOW_H_

#include "Category.h"

class OneJetLow: public Category {
public:
	OneJetLow(TString Name_, TString id_);
	virtual ~OneJetLow();

	// enumerator for the cuts in this category
	enum cuts_OneJetLow {
		NJet = CatCut1,
		NotVbf,
		TauPt,
		NCuts
	};

protected:
	void categoryConfiguration();
	bool categorySelection();

};

#endif /* ONEJETLOW_H_ */
