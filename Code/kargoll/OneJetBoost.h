/*
 * OneJetBoost.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ONEJETBOOST_H_
#define ONEJETBOOST_H_

#include "Category.h"

class OneJetBoost: public Category {
public:
	OneJetBoost(TString Name_, TString id_);
	virtual ~OneJetBoost();

	enum cuts_OneJetBoost {
		NJet = HToTaumuTauh::CatCut1,
		NotVbf,
		TauPt,
		HiggsPt,
		NCuts
	};


protected:
	void categoryConfiguration();
	bool categorySelection();
	void categoryPlotting();
};

#endif /* ONEJETBOOST_H_ */
