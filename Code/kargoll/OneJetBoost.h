/*
 * VBFTight.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ONEJETBOOST_H_
#define ONEJETBOOST_H_

#include "HToTaumuTauh.h"

class OneJetBoost: public HToTaumuTauh {
public:
	OneJetBoost(TString Name_, TString id_);
	virtual ~OneJetBoost();
};

#endif /* ONEJETBOOST_H_ */
