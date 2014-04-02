/*
 * ZeroJetHigh.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ZEROJETHIGH_H_
#define ZEROJETHIGH_H_

#include "HToTaumuTauh.h"

class ZeroJetHigh: public HToTaumuTauh {
public:
	ZeroJetHigh(TString Name_, TString id_);
	virtual ~ZeroJetHigh();
};

#endif /* ZEROJETHIGH_H_ */
