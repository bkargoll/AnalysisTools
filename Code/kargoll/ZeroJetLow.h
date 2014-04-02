/*
 * ZeroJetLow.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ZEROJETLOW_H_
#define ZEROJETLOW_H_

#include "HToTaumuTauh.h"

class ZeroJetLow: public HToTaumuTauh {
public:
	ZeroJetLow(TString Name_, TString id_);
	virtual ~ZeroJetLow();
};

#endif /* ZEROJETLOW_H_ */
