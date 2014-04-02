/*
 * OneJetLow.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef ONEJETLOW_H_
#define ONEJETLOW_H_

#include "HToTaumuTauh.h"

class OneJetLow: public HToTaumuTauh {
public:
	OneJetLow(TString Name_, TString id_);
	virtual ~OneJetLow();
};

#endif /* ONEJETLOW_H_ */
