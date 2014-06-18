/*
 * Inclusive.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef INCLUSIVE_H_
#define INCLUSIVE_H_

#include "HToTaumuTauh.h"

class Inclusive: public HToTaumuTauh {
public:
	Inclusive(TString Name_, TString id_);
	virtual ~Inclusive();
};

#endif /* INCLUSIVE_H_ */
