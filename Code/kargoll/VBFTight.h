/*
 * VBFTight.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef VBFTIGHT_H_
#define VBFTIGHT_H_

#include "HToTaumuTauh.h"

class VBFTight: public HToTaumuTauh {
public:
	VBFTight(TString Name_, TString id_);
	virtual ~VBFTight();
};

#endif /* VBFTIGHT_H_ */
