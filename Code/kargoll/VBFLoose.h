/*
 * VBFLoose.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef VBFLOOSE_H_
#define VBFLOOSE_H_

#include "HToTaumuTauh.h"

class VBFLoose: public HToTaumuTauh {
public:
	VBFLoose(TString Name_, TString id_);
	virtual ~VBFLoose();
};

#endif /* VBFLOOSE_H_ */
