/*
 * HToTaumuTauhBackgrounds.h
 *
 *  Created on: Jun 20, 2014
 *      Author: kargoll
 */

#ifndef HTOTAUMUTAUHBACKGROUNDS_H_
#define HTOTAUMUTAUHBACKGROUNDS_H_

#include "HToTaumuTauh.h"

class HToTaumuTauhBackgrounds: public HToTaumuTauh {
public:
	HToTaumuTauhBackgrounds(TString Name_, TString id_);
	virtual ~HToTaumuTauhBackgrounds();

protected:
	virtual void  Finish();
};

#endif /* HTOTAUMUTAUHBACKGROUNDS_H_ */
