/*
 * HToTaumuTauhSkim.h
 *
 *  Created on: May 14, 2014
 *      Author: kargoll
 */

#ifndef HTOTAUMUTAUHSKIM_H_
#define HTOTAUMUTAUHSKIM_H_

#include "HToTaumuTauh.h"

class HToTaumuTauhSkim: public HToTaumuTauh {
public:
	HToTaumuTauhSkim(TString Name_, TString id_);
	virtual ~HToTaumuTauhSkim();

	virtual void Configure();

protected:
	virtual void doEvent();
	virtual void Store_ExtraDist();

	// cuts to be relaxed for skimming
	double cSkim_Mu_relIso;
	double cSkim_Tau_rawIso;

};

#endif /* HTOTAUMUTAUHSKIM_H_ */
