/*
 * VBFLoose.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef VBFLOOSE_H_
#define VBFLOOSE_H_

#include "Category.h"
#include "HToTaumuTauh.h"

class VBFLoose: public Category {
public:
	VBFLoose(TString Name_, TString id_);
	virtual ~VBFLoose();

	// enumerator for the cuts in this category
	enum cuts_VBFLoose {
		NJet	= HToTaumuTauh::CatCut1,
		DeltaEta,
		NJetRapGap,
		JetInvM,
		NotVbfTight,
		NCuts
	};

protected:
	void categoryConfiguration();
	bool categorySelection();
	void categoryPlotting();
};

#endif /* VBFLOOSE_H_ */
