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

	virtual void Configure();

protected:
	virtual void Setup();
	virtual void doEvent();
	virtual void Store_ExtraDist();
	virtual void Finish();

	std::vector<TH1D> Cat0JetLowMt;
	std::vector<TH1D> Cat0JetLowMtSideband;
	std::vector<TH1D> Cat0JetLowMtExtrapolation;
	std::vector<TH1D> Cat0JetHighMt;
	std::vector<TH1D> Cat0JetHighMtSideband;
	std::vector<TH1D> Cat0JetHighMtExtrapolation;
	std::vector<TH1D> Cat1JetLowMt;
	std::vector<TH1D> Cat1JetLowMtSideband;
	std::vector<TH1D> Cat1JetLowMtExtrapolation;
	std::vector<TH1D> Cat1JetHighMt;
	std::vector<TH1D> Cat1JetHighMtSideband;
	std::vector<TH1D> Cat1JetHighMtExtrapolation;
	std::vector<TH1D> Cat1JetBoostMt;
	std::vector<TH1D> Cat1JetBoostMtSideband;
	std::vector<TH1D> Cat1JetBoostMtExtrapolation;
	std::vector<TH1D> CatVBFLooseMt;
	std::vector<TH1D> CatVBFLooseMtSideband;
	std::vector<TH1D> CatVBFLooseRelaxMt;
	std::vector<TH1D> CatVBFLooseRelaxMtExtrapolation;
	std::vector<TH1D> CatVBFTightMt;
	std::vector<TH1D> CatVBFTightMtSideband;
	std::vector<TH1D> CatVBFTightRelaxMt;
	std::vector<TH1D> CatVBFTightRelaxMtExtrapolation;
	std::vector<TH1D> CatInclusiveMt;
	std::vector<TH1D> CatInclusiveMtSideband;
	std::vector<TH1D> CatInclusiveMtExtrapolation;

	std::vector<TH1D> Cat0JetLowMtSS;
	std::vector<TH1D> Cat0JetLowMtSidebandSS;
	std::vector<TH1D> Cat0JetLowMtExtrapolationSS;
	std::vector<TH1D> Cat0JetHighMtSS;
	std::vector<TH1D> Cat0JetHighMtSidebandSS;
	std::vector<TH1D> Cat0JetHighMtExtrapolationSS;
	std::vector<TH1D> Cat1JetLowMtSS;
	std::vector<TH1D> Cat1JetLowMtSidebandSS;
	std::vector<TH1D> Cat1JetLowMtExtrapolationSS;
	std::vector<TH1D> Cat1JetHighMtSS;
	std::vector<TH1D> Cat1JetHighMtSidebandSS;
	std::vector<TH1D> Cat1JetHighMtExtrapolationSS;
	std::vector<TH1D> Cat1JetBoostMtSS;
	std::vector<TH1D> Cat1JetBoostMtSidebandSS;
	std::vector<TH1D> Cat1JetBoostMtExtrapolationSS;
	std::vector<TH1D> CatVBFLooseMtSS;
	std::vector<TH1D> CatVBFLooseMtSidebandSS;
	std::vector<TH1D> CatVBFTightMtSS;
	std::vector<TH1D> CatVBFTightMtSidebandSS;
	std::vector<TH1D> CatInclusiveMtSS;
	std::vector<TH1D> CatInclusiveMtSidebandSS;
	std::vector<TH1D> CatInclusiveMtExtrapolationSS;

	std::vector<TH1D> Cat0JetLowMtAntiIso;
	std::vector<TH1D> Cat0JetHighMtAntiIso;
	std::vector<TH1D> Cat1JetLowMtAntiIso;
	std::vector<TH1D> Cat1JetHighMtAntiIso;
	std::vector<TH1D> Cat1JetBoostMtAntiIso;
	std::vector<TH1D> CatVBFLooseMtAntiIso;
	std::vector<TH1D> CatVBFTightMtAntiIso;
	std::vector<TH1D> CatInclusiveMtAntiIso;
	std::vector<TH1D> Cat0JetLowMtAntiIsoSS;
	std::vector<TH1D> Cat0JetHighMtAntiIsoSS;
	std::vector<TH1D> Cat1JetLowMtAntiIsoSS;
	std::vector<TH1D> Cat1JetHighMtAntiIsoSS;
	std::vector<TH1D> Cat1JetBoostMtAntiIsoSS;
	std::vector<TH1D> CatVBFLooseMtAntiIsoSS;
	std::vector<TH1D> CatVBFTightMtAntiIsoSS;
	std::vector<TH1D> CatInclusiveMtAntiIsoSS;

	std::vector<TH1D> Cat0JetLowQcdAbcd;
	std::vector<TH1D> Cat0JetHighQcdAbcd;
	std::vector<TH1D> Cat1JetLowQcdAbcd;
	std::vector<TH1D> Cat1JetHighQcdAbcd;
	std::vector<TH1D> Cat1JetBoostQcdAbcd;
	std::vector<TH1D> CatVBFLooseQcdAbcd;
	std::vector<TH1D> CatVBFTightQcdAbcd;
	std::vector<TH1D> CatInclusiveQcdAbcd;

	std::vector<TH1D> Cat0JetLowQcdOSMuIso;
	std::vector<TH1D> Cat0JetLowQcdOSTauIso;
	std::vector<TH1D> Cat0JetLowQcdSSMuIso;
	std::vector<TH1D> Cat0JetLowQcdSSTauIso;
	std::vector<TH1D> Cat0JetHighQcdOSMuIso;
	std::vector<TH1D> Cat0JetHighQcdOSTauIso;
	std::vector<TH1D> Cat0JetHighQcdSSMuIso;
	std::vector<TH1D> Cat0JetHighQcdSSTauIso;
	std::vector<TH1D> Cat1JetLowQcdOSMuIso;
	std::vector<TH1D> Cat1JetLowQcdOSTauIso;
	std::vector<TH1D> Cat1JetLowQcdSSMuIso;
	std::vector<TH1D> Cat1JetLowQcdSSTauIso;
	std::vector<TH1D> Cat1JetHighQcdOSMuIso;
	std::vector<TH1D> Cat1JetHighQcdOSTauIso;
	std::vector<TH1D> Cat1JetHighQcdSSMuIso;
	std::vector<TH1D> Cat1JetHighQcdSSTauIso;
	std::vector<TH1D> Cat1JetBoostQcdOSMuIso;
	std::vector<TH1D> Cat1JetBoostQcdOSTauIso;
	std::vector<TH1D> Cat1JetBoostQcdSSMuIso;
	std::vector<TH1D> Cat1JetBoostQcdSSTauIso;
	std::vector<TH1D> CatVBFLooseQcdOSMuIso;
	std::vector<TH1D> CatVBFLooseQcdOSTauIso;
	std::vector<TH1D> CatVBFLooseQcdSSMuIso;
	std::vector<TH1D> CatVBFLooseQcdSSTauIso;
	std::vector<TH1D> CatVBFTightQcdOSMuIso;
	std::vector<TH1D> CatVBFTightQcdOSTauIso;
	std::vector<TH1D> CatVBFTightQcdSSMuIso;
	std::vector<TH1D> CatVBFTightQcdSSTauIso;
	std::vector<TH1D> CatInclusiveQcdOSMuIso;
	std::vector<TH1D> CatInclusiveQcdOSTauIso;
	std::vector<TH1D> CatInclusiveQcdSSMuIso;
	std::vector<TH1D> CatInclusiveQcdSSTauIso;

	std::vector<TH1D> CatInclusiveMtOSChargeSum;
	std::vector<TH1D> CatInclusiveMtSSChargeSum;

	std::vector<TH1D> CatVBFLooseQcdEff;
	std::vector<TH1D> CatVBFTightQcdEff;
	std::vector<TH1D> Cat1JetBoostQcdEff;
};

#endif /* HTOTAUMUTAUHBACKGROUNDS_H_ */
