/*
 * Category.h
 * Created on: Apr 20, 2015
 *     Author: kargoll
 *
 * Base class for Categories based on HToTaumuTauh analysis.
 * Categories (or general extensions of the main analysis) should derive from
 * this class, and implement the member functions "categoryXXX()". The code within
 * the categories is then executed after the main analysis code. All cuts inside
 * the main analysis and inside the category are evaluated before the plots are
 * created.
 */

#ifndef CATEGORY_H_
#define CATEGORY_H_

#include "HToTaumuTauh.h"

class Category: public HToTaumuTauh {
public:
	Category(TString Name_, TString id_);
	virtual ~Category();

	// override these functions to include category stuff
	void Configure();

protected:
	//*** implementation of the categories goes into these functions ***
	virtual void categoryConfiguration(){};
	virtual void categoryExtradist()	{};
	virtual bool categorySelection()	{return true;};
	virtual void categoryPlotting()		{};

	// override these functions to include category stuff
	virtual void doEvent();
	virtual void Store_ExtraDist();
};

#endif /* CATEGORY_H_ */
