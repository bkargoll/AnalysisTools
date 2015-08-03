/*
 * Inclusive.h
 *
 *  Created on: Mar 21, 2014
 *      Author: kargoll
 */

#ifndef INCLUSIVE_H_
#define INCLUSIVE_H_

#include "Category.h"

class Inclusive: public Category {
public:
	Inclusive(TString Name_, TString id_);
	virtual ~Inclusive();

protected:
	// no additional cuts are performed for Inclusive category
	// thus nothing needs to be set here except for migrateIntoMain
	bool categorySelection();
};

#endif /* INCLUSIVE_H_ */
