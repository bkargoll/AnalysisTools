/*
 * CategoryStandalone.cxx
 *
 *  Created on: May 8, 2015
 *      Author: kargoll
 */

#include "CategoryStandalone.h"

CategoryStandalone::CategoryStandalone() {
	// set values in constructor of derived class
}

CategoryStandalone::~CategoryStandalone() {
}

void CategoryStandalone::run(){
	setCutValues();
	execute();
}


bool CategoryStandalone::passed() {
	for (unsigned i = 0; i < nCuts_; i++) {
		if ( not passCut_.at(i)) return false;
	}
	return true;
}

