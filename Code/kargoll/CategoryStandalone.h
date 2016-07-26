/*
 * CategoryStandalone.h
 *
 *  Created on: May 8, 2015
 *      Author: kargoll
 *
 *	Base class similar to class "Category", but independent from HToTaumuTauh.
 *	It allows to implement the definition of a given category
 *	and run it, without actually discarding events, overriding plots, etc.
 *	This allows to access the result of the standalone category from within a
 *	different category.
 *
 *	With small adoptions, this class could be used to add additional plots, etc.,
 *	to any analysis class. Have a look into Category.h to see how this could be
 *	implemented.
 */

#ifndef CATEGORYSTANDALONE_H_
#define CATEGORYSTANDALONE_H_

#include <vector>
#include <cmath>

class CategoryStandalone {
public:
	CategoryStandalone();
	virtual ~CategoryStandalone();

	virtual void run();

	virtual bool passed();

	virtual unsigned get_nCuts() const {return nCuts_;}
	virtual const std::vector<float>& get_cutValues() {setCutValues(); return cutValues_;}
	virtual const std::vector<float>& get_eventValues() const {return eventValues_;}
	virtual const std::vector<bool>& get_passCut() const {return passCut_;}

protected:
	unsigned nCuts_;
	std::vector<float> cutValues_;
	std::vector<float> eventValues_;
	std::vector<bool> passCut_;

	virtual void execute() = 0;

	virtual void setCutValues() = 0;
};

#endif /* CATEGORYSTANDALONE_H_ */
