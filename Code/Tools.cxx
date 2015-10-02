#include "Tools.h"

namespace Tools {

	double phiUncertainty(double px, double py, double pxVar, double pyVar, double pxpyCov){
		double global = px / (px*px + py+py);
		double cx = py/px;
		double cy = 1;
		double cxy = cx * cy;

		return global * sqrt( cx*cx*pxVar + cy*cy*pyVar + 2*cxy*pxpyCov);
	}

}
