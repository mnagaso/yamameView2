/*
 * CorrFunc.h
 *
 *  Created on: Jul 16, 2013
 *      Author: nagaso
 */

#ifndef CORRFUNC_H_
#define CORRFUNC_H_

class CorrFunc
{
public:
	CorrFunc(double, double, double, int, double);
	bool readMap();
	void initialize();
	void getRatio();
	void searchNearPoint();
	void searchNearPointandLinearInterpolation();
	void widenMap();

	string filename = "refCorrectionMap2.csv";
	double angle;
	double distance;
	double reflect;
	const static int datapoints = 3248; //points of rotatemap .csv
	double map[datapoints][4]; // 1: distance, 2:angle, 3:ratio, 4:ref. strength
	double ratio;
	double refVal_rot;
	int slide_amount;
	double widen_times;
};




#endif /* CORRFUNC_H_ */
