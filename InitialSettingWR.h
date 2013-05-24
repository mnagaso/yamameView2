/*
 * InitialSettingWR.h
 *
 *  Created on: May 23, 2013
 *      Author: nagaso
 */

#ifndef INITIALSETTINGWR_H_
#define INITIALSETTINGWR_H_

#include <string>
using std::string;

#include <fstream>
using std::fstream;

class InitialSettingWR
{
public:
	int wrflag;
	int xs, xe, ys, ye;
	float zs, ze;
	int mounWid, peakNum;
	float tri_thre, tri_leng;
	int v_point[10000];
	InitialSettingWR(int , int , int , int , int, float, float, int, int, float, float, int , int[] , const string&);
private:


	//int filenum;

	void writeSetting( int, int[10000], string);
    void readSetting(string, int);
};



#endif /* INITIALSETTINGWR_H_ */
