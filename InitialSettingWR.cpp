/*
 * InitialSettingWR.cpp
 *
 *  Created on: May 23, 2013
 *      Author: nagaso
 */

#include <iostream>
using std::ios;
using std::cin;
using std::cout;
using std::endl;
#include <fstream>
using std::ofstream;
using std::ifstream;

#include <string.h>
using std::string;

#include <iostream>
#include <cstdlib>

#include "InitialSettingWR.h"

InitialSettingWR::InitialSettingWR(int flag, int xstart, int xend, int ystart, int yend,
									float zstart, float zend, int mountWidth, int mountPeakNum,
									float trithre, float trileng,  int fn, int vispoint[10000], const string& fishname )
{
	if(flag == 0)
	{
		wrflag = flag;
		xs = xstart;
		xe = xend;
		ys = ystart;
		ye = yend;
		zs = zstart;
		ze = zend;
		mounWid = mountWidth;
		peakNum = mountPeakNum;
		tri_thre = trithre;
		tri_leng = trileng;
		int filenum = fn;
		writeSetting( filenum, vispoint, fishname);
	}
	else
	{
		readSetting(fishname, fn);
	}


}

void InitialSettingWR::writeSetting( int fn, int vpoint[], string fname)
{
	string filename = fname + ".txt";
	const char *cstr = filename.c_str();
	ofstream ofs(cstr);
	ofs << xs << endl
			<< xe << endl
			<< ys << endl
			<< ye << endl
			<< zs << endl
			<< ze << endl
			<< mounWid << endl
			<< peakNum << endl
			<< tri_thre << endl
			<< tri_leng << endl;

	for( int i = 0; i < fn; i++)
		ofs << vpoint[fn] << endl;

}

void InitialSettingWR::readSetting(string fn, int fnum)
{
	string filename = fn + ".txt";
	const char *cstr = filename.c_str();
	ifstream ifs(cstr);
	string str;

	getline(ifs, str);
	xs = atoi(str.c_str());
	getline(ifs, str);
	xe = atoi(str.c_str());
	getline(ifs, str);
	ys = atoi(str.c_str());
	getline(ifs, str);
	ye = atoi(str.c_str());

	getline(ifs, str);
	zs = atof(str.c_str());
	getline(ifs, str);
	ze = atof(str.c_str());

	getline(ifs, str);
	mounWid = atoi(str.c_str());
	getline(ifs, str);
	peakNum = atoi(str.c_str());

	getline(ifs, str);
	tri_thre = atof(str.c_str());
	getline(ifs, str);
	tri_leng = atof(str.c_str());

	for(int i = 0; i < fnum; i++)
	{
		getline(ifs, str);
		v_point[i] = atoi( str.c_str());
	}
}
