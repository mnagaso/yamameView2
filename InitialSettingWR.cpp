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
									float zstart, float zend,
									float trithre, float trileng,  int fn, int vispoint[10000], const string& fishname,
									int mni, int siv, int stn, float sA, double cV, float bn, float gt, float ssv)
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
		tri_thre = trithre;
		tri_leng = trileng;
		int filenum = fn;
		mountnuminteg = mni;
		surfaceinterval = siv;
		surfthickness = stn;
		surfAlph = sA;
		correctVal = cV;
		brightness = bn;
		geta = gt;
		spinnerSurfVal = ssv;
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
			<< tri_thre << endl
			<< tri_leng << endl
			<< mountnuminteg << endl
			<< surfaceinterval << endl
			<< surfthickness << endl
			<< surfAlph << endl
			<< correctVal << endl
			<< brightness << endl
			<< geta << endl
			<< spinnerSurfVal << endl;

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
	tri_thre = atof(str.c_str());
	getline(ifs, str);
	tri_leng = atof(str.c_str());

	getline(ifs, str);
	mountnuminteg = atoi(str.c_str());
	getline(ifs, str);
	surfaceinterval = atoi(str.c_str());
	getline(ifs, str);
	surfthickness = atoi(str.c_str());
	getline(ifs, str);
	surfAlph = atof(str.c_str());
	getline(ifs, str);
	correctVal = atof(str.c_str());
	getline(ifs, str);
	brightness = atof(str.c_str());
	getline(ifs, str);
	geta = atof(str.c_str());
	getline(ifs, str);
	spinnerSurfVal = atof(str.c_str());
	for(int i = 0; i < fnum; i++)
	{
		getline(ifs, str);
		v_point[i] = atoi( str.c_str());
	}
}
