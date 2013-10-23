/*
 * ShishamoView.cpp
 *
 *  Created on: 2013/02/15
 *      Author: mnsaru
 */

using namespace std;

#include <stdio.h>
#include <float.h>
#include <iostream>
using std::ios;
using std::cin;
using std::cout;
using std::endl;

#include <math.h>

#include <vector>
#include <sys/types.h>
#include <dirent.h>
#include <sys/dir.h>
#include <sys/stat.h>

#include <GL/glui.h>
//#include <GL/glut.h>
//#include <GL/freeglut.h>

#include "ReadCSV.h"
#include "OffsetHilbert.h"
#include "ExtractSurface.h"
#include "InitialSettingWR.h"
#include "CorrFunc.h"
#include "Harmonic.h"

//#define filepath "./data/"
//#define filepath "/var/run/media/nagaso/sotoHD/iwana/2/"
//#define filepath "/var/run/media/nagaso/sotoHD/yamame/1/"
//#define filepath "/var/run/media/nagaso/sotoHD/rotate2/"

#define data_width 2
#define data_length 10000
#define sampling_rate 0.00000001


//iwana1
//#define fishname "iwana1"
//#define x_steps 158
//#define y_steps 41
//#define filepath  "/var/run/media/nagaso/sotoHD/iwana/1/"
//iwana2
//#define fishname "iwana2"
//#define x_steps 157
//#define y_steps 35
//#define filepath "/var/run/media/nagaso/sotoHD/iwana/2/"
//iwana3
//#define fishname "iwana3"
//#define x_steps 144
//#define y_steps 32
//#define filepath "/var/run/media/nagaso/sotoHD/iwana/3/"
//iwana4
//#define fishname "iwana4"
//#define x_steps 137
//#define y_steps 29
//#define filepath "/var/run/media/nagaso/sotoHD/iwana/4/"
//iwana5
//#define fishname "iwana5"
//#define x_steps 154
//#define y_steps 33
//#define filepath "/var/run/media/nagaso/sotoHD/iwana/5/"
//iwana11
//#define fishname "iwana11"
//#define x_steps 141
//#define y_steps 31
//#define filepath "/var/run/media/nagaso/sotoHD/iwana/11/"


//yamame1
#define fishname "yamame1"
#define x_steps 131
#define y_steps 32
#define filepath "/var/run/media/nagaso/sotoHD/yamame/1/"
//yamame2
//#define fishname "yamame2"
//#define x_steps 157
//#define y_steps 35
//#define filepath "/var/run/media/nagaso/sotoHD/yamame/2/"
//yamame3
//#define fishname "yamame3"
//#define x_steps 126
//#define y_steps 28
//#define filepath "/var/run/media/nagaso/sotoHD/yamame/3/"
//yamame4
//#define fishname "yamame4"
//#define x_steps 127
//#define y_steps 28
//#define filepath "/var/run/media/nagaso/sotoHD/yamame/4/"

//rotate2
//#define fishname "rotate2"
//#define x_steps 21
//#define y_steps 76
//rotate
//#define fishname "rotate"
//#define x_steps 21
//#define y_steps 155
//#define filepath "/var/run/media/nagaso/sotoHD/rotato/"

#define filenum x_steps * y_steps
#define sonicvelo 1500 //onsoku m/s
//#define ini_ignore 4000 //for coloring

vector<direct *> entries;

double dataview[filenum][data_length][data_width + 2];//time, envelope, law amp., harmonic law
//double dataview_env[filenum][data_length][data_width];
int point_status_3D[filenum][10000];// secondpeak( == 3), peak( ==2), mountain(==1), under thre(==0)
int integrated_flag[filenum][10000];// integrated( == 1) or not( == 0)
//double integlatedpeakval[filenum];
int visualize_point[filenum] = {0};// numbers of peaks which are indicated
int visualize_point_second[filenum] = {0};// numbers of second surface peaks which are indicated
int mountain[filenum][10000][2];
int secondSurfStart[filenum] = {0};// memo starting point of second mount begining
int viewchanged_flag = 0;
double colorlange;
double colorlange_corrected;
double noiseLev;
int reftimebottom;
//#define headerchangenum 1 //the x axe number one before header name changes
#define thresholdVal spinnerValf

//-------varieties for Delaunay view----
#include "Delaunay2D.h"
std::set<Tercel::Vector> vertices;
std::set<Tercel::Triangle> triangles;
std::set<Tercel::Vector> vertices2;
std::set<Tercel::Triangle> triangles2;
double zpoint[filenum] = {0.0};// memo z cord. of each surface point
double zamp[filenum] = {0.0};//memo z amp at first peak
double zamp_norm[filenum] = {0.0};
double zamp_norm_corrected[filenum] = {0.0};
double ratio[filenum] = {0.0}; //memo correction func
double angle[filenum][4] = {0.0};//0~2: elemonts of normalized vertical vector, 3:the number of patch relating on the nod(filenum)
double zpoint2[filenum] = {0.0};// memo z cord. of each surface point of second surf.
double zamp2[filenum] = {0.0};//memo z amp at first peak of second..
double zamp_norm2[filenum] = {0.0};
double zamp_norm2_corrected[filenum] = {0.0};
double ratio2[filenum] = {0.0};
double angle2[filenum][4] = { 0.0};
double angle_final[filenum] = {0.0};
double calced_mountainWidth[filenum] = {0.0};

//-------varieties for glut-----------
int WinID[3]; //ウィンドウID
int WindowNum = 0;
int WinFlag[2];
const char *WindowName[]={"main_window", "sub_window", "sub2_window"};

//-------varieties for glui-----------
float rotary[16] = {
	1.0, 0.0, 0.0, 0.0,
	0.0, 1.0, 0.0, 0.0,
	0.0, 0.0, 1.0, 0.0,
	0.0, 0.0, 0.0, 1.0
};

float spinnerValf;
float spinnerSurfVal;
int x_begin;
int x_end;
int y_begin;
int y_end;
int y_begin2;
int y_end2;
double z_begin;
double z_end;
float z_front;
float z_back;
float trans_ary[2] = {0.0, 0.0};
float trans_aryZ[] = {0.0};
int rot_deg = 0;
float rot_axs_x = 0;
float rot_axs_y = 0;
float rot_axs_z = 0;
int sub_x;
int sub_y;
float trans_ary_sub[] = {0.0};
float trans_aryZ_sub[] = {0.0};
float trans_ary_sub2[] = {0.0};
float trans_aryZ_sub2[] = {0.0};
float brightness;
float geta;
int pointsize = 1;
int radiobuttonVal;
double peak_threshold;
float tri_threshold;
float tri_length;
int mountainWidth;
int mountainPeakNum;
int addPeakNum;
int mountnuminteg;
int layerflag;
int surfaceinterval;
int surfthickness;
int radiocolorVal;
float surfAlph;
double correctVal;
int correctionfuncOFflag = 0;
int correctionfuncISALREADYGOTTENflag = 0;
int corrpress;
int corrangle;
float corrwiden;
float mount_interval_thre;
int mountHeightAverageFlag = 0;

double getmaxvol()
{
	double volmax = 0;
	for(int j = 0; j < data_length; j++)
	{
		if( volmax < dataview[0][j][1] )
			volmax = dataview[0][j][1];
	}
	cout << "max vol is " << volmax << endl;
	return volmax;
}
double getmaxvol_corrected()
{
	double volmax = 0;
	for(int i = 0; i < filenum; i++)
		for( int j = 0; j < data_length; j++)
		{
			if( volmax < dataview[i][j][1] / ratio[i] )
				volmax = dataview[i][j][1] / ratio[i];
		}
	cout << "max vol is " << volmax << endl;
	return volmax;
}
double gettimemax()
{
	double volmax = 0;
	for(int i = 0; i < filenum; i++)
		for( int j = 0; j < data_length; j++)
		{
			if( volmax < dataview[i][j][0] )
				volmax = dataview[i][j][0];
		}
	return volmax;
}
double gettimemin()
{
	double volmin = 0;
	for(int i = 0; i < filenum; i++)
		for( int j = 0; j < data_length; j++)
		{
			if( volmin > dataview[i][j][1] )
				volmin = dataview[i][j][1];
		}
	return volmin;
}
int gettimebottom()
{
	double valmax = 0;
	int timebottom = 0;
	for( int i = 0; i < data_length; i++)
		if( dataview[0][i][1] > valmax)
		{
			valmax = dataview[0][i][1];
			timebottom = i;
		}
	return timebottom;
}
double makeDistance( double time )
{
	return time * sonicvelo * 1000 / 2;
}
double makeDistance_waveraise( int pointnum)
{
	int peak_point = 0;
	int peakcount = 0;

	for(int i = 0; i < 10000; i++)// loop among data points: to avoid the front noise peak
		if(point_status_3D[pointnum][i] == 2)
		{
			if(peakcount == visualize_point[pointnum]){
				peak_point = i;
				goto OUT;
			}
			peakcount++;
		}
	OUT:

	double time = 0;
	for(int i = 0; i < 10000; i++) //loop par mountain num
		if( mountain[pointnum][i][0] <= peak_point && mountain[pointnum][i][1] >= peak_point)
		{
			time = dataview[pointnum][mountain[pointnum][i][0]][0];
			goto OUT2;
		}
	OUT2:

	return time * sonicvelo * 1000 / 2;
}
void readOnly(int xs, int ys, int count)
{
	bool headerflag = 0;
	string filename = entries[count]->d_name;
	ReadCSV test(filename, xs, ys, headerflag);

	for(int k = 0; k < data_length; k++){
		dataview[count][k][0] = test.data[k][0];
		dataview[count][k][2] = test.data[k][1];
	}
}
void HilbertOnly(int xs, int ys, int count)
{
	double temp[data_length];
	for(int k = 0; k < data_length; k++)
		temp[k] = dataview[count][k][2];

	OffsetHilbert test2(temp);

	for(int k = 0; k < data_length; k++){
		//dataview[count][k][0] = dataview[count][k][0];
		dataview[count][k][1] = test2.data2[k][1];
		//dataview[count][k][1] = test2.data2[k][1] / ratio;
		if(xs == 0 && ys == 0)
			cout << dataview[count][k][1] << endl;
	}

	if(xs == 0 && ys == 0){
			cout << "testing OffsetHilbert..." << endl;
			reftimebottom = gettimebottom(); //get the distance to the mounting board from 0,0 data signal
			noiseLev = test2.average;
	}
}
void HarmonicOnly(int xs, int ys, int count)
{
	double temp[data_length];
	for(int k = 0; k < data_length; k++)
		temp[k] = dataview[count][k][2];

	Harmonic test2(temp);

	for(int k = 0; k < data_length; k++){
		dataview[count][k][3] = test2.data2[k][1];
	}

	if(xs == sub_x && ys == sub_y){
		cout << "test Harmonic..." << endl;
		for(int k = 0; k < data_length; k++)
			cout << dataview[count][k][2] << "," << dataview[count][k][3] << endl;
	}

}
void readAndHilbert(int xs, int ys, int count)
{
	bool headerflag = 0;
	string filename = entries[count]->d_name;
	ReadCSV test(filename, xs, ys, headerflag);

	double temp[data_length];
	for(int k = 0; k < data_length; k++){
		temp[k] = test.data[k][1];
		dataview[count][k][2] = temp[k];
	}

	OffsetHilbert test2(temp);

	for(int k = 0; k < data_length; k++){
		dataview[count][k][0] = test.data[k][0];
		dataview[count][k][1] = test2.data2[k][1];
		//dataview[count][k][1] = test2.data2[k][1] / ratio;
		if(xs == 0 && ys == 0)
			cout << dataview[count][k][1] << endl;
	}
	if(xs == 0 && ys == 0){
			cout << "testing OffsetHilbert..." << endl;
			reftimebottom = gettimebottom(); //get the distance to the mounting board from 0,0 data signal
			noiseLev = test2.average;
	}
}
double matchingAmpFunc(double data0[data_length], double (*data1)[data_width])
{
	double max_0, max_1;
	double ratio;

	for(int i = 0; i < data_length; i++)
	{
		if(max_0 < fabs(data0[i]))
			max_0 = fabs(data0[i]);
		if(max_1 < fabs(data1[i][1]))
			max_1 = fabs(data1[i][1]);
	}

	ratio = max_1 / max_0;
	return ratio;
}

void extractSurface(int x, int y, int count, int point_status[])
{
	double temp[data_length];
	for( int i = 0; i < data_length; i++)
		temp[i] = dataview[count][i][1];

	if(fishname == "rotate")
		reftimebottom = 5000;
	ExtractSurface ES(reftimebottom, temp, spinnerSurfVal, mountainWidth, mountainPeakNum);

	for( int i = 0; i < data_length; i++)
		point_status[i] = 0;

	for( int i = 0; i <data_length; i++)
		point_status[i] = ES.point_status[i];

	if(x == sub_x && y == sub_y){
		peak_threshold = ES.threshold;
	}

	for(int i = 0; i < 10000; i++) //memo mountain_place
	{
		for(int j = 0; j < 2; j++)
		{
			if(i < ES.mountain_num)
				mountain[count][i][j] = ES.mountain_place[i][j];
			else
				mountain[count][i][j] = 0;
		}
	}
}
bool triDistinguish(double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
{
	double crossProduct = sqrt(((y1 - y3)*(z2 - z3) - (z1 - z3)*(y2 - y3)) * ((y1 - y3)*(z2 - z3) - (z1 - z3)*(y2 - y3))
							 + ((z1 - z3)*(x2 - x3) - (x1 - x3)*(z2 - z3)) * ((z1 - z3)*(x2 - x3) - (x1 - x3)*(z2 - z3))
							 + ((x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3)) * ((x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3)));
	/*double crossProduct2D = sqrt(((y1-y3)*(x2-x3)-(x1-x3)*(y2-y3)) * ((y1-y3)*(x2-x3)-(x1-x3)*(y2-y3))
								+ ((x1-x3)*(y2-y3)-(y1-y3)*(x2-x3)) * ((x1-x3)*(y2-y3)-(y1-y3)*(x2-x3)));*/
	double length1 = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
	double length2 = sqrt((x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2));
	double length3 = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) + (z1-z3)*(z1-z3));

	if( crossProduct <= tri_threshold
		&& length1 < tri_length
		&& length2 < tri_length
		&& length3 < tri_length)
		return 0;
	else
		return 1;
}
double calcAngle(double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
{
	double crossProduct1 = (y1 - y3)*(z2 - z3) - (z1 - z3)*(y2 - y3); //x
	double crossProduct2 = (z1 - z3)*(x2 - x3) - (x1 - x3)*(z2 - z3); //y
	double crossProduct3 = (x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3); //z
	double crossAbs = sqrt(crossProduct1 * crossProduct1 + crossProduct2 * crossProduct2 + crossProduct3 * crossProduct3);
	double angle = acos(crossProduct3 / crossAbs) / M_PI * 180;
	return angle;
}
class calcAngle2
{
public:
	calcAngle2(double x1, double x2, double x3, double y1, double y2, double y3, double z1, double z2, double z3)
	{
		cP1 = (y1 - y3)*(z2 - z3) - (z1 - z3)*(y2 - y3); //x
		cP2 = (z1 - z3)*(x2 - x3) - (x1 - x3)*(z2 - z3); //y
		cP3 = (x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3); //z
		crossAbs = sqrt(cP1 * cP1 + cP2 * cP2 + cP3 * cP3);
		angle = acos(cP3 / crossAbs) / M_PI * 180;
		e1 = cP1 / crossAbs;
		e2 = cP2 / crossAbs;
		e3 = cP3 / crossAbs;
	}
	double cP1;
	double cP2;
	double cP3;
	double crossAbs;
	double angle;
	double e1, e2, e3;//normalized vertical vecter
};
void outlierOmit(int points[filenum][10000])
{
	for(int i = 0; i < x_steps; i++)
	{
		double dump_y[y_steps][2]= {0.00000000000000, 0.0000000000000000}; //[0]:,[1]:

		for(int j = 0; j < y_steps; j++){
			int peak_count = 0;
			for( int k = 0; k < 10000; k++)
			{
				if(points[i + j * x_steps][k] == 2  && peak_count == 0)
				{
					dump_y[j][0] = i + j * x_steps;
					dump_y[j][1] = k;
					peak_count++;
				}
			}
		}

		int count = 0;
		double dump[3][2] = {0.0, 0.0};
		for(int j = 0; j < y_steps - 3; j++)
		{
			double dumptemp = 0.0000;
			if( count < 3 && dump_y[j][0] != 0)
			{
				dump[count][0] = dump_y[j][0];
				dump[count][1] = dump_y[j][1];
 				count++;
			}
			else if( count == 3 && dump_y[j][0] != 0)
			{
					dump[0][0] = dump[1][0]; dump[0][1] = dump[1][1];
					dump[1][0] = dump[2][0]; dump[1][1] = dump[2][1];
					dump[2][0] = dump_y[j][0]; dump[2][1] = dump_y[j][1];

					if(abs((dump[0][1] + dump[2][1]) / 2) * 0.8 > dump[1][1]
					    && (dump[1][1] - dump[0][1]) * (dump[1][1] - dump[2][1]) < 0)
					{
						points[int(dump[1][0])][int((dump[0][1] + dump[2][1]) / 2)] = 2;
						int k = 0;
						while(k == int(dump[0][1] + dump[2][1]) / 2)
						{
							points[int(dump[1][0])][k] = 0;
							k++;
						}
					}
			}
		}
	}
}
double calcIntegration(int pointnum)
{
	int peak_point = 0;
	int peakcount = 0;
	for(int i = 0; i < 10000; i++)
		integrated_flag[pointnum][i] = 0;

	for(int i = 0; i < 10000; i++)// loop among data points: to avoid the front noise peak
		if(point_status_3D[pointnum][i] == 2)
		{
			if(peakcount == visualize_point[pointnum]){
				peak_point = i;
				goto OUT;
			}
			peakcount++;
		}
	OUT:

	double integratedVal = 0;
	int integrated_mountWidth = 0;
	for(int i = 0; i < 10000; i++) //loop par mountain num
		if( mountain[pointnum][i][0] <= peak_point && mountain[pointnum][i][1] >= peak_point)
		{
			//integration for first mountain
			for( int j = mountain[pointnum][i][0]; j < mountain[pointnum][i][1]; j++)
			{
				integratedVal += dataview[pointnum][j][1];
				integrated_mountWidth++;
				integrated_flag[pointnum][j] = 1;
			}

			int k = 1;
			//threshold the distance from the mountain end point to next mountain beginning point

			//cout <<  mountain[pointnum][i + k][0] - mountain[pointnum][i + k - 1][1];
			while( mountain[pointnum][i + k][0] - mountain[pointnum][i + k - 1][1] < mount_interval_thre )
			{
				for( int j = mountain[pointnum][i + k][0]; j < mountain[pointnum][i + k][1]; j++)
				{
					integratedVal += dataview[pointnum][j][1];
					integrated_mountWidth++;
					integrated_flag[pointnum][j] = 1;
				}
				k++;
				if( mountain[pointnum][i + k][0] == 0)
					goto OUT2;
			}

			OUT2:
			secondSurfStart[pointnum] = mountain[pointnum][i + k - 1][1];

			goto OUT3;
		}
	OUT3:
	calced_mountainWidth[pointnum] = integrated_mountWidth;
	return integratedVal;
}
double calcIntegSecSurf(int pointnum)
{
	double integVal = 0;
	for(int i = visualize_point_second[pointnum]; i < 10000; i++)
		if(point_status_3D[pointnum][i] == 3 || point_status_3D[pointnum][i] == 4 )
			integVal += dataview[pointnum][i][1];

	return integVal;
}
void getSecondSurf()
{
	for(int i = 0; i < filenum; i++)
	{
		int peakcount = 0;
		if( secondSurfStart[i] != 0)
		{
			int j = 0;
			for(j = secondSurfStart[i] + surfaceinterval; j < 10000; j++)
				if(point_status_3D[i][j] == 2)
				{
					//if( j < 10000 && peakcount != 0)
					if( j < 10000 && peakcount == 0)
					{
						visualize_point_second[i] = j;
						point_status_3D[i][j] = 4; // flag 4 for first peak of second surface
						peakcount++;

						//break;
					}
					else if( j < 10000 && peakcount != 0)
					{
						point_status_3D[i][j] = 3; // flag 3 for the other peaks of second surf.
						peakcount++;
					}
				}
		}
	}
}
void makeArray()
{
	int count = 0;

	for( int i = 0; i < y_steps; i++)
		for( int j = 0; j < x_steps; j++)
		{
			readOnly(i, j, count);
			HilbertOnly(i, j, count);
			//HarmonicOnly(i, j, count);
			//readAndHilbert(j, i, count);
			//extractSurface(i, j, count);
			count++;
		}
	z_begin = makeDistance(0);
	z_end =  makeDistance(gettimemax());
	//calcImpedance();

	cout << "z min/max..." << z_begin << "/" << z_end << endl;
	cout << "last count filenum = " << count << endl;

	colorlange = getmaxvol();
	cout << "max_vol..." << colorlange << endl;
}
//----------------------------functions for indicate main_window-----------------------------------------
void setColor(double incol, double lancol)
{
	float R, G, B;
	//double geta = 0.3;
	//monocolor
	//lancol = 0.02;
	//if( correctionfuncOFflag == 0)
		lancol *= brightness;
	//else
		//lancol *= (brightness * 1000);

	double noise = noiseLev / lancol;
	double unit_width = (1 - noise) / 256;
	double geta_grad = (geta - 1) / 256; //calc. gradient

	double val = (incol / lancol - noise) +  (1 - (incol / lancol)) / unit_width * geta_grad;

	if(radiocolorVal == 1)
	{
		if (val > 0.75)
		{
			R = 1.0;
			G = (1 - val) * 4;
			B = 0.0;
		}
		else if( val > 0.5)
		{
			R = (val - 0.5) * 4;
			G = 1.0;
			B = 0.0;
		}
		else if( val > 0.25)
		{
			R = 0.0;
			G = 1.0;
			B = 0.5 - val;
		}
		else
		{
			R = 0.0;
			G = val * 4;
			B = 1.0;
		}

		if(R >= 1.0)
			R = 1.0;
		if(G >= 1.0)
			G = 1.0;
		if(B >= 1.0)
			B = 1.0;

		glColor4f(R, G, B, 1.0);
	}
	else if( radiocolorVal == 0 && radiobuttonVal == 6)
	{
		R = val;
		G = val;
		B = val;
		if(R >= 1.0)
			R = 1.0;
		if(G >= 1.0)
			G = 1.0;
		if(B >= 1.0)
			B = 1.0;

		glColor4f(R, G, B, surfAlph);
	}
	else
	{
		R = val;
		G = val;
		B = val;
		if(R >= 1.0)
			R = 1.0;
		if(G >= 1.0)
			G = 1.0;
		if(B >= 1.0)
			B = 1.0;
		glColor4f(R, G, B, 1.0);
	}
}
void tri_color(GLdouble a1, GLdouble a2, GLdouble a3)
{
	float R,G,B;
	double ave = (a1 + a2 + a3) / 3;
	R = ave / colorlange;
	G = ave / colorlange;
	B = ave / colorlange;
	//cout << "R G B " << R << " "<< G << " "<< B << " "<< ave << " "<< colorlange << endl;
	if(R >= 1.0)
		R = 1.0;
	if(G >= 1.0)
		G = 1.0;
	if(B >= 1.0)
		B = 1.0;

	glColor3f(R, G, B);
}
void DrawLaw()
{
	int multnum_y_axis = 100;
		int mountain_count = 0;
		bool mountain_inout_flag = 0;
		double y_min, y_max;
		int filecount = sub_y * x_steps + sub_x;

		glPointSize(1.0);
		glBegin(GL_POINTS);

		for(int i = 0; i < data_length; i++){
			GLdouble x_point = makeDistance(dataview[filecount][i][0]);
			GLdouble y_point = dataview[filecount][i][3] / colorlange * multnum_y_axis;
			//double distance = makeDistance(dataview[filecount][i][0]);

			if(	x_point > z_front && x_point < z_back )
			{
				glColor3f(1.0, 1.0, 1.0);
				glVertex2d(x_point, y_point);

				if(y_point > y_max)
					y_max = x_point;
				if(y_point < y_min)
					y_min = x_point;
			}
		}
		glEnd();
}
void DrawEnvelope()
{
	int multnum_y_axis = 100;
	int mountain_count = 0;
	bool mountain_inout_flag = 0;
	double y_min, y_max;
	int filecount = sub_y * x_steps + sub_x;

	glPointSize(1.0);
	glBegin(GL_POINTS);
	if(viewchanged_flag == 0)// drawing for pointview
	{
		for(int i = 0; i < data_length; i++){
			GLdouble x_point = makeDistance(dataview[filecount][i][0]);
			GLdouble y_point = dataview[filecount][i][1] / colorlange * multnum_y_axis;
			double distance = makeDistance(dataview[filecount][i][0]);

			if(	distance > z_front && distance < z_back )
			{
				glColor3f(1.0, 1.0, 1.0);
				glVertex2d(x_point, y_point);

				if(y_point > y_max)
					y_max = distance;
				if(y_point < y_min)
					y_min = distance;
			}
		}
	}
	else if( viewchanged_flag == 3)
	{
		for( int i = 0; i < data_length; i++)
		{
			GLdouble x_point = makeDistance(dataview[filecount][i][0]);
			GLdouble y_point = dataview[filecount][i][1] / colorlange * multnum_y_axis;
			double distance = makeDistance(dataview[filecount][i][0]);

			if(	distance > z_front && distance < z_back )
			{
				glColor3f(1.0, 1.0, 1.0);
				if( i >= secondSurfStart[filecount])
					glColor3f(1.0, 0.0, 0.0);

				glVertex2d(x_point, y_point);

				if(y_point > y_max)
					y_max = distance;
				if(y_point < y_min)
					y_min = distance;
			}
		}
	}
	else //for peak point view and one point view
	{
		int peak_count = 0;
		for(int i = 0; i < data_length; i++){
			GLdouble x_point = makeDistance(dataview[filecount][i][0]);
			GLdouble y_point = dataview[filecount][i][1] / colorlange * multnum_y_axis;
			double distance = makeDistance(dataview[filecount][i][0]);

			if(	distance > z_front && distance < z_back )
			{
				if(point_status_3D[filecount][i] == 0)
					glColor3f(1.0, 1.0, 1.0);
				else if( integrated_flag[filecount][i] == 1)
					glColor3f(1.0, 1.0, 0.0);
				else if( point_status_3D[filecount][i] == 1)
					glColor3f(0.0, 1.0, 0.0);
				else
				{
					if(peak_count == visualize_point[filecount])
					{
						glColor3f(0.0, 0.0, 1.0);
						peak_count++;
					}
					else
					{
						glColor3f(1.0, 0.0, 0.0);
						peak_count++;
					}
				}

				glVertex2d(x_point, y_point);

				if(y_point > y_max)
					y_max = distance;
				if(y_point < y_min)
					y_min = distance;
			}
		}
	}

	glEnd();
	//draw threshold line of point view
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex2d(z_front, thresholdVal * multnum_y_axis);
	glVertex2d(z_back, thresholdVal * multnum_y_axis);
	glEnd();
	//draw threshold line of peak detection
	glBegin(GL_LINES);
	glColor3f(0.0, 0.0, 1.0);
	glVertex2d(z_front, peak_threshold * multnum_y_axis);
	glVertex2d(z_back, peak_threshold * multnum_y_axis);
	glEnd();
}
void strengthNormalization()
{
	if(correctionfuncOFflag == 0)
	{
		double tempmax = 0;
		for(int i = 0; i < filenum; i++)
		{
			if(tempmax < zamp[i])
				tempmax = zamp[i];
		}
		for(int i = 0; i < filenum; i++)
			zamp_norm[i] = zamp[i] / tempmax;
			//zamp[i] /= tempmax;

		if(radiobuttonVal == 6)
		{
			double tempmax = 0;
			for(int i = 0; i < filenum; i++)
			{
				if(tempmax < zamp2[i])
					tempmax = zamp2[i];
			}
			for(int i = 0; i < filenum; i++)
				zamp_norm2[i] = zamp2[i] / tempmax;
		}
	}
	else
	{
		double tempmax = 0;
		for(int i = 0; i < filenum; i++)
		{
			if(tempmax < zamp[i] / ratio[i] && finite(zamp[i] / ratio[i]) != 0 )
				tempmax = zamp[i] / ratio[i];
		}
		for(int i = 0; i < filenum; i++)
			zamp_norm_corrected[i] = zamp[i] / ratio[i] / tempmax;
			//zamp[i] /= tempmax;

		if(radiobuttonVal == 6)
		{
			double tempmax = 0;
			for(int i = 0; i < filenum; i++)
			{
				if(tempmax < zamp2[i])
					tempmax = zamp2[i] / ratio2[i];
			}
			for(int i = 0; i < filenum; i++)
				zamp_norm2_corrected[i] = zamp2[i] / ratio2[i] / tempmax;
		}
	}
}
void calcCorrectionfunc()
{
	typedef std::set<Tercel::Triangle> TriangleSet;
	for( TriangleSet::iterator it = triangles.begin(); it != triangles.end(); ++it) //draw triangulated triangles
	{
		Tercel::Triangle t = *it;
		int x1 = t.p1->x;int x2 = t.p2->x;int x3 = t.p3->x;
		int y1 = t.p1->y;int y2 = t.p2->y;int y3 = t.p3->y;
		int filecount1 = y1 * x_steps + x1;
		int filecount2 = y2 * x_steps + x2;
		int filecount3 = y3 * x_steps + x3;
		GLdouble z1 = zpoint[filecount1];
		GLdouble z2 = zpoint[filecount2];
		GLdouble z3 = zpoint[filecount3];

		if(z1 != 0 && z2 != 0 && z3 != 0
			&& z1 > z_front && z2 > z_front && z3 > z_front
			&& z1 < z_back && z2 < z_back && z3 < z_back)
		{
			if(triDistinguish(x1,x2,x3,y1,y2,y3,z1,z2,z3) == 0){
				calcAngle2 cA(x1,x2,x3,y1,y2,y3,z1,z2,z3);
				angle[filecount1][0] += cA.e1; angle[filecount1][1] += cA.e2; angle[filecount1][2] += cA.e3; angle[filecount1][3] +=1;
				angle[filecount2][0] += cA.e1; angle[filecount2][1] += cA.e2; angle[filecount2][2] += cA.e3; angle[filecount2][3] +=1;
				angle[filecount3][0] += cA.e1; angle[filecount3][1] += cA.e2; angle[filecount3][2] += cA.e3; angle[filecount3][3] +=1;
			}
		}
	}

	for(int i = 0; i < filenum; i++)
	{
		if(angle[i][3] != 0)
		{
			double Abs = sqrt( angle[i][0]*angle[i][0] + angle[i][1]*angle[i][1] + angle[i][2]*angle[i][2]);
			double Ag = acos( angle[i][2] / Abs ) / M_PI * 180;
			/*if(Ag < temp_smallest_angle)
			{
				temp_smallest_angle = Ag;
				smallest_filenum = i;
			}*/
			CorrFunc CF(Ag, zpoint[i], zamp[i], corrangle, corrwiden);
			//ratio[i] = CF.ratio;
			ratio[i] = CF.ratio / corrpress + (corrpress - 1) / corrpress;
			angle_final[i] = Ag;
			//cout << zpoint[i] << "," << Ag  <<  "," << zamp[i] << endl;
		}
		else
			ratio[i] = 99999;
	}

	//CorrFunc CFF(temp_smallest_angle, zpoint[smallest_filenum], zamp[smallest_filenum]);
	//cout << "smallest angle" << temp_smallest_angle << "ref strength rot" << CFF.refVal_rot << endl;

	if( radiobuttonVal == 6)
	{

		typedef std::set<Tercel::Triangle> TriangleSet2;
		for( TriangleSet2::iterator it2 = triangles2.begin(); it2 != triangles2.end(); ++it2) //draw triangulated triangles
		{
			Tercel::Triangle t2 = *it2;
			int x1 = t2.p1->x;int x2 = t2.p2->x;int x3 = t2.p3->x;
			int y1 = t2.p1->y;int y2 = t2.p2->y;int y3 = t2.p3->y;
			int filecount1 = y1 * x_steps + x1;
			int filecount2 = y2 * x_steps + x2;
			int filecount3 = y3 * x_steps + x3;
			GLdouble z1 = zpoint2[filecount1];
			GLdouble z2 = zpoint2[filecount2];
			GLdouble z3 = zpoint2[filecount3];

			if(z1 != 0 && z2 != 0 && z3 != 0
				&& z1 > z_front && z2 > z_front && z3 > z_front
				&& z1 < z_back && z2 < z_back && z3 < z_back
				&& y1 > y_begin2 && y2 > y_begin2 && y3 > y_begin2
				&& y1 < y_end2 && y2 < y_end2 && y3 < y_end2)
			{
				if(triDistinguish(x1,x2,x3,y1,y2,y3,z1,z2,z3) == 0){
					calcAngle2 cA(x1,x2,x3,y1,y2,y3,z1,z2,z3);
					angle2[filecount1][0] += cA.e1; angle2[filecount1][1] += cA.e2; angle2[filecount1][2] += cA.e3; angle2[filecount1][3] +=1;
					angle2[filecount2][0] += cA.e1; angle2[filecount2][1] += cA.e2; angle2[filecount2][2] += cA.e3; angle2[filecount2][3] +=1;
					angle2[filecount3][0] += cA.e1; angle2[filecount3][1] += cA.e2; angle2[filecount3][2] += cA.e3; angle2[filecount3][3] +=1;
				}
			}
		}

		for(int i = 0; i < filenum; i++)
		{
			if(angle2[i][3] != 0)
			{
				double Abs = sqrt( angle2[i][0]*angle2[i][0] + angle2[i][1]*angle2[i][1] + angle2[i][2]*angle2[i][2]);
				double Ag = acos( angle2[i][2] / Abs ) / M_PI * 180;
				CorrFunc CF(Ag, zpoint2[i], zamp2[i], corrangle, corrwiden);
				//ratio2[i] = CF.ratio;
				ratio2[i] = CF.ratio / corrpress + (corrpress - 1) / corrpress;
			}
			else
				ratio2[i] = 99999;
		}
	}
	strengthNormalization();
}

void invokeDelaunay()
{
	for( int i = 0; i < filenum; i++)
	{
		if( visualize_point[i] != 99)
		{
			GLdouble y_point = int( i / x_steps);
			GLdouble x_point = i % x_steps;
			int top_point_count = 0;

			if( x_point > x_begin
					&& x_point < x_end
					&& y_point > y_begin
					&& y_point < y_end
				)
			{

				for( int j  = 0; j < data_length; j++)
				{
					if( top_point_count != visualize_point[i] && point_status_3D[i][j] == 2)
						top_point_count++;
					else if( top_point_count == visualize_point[i] && point_status_3D[i][j] == 2)
					{
						Tercel::Vector v;
						v.x = x_point;
						v.y = y_point;
						vertices.insert(v);

						zpoint[i] = makeDistance_waveraise(i);
						zamp[i] = calcIntegration(i);
						top_point_count++;
					}
				}
				if( radiobuttonVal == 6 && visualize_point_second[i] != 0)
				{
					Tercel::Vector v2;
					v2.x = x_point;
					v2.y = y_point;
					vertices2.insert(v2);

					zpoint2[i] = makeDistance(dataview[i][visualize_point_second[i]][0]);
					zamp2[i] = calcIntegSecSurf(i);
				}
			}
		}
	}
	strengthNormalization();
	if(radiobuttonVal == 6) // inner delaunay
	{
		Tercel::Delaunay2d::getDelaunayTriangles( vertices2, &triangles2); // Delaunay triangulation invoke
	}
	Tercel::Delaunay2d::getDelaunayTriangles( vertices, &triangles); // Delaunay triangulation invoke
}
void Draw3d()
{
	if(radiobuttonVal != 4 && radiobuttonVal != 5 && radiobuttonVal != 6 && radiobuttonVal != 7)//drawing not for delaunay view
	{
		glPointSize(pointsize);
		glBegin(GL_POINTS);
		for( int i = 0; i < filenum; i++)
		{
			double y_measure;
			if(fishname != "rotate" && fishname != "rotate2")
				y_measure = 1.0;
			else if( fishname == "rotate")
				y_measure = 0.1;
			else
				y_measure = 0.2;
			GLdouble y_point = int( i / x_steps) * y_measure;
			GLdouble x_point = i % x_steps;

			if( radiobuttonVal == 0) //point view
			{
				for(int j = 0; j < data_length; j++)
				{
					GLdouble z_point = makeDistance(dataview[i][j][0]);

					if( fabs(dataview[i][j][1] / colorlange) > thresholdVal
							&& x_point > x_begin
							&& x_point <= x_end
							&& y_point > y_begin
							&& y_point <= y_end
							&& z_point > z_front
							&& z_point <= z_back
						)
					{
						setColor( dataview[i][j][1], colorlange);
						if(x_point == sub_x && y_point == sub_y)
							glColor3f(1.0, 0.0, 0.0);

						glVertex3d(x_point, y_point, z_point);
					}
				}
			}
			else if( radiobuttonVal == 1) //peak point view
			{
				for(int j = 0; j < data_length; j++)
				{
					if(point_status_3D[i][j] == 2)// indicate peak only
					{
						GLdouble z_point = makeDistance(dataview[i][j][0]);
						if( x_point > x_begin
							&& x_point <= x_end
							&& y_point > y_begin
							&& y_point <= y_end
							&& z_point > z_front
							&& z_point <= z_back
							)
						{
							setColor( dataview[i][j][1], colorlange);
							if(x_point == sub_x && y_point == sub_y)
								glColor3f(1.0, 0.0, 0.0);

							glVertex3d(x_point, y_point, z_point);
						}
					}

				}
			}
			else if( radiobuttonVal == 2) //one point view
			{
				int top_point_count = 0;
				for(int j = 0; j < data_length; j++)
				{
					if(top_point_count != visualize_point[i] && point_status_3D[i][j] == 2 )
						top_point_count++;
					else if(top_point_count == visualize_point[i] && point_status_3D[i][j] == 2 )
					{
						GLdouble z_point = makeDistance_waveraise(i);
						if( x_point > x_begin
							&& x_point <= x_end
							&& y_point > y_begin
							&& y_point <= y_end
							&& z_point > z_front
							&& z_point <= z_back
							)
						{
							double integCol = calcIntegration(i);
							setColor( integCol, colorlange);
							if(x_point == sub_x && y_point == sub_y)
								glColor3f(1.0, 0.0, 0.0);

							glVertex3d(x_point, y_point, z_point);
						}
						top_point_count++;
					}
				}
			}
			else if( radiobuttonVal == 3) //inner view
			{
				int top_point_count = 0;
				for(int j = 0; j < data_length; j++)
				{
					if(visualize_point[i] != top_point_count && point_status_3D[i][j] == 2 )
						top_point_count++;
					else if(point_status_3D[i][j] == 2 && top_point_count == visualize_point[i])
					{
						GLdouble z_point = makeDistance(dataview[i][j][0]);
						if( x_point > x_begin
							&& x_point <= x_end
							&& y_point > y_begin
							&& y_point <= y_end
							&& z_point > z_front
							&& z_point <= z_back
							)
						{
							double integCol = calcIntegration(i);
							radiocolorVal = 0;
							setColor( integCol, colorlange);
							if(x_point == sub_x && y_point == sub_y)
								glColor3f(1.0, 0.0, 0.0);

							glVertex3d(x_point, y_point, z_point);
						}
						top_point_count++;
					}
					else if(point_status_3D[i][j] == 4)
					{
						GLdouble z_point = makeDistance(dataview[i][j][0]);
						if( x_point > x_begin
							&& x_point <= x_end
							&& y_point > y_begin
							&& y_point <= y_end
							&& z_point > z_front
							&& z_point <= z_back
							)
						{
							double integCol = calcIntegSecSurf(i);
							radiocolorVal = 1;
							setColor( integCol, colorlange );
							if(x_point == sub_x && y_point == sub_y)
								glColor3f(1.0, 0.0, 0.0);

							glVertex3d(x_point, y_point, z_point);
						}
						//top_point_count++;
					}
				}
			}
		}
		glEnd();
	}
	else if( radiobuttonVal == 4 || radiobuttonVal == 5 || radiobuttonVal == 6) // for delaunay view
	{
		if(radiobuttonVal == 6){
			typedef std::set<Tercel::Triangle> TriangleSet2;
			for( TriangleSet2::iterator it2 = triangles2.begin(); it2 != triangles2.end(); ++it2) //draw triangulated triangles
			{
				Tercel::Triangle t2 = *it2;
				int x1 = t2.p1->x;int x2 = t2.p2->x;int x3 = t2.p3->x;
				int y1 = t2.p1->y;int y2 = t2.p2->y;int y3 = t2.p3->y;
				int filecount1 = y1 * x_steps + x1;
				int filecount2 = y2 * x_steps + x2;
				int filecount3 = y3 * x_steps + x3;
				GLdouble z1 = zpoint2[filecount1];
				GLdouble z2 = zpoint2[filecount2];
				GLdouble z3 = zpoint2[filecount3];
				GLdouble amp1;
				GLdouble amp2;
				GLdouble amp3;
				if(correctionfuncOFflag == 0)
				{
					amp1 = zamp_norm2[filecount1];
					amp2 = zamp_norm2[filecount2];
					amp3 = zamp_norm2[filecount3];
				} else
				{
					amp1 = zamp_norm2_corrected[filecount1];
					amp2 = zamp_norm2_corrected[filecount2];
					amp3 = zamp_norm2_corrected[filecount3];
				}

				if(z1 != 0 && z2 != 0 && z3 != 0
						&& z1 > z_front && z2 > z_front && z3 > z_front
						&& z1 <= z_back && z2 <= z_back && z3 <= z_back
						&& y1 > y_begin2 && y2 > y_begin2 && y3 > y_begin2
						&& y1 <= y_end2 && y2 <= y_end2 && y3 <= y_end2)
				{
					if(triDistinguish(x1,x2,x3,y1,y2,y3,z1,z2,z3) == 0){
						//double patch_angle = calcAngle(x1,x2,x3,y1,y2,y3,z1,z2,z3);
						//cout << patch_angle << endl;

						glBegin(GL_TRIANGLES);
						radiocolorVal = 1;
						setColor(amp1, 1);
						glVertex3d(x1, y1, z1);
						radiocolorVal = 1;
						setColor(amp2, 1);
						glVertex3d(x2, y2, z2);
						radiocolorVal = 1;
						setColor(amp3, 1);
						glVertex3d(x3, y3, z3);
						glEnd();
					}
				}
			}
		}

		typedef std::set<Tercel::Triangle> TriangleSet;
		for( TriangleSet::iterator it = triangles.begin(); it != triangles.end(); ++it) //draw triangulated triangles
		{
			Tercel::Triangle t = *it;
			int x1 = t.p1->x;int x2 = t.p2->x;int x3 = t.p3->x;
			int y1 = t.p1->y;int y2 = t.p2->y;int y3 = t.p3->y;
			int filecount1 = y1 * x_steps + x1;
			int filecount2 = y2 * x_steps + x2;
			int filecount3 = y3 * x_steps + x3;
			GLdouble z1 = zpoint[filecount1];
			GLdouble z2 = zpoint[filecount2];
			GLdouble z3 = zpoint[filecount3];

			GLdouble amp1;
			GLdouble amp2;
			GLdouble amp3;
			if(correctionfuncOFflag == 0)
			{
				//amp1 = zamp_norm[filecount1];
				//amp2 = zamp_norm[filecount2];
				//amp3 = zamp_norm[filecount3];
				amp1 = zamp[filecount1] / colorlange;
				amp2 = zamp[filecount2] / colorlange;
				amp3 = zamp[filecount3] / colorlange;
			} else
			{
				/*
				amp1 = zamp_norm_corrected[filecount1];
				amp2 = zamp_norm_corrected[filecount2];
				amp3 = zamp_norm_corrected[filecount3];
				*/
				amp1 = zamp[filecount1] / ratio[filecount1];
				amp2 = zamp[filecount2] / ratio[filecount2];
				amp3 = zamp[filecount3] / ratio[filecount3];
			}

			if(z1 != 0 && z2 != 0 && z3 != 0
				&& z1 > z_front && z2 > z_front && z3 > z_front
				&& z1 <= z_back && z2 <= z_back && z3 <= z_back
				&& x1 > x_begin
				&& x1 <= x_end
				&& y1 > y_begin
				&& y1 <= y_end)
			{
				if(triDistinguish(x1,x2,x3,y1,y2,y3,z1,z2,z3) == 0){
					if( radiobuttonVal == 4)
					{
						glBegin(GL_TRIANGLES);
						setColor(amp1, 1);
						glVertex3d(x1, y1, z1);
						setColor(amp2, 1);
						glVertex3d(x2, y2, z2);
						setColor(amp3, 1);
						glVertex3d(x3, y3, z3);
						glEnd();
					}
					else if(radiobuttonVal == 5)
					{
						glBegin(GL_LINE_STRIP);
						glVertex3d(x1, y1, z1);
						glVertex3d(x2, y2, z2);
						glVertex3d(x3, y3, z3);
						glVertex3d(x1, y1, z1);
						glEnd();
					}
					else if( radiobuttonVal == 6)
					{
						glBegin(GL_TRIANGLES);
						radiocolorVal = 0;
						setColor(amp1, 1);
						glVertex3d(x1, y1, z1);
						radiocolorVal = 0;
						setColor(amp2, 1);
						glVertex3d(x2, y2, z2);
						radiocolorVal = 0;
						setColor(amp3, 1);
						glVertex3d(x3, y3, z3);
						glEnd();
					}
				}
			}
		}
	}
	else if( radiobuttonVal == 7)
	{
		typedef std::set<Tercel::Triangle> TriangleSet;
		for( TriangleSet::iterator it = triangles.begin(); it != triangles.end(); ++it) //draw triangulated triangles
		{
			Tercel::Triangle t = *it;
			int x1 = t.p1->x;int x2 = t.p2->x;int x3 = t.p3->x;
			int y1 = t.p1->y;int y2 = t.p2->y;int y3 = t.p3->y;
			int filecount1 = y1 * x_steps + x1;
			int filecount2 = y2 * x_steps + x2;
			int filecount3 = y3 * x_steps + x3;
			GLdouble z1 = zpoint[filecount1];
			GLdouble z2 = zpoint[filecount2];
			GLdouble z3 = zpoint[filecount3];
			GLdouble amp1 = zamp[filecount1];
			GLdouble amp2 = zamp[filecount2];
			GLdouble amp3 = zamp[filecount3];
			if(correctionfuncOFflag == 1)
			{
				amp1 = ratio[filecount1];
				amp2 = ratio[filecount2];
				amp3 = ratio[filecount3];
			}

			if(z1 != 0 && z2 != 0 && z3 != 0
				&& z1 > z_front && z2 > z_front && z3 > z_front
				&& z1 <= z_back && z2 <= z_back && z3 <= z_back)
			{
				if(triDistinguish(x1,x2,x3,y1,y2,y3,z1,z2,z3) == 0
					&& x1 > x_begin
					&& x1 <= x_end
					&& y1 > y_begin
					&& y1 <= y_end
					&& z1 > z_front
					&& z1 <= z_back
					)
				{
					glBegin(GL_TRIANGLES);
					setColor(amp1, colorlange);
					glVertex3d(x1, y1, z1);
					setColor(amp2, colorlange);
					glVertex3d(x2, y2, z2);
					setColor(amp3, colorlange);
					glVertex3d(x3, y3, z3);
					glEnd();
					//draw normalized vector
					glBegin(GL_LINES);

					glColor3d(1.0,0.0,0.0);
					glVertex3d(x1,y1,z1);
					double Abs1 = sqrt( angle[filecount1][0]*angle[filecount1][0] + angle[filecount1][1]*angle[filecount1][1] + angle[filecount1][2]*angle[filecount1][2]);
					glVertex3d(angle[filecount1][0]/Abs1 + x1, angle[filecount1][1]/Abs1 + y1, angle[filecount1][2]/Abs1 + z1);
					glColor3d(1.0,1.0,0.0);
					glVertex3d(x1,y1,z1);
					glVertex3d(x1,y1,z1+1);

					glColor3d(1.0,0.0,0.0);
					glVertex3d(x2,y2,z2);
					double Abs2 = sqrt( angle[filecount2][0]*angle[filecount2][0] + angle[filecount2][1]*angle[filecount2][1] + angle[filecount2][2]*angle[filecount2][2]);
					glVertex3d(angle[filecount2][0]/Abs2 + x2, angle[filecount2][1]/Abs2 + y2, angle[filecount2][2]/Abs2 + z2);
					glColor3d(1.0,1.0,0.0);
					glVertex3d(x2,y2,z2);
					glVertex3d(x2,y2,z2+1);

					glColor3d(1.0,0.0,0.0);
					glVertex3d(x3,y3,z3);
					double Abs3 = sqrt( angle[filecount3][0]*angle[filecount3][0] + angle[filecount3][1]*angle[filecount3][1] + angle[filecount3][2]*angle[filecount3][2]);
					glVertex3d(angle[filecount3][0]/Abs3 + x3, angle[filecount3][1]/Abs3 + y3, angle[filecount3][2]/Abs3 + z3);
					glColor3d(1.0,1.0,0.0);
					glVertex3d(x3,y3,z3);
					glVertex3d(x3,y3,z3+1);

					glEnd();
				}
			}
		}
	}
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();
	gluLookAt( 0.0, 0.0, 300.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	//glRotatef( 90.0, 0.0, 0.0, 1.0);
	glRotatef( rot_deg, rot_axs_x, rot_axs_y, rot_axs_z);
	glTranslated( -1 * x_steps / 2, -1 * y_steps / 2, 0.0);

	glTranslatef(trans_ary[0],trans_ary[1], trans_aryZ[0]);
	glMultMatrixf( rotary );

	//Draw2d();

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	Draw3d();
	glPopMatrix();

	glutSwapBuffers();
	glutPostRedisplay();
}


void gluiCallback(int num)
{
	exit(0);
}
void funcReset(int num)
{
	for(int i = 0; i < 16; i++){
		rotary[i] = 0.0;
		if(i == 0 || i == 5 || i == 10 || i == 15)
			rotary[i] = 1.0;
	}
	GLUI_Master.sync_live_all();
}
void funcResetTrans(int num)
{
	trans_ary[0] = 0.0; trans_ary[1] = 0.0;
	trans_aryZ[0] = 0.0;
}
void funcResetSub(int num)
{
	sub_x = 99;
	sub_y = 99;
}

void init()
{
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glEnable(GL_DEPTH_TEST | GL_CULL_FACE );
	glCullFace(GL_BACK);
}

void idle()
{
	for(int loop = 0; loop < WindowNum; ++loop){
		if(WinFlag[loop]>0)
		{
			//printf("idle, loop=%d, WinFlag=%d\n", loop, WinFlag[loop]);
			glutSetWindow(WinID[loop]);
			glutPostRedisplay(); //redraw (calls display() fnc)
		}
	}
}


void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 'q':
	case 'Q':
	case '\033': // esc
		exit(0);
	default:
		break;
	}
}

void resize( GLsizei w, GLsizei h)
{
	GLfloat fAspect;
	//configure region for drawing picture in window
	glViewport(0, 0, w, h);

	fAspect = (GLfloat)w / (GLfloat)h;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective( 30.0f, fAspect, 1.0, 600.0);//sawaranai

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Init()
{

	glutInitWindowPosition(0,0);
	glutInitWindowSize(1000,1000);
}

void callBacks()
{
	//call back func
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(resize);
	init();
}


//-----------------------functions to set up sub window-------------------------
void display2()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glPushMatrix();
	gluLookAt( 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	//glRotatef( 90.0, 0.0, 0.0, 1.0);
	//glRotatef( 180.0, 0.0, 1.0, 0.0);
	glTranslated(0.0, 0.0, 0.0);

	glTranslatef(trans_ary_sub[0], 0.0, trans_aryZ_sub[0]);

	DrawEnvelope();

	glPopMatrix();

	glutSwapBuffers();
	glutPostRedisplay();
}

void resize2( GLsizei w, GLsizei h)
{
	GLfloat fAspect;
	//configure region for drawing picture in window
	glViewport(0, 0, w, h);

	fAspect = (GLfloat)w / (GLfloat)h;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective( 30.0f, fAspect, 1.0, 600.0);//sawaranai

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Init2()
{
	glutInitWindowPosition(1000,0);
	glutInitWindowSize(800,500);
}
void callBacks2()
{
	glutDisplayFunc(display2);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(resize2);

	//init();
}
//-----------------------functions to set up sub2 window-------------------------
void display3()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glPushMatrix();
	gluLookAt( 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	//glRotatef( 90.0, 0.0, 0.0, 1.0);
	//glRotatef( 180.0, 0.0, 1.0, 0.0);
	glTranslated(0.0, 0.0, 0.0);

	glTranslatef(trans_ary_sub2[0], 0.0, trans_aryZ_sub2[0]);

	DrawLaw();

	glPopMatrix();

	glutSwapBuffers();
	glutPostRedisplay();
}

void resize3( GLsizei w, GLsizei h)
{
	GLfloat fAspect;
	//configure region for drawing picture in window
	glViewport(0, 0, w, h);

	fAspect = (GLfloat)w / (GLfloat)h;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective( 30.0f, fAspect, 1.0, 600.0);//sawaranai

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Init3()
{
	glutInitWindowPosition(1000,0);
	glutInitWindowSize(800,500);
}
void callBacks3()
{
	glutDisplayFunc(display3);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(resize3);

	//init();
}



char* targetdir;

int dotOmit(struct dirent *darray)
{
	if( strcmp( darray->d_name,".") && strcmp( darray->d_name, ".."))
	{
		return 1;
	}

	return 0;
}

int timesort(const dirent **v1, const dirent **v2)
{
	struct dirent *d1, *d2;
	struct stat s1, s2;
	char fnm[256];

	d1 = *(struct dirent **)v1;
	d2 = *(struct dirent **)v2;

	//cout << targetdir << endl;
	//cout << fnm << endl;

	strcpy(fnm, targetdir);
	strcat(fnm, d1->d_name);
	stat(fnm, &s1);

	strcpy(fnm, targetdir);
	strcat(fnm, d2->d_name);
	stat(fnm, &s2);

	return (int)(s1.st_mtime - s2.st_mtime);
}

int fileFilter(const direct *dir)
{
	const char * s = dir->d_name;
	int csvlen = strlen(s) - 4;
	if ( csvlen >= 0 )
		if ( strncmp(s+csvlen,".csv",4) == 0)
			return 1;
	return 0;
}

void getDirEntries(const string& dirname, vector<direct *> &entries)
{
	direct** darray;

	targetdir = const_cast<char*>(dirname.c_str());
	//cout << dirname << "  " << targetdir << endl;
	int entryCount = scandir(const_cast<char*>(dirname.c_str()), &darray, fileFilter, timesort);

	for(int k = 0; k < entryCount; k++)
	{
		entries.push_back(darray[k]);
	}
}

void getFileList()
{
	string filename = filepath;
	//vector<direct *> entries;
	getDirEntries(filename, entries);

	cout << "making file list... " << endl;

	for(int k = 0; k < entries.size(); k++){
		cout << k << ".\t" << entries[k]->d_name << endl;
	}

}

void viewchange_callback(int val)
{
	//cout << "val is" << val << " " << radiobuttonVal << endl;
	if( radiobuttonVal == 0){
		cout << "point view" << endl;
		viewchanged_flag = 0;
	}
	else if( radiobuttonVal == 1)
	{
		cout << "peak point view" << endl;
		viewchanged_flag = 1;
	}
	else if ( radiobuttonVal == 2)
	{
		cout << "one point view" << endl;
		viewchanged_flag = 2;
	}
	else if( radiobuttonVal == 3)
	{
		cout << "inner view" << endl;
		viewchanged_flag = 3;
	}
	else if( radiobuttonVal == 4)
	{
		cout << "delaunay view" << endl;
		viewchanged_flag = 4;
	}
	else if( radiobuttonVal == 5)
	{
		cout << "delaunay wire" << endl;
		viewchanged_flag = 5;
	}
	else if( radiobuttonVal ==6)
	{
		cout << "mult layer delaunay" << endl;
		viewchanged_flag = 6;
	}
	else
	{
		cout << "correction ratio view" << endl;
		viewchanged_flag = 7;
	}


	for(int i = 0; i < filenum; i++)
	{
		int y_point = int( i / x_steps);
		int x_point = i % x_steps;
		int point_status[10000];
		extractSurface( x_point, y_point, i, point_status);
		for(int j = 0; j < 10000; j++)
			point_status_3D[i][j] = point_status[j];
	}

	if( radiobuttonVal == 3 || radiobuttonVal == 6)// extract second surface
	{
		getSecondSurf();
	}

	if( radiobuttonVal == 4 || radiobuttonVal == 5 || radiobuttonVal == 6 || radiobuttonVal == 7)
	{
		outlierOmit(point_status_3D);
		invokeDelaunay();
	}
	//viewchanged_flag = 1;
	cout << "peak_place list has changed" << endl;

	GLUI_Master.sync_live_all();
}
void correctionOF_callback(int val)
{
	if(correctionfuncOFflag == 1 && correctionfuncISALREADYGOTTENflag == 0){
		calcCorrectionfunc();
		correctionfuncISALREADYGOTTENflag = 1;
	}

}
void correctionPress_callback(int val)
{
	calcCorrectionfunc();
}

void modifyPeak_callback(int val)
{
	int modfilenum = sub_x + sub_y * x_steps;

	if(visualize_point[modfilenum] != addPeakNum)
		visualize_point[modfilenum] = addPeakNum;


	viewchange_callback(radiobuttonVal);

}


void settingWR_callback( int val)
{
	if(val == 0)
	{
		//InitialSettingWR W(val, x_begin, x_end, y_begin, y_end, z_front, z_back, mountainWidth, mountainPeakNum, tri_threshold, tri_length, filenum, visualize_point, fishname);
		InitialSettingWR W(val, x_begin, x_end, y_begin, y_end, z_front, z_back, tri_threshold, tri_length, filenum, visualize_point, fishname, mountnuminteg, surfaceinterval, surfthickness, surfAlph, correctVal, brightness, geta, spinnerSurfVal);
	}
	else
	{
		//InitialSettingWR R(val, x_begin, x_end, y_begin, y_end, z_front, z_back, mountainWidth, mountainPeakNum, tri_threshold, tri_length, filenum, visualize_point, fishname);
		InitialSettingWR R(val, x_begin, x_end, y_begin, y_end, z_front, z_back, tri_threshold, tri_length, filenum, visualize_point, fishname, mountnuminteg, surfaceinterval, surfthickness, surfAlph, correctVal, brightness, geta, spinnerSurfVal);
		x_begin = R.xs;
		x_end = R.xe;
		y_begin = R.ys;
		y_end = R.ye;
		z_front = R.zs;
		z_back = R.ze;
		tri_threshold = R.tri_thre;
		tri_length = R.tri_leng;
		mountnuminteg = R.mountnuminteg;
		surfaceinterval = R.surfaceinterval;
		surfthickness = R.surfthickness;
		surfAlph = R.surfAlph;
		correctVal = R.correctVal;
		brightness = R.brightness;
		geta = R.geta;
		spinnerSurfVal = R.spinnerSurfVal;

		for(int i = 0; i < filenum; i++)
			visualize_point[i] = R.v_point[i];

		//viewchange_callback(radiobuttonVal);
		GLUI_Master.sync_live_all();
	}
}
void valWrite_callback(int val)
{
	/*
	int filecount = sub_y * x_steps + sub_x;
	int multnum_y_axis = 100;
	for(int i = 0; i < data_length; i++){
		GLdouble x_point = makeDistance(dataview[filecount][i][0]);
		GLdouble y_point = dataview[filecount][i][1] / colorlange * multnum_y_axis;

		cout << x_point << "," << y_point << endl;
	}
*/


	//cout << "x, " << "y, " << "z, " << "ref. str., " << "corr. func(ratio), " << "angle, " << endl;


	cout << "***" << endl;
	cout << "***" << endl;
	cout << "***" << endl;
	cout << "***" << endl;
	cout << "***" << endl;
	//int count_strength[10] = {0};
	for(int i = 0; i < filenum; i++)
	{
		int y_point = int( i / x_steps);
		int x_point = i % x_steps;

		if(correctionfuncOFflag == 0)
		{
			if(mountHeightAverageFlag == 0){
				if(i == 0)
					cout << zamp[i] << ",";
				else if(i != 0 && int( (i + 1) / x_steps) == y_point)
					cout << zamp[i] << ",";
				else if( i != 0 && int( (i + 1) / x_steps) != y_point)
					cout << zamp[i] << "," << endl;
			} else
			{
				double value;
				if(calced_mountainWidth[i] != 0)
					value = zamp[i] / calced_mountainWidth[i];
				else
					value = 0;

				if(i == 0)
					cout << value << ",";
				else if(i != 0 && int( (i + 1) / x_steps) == y_point)
					cout << value << ",";
				else if( i != 0 && int( (i + 1) / x_steps) != y_point)
					cout << value << "," << endl;
			}
		}
		else if(correctionfuncOFflag == 1)
		{
			if(i == 0)
				cout << zamp[i] / ratio[i] << ",";
			else if(i != 0 && int( (i + 1) / x_steps) == y_point)
				cout << zamp[i] / ratio[i] << ",";
			else if( i != 0 && int( (i + 1) / x_steps) != y_point)
				cout << zamp[i] / ratio[i] << "," << endl;
		}
	}
}
void angleValWrite_callback( int val)
{
	cout << "***" << endl;
	cout << "***" << endl;
	cout << "***" << endl;
	cout << "***" << endl;
	cout << "***" << endl;
	cout << "dist, angle, strength " << endl;
	//int count_strength[10] = {0};
	for(int i = 0; i < filenum; i++)
	{
		int y_point = int( i / x_steps);
		int x_point = i % x_steps;
		//double z_point = makeDistance_waveraise(i);

		if( x_point > x_begin
			&& x_point <= x_end
			&& y_point > y_begin
			&& y_point <= y_end
			&& zpoint[i] > z_front
			&& zpoint[i] <= z_back
			)
		{
			cout << zpoint[i] << "," << angle_final[i] << "," << zamp[i] << endl;
		}
	}
}


int main(int argc, char *argv[])
{
	getFileList();

	//3dview---------------------------------------
	makeArray();




	//initiation GLUT
	glutInit(&argc, argv);//initialize OpenGL environment
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );//double buffering

	//funcs for main_window
	Init();
	WinID[WindowNum] = glutCreateWindow(WindowName[WindowNum]);
	callBacks();
	WinFlag[WindowNum] = 1;
	WindowNum = WindowNum + 1;

	//funcs for sub_window
	Init2();
	WinID[WindowNum] = glutCreateWindow(WindowName[WindowNum]);
	callBacks2();
	WinFlag[WindowNum] = 1;
	WindowNum=WindowNum+1;

	//funcs for sub2_window
	Init3();
	WinID[WindowNum] = glutCreateWindow(WindowName[WindowNum]);
	callBacks3();
	WinFlag[WindowNum] = 1;
	WindowNum=WindowNum+1;

	glutIdleFunc(idle);

	//GLUI window and add controls

	//controls for 3D view
	GLUI *glui = GLUI_Master.create_glui("control", 0);
	GLUI *glui2 = GLUI_Master.create_glui("control_subwindow", 0);
	GLUI *glui3 = GLUI_Master.create_glui("control_color", 0);
	GLUI_Rotation *view_rot = glui->add_rotation("Rotation", rotary);
		glui->add_button("rot_reset", 0, funcReset);
	GLUI_EditText *segment_edittext_rot_deg =
			glui->add_edittext( "rot_deg", GLUI_EDITTEXT_INT, &rot_deg);
	GLUI_EditText *segment_edittext_rot_axs_x =
				glui->add_edittext( "rot_axs_x", GLUI_EDITTEXT_FLOAT, &rot_axs_x);
	GLUI_EditText *segment_edittext_rot_axs_y =
				glui->add_edittext( "rot_axs_y", GLUI_EDITTEXT_FLOAT, &rot_axs_y);
	GLUI_EditText *segment_edittext_rot_axs_z =
				glui->add_edittext( "rot_axs_z", GLUI_EDITTEXT_FLOAT, &rot_axs_z);

	GLUI_Translation *translation_xy = glui->add_translation( "TranslationXY", GLUI_TRANSLATION_XY, trans_ary);
	GLUI_Translation *translation_z = glui->add_translation( "TranslationZ", GLUI_TRANSLATION_Z, trans_aryZ);
	glui->add_button("trans_reset", 0, funcResetTrans);


	GLUI_Spinner *segment_spinner =
			glui->add_spinner( "threshold", GLUI_SPINNER_FLOAT, &spinnerValf);
	segment_spinner->set_float_limits(0.0, 1.0, GLUI_LIMIT_CLAMP);
	segment_spinner->set_float_val(0.006);

	GLUI_EditText *segment_edittext_xmin =
			glui->add_edittext( "x_min", GLUI_EDITTEXT_INT, &x_begin);
	segment_edittext_xmin->set_int_limits( 0, x_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_xmin->set_int_val(0);
	GLUI_EditText *segment_edittext_xmax =
			glui->add_edittext( "x_max", GLUI_EDITTEXT_INT, &x_end);
	segment_edittext_xmax->set_int_limits( 0, x_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_xmax->set_int_val(x_steps);
	GLUI_EditText *segment_edittext_ymin =
			glui->add_edittext( "y_min", GLUI_EDITTEXT_INT, &y_begin);
	segment_edittext_ymin->set_int_limits( 0, y_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_ymin->set_int_val(0);
	GLUI_EditText *segment_edittext_ymax =
			glui->add_edittext( "y_max", GLUI_EDITTEXT_INT, &y_end);
	segment_edittext_ymax->set_int_limits( 0, y_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_ymax->set_int_val(y_steps);
	GLUI_EditText *segment_edittext_zmin =
			glui->add_edittext( "z_min", GLUI_EDITTEXT_FLOAT, &z_front);
	segment_edittext_zmin->set_float_limits( 0.0, 100.0, GLUI_LIMIT_CLAMP);
	segment_edittext_zmin->set_float_val(10.0);
	//z_begin = segment_edittext_zmin->get_float_val();
	GLUI_EditText *segment_edittext_zmax =
			glui->add_edittext( "z_max", GLUI_EDITTEXT_FLOAT, &z_back);
	segment_edittext_zmax->set_float_limits( 0.0, 100.0, GLUI_LIMIT_CLAMP);
	segment_edittext_zmax->set_float_val(35);
	GLUI_EditText *segment_edittext_brightness =
				glui3->add_edittext( "brightness", GLUI_EDITTEXT_FLOAT, &brightness);
		segment_edittext_brightness->set_float_limits( 0.0, 10000.0, GLUI_LIMIT_CLAMP);
		segment_edittext_brightness->set_float_val(1.0);
	GLUI_EditText *segment_edittext_geta =
				glui3->add_edittext( "geta", GLUI_EDITTEXT_FLOAT, &geta);
		segment_edittext_geta->set_float_limits( 1, 100.0, GLUI_LIMIT_CLAMP);
		segment_edittext_geta->set_float_val(1.05);
	GLUI_EditText *segment_edittext_point =
				glui->add_edittext( "point size", GLUI_EDITTEXT_INT, &pointsize);
		segment_edittext_point->set_int_limits( 0, 10000, GLUI_LIMIT_CLAMP);
		segment_edittext_point->set_int_val(1);
	GLUI_Spinner *segment_spinner_surface =
			glui->add_spinner( "surface_extraction", GLUI_SPINNER_FLOAT, &spinnerSurfVal, radiobuttonVal, viewchange_callback);
		segment_spinner_surface->set_float_limits(0.0, 10000.0, GLUI_LIMIT_CLAMP);
		segment_spinner_surface->set_float_val(6.0);
	GLUI_EditText *segment_edittext_mountwidth =
			glui->add_edittext( "**useless**mount w threshold", GLUI_EDITTEXT_INT, &mountainWidth, radiobuttonVal, viewchange_callback);
		segment_edittext_mountwidth->set_int_limits(0, 10000, GLUI_LIMIT_CLAMP);
		segment_edittext_mountwidth->set_int_val(0);
	GLUI_EditText *segment_edittext_mountPeakNum =
			glui->add_edittext( "**useless**mount peak num threshold", GLUI_EDITTEXT_INT, &mountainPeakNum, radiobuttonVal, viewchange_callback);
		segment_edittext_mountPeakNum->set_int_limits(0, 10000, GLUI_LIMIT_CLAMP);
		segment_edittext_mountPeakNum->set_int_val(0);
	GLUI_EditText *segment_tri_threshold =
				glui->add_edittext( "tri_threshold", GLUI_EDITTEXT_FLOAT, &tri_threshold);
		segment_tri_threshold->set_float_limits( 0.0, 100.0, GLUI_LIMIT_CLAMP);
		segment_tri_threshold->set_float_val(100.0);
	GLUI_EditText *segment_tri_length =
				glui->add_edittext( "tri_string_length", GLUI_EDITTEXT_FLOAT, &tri_length);
		segment_tri_length->set_float_limits( 0.0, 15.0, GLUI_LIMIT_CLAMP);
		segment_tri_length->set_float_val(5.0);
	/*
	//set a number for calc. integration val. of peak and its neighbor
	GLUI_EditText *segment_edittext_integmountnum =
				glui->add_edittext( "mount num integ", GLUI_EDITTEXT_INT, &mountnuminteg, radiobuttonVal, viewchange_callback);
			segment_edittext_integmountnum->set_int_limits(1, 100, GLUI_LIMIT_CLAMP);
			segment_edittext_integmountnum->set_int_val(10);
	*/
	GLUI_EditText *segment_edittext_mount_int_thre =
			glui->add_edittext( "mount interval setting", GLUI_EDITTEXT_FLOAT, &mount_interval_thre, radiobuttonVal, viewchange_callback);
		segment_edittext_mount_int_thre->set_float_limits(0.0, 100.0, GLUI_LIMIT_CLAMP);
		segment_edittext_mount_int_thre->set_float_val(40);
	GLUI_EditText *segment_edittext_surfaceinterval =
				glui->add_edittext( "interval for second surface", GLUI_EDITTEXT_INT, &surfaceinterval, radiobuttonVal, viewchange_callback);
			segment_edittext_surfaceinterval->set_int_limits(0, 5000, GLUI_LIMIT_CLAMP);
			segment_edittext_surfaceinterval->set_int_val(300);
			/*
	GLUI_EditText *segment_edittext_surface_thickness =
					glui->add_edittext( "surface thickness", GLUI_EDITTEXT_INT, &surfthickness, radiobuttonVal, viewchange_callback);
				segment_edittext_surface_thickness->set_int_limits(0, 5000, GLUI_LIMIT_CLAMP);
				segment_edittext_surface_thickness->set_int_val(300);
*/
	GLUI_RadioGroup *radio_button =
			glui->add_radiogroup(&radiobuttonVal, radiobuttonVal, viewchange_callback);
		glui->add_radiobutton_to_group(radio_button, "point view");
		glui->add_radiobutton_to_group(radio_button, "peak view");
		glui->add_radiobutton_to_group(radio_button, "one point view");
		glui->add_radiobutton_to_group(radio_button, "inner view");
		glui->add_radiobutton_to_group(radio_button, "delaunay view");
		glui->add_radiobutton_to_group(radio_button, "delaunay wire");
		glui->add_radiobutton_to_group(radio_button, "mult. delaunay");
		glui->add_radiobutton_to_group(radio_button, "correction ratio view");


		//-------sub window-----------
	glui3->add_separator();
	GLUI_RadioGroup *radio_color =
			glui3->add_radiogroup(&radiocolorVal, radiobuttonVal, viewchange_callback);
		glui3->add_radiobutton_to_group(radio_color, "mono");
		glui3->add_radiobutton_to_group(radio_color, "color");
	GLUI_EditText *segment_alpha =
					glui3->add_edittext( "surface alpha val ", GLUI_EDITTEXT_FLOAT, &surfAlph, radiobuttonVal, viewchange_callback);
			segment_alpha->set_float_limits( 0.0, 1.0, GLUI_LIMIT_CLAMP);
			segment_alpha->set_float_val(0.5);
	GLUI_EditText *segment_surfcorrect =
					glui3->add_edittext( "surface correction", GLUI_EDITTEXT_FLOAT, &correctVal, radiobuttonVal, viewchange_callback);
			segment_surfcorrect->set_float_limits( 0.001, 100.0, GLUI_LIMIT_CLAMP);
			segment_surfcorrect->set_float_val(1.0);
	GLUI_EditText *segment_edittext_ymin2 =
				glui2->add_edittext( "y_min", GLUI_EDITTEXT_INT, &y_begin2);
		segment_edittext_ymin2->set_int_limits( 0, y_steps, GLUI_LIMIT_CLAMP);
		segment_edittext_ymin2->set_int_val(0);
	GLUI_EditText *segment_edittext_ymax2 =
			glui2->add_edittext( "y_max", GLUI_EDITTEXT_INT, &y_end2);
		segment_edittext_ymax2->set_int_limits( 0, y_steps, GLUI_LIMIT_CLAMP);
		segment_edittext_ymax2->set_int_val(y_steps);

	glui2->add_separator();

	//controls for 2D view
	GLUI_EditText *segment_edittext_sub_x =
			glui2->add_edittext( "sub_x", GLUI_EDITTEXT_INT, &sub_x);
	segment_edittext_sub_x->set_int_limits(0.0, x_steps, GLUI_LIMIT_CLAMP);
	segment_edittext_sub_x->set_int_val(30);
	GLUI_EditText *segment_edittext_sub_y =
				glui2->add_edittext( "sub_y", GLUI_EDITTEXT_INT, &sub_y);
		segment_edittext_sub_y->set_int_limits(0.0, y_steps, GLUI_LIMIT_CLAMP);
		segment_edittext_sub_y->set_int_val(20);
	GLUI_Translation *translation_x_sub = glui2->add_translation( "TranslationX", GLUI_TRANSLATION_X, trans_ary_sub);
		translation_x_sub->set_speed( 0.5 );
	GLUI_Translation *translation_z_sub = glui2->add_translation( "TranslationZ", GLUI_TRANSLATION_Z, trans_aryZ_sub);
		translation_z_sub->set_speed( 0.5 );
	GLUI_Translation *translation_x_sub2 = glui2->add_translation( "TranslationX2", GLUI_TRANSLATION_X, trans_ary_sub2);
		translation_x_sub2->set_speed( 0.5 );
	GLUI_Translation *translation_z_sub2 = glui2->add_translation( "TranslationZ2", GLUI_TRANSLATION_Z, trans_aryZ_sub2);
		translation_z_sub2->set_speed( 0.5 );
	glui2->add_button("reset_sub", 0, funcResetSub);

	glui->add_separator();
	//modification panels
	GLUI_EditText *segment_edittext_addpoint =
			glui2->add_edittext( "add peak num for visualizing", GLUI_EDITTEXT_INT, &addPeakNum, 0, modifyPeak_callback);
		segment_edittext_addpoint->set_int_limits(0, 100, GLUI_LIMIT_CLAMP);
		segment_edittext_addpoint->set_int_val(0);


	//glui2->add_button("change", 0, modifyPeak_callback);
	glui2->add_button("setting write", 0, settingWR_callback);
	glui2->add_button("setting read", 1, settingWR_callback);
	glui2->add_separator();
	GLUI_RadioGroup *radio_correction =
		glui2->add_radiogroup(&correctionfuncOFflag, radiobuttonVal, correctionOF_callback);
	glui2->add_radiobutton_to_group(radio_correction, "correction off");
	glui2->add_radiobutton_to_group(radio_correction, "correction on");
	GLUI_EditText *segment_edittext_corrmap =
			glui2->add_edittext( "corr. press", GLUI_EDITTEXT_INT, &corrpress, 0, correctionPress_callback);
		segment_edittext_corrmap->set_int_limits(1, 10000, GLUI_LIMIT_CLAMP);
		segment_edittext_corrmap->set_int_val(1);
	GLUI_EditText *segment_edittext_corrmap_angle =
			glui2->add_edittext( "corr. angle slide", GLUI_EDITTEXT_INT, &corrangle, 0, correctionPress_callback);
		segment_edittext_corrmap_angle->set_int_limits(0, 100, GLUI_LIMIT_CLAMP);
		segment_edittext_corrmap_angle->set_int_val(0);
	GLUI_EditText *segment_edittext_corrmap_widen =
				glui2->add_edittext( "widen angle corr.", GLUI_EDITTEXT_FLOAT, &corrwiden, 0, correctionPress_callback);
			segment_edittext_corrmap_widen->set_float_limits(1.0, 100.0, GLUI_LIMIT_CLAMP);
			segment_edittext_corrmap_widen->set_float_val(1.0);
	glui3->add_button("surface val write", 0, valWrite_callback);
	glui3->add_button("surface angle'n'distance'n'val write", 0, angleValWrite_callback);
	GLUI_RadioGroup *radio_average =
		glui3->add_radiogroup(&mountHeightAverageFlag, radiobuttonVal, viewchange_callback);
	glui3->add_radiobutton_to_group(radio_average, "high average off");
	glui3->add_radiobutton_to_group(radio_average, "high average on");

	glui->add_button("Exit", 0, gluiCallback);


	glutMainLoop();

	return 0;
}


