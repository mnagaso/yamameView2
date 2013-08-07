/*
 * CorrFunc.cpp
 *
 *  Created on: Jul 16, 2013
 *      Author: nagaso
 */

#include <stdio.h>
#include<math.h>
#include <iostream>
using std::ios;
using std::cin;
using std::cout;
using std::endl;

#include <fstream>
using std::fstream;

#include <cstdlib>

#include <string.h>
using std::string;

#include <sstream>
#include <stdio.h>
#include <glob.h>


#include "CorrFunc.h"
#define filepath2 "/home/nagaso/workspace/yamameView2/"

CorrFunc::CorrFunc(double ang, double dist, double ref, int sl, double wt)
: angle( ang ), distance( dist ), slide_amount( sl ), widen_times( wt )
{
	//cout << dist << "," << ang << "," << ref << endl;
	initialize();
	readMap();
	widenMap();
	//searchNearPoint();
	searchNearPointandLinearInterpolation();
}

void CorrFunc::searchNearPoint()
{
	double minanglediff = 999;
	double tempnearestangle = 0;
	double second_nearestangle;
	angle -= slide_amount;

	for(int i = 0; i < datapoints; i++)
	{

		if(angle < 0)
			angle = fabs(angle);
		if(minanglediff > fabs(angle - map[i][1]))
		{
			minanglediff = fabs(angle - map[i][1]);
			tempnearestangle = map[i][1];
		}
	}

	double mindistdiff = 999;
	int nearest_filenum;

	for(int i = 0; i < datapoints; i++){
		if( map[i][1] == tempnearestangle)
		{
			double tempdiff = fabs(distance - map[i][0]);
			if(tempdiff < mindistdiff)
			{
				nearest_filenum = i;
				//tempnearestdist = map[i][0];
				mindistdiff = tempdiff;

			}
		}

	}
	ratio = map[nearest_filenum][2];
	refVal_rot = map[nearest_filenum][3];

}

void CorrFunc::searchNearPointandLinearInterpolation()
{
	double minanglediff = 999;
	int tempanglepoint;
	double minanglediff2 = 999;
	int tempanglepoint2;
	angle -= slide_amount;

	if(angle < 0)
		angle = fabs(angle);

	for( int i = 0; i < datapoints; i++ )
	{
		//cout << map[i][0] << ", " << map[i][1] << ", " << map[i][2] << ", " << map[i][3] << endl;
		if(fabs(angle - map[i][1]) < minanglediff)
		{
			minanglediff = fabs(angle - map[i][1]);
			tempanglepoint = i;
		}
	}

	for(int i = 0; i < datapoints; i++)
	{
		if(fabs(angle - map[i][1]) < minanglediff2 && map[i][1] != map[tempanglepoint][1])
		{
			minanglediff2 = fabs(angle - map[i][1]);
			tempanglepoint2 = i;
		}
	}

	if(map[tempanglepoint][1] > map[tempanglepoint2][1])
	{
		int temp = tempanglepoint;
		tempanglepoint = tempanglepoint2;
		tempanglepoint2 = temp;
	}

	int tempdistpoint11 = 0;
	int tempdistpoint12 = 0;
	int tempdistpoint21 = 0;
	int tempdistpoint22 = 0;

	double mindistdiff11 = 999;
	double mindistdiff12 = 999;
	double mindistdiff21 = 999;
	double mindistdiff22 = 999;

	for( int i = 0; i < datapoints; i++ )
	{
		if( map[i][1] == map[tempanglepoint][1] )
		{
			if( fabs( distance - map[i][0] ) <= mindistdiff11 )
			{
				mindistdiff12 = mindistdiff11;
				tempdistpoint12 = tempdistpoint11;
				mindistdiff11 = fabs( distance - map[i][0]);
				tempdistpoint11 = i;

			} else if( fabs( distance - map[i][0] ) <= mindistdiff12 )
			{
				mindistdiff12 = fabs( distance - map[i][0] );
				tempdistpoint12 = i;
			}
		} else if( map[i][1] == map[tempanglepoint2][1] )
		{
			if( fabs( distance - map[i][0] ) <= mindistdiff21 )
			{
				mindistdiff22 = mindistdiff21;
				tempdistpoint22 = tempdistpoint21;
				mindistdiff21 = fabs( distance - map[i][0]);
				tempdistpoint21 = i;
			} else if( fabs( distance - map[i][0] ) <= mindistdiff22 )
			{
				mindistdiff22 = fabs( distance - map[i][0] );
				tempdistpoint22 = i;
			}
		}
	}

	if(map[tempdistpoint11][0] > map[tempdistpoint12][0])
	{
		int temp = tempdistpoint11;
		tempdistpoint11 = tempdistpoint12;
		tempdistpoint12 = temp;
	}
	if(map[tempdistpoint21][0] > map[tempdistpoint22][0])
	{
		int temp = tempdistpoint21;
		tempdistpoint21 = tempdistpoint22;
		tempdistpoint22 = temp;
	}


	double ratioangle1;
	double ratioangle2;
	if(map[tempdistpoint12][0] < distance)
	{
		ratioangle1 = map[tempdistpoint12][2];
	}
	if(map[tempdistpoint22][0] < distance)
	{
		ratioangle2 = map[tempdistpoint22][2];
	}
	if(map[tempdistpoint11][0] > distance)
	{
		ratioangle1 = map[tempdistpoint11][2];
	}
	if(map[tempdistpoint21][0] > distance)
	{
		ratioangle2 = map[tempdistpoint21][2];
	}
	if(map[tempdistpoint11][0] < distance && map[tempdistpoint12][0] > distance && map[tempdistpoint12][0] != map[tempdistpoint11][0])
	{
		ratioangle1 = map[tempdistpoint11][2] + (map[tempdistpoint12][2] - map[tempdistpoint11][2]) * (distance - map[tempdistpoint11][0]) / (map[tempdistpoint12][0] - map[tempdistpoint11][0]);
	}
	if(map[tempdistpoint12][0] == map[tempdistpoint11][0])
	{
		ratioangle1 = map[tempdistpoint11][2];
	}
	if(map[tempdistpoint21][0] < distance && map[tempdistpoint22][0] > distance && map[tempdistpoint22][0] != map[tempdistpoint21][0])
	{
		ratioangle2 = map[tempdistpoint21][2] + (map[tempdistpoint22][2] - map[tempdistpoint21][2]) * (distance - map[tempdistpoint21][0]) / (map[tempdistpoint22][0] - map[tempdistpoint21][0]);
	}
	if(map[tempdistpoint22][0] == map[tempdistpoint21][0])
	{
		ratioangle2 = map[tempdistpoint21][2];
	}

	if(map[tempanglepoint][1] != map[tempanglepoint2][1])
		ratio = ratioangle1 + (ratioangle2 - ratioangle1) * (angle - map[tempanglepoint][1]) / (map[tempanglepoint2][1] - map[tempanglepoint][1]);
	else
		ratio = ratioangle1;
	//cout << "ratio:  " << ratio << endl;
}

void CorrFunc::initialize()
{
	for(int i = 0; i < datapoints; i++)
		for(int j = 0; j < 3; j++)
			map[i][j] = 0.00000000000000000;
}

void CorrFunc::widenMap()
{
	for(int i = 0; i < datapoints; i++)
	{
		map[i][1] = map[i][1] * widen_times;
	}

}

bool CorrFunc::readMap()
{
	FILE *fp;
	int separator = ',';
	int len, i, j, k;
	string fn = filepath2;
	fn += filename;
	const char *fname = fn.c_str();
	char s[100],word[100], *p, *pos;
	bool file_flag;

	i = 0;
	fp = fopen( fname, "r");
	if( fp == NULL ) {
		cout << "no file existing..." << endl;
		return -1;
	}

	while( fgets( s, 100, fp ) != NULL){
		len = strlen( s );
		if( len > 0 && s[ len - 1 ] == '\n')
			s[ len - 1 ] = '\0';
		j = 0, p = s, k = 0;
		file_flag = 0;
		while ( *p != '\0' ){
			pos = strchr( p, separator );
			if( pos != NULL ){
				strncpy( word, p, pos - p );
				word[pos - p] = '\0';
				p = ++pos;
			}
			else{
				strcpy( word, p);
				*p = '\0';
			}

			map[i][j] = atof( word );
			j++;
		}
		//cout << "test read file... " <<map[i][0] << "  " << map[i][1] << endl;

		i++;
	}
	fclose( fp );
}
