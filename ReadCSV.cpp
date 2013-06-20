/*
 * ReadCSV.cpp
 *
 *  Created on: 2013/02/15
 *      Author: mnsaru
 */

#include <stdio.h>
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

#include <stdio.h>

#include "ReadCSV.h"

//#define headerchangenum 55 //the x axe number one before header name changes
#define eleNum 10000 //numbers of value in each file
#define ignore 2000 //points ignored
#define sonicvelo 1500 //onsoku m/s
//#define filepath "/var/run/media/nagaso/sotoHD/iwana/2/" //"./data/"
//#define filepath "/var/run/media/nagaso/sotoHD/yamame/1/"
#define filepath "/var/run/media/nagaso/sotoHD/rotate2/"
//#define filepath "/var/run/media/nagaso/sotoHD/rotato/"
//#define filepath "./data/"

ReadCSV::ReadCSV( const string &fn, int xstep, int ystep, bool flag)
	: flag_headerchange( flag )
{
	filename = fn;
	//fileOpen();
	//setData();
	fileOpenfgets();
	//makeDistance();
	//writeData2();

}

bool ReadCSV::fileOpenfgets()
{
	FILE *fp;
	int separator = ',';
	int len, i, j;
	string fn = filepath;
	fn += filename;
	const char *fname = fn.c_str();
	char s[100],word[100], *p, *pos;

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
		j = 0, p = s;
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
			data[i][j] = atof( word );
			j++;
		}
		//cout << "test read file... " << data[i][0] << "  " << data[i][1] << endl;

		i++;
	}
	fclose( fp );
}


