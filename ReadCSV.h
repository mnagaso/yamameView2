/*
 * ReadCSV.h
 *
 *  Created on: 2013/02/14
 *      Author: mnsaru
 */

#ifndef READCSV_H_
#define READCSV_H_

#include <string>
using std::string;

#include <fstream>
using std::fstream;

class ReadCSV
{
public:
	const static int data_width = 2;
	const static int data_length = 10000;
	//const static int file_num = 4192;
	//const static int fake_file_num = 1122;
	double data[data_length][data_width];

	ReadCSV( const string&, int, int, bool);

	bool fileOpenfgets();

private:
	int flag_headerchange;
	string filehead;
	string filename;
	fstream file;

};


#endif /* READCSV_H_ */
