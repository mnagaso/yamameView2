/*
 * harmonic.h
 *
 *  Created on: Sep 10, 2013
 *      Author: nagaso
 */

#ifndef HARMONIC_H_
#define HARMONIC_H_

#include <fftw3.h>

class Harmonic
{
public:
	const static int data_width = 2;
	const static int data_length = 10000;
	//const static int file_num = 4192;
	//const static int offset_sample_start = 5500;
	//const static int offset_sample_end = 6500;
	//const static int offset_sample_start = 9000;//for rotate
	//const static int offset_sample_end = 10000;
	double data[data_length];
	//double data[data_length][data_width];
	double data2[data_length][data_width];
	double average;
	//#define  extern fishname //from iwashiView.cpp

	Harmonic( double data_in[data_length]);
	void readData( double dt_in[data_length]);
	double calcAverage();
	void offset();
	void FFT( int );
	void windowFunc();
	void getEnvelope();
	void resetData();

private:

};


#endif /* HARMONIC_H_ */
