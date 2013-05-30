/*
 * OffsetHilbert.h
 *
 *  Created on: 2013/03/02
 *      Author: mnsaru
 */

#ifndef OFFSETHILBERT_H_
#define OFFSETHILBERT_H_


#include <fftw3.h>

class OffsetHilbert
{
public:
	const static int data_width = 2;
	const static int data_length = 10000;
	//const static int file_num = 4192;
	const static int offset_sample_start = 5500;
	const static int offset_sample_end = 6500;
	double data[data_length];
	//double data[data_length][data_width];
	double data2[data_length][data_width];
	double average;

	OffsetHilbert( double data_in[data_length]);
	void readData( double dt_in[data_length]);
	double calcAverage();
	void offset();
	void FFT( int );
	void windowFunc();
	void getEnvelope();
	void resetData();

private:

};


#endif /* OFFSETHILBERT_H_ */
