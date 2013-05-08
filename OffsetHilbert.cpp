/*
 * OffsetHilbert.cpp
 *
 *  Created on: 2013/03/02
 *      Author: mnsaru
 */

#include <iostream>
using std::ios;
using std::cin;
using std::cout;
using std::endl;

#include <fftw3.h>
#include <string>
using std::string;
#include <math.h>
#include "OffsetHilbert.h"

OffsetHilbert::OffsetHilbert(double data_in[data_length] )
{
	resetData();
	//cout << data2[100][0] << endl;

	readData( data_in);
	average = calcAverage();
	//offset();

	int flag = -1;
	FFT(flag);

	windowFunc();
	flag = 1;
	FFT(flag);

	getEnvelope();


	//resetData();
}

void OffsetHilbert::readData(double dt_in[data_length])
{
	for(int i = 0; i < data_length; i++)
		data[i] = dt_in[i];
}

void OffsetHilbert::resetData()
{
	for(int i = 0; i < data_length; i++)
		for(int j = 0; j < data_width; j++)
		data2[i][j] = 0.0;

}

double OffsetHilbert::calcAverage()
{
	double temp_sum = 0;

	for( int i = offset_sample_start; i < offset_sample_end ; i++)
			temp_sum += data[i];

	return temp_sum / (offset_sample_end - offset_sample_start);
}

void OffsetHilbert::offset()
{
	if( average >= 0)
		for( int i = 0; i < data_length; i++)
			data[i] -= average;
	else
		for( int i = 0; i < data_length; i++)
			data[i] += average;
}

void OffsetHilbert::FFT( int flag)
{
	fftw_complex *in, *out;
	fftw_plan p;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data_length);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data_length);

	if( flag == -1)
		for( int i = 0; i < data_length; i++){
			in[i][0] = data[i];
			in[i][1] = 0;
		}
	else if( flag == 1)
		for( int i = 0; i < data_length; i++){
			in[i][0] = data2[i][0];
			in[i][1] = data2[i][1];
		}


	p = fftw_plan_dft_1d(data_length, in, out, flag, FFTW_ESTIMATE); //flag == -1 FFTW_FORWARD  flag = 1 FFTW_BACK

	fftw_execute(p);
	//windowFunc(out);

	for( int i = 0; i < data_length; i++){
		data2[i][0] = out[i][0];
		data2[i][1] = out[i][1];
	}

	fftw_destroy_plan(p);
	fftw_free(in);fftw_free(out);

}

void OffsetHilbert::windowFunc()
{
	/*
	double buffer;
	for(int i = 0; i < data_length; i++){
		buffer = data2[i][0];
		data2[i][0] = data2[i][1];
		data2[i][1] = buffer;
	}
	*/

	int midpoint = data_length / 2;
	for( int i = 0; i< data_length; i++)
	{
		if(i < midpoint){
			data2[i][0] *= 2;
			data2[i][1] *= 2;
		}
		else if( i == midpoint){
			data2[i][0] *= 1;
			data2[i][1] *= 1;
		}
		else{//how should if data_length has odd num
			data2[i][0] *= 0;
			data2[i][1] *= 0;
		}
	}
}

void OffsetHilbert::getEnvelope()
{
	for(int i = 0; i < data_length; i++)
		data2[i][1] = sqrt( data2[i][0] * data2[i][0] + data2[i][1] * data2[i][1]) / data_length; //normalize
		//data[i][1] = sqrt( data[i][1] * data[i][1] + data2[i][1] * data2[i][1]);
		//data2[i][1] = data[i];
	//cout << data2[1][1] << endl;

}
