/*
 * ExtractSurface.cpp
 *
 *  Created on: May 8, 2013
 *      Author: nagaso
 */
#include <iostream>
using std::ios;
using std::cin;
using std::cout;
using std::endl;
#include "ExtractSurface.h"


ExtractSurface::ExtractSurface( int timebottom, double data[data_length], int mult_times)
: sample_end( timebottom )
{
	resetData(data);

	average = calcAverage();
	times = mult_times;
	threshold = average * times;

	init();
	findMountain();
	extractPeak();
	//cout << "";
}

void ExtractSurface::resetData(double data_in[data_length])
{

	for(int i = 0; i < data_length; i++)
		data_sample[i] = data_in[i];
}

void ExtractSurface::findMountain()
{
	for(int i = sample_start; i < sample_end; i++)
	{
		if( data_sample[i] >= threshold)
			point_status[i] = 1;
		else
			point_status[i] = 0;
	}

	int mountain_count = 0;
	for( int i = sample_start; i < sample_end; i++)
	{
		if(point_status[i] == 1 && point_status[i - 1] == 0)
			mountain_place[mountain_count][0] = i;
		else if( point_status[i] == 0 && point_status[i - 1] == 1)
		{
			mountain_place[mountain_count][1] = i;
			mountain_count++;
		}
		else if(i + 1 == sample_end && point_status[i] == 1)
		{
			mountain_place[mountain_count][1] = i;
			mountain_count++;
		}
	}
	mountain_num = mountain_count;
}

void ExtractSurface::extractPeak()
{
	int peak_count = 0;
	if( mountain_num != 0)
		for( int i = 0; i < mountain_num; i++)
			for( int j = mountain_place[i][0]; j < mountain_place[i][1]; j++)
			{
				if(data_sample[j] - data_sample[j - 1] >= 0 && data_sample[j + 1] - data_sample[j] < 0)
				{
					peak_place[peak_count] = j;
					point_status[j] = 2;
					peak_count++;
				}
			}

	peak_num = peak_count;
}

double ExtractSurface::calcAverage()
{
	//double temp_sum = 0;
	int calc_width = 1000;
	int calc_times = (int)data_length/calc_width;
	double temp_average[ calc_times];
	double temp_amount;
	int temp_calcstart = 0;

	for( int i = 0; i < calc_times; i++)
	{
		for(int j = temp_calcstart; j < temp_calcstart + calc_width; j++)
			temp_amount += data_sample[j];

		temp_average[i] = temp_amount / calc_width;
		temp_amount = 0;
		temp_calcstart += calc_width;
 	}

	double temp_small = 10000000;
	for( int i = 0; i < calc_times; i++)
		if( temp_average[i] < temp_small)
			temp_small = temp_average[i];

	return temp_small;
}

void ExtractSurface::init()
{
	for(int i = 0; i < 10000; i++){
		point_status[i] = 0;
		peak_place[i] = 0;
		for(int j = 0; j < 2; j++)
			mountain_place[i][j] = 0;
	}
}
