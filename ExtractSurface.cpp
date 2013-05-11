/*
 * ExtractSurface.cpp
 *
 *  Created on: May 8, 2013
 *      Author: nagaso
 */

#include "ExtractSurface.h"


ExtractSurface::ExtractSurface( int timebottom, double data[data_length])
: sample_end( timebottom )
{
	resetData(data);

	average = calcAverage();
	threshold = average * times;

	findMountain();
	extractPeak();
}

void ExtractSurface::resetData(double data_in[data_length])
{

	for(int i = 0; i < data_length; i++)
		data_sample[i] = data_in[i];
}

void ExtractSurface::findMountain()
{
	int mountain_count = 0;
	bool mountain_inout_flag = 0;
	for(int i = sample_start; i < sample_end; i++)
	{
		if( mountain_inout_flag == 0 && data_sample[i] > threshold)
		{
			mountain_inout_flag = 1;
			mountain_place[mountain_count][0] = i;
		}
		else if(mountain_inout_flag == 1 && data_sample[i] < threshold)
		{
			mountain_inout_flag = 0;
			mountain_place[mountain_count][1] = i;
			mountain_count++;
		}
	}

	mountain_num = mountain_count;
}

void ExtractSurface::extractPeak()
{
	bool flag_overthreshold = 0;
	int peak_count = 0;

	for(int i = 0; i < mountain_num; i++)
		for(int j = mountain_place[i][0]; j < mountain_place[i][1]; j++ )
		{
			double difference_forward = ( data_sample[j] - data_sample[j - 1]);
			double difference_behind = ( data_sample[j + 1] - data_sample[j]);

			if( difference_forward > 0 && difference_behind < 0)
			{
				peak_place[peak_count] = j;
				peak_count++;
			}

		}

	peak_num = peak_count;
}

double ExtractSurface::calcAverage()
{
	double temp_sum = 0;

	for( int i = offset_sample_start; i < offset_sample_end ; i++)
			temp_sum += data_sample[i];

	return temp_sum / (offset_sample_end - offset_sample_start);
}

