/*
 * ExtractInnerSurface.h
 *
 *  Created on: May 27, 2013
 *      Author: nagaso
 */

#ifndef EXTRACTINNERSURFACE_H_
#define EXTRACTINNERSURFACE_H_

class ExtractInnerSurface
{
public:
	const static int data_length = 10000;
	int sample_start = 0;
	double data_sample[data_length];
	int mountain_place[10000][2];//[mountain_num][in or out]
	int point_status[10000];
	int mountain_num = 0;
	int peak_num = 0;
	int peak_place[10000];
	int sample_end = 0;
	double threshold = 0.0;
	const static int offset_sample_start = 5500;
	const static int offset_sample_end = 6500;
	int times = 0; //multiple num of noise lever(average)for make threshold

	ExtractInnerSurface(int time_start, int timebottom, double data[data_length], int mult_times, int mountainWidth, int peakNum);
	void resetData(double data_in[data_length]);
	void findMountain();
	void mountainWidthThreshold();
	void mountainNumbering();
	void extractPeak();
	void eraceNoisePeak();
	double calcAverage();
	void init();
private:
	double average;
	int mountainWidthThre = 0;
	int peakNum = 0;
};



#endif /* EXTRACTINNERSURFACE_H_ */
