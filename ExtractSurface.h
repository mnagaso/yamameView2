/*
 * ExtractSurface.h
 *
 *  Created on: May 8, 2013
 *      Author: nagaso
 */

#ifndef EXTRACTSURFACE_H_
#define EXTRACTSURFACE_H_

class ExtractSurface
{
public:
	const static int data_length = 10000;
	const static int sample_start = 1138;
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
	double times = 0; //multiple num of noise lever(average)for make threshold

	ExtractSurface(int timebottom, double data[data_length], double mult_times, int mountainWidth, int peakNum);
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


#endif /* EXTRACTSURFACE_H_ */
