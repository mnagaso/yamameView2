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
	int mountain_place[100][2];//[mountain_num][in or out]
	int mountain_num;
	int peak_num;
	int peak_place[100];
	int sample_end;
	double threshold;
	const static int offset_sample_start = 5500;
	const static int offset_sample_end = 6500;
	int times = 3; //multiple num of noise lever(average)for make threshold

	ExtractSurface(int timebottom, double data[data_length]);
	void resetData(double data_in[data_length]);
	void findMountain();
	void extractPeak();
	double calcAverage();
private:
	double average;
};


#endif /* EXTRACTSURFACE_H_ */
