/*
 * Utils.cpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */

#include "Utils.hpp"

using namespace std;

//summation over data
double sum(double * data, int size)
{
	double tmp = 0;
	for (int i = 0; i < size; ++i)
	{
		tmp += data[i];
	}
	return tmp;
}

inline double max(double a, double b)
{
	return a > b ? a : b;
}

