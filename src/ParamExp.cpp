/*
 * ParamExp.cpp
 *
 *  Created on: May 26, 2015
 *      Author: young
 */

#include <ctime>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include "ParamExp.hpp"

using namespace std;

//constructor
ParamExp::ParamExp(int K) :
		K(K)
{
	lambda = new double[K];
	for (int j = 0; j < K; ++j)
	{
		lambda[j] = 1.0 * rand() / RAND_MAX;
	}
	prop = new double[K];
	double tmp = 0;
	for (int j = 0; j < (K - 1); ++j)
	{
		prop[j] = rand() / (1.0 * K * RAND_MAX);
		tmp += prop[j];
	}
	prop[K - 1] = 1 - tmp;
}

//another constructor
ParamExp::ParamExp(int K, double *l, double *p) :
		K(K)
{
	lambda = new double[K];
	prop = new double[K];
	for (int j = 0; j < K; ++j)
	{
		lambda[j] = l[j];
		prop[j] = p[j];
	}
}

//copy constructor
ParamExp::ParamExp(const ParamExp &rhs)
{
	K = rhs.K;
	lambda = new double[K];
	prop = new double[K];
	for (int j = 0; j < K; ++j)
	{
		lambda[j] = rhs.lambda[j];
		prop[j] = rhs.prop[j];
	}
}

//overloading operator assignment
ParamExp & ParamExp::operator=(const ParamExp &rhs)
{
	if (this != &rhs)
	{
		K = rhs.K;
		if (lambda != NULL)
			delete[] lambda;
		if (prop != NULL)
			delete[] prop;
		lambda = new double[K];
		prop = new double[K];
		for (int j = 0; j < K; ++j)
		{
			lambda[j] = rhs.lambda[j];
			prop[j] = rhs.prop[j];
		}
	}
	return *this;
}

int ParamExp::getK() const
{
	return K;
}

double ParamExp::getLambda(int index) const
{
	return lambda[index];
}

double ParamExp::getProp(int index) const
{
	return prop[index];
}

//test convergence
bool ParamExp::isConverge(const ParamExp & par)
{
	bool converge = 1;
	for (int i = 0; i < K; ++i)
	{
		if (abs(lambda[i] - par.getLambda(i)) > kDelta || abs(prop[i] - par.getProp(i)) > kDelta)
		{
			converge = 0;
			break;
		}
	}
	return converge;
}

void ParamExp::sortByLambda()
{
	map<double, double> temp;
	int i;
	for (i = 0; i < K; ++i)
	{
		temp[lambda[i]] = prop[i];
	}
	i = 0;
	for (map<double, double>::iterator it = temp.begin(); it != temp.end(); it++)
	{
		lambda[i] = it->first;
		prop[i] = it->second;
		i++;
	}
}

void ParamExp::print()
{
	cout << "par=(";
	for (int j = 0; j < K; ++j)
	{
		cout << lambda[j] << ", " << prop[j] << "; ";
	}
	cout << ")" << endl;
}

//
ParamExp::~ParamExp()
{
	if (lambda != NULL)
		delete[] lambda;
	if (prop != NULL)
		delete[] prop;
}
