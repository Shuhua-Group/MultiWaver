/*
 * ParamExp.cpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */
#include <cmath>
#include <cstdlib>
#include <map>
#include <iostream>

#include "ParamExp.hpp"

using namespace std;

/* constructor
 * Initial with K exponential distributions
 * parameters lambda are randomly initial with value between 0 and 1
 * proportions are randomly initial with value between 0 and 1, and
 * the summation equals 1
 */
ParamExp::ParamExp(int K) :
		K(K)
{
	lambda = new double[K];
	for (int i = 0; i < K; ++i)
	{
		lambda[i] = 1.0 * rand() / RAND_MAX;
	}
	prop = new double[K];
	double tmp = 0;
	for (int i = 0; i < (K - 1); ++i)
	{
		prop[i] = rand() / (1.0 * K * RAND_MAX);
		tmp += prop[i];
	}
	//ensure sum to one
	prop[K - 1] = 1 - tmp;
}

/* another constructor
 * initialize with K, lambdas and proportions
 */
ParamExp::ParamExp(int K, double *l, double *p) :
		K(K)
{
	lambda = new double[K];
	prop = new double[K];
	for (int i = 0; i < K; ++i)
	{
		lambda[i] = l[i];
		prop[i] = p[i];
	}
}

//copy constructor
ParamExp::ParamExp(const ParamExp &rhs)
{
	K = rhs.K;
	lambda = new double[K];
	prop = new double[K];
	for (int i = 0; i < K; ++i)
	{
		lambda[i] = rhs.lambda[i];
		prop[i] = rhs.prop[i];
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
		for (int i = 0; i < K; ++i)
		{
			lambda[i] = rhs.lambda[i];
			prop[i] = rhs.prop[i];
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
bool ParamExp::isConverge(const ParamExp & par, double epsilon)
{
	bool converge = 1;
	for (int i = 0; i < K; ++i)
	{
		if (abs(lambda[i] - par.getLambda(i)) > epsilon || abs(prop[i] - par.getProp(i)) > epsilon)
		{
			converge = 0;
			break;
		}
	}
	return converge;
}

void ParamExp::sortByLambda()
{
	/*
	 * the key of map are automatically sorted, therefore can be used
	 * to sort lambda and corresponding proportion accordingly
	 */
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
	for (int i = 0; i < K; ++i)
	{
		cout << lambda[i] << ", " << prop[i] << "; ";
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
