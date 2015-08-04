/*
 * EMExp.cpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */

#include <cmath>
#include <iostream>
#include "EMExp.hpp"

double sum(double * data, int size)
{
	double tmp = 0;
	for (int i = 0; i < size; ++i)
	{
		tmp += data[i];
	}
	return tmp;
}

using namespace std;

//constructor
EMExp::EMExp(const ParamExp &par, const vector<double> &observ) :
		par(par), observ(observ)
{
	updateLik();
}

double EMExp::getLik() const
{
	return lik;
}

ParamExp & EMExp::getPar()
{
	return par;
}

void EMExp::setPar(const ParamExp &par)
{
	this->par = par;
}
//update likelihood
/*
 * Log-lik = sum(log(sum(m_j*l_j*exp(-l_j*x_i))))
 */
void EMExp::updateLik()
{
	double tmp = 0;
	for (size_t i = 0; i < observ.size(); ++i)
	{
		int k = par.getK();
		double fval = 0;
		for (int j = 0; j < k; ++j)
		{
			double m = par.getProp(j);
			double l = par.getLambda(j);
			fval += m * l * exp(-l * observ.at(i));
		}
		tmp += log(fval);
	}
	lik = tmp;
}

//EM iteration, converge or reach max iteration time will terminate
void EMExp::iterate(int maxIter, double epsilon)
{
	int it = 0; //iteration number
	int kval = par.getK();
	int size = observ.size();
	double *nlambda; //new lambda
	double *nprop; //new proportion
	double **pval;
	double **pxval;
	nlambda = new double[kval];
	nprop = new double[kval];
	pval = new double *[kval];
	pxval = new double *[kval];
	for (int i = 0; i < kval; ++i)
	{
		pval[i] = new double[size];
		pxval[i] = new double[size];
	}
	while (++it < maxIter)
	{
		//E-step
		for (int i = 0; i < size; ++i)
		{
			double denormator = 0;
			for (int j = 0; j < kval; ++j)
			{
				double mj = par.getProp(j);
				double lj = par.getLambda(j);
				pval[j][i] = mj * lj * exp(-lj * observ.at(i));
				denormator += pval[j][i];
			}
			for (int j = 0; j < kval; ++j)
			{
				pval[j][i] /= denormator;
				pxval[j][i] = pval[j][i] * observ.at(i);
			}
		}
		//M-step
		for (int i = 0; i < kval; ++i)
		{
			double sump = sum(pval[i], size);
			nprop[i] = sump / size;
			nlambda[i] = sump / sum(pxval[i], size);
		}
		ParamExp updatedPar(kval, nlambda, nprop);
		/*
		 * check converge or not, if converge, jump out of loop
		 */
		if (par.isConverge(updatedPar, epsilon))
		{
			par = updatedPar;
			updateLik();
			break;
		}
		else
		{
			par = updatedPar;
			updateLik();
		}
		//cout << "Iteration " << ++it << " --> llk: " << getLik() << "; ";
		//getPar().print();
	}
	if (it >= maxIter)
	{
		cerr << "Warning: Max interation reached before convergence"
	}
	//clean stuff
	delete[] nlambda;
	delete[] nprop;
	for (int i = 0; i < kval; ++i)
	{
		delete[] pval[i];
		delete[] pxval[i];
	}
	delete[] pval;
	delete[] pxval;
}

EMExp::~EMExp()
{
}

