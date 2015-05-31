/*
 * EMExp.cpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */

#include <cmath>
#include "EMExp.hpp"
#include "Utils.hpp"

using namespace std;

//constructor
EMExp::EMExp(const ParamExp &par, const std::vector<double> &observ) :
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
void EMExp::iterate(int maxIter)
{
	int it = 0;
	int kval = par.getK();
	int size = observ.size();
	double *nlambda;
	double *nprop;
	double **pval;
	double **pxval;
	nlambda = new double[kval];
	nprop = new double[kval];
	pval = new double *[kval];
	pxval = new double *[kval];
	for (int j = 0; j < kval; ++j)
	{
		pval[j] = new double[size];
		pxval[j] = new double[size];
	}
	while (it < maxIter)
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
		for (int j = 0; j < kval; ++j)
		{
			double sump = sum(pval[j], size);
			nprop[j] = sump / size;
			nlambda[j] = sump / sum(pxval[j], size);
		}
		ParamExp updatedPar(kval, nlambda, nprop);
		if (par.isConverge(updatedPar))
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

	//clean stuff
	delete[] nlambda;
	delete[] nprop;
	for (int j = 0; j < kval; ++j)
	{
		delete[] pval[j];
		delete[] pxval[j];
	}
	delete[] pval;
	delete[] pxval;
}

EMExp::~EMExp()
{
}

