/*
 * EMExp.hpp
 *
 *  Created on: May 26, 2015
 *      Author: young
 */
/******************************************************************************
 * @brief A class perform EM
 * Class EMExp to perform EM algorithm in parameters estimation with hidden variable
 * The EMExp class has three attributes:
 * 1) parameters to be estimated;
 * 2) a sequence of observations, which used to estimate the parameters;
 * 3) and the likelihood corresponding to the observations and parameters.

 * The method iterate is used to perform EM algorithm, including two steps:
 * E-Step:
 * 1) calculate p_j=prob(z_i=j|x_i, theta_t)
 * =prob(x_i|z_i=j,theta_t)*prob(z_i=j,theta_t)/(sum_j{from 1 to K}(.)
 * denotes the numerator
 * 2) calculate p_j*xi ; j=1,...,K
 * theta_t: parameters at the t-th iteration
 * M-Step:
 * 1) update prop, refer as m, m_j=sum(p_j)/sum(p_1+p_2+...+p_K)=sum(p_j)/n
 * 2) update lambda, lambda_j=sum(p_j)/sum(p_j*x_i)
 * note: update the parameters for t+1 times iteration
 ******************************************************************************/
#ifndef EMEXP_HPP_
#define EMEXP_HPP_

#include <vector>
#include <string>
#include "ParamExp.hpp"

class EMExp
{
public:
	EMExp(const ParamExp &par, const std::vector<double> &observ);
	double getLik() const;
	ParamExp & getPar();
	void setPar(const ParamExp &par);
	void updateLik();
	void iterate(int maxIter);
	virtual ~EMExp();
private:
	double lik;
	ParamExp par;
	std::vector<double> observ;
};

#endif /* EMEXP_HPP_ */
