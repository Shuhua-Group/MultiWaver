/*
 * MultiWaveInfer.hpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */

#ifndef MULTIWAVEINFER_HPP_
#define MULTIWAVEINFER_HPP_

#include <vector>
#include <string>
#include "ParamExp.hpp"

const double kMinP = 0.05; //minimum proportion for a wave
//const double kCriticalValue = 5.991; //chi-squared value for p-value 0.05 with degree of freedom 2
const std::string kVersion = "1.0.1";
const std::string kProgramName = "MultiWaveInfer";

/*
 * perform permutation over a sequence
 */
template<class T>
std::vector<std::vector<T> > perm(std::vector<T> &seq);

/*
 * calculate the critical value of chi-squared distribution with significance level alpha
 */
double cv_chisq(int df, double alpha = 0.05);

/*
 * perform EM to find the best number of exponential distribution and estimate
 * corresponding parameters
 */
ParamExp findOptPar(const std::vector<double> &observ, int maxIter, double ancestryProp, double criticalValue = 5.991);

/*
 * print help information
 */
void help();

#endif /* MULTIWAVEINFER_HPP_ */
