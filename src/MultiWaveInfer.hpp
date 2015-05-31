/*
 * MultiWaveInfer.hpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */

#ifndef MULTIWAVEINFER_HPP_
#define MULTIWAVEINFER_HPP_

#include <vector>
#include "ParamExp.hpp"

const double kminP = 0.01; //minimum proportion for a wave
const std::string kVersion = "0.0.1";
template <class T>
std::vector<std::vector<T> > perm(std::vector<T> &seq);
ParamExp findOptPar(const std::vector<double> &observ, int maxIter);
void help();

#endif /* MULTIWAVEINFER_HPP_ */
