/*
 * MultiWaveInfer.cpp
 *
 *  Created on: May 26, 2015
 *  Author: young
 */
#include <map>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "EMExp.hpp"
#include "MultiWaveInfer.hpp"

using namespace std;

const double kCriticalValue = 5.991;

template<class T>
vector<vector<T> > perm(vector<T> &seq)
{
	vector<vector<T> > result;
	sort(seq.begin(), seq.end());
	do
	{
		if (seq.at(0) != seq.at(1))
		{
			result.push_back(seq);
		}
	} while (next_permutation(seq.begin(), seq.end()));
	return result;
}

ParamExp findOptPar(const std::vector<double> &observ, int maxIter)
{
	bool findBest = false;
	int k = 1;
	ParamExp parPrev(k);
	EMExp em(parPrev, observ);
	em.iterate(maxIter);
	double llkPrev = em.getLik();
	parPrev = em.getPar();
	while (!findBest)
	{
		k++;
		ParamExp parCur(k);
		em.setPar(parCur);
		//em.updateLik();
		em.iterate(maxIter);
		//cout << "K=" << k << "; Likelihood=" << setprecision(8) << em.getLik() << "; ";
		//em.getPar().print();
		double llkCur = em.getLik();
		if (2 * (llkCur - llkPrev) < kCriticalValue)
		{
			//cout << "Optimal K " << k - 1 << endl;
			findBest = true;
		}
		else
		{
			llkPrev = llkCur;
			parPrev = em.getPar();
		}
	}
	//parPrev.print();
	parPrev.sortByLambda();
	return parPrev;
}

void help()
{
	cout << "MultiWaveInfer version " << kVersion << endl;
	cout << "MultiWaveInfer is used to estimated parameters of mixture exponential distributions" << endl;
	cout << "Arguments:" << endl;
	cout << "	-h/--help	print help message [optional]" << endl;
	cout << "	-i/--input	file name of observations [required]" << endl;
	cout << "	-l/--lower	number of exponential distributions [optional,0.0]" << endl;
	cout << "	-m/--maxIt	max number of iteration to perform EM [optional, 10000]" << endl;
}

int main(int argc, char **argv)
{
	if (argc < 2)
	{
		cerr << "Need more arguments than provided, please use -h/--help to get help" << endl;
		exit(1);
	}
	string filename = "";
	double lower = 0;
	int maxIter = 10000;
	for (int i = 1; i < argc; ++i)
	{
		string arg(argv[i]);
		if (arg == "-h" || arg == "--help")
		{
			help();
			exit(0);
		}
		else if (arg == "-i" || arg == "--input")
		{
			filename = string(argv[++i]);
		}
		else if (arg == "-l" || arg == "--lower")
		{
			lower = atof(argv[++i]);
		}
		else if (arg == "-m" || arg == "--maxIt")
		{
			maxIter = atoi(argv[++i]);
		}
	}
	if (filename.size() == 0)
	{
		cerr << "File name required, please check help" << endl;
		exit(1);
	}
	srand(time(NULL));
	vector<double> observ;
	//Reading test data from file named as test.dat
	ifstream fin(filename.c_str());
	if (!fin.is_open())
	{
		cout << "Can't open file " << filename << endl;
		exit(1);
	}
	map<string, vector<double> > segs;
	map<string, double> sumLens;
	vector<string> labels;
	double tLen = 0;
	string line;
	while (getline(fin, line))
	{
		istringstream ss(line);
		double start, end;
		string label;
		ss >> start;
		ss >> end;
		ss >> label;
		double len = end - start;
		//no key not found, then create an new one
		if (segs.find(label) == segs.end())
		{
			vector<double> tmp;
			segs[label] = tmp;
			sumLens[label] = 0.0;
			labels.push_back(label);
		}
		sumLens.at(label) += len;
		tLen += len;
		//greater than or equal to cutoff
		if (len > lower)
		{
			segs.at(label).push_back(len);
		}
	}
	fin.close();
	int numLabel = static_cast<int>(labels.size());
	map<string, double> props;
	map<string, ParamExp> optPars;
	for (int i = 0; i < numLabel; ++i)
	{
		string label = labels.at(i);
		props[label] = sumLens.at(label) / tLen;
		optPars[label] = findOptPar(segs.at(label), maxIter);
	}
	for (int i = 0; i < numLabel; ++i)
	{
		string label = labels.at(i);
		cout << "Population: " << label << "; Mix proportion: " << props.at(label) << endl;
		cout << "Optimal Parameters: ";
		optPars.at(label).print();
	}

	//gw begin
	//gw
	// get m
//	map<int, vector<double> > optLabelM;
//	for (int i = 0; i < numLabel; ++i)
//	{
//		string label = labels.at(i);
//		double SurProp = props[label]; // total contribute for a pop
//		ParamExp tmpKpara = optPars[label];
//		int k = tmpKpara.getK();
//		vector<double> m; //mk_j
//		double r[k];
//		double rSum = 0;
//		for (int j = 0; j < k; ++j)
//		{
//			r[j] = tmpKpara.getProp(j) / tmpKpara.getLambda(j); // the ratio of mk_j
//			rSum += r[j] / r[0];
//		}
//
//		double temp = SurProp / rSum;
//
//		for (int j = 0; j < k; ++j)
//		{
//			m.push_back(temp * r[j] / r[0]);
//		}
//
//		optLabelM[i] = m;
//	}
//
//	//get an order
//	vector<int> EasyOrder; // start from pop 0
//	for (int i = 0; i < numLabel; ++i)
//	{
//		string label = labels.at(i);
//		ParamExp tmpKpara = optPars[label];
//		int k = tmpKpara.getK();
//		for (int j = 0; j < k; ++j)
//		{
//			EasyOrder.push_back(i);
//		}
//	}
//
//	//calculate the sum of admiture events
//	int kSum = 0;
//	for (int i = 0; i < numLabel; ++i)
//	{
//		string label = labels.at(i);
//		ParamExp tmpKpara = optPars[label];
//		int k = tmpKpara.getK();
//		kSum += k;
//	}
//
//	//get all order
//	vector<vector<int> > AllOrder = perm(EasyOrder);
//	int order[kSum];
//	// get m in order
//	double MinOrder[kSum]; // m in order
//	for (int i = 0; i < kSum; ++i)
//	{
//		for (int j = 0; j < numLabel; ++j)
//		{
//			vector<double> m_j = optLabelM[j];
//			vector<double>::iterator iter = m_j.begin();
//			if (order[i] == j + 1)
//			{
//				MinOrder[i] = *iter;
//				++iter;
//			}
//		}
//
//	}
//
//	//gw
//	//get a
//	double a[kSum];
//	a[kSum - 1] = MinOrder[kSum - 1];
//	for (int i = 0; i < kSum; ++i)
//	{
//		int MutiA;
//		MutiA *= (1 - a[kSum - i - 1]);
//		a[kSum - i - 2] = MinOrder[kSum - i - 2] / MutiA;
//	}
//
//	//gw
//	//get H
//	double H[numLabel][kSum];
//	for (int k = 0; k < numLabel; ++k)
//	{
//		for (int h = 0; h < kSum; ++h)
//		{
//			//get w
//			int p = h;
//			int w = 0;
//			if (order[p] != k)
//			{
//				p--;
//			}
//			for (int c = 0; c < p; ++p)
//			{
//				if (order[c] == k)
//					w++;
//			}
//
//			for (int j = 0; j < w; ++j)
//			{
//				//get Ikj
//				int J = 0;
//				int I_j = 0;
//				for (int c = 0; c < w; ++c)
//				{
//					if (order[c] == k)
//					{
//						J++;
//						if (J - 1 == j)
//							I_j = c;
//					}
//				}
//
//				//get multiply
//				double multip = 0;
//				for (int t = I_j; t < h; ++t)
//				{
//					multip *= (1 - a[t]);
//				}
//				H[k][h] += a[I_j] * multip;
//			}
//
//		}
//	}
//
//	//gw
//	//get T
//	double Tdelta[kSum];
//	double laststep;
//	double TotT = 0;
//	for (int i = kSum; i > 0; --i)
//	{
//		for (int k = 0; k < numLabel; ++k)
//		{
//
//			if (order[i] == k)
//			{
//				string label = labels.at(k);
//				ParamExp tmpKpara = optPars[label];
//				int num = 0;
//				for (int j = 0; j < i; ++j)
//				{
//					if (order[j] == k)
//						num++;
//				}
//				double uk_j = tmpKpara.getLambda(num);
//
//				Tdelta[i] = (uk_j - laststep) / (1 - H[k][i]);
//				laststep = (1 - H[k][i]) * Tdelta[i];
//			}
//		}
//		TotT += Tdelta[i];
//	}

	return 0;
}
