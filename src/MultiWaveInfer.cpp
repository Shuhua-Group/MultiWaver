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
#include <boost/math/distributions/chi_squared.hpp>

using namespace std;

template<class T>
vector<vector<T> > perm(vector<T> &seq)
{
	vector<vector<T> > result;
	int size = seq.size();
	sort(seq.begin(), seq.end());
	do
	{
		if (seq.at(size - 1) == seq.at(size - 2))
			continue;
		int rsize = result.size();
		if (rsize > 0 && result.at(rsize - 1).at(size - 2) == seq.at(size - 1))
			continue;
		result.push_back(seq);
	} while (next_permutation(seq.begin(), seq.end()));
	return result;
}

double cv_chisq(int df, double alpha)
{
	boost::math::chi_squared dist(df);
	return boost::math::quantile(dist, 1 - alpha);
}

ParamExp findOptPar(const std::vector<double> &observ, int maxIter, double ancestryProp, double criticalValue, double epsilon)
{
	bool findOpt = false;
	int k = 1;
	ParamExp parPrev(k);
	EMExp em(parPrev, observ);
	em.iterate(maxIter, epsilon);
	double llkPrev = em.getLik();
	parPrev = em.getPar();
	while (!findOpt)
	{
		k++;
		ParamExp parCur(k);
		em.setPar(parCur);
		//em.updateLik();
		em.iterate(maxIter, epsilon);
		//cout << "K=" << k << "; Likelihood=" << setprecision(8) << em.getLik() << "; ";
		//em.getPar().print();
		double llkCur = em.getLik();
		if (2 * (llkCur - llkPrev) < criticalValue)
		{
			//cout << "Optimal K " << k - 1 << endl;
			findOpt = true;
		}
		else
		{
			llkPrev = llkCur;
			parPrev = em.getPar();
		}
		for (int i = 0; i < parPrev.getK(); ++i)
		{
			if (parPrev.getProp(i) * ancestryProp < kMinP)
			{
				findOpt = true;
				break;
			}
		}
	}
	//parPrev.print();
	parPrev.sortByLambda();
	return parPrev;
}

void help()
{
	cout << kProgramName << " v" << kVersion << endl;
	cout << kProgramName << " is designed to scan for multiple waves of admixture events and estimated corresponding parameters." << endl;
	cout << "General usage: " << kProgramName << " <arguments>" << endl;
	cout << "Arguments" << endl;
	cout << "\t-i/--input\t<string>\tInput of the ancestral tracks [required]" << endl;
	cout << "\t-a/--alpha\t[double]\tSignificance level to reject null hypothesis in LRT [optional, default 0.05]" << endl;
	cout << "\t-e/--epsilon\t[double]\tEpsilon to check whether a parameter converge or not [optional, default 0.000001] " << endl;
	cout << "\t-l/--lower\t[double]\tLower bound to discard short tracks [optional, default 0]" << endl;
	cout << "\t-m/--maxIt\t[integer]\tMax number of iterations to perform EM [optional, default 10000]" << endl;
	cout << "Option" << endl;
	cout << "\t-h/--help\tPrint help message." << endl;
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
	double alpha = 0.05;
	double epsilon = 0.000001;
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
		else if (arg == "-a" || arg == "--alpha")
		{
			alpha = atof(argv[++i]);
		}
		else if (arg == "-e" || arg == "--epsilon")
		{
			epsilon = atof(argv[++i]);
		}
	}
	if (filename.size() == 0)
	{
		cerr << "Input file name required, please check help" << endl;
		exit(1);
	}
	//echo command entered
	cout << endl << "// COMMAND ";
	for (int i = 0; i < argc; ++i)
	{
		cout << argv[i] << " ";
	}
	cout << endl << endl;

	srand(time(NULL));
	vector<double> observ;
	//Reading test data from file named as test.dat
	ifstream fin(filename.c_str());
	if (!fin.is_open())
	{
		cout << "Can't open file " << filename << endl;
		exit(1);
	}
	cout << "Reading data from " << filename << "..." << endl;
	map<string, vector<double> > segs;
	map<string, double> sumLengths;
	vector<string> labels;
	double totalLength = 0;
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
			sumLengths[label] = 0.0;
			labels.push_back(label);
		}
		sumLengths.at(label) += len;
		totalLength += len;
		//greater than or equal to cutoff
		if (len > lower)
		{
			segs.at(label).push_back(len - lower);
		}
	}
	fin.close();

	cout << "Start scan for admixture waves... " << endl;
	int numLabel = static_cast<int>(labels.size());
	map<string, double> mixtureProps; //S_k
	map<string, ParamExp> optPars;
	double criticalValue = cv_chisq(2, alpha);
	for (int i = 0; i < numLabel; ++i)
	{
		string label = labels.at(i);
		cout << "Perform EM scan for waves of population " << label << "..." << endl;
		mixtureProps[label] = sumLengths.at(label) / totalLength;
		optPars[label] = findOptPar(segs.at(label), maxIter, mixtureProps.at(label), criticalValue, epsilon);
	}
	cout << "Finished scanning for admixture waves." << endl << endl;
//	for (int i = 0; i < numLabel; ++i)
//	{
//		string label = labels.at(i);
//		cout << "Population: " << label << "; Mix proportion: " << mixtureProps.at(label) << endl;
//		cout << "Optimal Parameters: ";
//		optPars.at(label).print();
//	}

	int totalNumOfWaves = 0;
	vector<int> popOrder;
	/* calculate survivalProportion m_Ik(j) of each wave */
	map<int, vector<double> > survivalProps; //m_Ik(j)
	for (int i = 0; i < numLabel; ++i)
	{
		string label = labels.at(i);
		int numOfExp = optPars.at(label).getK();
		totalNumOfWaves += numOfExp;
		vector<double> tempMIK;
		double tempSum = 0;
		double temp[numOfExp];
		for (int j = 0; j < numOfExp; ++j)
		{
			popOrder.push_back(i);
			temp[j] = optPars.at(label).getProp(j) / optPars.at(label).getLambda(j);
			tempSum += temp[j] / temp[0];
		}
		tempSum = mixtureProps.at(label) / tempSum;
		for (int j = 0; j < numOfExp; ++j)
		{
			tempMIK.push_back(tempSum * temp[j] / temp[0]);
		}
		survivalProps[i] = tempMIK;
	}

	bool hasSolution = false;
	cout << "There is(are) " << totalNumOfWaves - 1 << " wave(s) of admixture event(s) detected" << endl;
	cout << "-----------------------------------------------------------------------------" << endl;
	cout << setw(44) << "Results summary" << endl << endl;
	cout << setw(32) << "Parental population" << setw(32) << "Admixture proportion" << endl;
	for (map<string, double>::iterator iter = mixtureProps.begin(); iter != mixtureProps.end(); iter++)
	{
		cout << setw(32) << iter->first << setw(32) << iter->second << endl;
	}
	cout << endl;

	int scenarioCount = 0;
	vector<vector<int> > allOrder = perm(popOrder);
	for (vector<vector<int> >::iterator iter = allOrder.begin(); iter != allOrder.end(); ++iter)
	{

		/*for a given order, calculate alpha, H and time T */
		double mInOrder[totalNumOfWaves];
		for (int i = 0; i < numLabel; ++i)
		{
			int index = 0;
			for (int j = 0; j < totalNumOfWaves; ++j)
			{
				if (iter->at(j) == i)
				{
					mInOrder[j] = survivalProps[i].at(index);
					index++;
				}
			}
		}

		/* calculate alpha */
		double alphaInOrder[totalNumOfWaves];
		double multiplier = 1.0;
		alphaInOrder[0] = mInOrder[0];
		for (int i = 1; i < totalNumOfWaves - 1; ++i)
		{
			multiplier *= (1 - alphaInOrder[i - 1]);
			alphaInOrder[i] = mInOrder[i] / multiplier;
		}
		alphaInOrder[totalNumOfWaves - 1] = 1.0 - alphaInOrder[totalNumOfWaves - 2];

		/* calculate total ancestry proportion of kth ancestral population at t generation H_k(t)*/
		double hInOrder[numLabel][totalNumOfWaves - 1];
		/*
		 * Initial last element of H for each population.
		 * look at the last two positions of population order: if the second last
		 * position is from population k, then set the last value of H for population
		 * k as the alpha in order of the second last value; if the last position is
		 * from population k, then set the last value of H for population k as the alpha
		 * in order of the last value; otherwise, set the last value of H for population
		 * k as 0.0
		 */
		for (int i = 0; i < numLabel; ++i)
		{
			int lastIndexOfH = totalNumOfWaves - 2;
			for (int j = 0; j < numLabel; ++j)
			{
				if (iter->at(lastIndexOfH) == j)
				{
					hInOrder[j][lastIndexOfH] = alphaInOrder[lastIndexOfH];
				}
				else if (iter->at(lastIndexOfH + 1) == j)
				{
					hInOrder[j][lastIndexOfH] = alphaInOrder[lastIndexOfH + 1];
				}
				else
				{
					hInOrder[j][lastIndexOfH] = 0;
				}
			}
			/*
			 * Recursively update H
			 * if alpha in order from population k, then update H as H*(1-alpha)+alpha;
			 * else update H as H*(1-alpha)
			 */
			for (int j = totalNumOfWaves - 3; j >= 0; --j)
			{
				if (iter->at(j) == i)
				{
					hInOrder[i][j] = hInOrder[i][j + 1] * (1 - alphaInOrder[j]) + alphaInOrder[j];
				}
				else
				{
					hInOrder[i][j] = hInOrder[i][j + 1] * (1 - alphaInOrder[j]);
				}
			}
		}

		double admixTime[totalNumOfWaves];
		int indexes[numLabel];
		for (int i = 0; i < numLabel; ++i)
		{
			indexes[i] = 0;
		}

		for (int i = 0; i < totalNumOfWaves; ++i)
		{
			for (int j = 0; j < numLabel; ++j)
			{
				if (iter->at(i) == j)
				{
					double tempSum = 0;
					for (int k = 0; k < i; ++k)
					{
						tempSum += (1 - hInOrder[j][k]) * admixTime[k];
					}
					double rate = optPars.at(labels.at(j)).getLambda(indexes[j]) - tempSum;
					admixTime[i] = rate / (1 - hInOrder[j][i]);
					indexes[j]++;
					break;
				}
			}
		}

		//check whether the results is reasonable or not
		bool isReasonable = true;
		for (int i = 0; i < totalNumOfWaves - 1; ++i)
		{
			if (admixTime[i] < 0)
			{
				isReasonable = false;
				break;
			}
		}
		for (int i = 1; i < totalNumOfWaves; ++i)
		{
			admixTime[i] += admixTime[i - 1];
		}
		if (abs(admixTime[totalNumOfWaves - 1] - admixTime[totalNumOfWaves - 2]) / admixTime[totalNumOfWaves - 1] > 0.05)
		{
			isReasonable = false;
		}
		if (isReasonable)
		{
			if (!hasSolution)
			{
				hasSolution = true;
			}
			cout << endl;
			cout << "Possible scenario: #" << ++scenarioCount << endl;
			cout << setw(10) << admixTime[totalNumOfWaves - 1]; //time
			cout << ": (" << setw(1) << iter->at(totalNumOfWaves - 1); //population
			cout << ", " << setw(10) << alphaInOrder[totalNumOfWaves - 1]; //proportion
			//cout << ") -----------|----------- (";
			cout << ") =========>||<========= (";
			cout << setw(1) << iter->at(totalNumOfWaves - 2) << ", "; //population
			cout << setw(10) << alphaInOrder[totalNumOfWaves - 2] << ") :"; //proportion
			cout << setw(10) << admixTime[totalNumOfWaves - 2] << endl; //time
			for (int i = totalNumOfWaves - 3; i >= 0; --i)
			{
				//cout << setw(40) << "|" << endl << setw(40) << "|" << endl << setw(40) << "|" << endl;
				cout << setw(40) << "||" << endl << setw(40) << "||" << endl << setw(40) << "||" << endl;
				cout << setw(10) << admixTime[i] << ": (" << setw(1) << iter->at(i);
				//cout << ", " << setw(10) << alphaInOrder[i] << ") -----------|" << endl;
				cout << ", " << setw(10) << alphaInOrder[i] << ") =========>||" << endl;
			}
			//cout << setw(40) << "|" << endl << setw(40) << "|" << endl << setw(40) << "|" << endl << endl;
			cout << setw(40) << "||" << endl << setw(40) << "||" << endl << setw(40) << "||" << endl << endl;
		}
	}
	if (!hasSolution)
	{
		cout << "No possible scenario for current settings, please change settings and retry!" << endl;
	}
	else
	{
		cout << "Hint: " << endl;
		for (int i = 0; i < numLabel; ++i)
		{
			cout << i << ": population-" << labels.at(i) << "; ";
		}
		cout << endl;
	}
	cout << "-----------------------------------------------------------------------------" << endl;

	return 0;
}
