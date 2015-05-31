/*
 * Test.cpp
 *
 *  Created on: May 31, 2015
 *      Author: young
 */
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
template<class T>
vector<vector<T> > perm(vector<T> &seq)
{
	vector<vector<T> > result;
	sort(seq.begin(), seq.end());
	do
	{
		result.push_back(seq);
	} while (next_permutation(seq.begin(), seq.end()));
	return result;
}

void print(vector<int> &data)
{
	for (int i = 0; i < data.size(); ++i)
	{
		cout << data.at(i) << " ";
	}
	cout << endl;
}
int main2()
{
	vector<int> test;
	test.push_back(0);
	test.push_back(0);
	test.push_back(0);
	test.push_back(1);
	test.push_back(1);
	vector<vector<int> > pm = perm(test);
	for (int i = 0; i < pm.size(); ++i)
	{
		print(pm.at(i));
	}
	return 0;
}

