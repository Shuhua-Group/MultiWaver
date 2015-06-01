/*
 * Test.cpp
 *
 *  Created on: May 31, 2015
 *      Author: young
 */
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

using namespace std;
template<class T>
vector<vector<T> > perm(vector<T> &seq, int kSum)
{
	vector<vector<T> > result;
	sort(seq.begin(), seq.end());
	do
	{
		if (seq.at(0) != seq.at(1))
			result.push_back(seq);
	} while (next_permutation(seq.begin(), seq.end()));
//	bool flag = 1;
//	while (flag)
//	{
//		vector<T> temp;
//		for (typename vector<T>::reverse_iterator it = seq.rbegin(); it != seq.rbegin() + kSum; ++it)
//		{
//			temp.push_back(*it);
//			//cout << "t" << *it;
//		}
//
//		sort(temp.begin(), temp.end());
//		temp.erase(unique(temp.begin(), temp.end()), temp.end());
//		//cout << "temp" << temp.size() << endl;
//
//		//confirm no duplicate
//		if (temp.size() == kSum)
//			result.push_back(seq);
//
//		flag = next_permutation(seq.begin(), seq.end());
//	}

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
int main()
{
	vector<int> test;
	test.push_back(0);
	test.push_back(0);
	test.push_back(0);
	test.push_back(2);
	test.push_back(1);
	test.push_back(1);
	vector<vector<int> > pm = perm(test, 2);
	for (int i = 0; i < pm.size(); ++i)
	{
		print(pm.at(i));
	}
	map<int, int> mp;
	mp[7] = 1;
	mp[5] = 2;
	mp[2] = 0;
	mp[9] = 7;
	for (map<int, int>::iterator it = mp.begin(); it != mp.end(); it++)
	{
		cout << "key=" << it->first << "; value=" << it->second << endl;
	}
	return 0;
}

