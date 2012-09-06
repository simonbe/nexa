/// Collection of data generating functions and help functions related to modifying data

#pragma once
#ifndef DATASRCS_H
#define DATASRCS_H

#include <vector>
#include <climits>

using namespace std;

class DataSources
{
public:

	// generates overlapping bars
	vector<vector<float> > GetBars(int sizeX,int sizeY, int nrItems);
	
	// binary random data of activity probabilitt to be 1
	vector<vector<float> > GetRandomBinary(int size, float activity, int nrItems);

	// random data
	vector<vector<float> > GetRandom(int size, float valMin, float valMax, int nrItems);

	// random data divided into hypercolumns
	vector<vector<float> > GetRandomHCs(int nrHypercolumns, int nrRateUnits, int nrItems);

	// orthogonal random data divided into hypercolumns
	vector<vector<float> > GetRandomHCsOrthogonal(int nrHypercolumns, int nrRateUnits, int nrItems);

	// orthogonal random data divided into hypercolumns depending on input structure localIds
	vector<vector<float> > GetRandomHCsOrthogonalLocal(int nrHypercolumns, int nrRateUnits, vector<int> localIds, int nrItems);

	// reshapes matrix
	vector<vector<float> > Reshape(vector<vector<float> > data, int rows, int cols);

	// float data to a binary representation
	vector<vector<float> > toBinary(vector<vector<float> > data, int nrHc, int nrMc);

	// float data to a binary representation
	vector<vector<float> > GetSubset(vector<vector<float> > data, int nrItems);

	// if delimiter NULL, split all into their characters
	vector<string> SplitString(string s, char delimiter = NULL, bool removeWhitespaces = true);

	// returns 0-based index of max value in data
	int MaxIndex(vector<float> data);

private:

};

#endif