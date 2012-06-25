#pragma once
#ifndef DATASRCS_H
#define DATASRCS_H

#include <vector>
#include <climits>

using namespace std;

class DataSources
{
public:

	vector<vector<float> > GetBars(int sizeX,int sizeY, int nrItems);
	vector<vector<float> > GetRandomBinary(int size, float activity, int nrItems);
	vector<vector<float> > GetRandom(int size, float valMin, float valMax, int nrItems);
	vector<vector<float> > GetRandomHCs(int nrHypercolumns, int nrRateUnits, int nrItems);
	vector<vector<float> > GetRandomHCsOrthogonal(int nrHypercolumns, int nrRateUnits, int nrItems);
	vector<vector<float> > GetRandomHCsOrthogonalLocal(int nrHypercolumns, int nrRateUnits, vector<int> localIds, int nrItems);

	vector<vector<float> > Reshape(vector<vector<float> > data, int rows, int cols);
	vector<vector<float> > toBinary(vector<vector<float> > data, int nrHc, int nrMc);
	vector<vector<float> > GetSubset(vector<vector<float> > data, int nrItems);
	int MaxIndex(vector<float> data);

private:

};

#endif