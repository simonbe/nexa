#include <iostream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <algorithm>

#include "DataSources.h"


vector<vector<float> > DataSources::GetRandomBinary(int size, float activity, int nrItems)
{
	//
	
	vector<vector<float> > data(nrItems, vector<float>(size,0.0));

	for(int i=0;i<nrItems;i++)
	{
		for(int j=0;j<size;j++)
		{
			float r = (float)rand()/(float)RAND_MAX;

			if(r<activity)
				data[i][j] = 1.0;
		}
	}

	return data;
}

vector<vector<float> > DataSources::GetRandom(int size, float valMin, float valMax, int nrItems)
{
	vector<vector<float> > data(nrItems, vector<float>(size,0.0));

	for(int i=0;i<nrItems;i++)
	{
		for(int j=0;j<size;j++)
		{
			float r = rand()/(float(RAND_MAX)+1);
			data[i][j]   = (valMax-valMin)*r + valMin;
		}
	}

	return data;
}

vector<vector<float> > DataSources::GetRandomHCs(int nrHypercolumns, int nrRateUnits, int nrItems)
{
	int size = nrHypercolumns*nrRateUnits;
	vector<vector<float> > data(nrItems, vector<float>(size,0.0));

	for(int i=0;i<nrItems;i++)
	{
		int loc = 0;
		for(int j=0;j<nrHypercolumns;j++)
		{
			float r = nrRateUnits*((float)rand()/(float)RAND_MAX);
			
			for(int k=0;k<nrRateUnits;k++)
			{
				if(r>k && r<k+1)
				{
					loc = k;
					break;
				}
			}

			data[i][nrRateUnits*j + loc] = 1.0;
		}
		
	}

	return data;
}

vector<vector<float> > DataSources::GetRandomHCsOrthogonal(int nrHypercolumns, int nrRateUnits, int nrItems)
{
	int size = nrHypercolumns*nrRateUnits;
	vector<vector<float> > data(nrItems, vector<float>(size,0.0));

	for(int i=0;i<nrItems;i++)
	{
		int loc = 0;
		for(int j=0;j<nrHypercolumns;j++)
		{
			data[i][nrRateUnits*j + i] = 1.0;
		}
	}

	return data;
}

vector<vector<float> > DataSources::GetRandomHCsOrthogonalLocal(int nrHypercolumns, int nrRateUnits, vector<int> localIds, int nrItems)
{
	int size = nrHypercolumns*nrRateUnits;
	vector<vector<float> > data(nrItems, vector<float>(localIds.size(), 0.0));
	int start = localIds[0];
	int end = localIds[localIds.size()-1]; // assumed distribution

	int loc;
	for(int i=0;i<nrItems;i++)
	{
		for(int j=0;j<nrHypercolumns;j++)
		{
			loc = nrRateUnits*j + i;
			if(loc >= start && loc<=end)
				data[i][loc - start] = 1.0;
		}
	}

	return data;
}


vector<vector<float> > DataSources::GetBars(int sizeX,int sizeY, int nrItems)
{
	vector<vector<float> > out;

	for(int i=0;i<nrItems;i++)
	{
		vector<vector<float> > item(sizeY,vector<float>(sizeX));

		// 10% chance of every bar occuring

		float r;
		bool onefound = false;

		for(int k=0;k<sizeY;k++)
		{	
			r = (float)rand()/(float)RAND_MAX;
			if(r < 1.0/(float)sizeY)//0.1)
			{
				for(int j=0;j<sizeX;j++)
				{
					item[k][j] = 1.0;
				}
				onefound = true;
			}

			if(onefound == true)
				break;

			if(k == sizeY-1)
				k=-1;
		}

		onefound = false;

		for(int k=0;k<sizeX;k++)
		{	
			r = (float)rand()/(float)RAND_MAX;
			if(r < 1.0/(float)sizeX)
			{
				for(int j=0;j<sizeY;j++)
				{
					item[j][k] = 1.0;
				}
				onefound = true;
			}

			if(onefound == true)
				break;

			if(k == sizeY-1)
				k=-1;
		}

		vector<float> o2(sizeX*sizeY);

		int index = 0;
		for(int x=0;x<sizeX;x++)
		{
			for(int y=0;y<sizeY;y++)
			{
				o2[index] = item[y][x];
				index++;
			}
		}

		out.push_back(o2);
	}

	return out;
}


vector<vector<float> > DataSources::Reshape(vector<vector<float> > data, int rows, int cols)
{
	vector<vector<float> > out(rows,vector<float>(cols));

	if(data.size()*data[0].size() != rows*cols)
	{
		cout<<"WARNING, not same sizes in reshape. ("<<data.size()<<"/"<<data[0].size()<<", "<<rows<<"/"<<cols<<"]\n";
		if(data.size()*data[0].size() < rows*cols)
		{
			cout<<"ERROR, data size smaller than new matrix.\n";
			cout.flush();
		}
	}

	int inew = 0, jnew = 0;
	for(int i=0;i<data.size();i++)
	{
		for(int j=0;j<data[i].size();j++)
		{
			out[inew][jnew] = data[i][j];
			jnew++;
			if(jnew >= cols)
			{
				inew++;
				jnew = 0;
			}
		}
	}

	return out;
}

vector<vector<float> > DataSources::toBinary(vector<vector<float> > data, int nrHc, int nrMc)
{
	vector<vector<float> > out(data.size(),vector<float>(nrHc*nrMc));

	int minData = 1e6;
	for(int i=0;i<data.size();i++)
	{
		for(int j=0;j<data[i].size();j++)
		{
			if(data[i][j]<minData)
				minData = data[i][j];
		}
	}

	for(int j=0;j<data.size();j++)
	{
		vector<float> f(nrHc*nrMc);
		int currentI = 0;
		for(int i=0;i<data[j].size();i++)
		{
			out[j][currentI+data[j][i]-minData] = 1.0;
			currentI+=nrMc;		
		}
	}

	return out;
}

vector<vector<float> > DataSources::GetSubset(vector<vector<float> > data, int nrItems)
{
	vector<vector<float> > out(nrItems);
	for(int i=0;i<nrItems;i++)
		out[i] = data[i];

	return out;
}

int DataSources::MaxIndex(vector<float> data)
{
	int maxIndex = -1;
	float maxValue = -1e8;
	for(int i=0;i<data.size();i++)
	{
		if(data[i]>maxValue)
		{
			maxIndex = i;
			maxValue = data[i];
		}
	}

	return maxIndex;
}

vector<string> DataSources::SplitString(string s, char delimiter, bool removeWhitespaces)
{
	vector<string> out;

	if(removeWhitespaces)
			s.erase(remove(s.begin(),s.end(),' '),s.end()); 

	if(delimiter == NULL)
	{
		for(int i=0;i<s.length();i++)
			out.push_back(string(1,s[i]));
	}
	else
	{
		string spart;
		stringstream stream(s);
		while( getline(stream, spart, delimiter) )
			out.push_back(spart);
	}

	return out;
}
