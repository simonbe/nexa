#pragma once
#include <vector>
#include <iostream>

class Population;

using namespace std;

/// <summary>	Contains functions for handling the distribution of units according
/// 			to a parallelization scheme across a machine. </summary>

class MPIDistribution
{

public:

	enum ParallelizationSchemeLayer // necessary, move (?) 
	{
		ParallelizationSplitPopulations,
		ParallelizationRateUnits,
		ParallelizationHypercolumns,
		ParallelizationDefault // = ParallelizationSplitPopulations
	};

	MPIDistribution(Population* layer)
	{
		m_population = layer;
		m_commsHCsCreated = false;
		m_commsLayersCreated = false;
	}

	/// <summary>	Divide a fixed number of units equally. </summary>
	///
	/// <param name="nrUnits">	The nr units. </param>
	/// <param name="mpiRank">	This mpi process. </param>
	/// <param name="mpiSize">	Nr total processes. </param>
	///
	/// <returns>	Distribution. </returns>

	vector<vector<long> > DivideEqualByUnits(unsigned long nrUnits,int mpiRank,int mpiSize)
	{
		float part = (float)nrUnits/(float)mpiSize;

		if(mpiRank == 0)
		if(part<1.0)
			cout<<"WARNING in DivideEqualByUnits.";

		vector<vector<long> > indexesInterval;

		for(int i=0;i<mpiSize;i++)
		{
			long lb = (long)(part * (float)i);
			long ub = (long)(part * (float)(i+1));

			if(i == mpiSize-1 )
				ub = nrUnits;

			vector<long> v;
			v.push_back(lb);
			v.push_back(ub);
			indexesInterval.push_back(v);
		}

		return indexesInterval;
	}

	/// <summary>	Divides a single hypercolumn across processes.
	/// 			Used in DivideEquaByHypercolumn
	///				(TODO: check dependencies, make private) </summary>
	///
	/// <param name="nrRateUnits">	Nr of rate units in each column. </param>
	/// <param name="mpiRank">	  	This mpi process. </param>
	/// <param name="mpiSize">	  	Nr total processes. </param>
	///
	/// <returns>	Distribution. </returns>

	vector<vector<long> > DivideEqualByUnitsSingleHC(vector<int> nrRateUnits,int mpiRank,int mpiSize)
	{
		float part = (float)nrRateUnits.size()/(float)mpiSize;

		if(part > 1) // here because of specific parallelization strategy
		{
			vector<vector<long> > indexesInterval(mpiSize);

			int totUnits = 0;
			for(int i=0;i<nrRateUnits.size();i++)
				totUnits+=nrRateUnits[i];

			float frac = (float)totUnits/(float)mpiSize;

			float ub=frac,lb=0;
			for(int i=0;i<mpiSize;i++)
			{
				vector<long> v(2);
				v[0] = (long)lb; v[1] = (long)ub;
				ub+=frac; lb+=frac;
				if(i==mpiSize-1)
					v[1] = totUnits;
				indexesInterval[i] = v;
			}
			
			return indexesInterval;
		}

		vector<vector<int> > nodeHC(nrRateUnits.size());
		vector<int> hcNode(mpiSize);

		for(int i=0;i<mpiSize;i++)
		{
			int k = i*part;
			nodeHC[k].push_back(i);
			hcNode[i] = k;
		}

		if(mpiRank == 0)
		if(part<1.0)
			cout<<"WARNING in DivideEqualByUnits.";

		vector<vector<long> > indexesInterval;

		//int i=mpiRank;

		// may not scale perfectly
		for(int i=0;i<mpiSize;i++)
		{
			int k=hcNode[i];
			int nrUnits = 0;
			for(int j=0;j<k;j++)
				nrUnits+=nrRateUnits[j];

			int loc = 0;
			for(int j=0;j<nodeHC[k].size();j++)
				if(nodeHC[k][j] == i)
					loc = j;
			
			part = 1.0/(float)nodeHC[k].size();
			long lb = long((float)nrRateUnits[k]*(part * (float)loc)+nrUnits);
			long ub = long((float)nrRateUnits[k]*(part * (float)(loc+1))+nrUnits);

			if(loc == nodeHC[k][nodeHC[k].size()-1])
				ub = nrUnits+nrRateUnits[k];

			vector<long> v;
			v.push_back(lb);
			v.push_back(ub);
			indexesInterval.push_back(v);
		}

		// remove info about non-local related hcs
		//for(int i=0;i<indexesInter

		return indexesInterval;
	}

	///
	// returns [0] as process indexes for each hypercolumn
	// returns [1] as unit ids division for each process
	vector<vector<vector<long> > > DivideEqualByHypercolumn(vector<int> nrRateUnits, ParallelizationSchemeLayer parallelizationScheme, int mpiRank,int mpiSize);

	//Default division, puts different populations on different processes if possible
	vector<vector<vector<long> > > DivideEqualByTotalUnits(vector<int> nrRateUnits, long currentUnit, long totalNrUnits, int mpiRank,int mpiSize);

	//Default division, puts different populations on different processes if possible
	vector<vector<vector<long> > > DivideEqualByPopulations(vector<int> nrRateUnits, long currentUnit, long totalNrUnits, int mpiRank,int mpiSize);


	void MPICreateCommLayer();

	vector<MPI_Comm*> MPIGetCommHCs() // MPI communicators for each hypercolumn
	{
		return m_mpiCommHCs;
	}

	MPI_Comm* MPIGetCommLayer() // MPI communicator for entire layer
	{
		return m_mpiCommLayer;
	}

	void MPICreateCommsHypercolumns();

	// Distribution of values across nodes

	void MPIMakeHypercolumnsValuesLocal();
	void MPIMakeLayerValuesLocal();
	
	// Fetching operations

	//vector<long> MPIGetMaxIdInHypercolumns();

	vector<float> MPILayerReduceVariable(vector<float> data);
	vector<vector<float> > MPILayerReduceVariable(vector<vector<float> > data);
	vector<vector<int> > MPILayerGatherVariables(vector<vector<int> > data);

	void MakeActivityValuesLocal(bool fullLayer = false);

private:

	
	MPI_Comm* m_mpiCommLayer;
	int m_mpiSizeLocal;
	int m_mpiRankLocal;

	bool m_commsHCsCreated;
	bool m_commsLayersCreated;

	Population* m_population;
	vector<MPI_Comm*> m_mpiCommHCs;

};