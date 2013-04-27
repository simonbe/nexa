#include <mpi.h>
#include "MPIDistribution.h"
#include "NetworkPopulation.h"

/// <summary>	Creates communicator for population. 
/// 			TODO: Check if necessary any longer to call in network initialization. (should be able to remove)</summary>
 
void MPIDistribution::MPICreateCommLayer()
{
	if(m_population->network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"Creating MPI Communicator for layer...";cout.flush();
	}

	vector<int> localHypercolumns = ((PopulationColumns*)m_population)->GetLocalHypercolumnIndexes();
	MPI_Group orig_group, new_group; 
	MPI_Comm_group(NETWORK_COMM_WORLD, &orig_group);

	MPI_Comm* new_comm = new MPI_Comm();
	//debug
	//*new_comm = NETWORK_COMM_WORLD;
	//end debug

	vector<int> nodeLayerIndexes =  m_population->GetNodeLayerIndexes();
	vector<int> mpiProcsUsed = m_population->MPIGetProcessesUsed();

	// currently not set here in this way but may be changed
	/*
	MPI_Group_incl(orig_group, nodeLayerIndexes.size(), &(nodeLayerIndexes[0]), &new_group);
	MPI_Comm_create(NETWORK_COMM_WORLD, new_group, new_comm);
	*/
	
	//m_mpiCommLayer = new_comm;
	m_mpiCommLayer = new MPI_Comm();
	
	if(mpiProcsUsed.size()==0) // all used
	{
		*m_mpiCommLayer = NETWORK_COMM_WORLD; // currently may result in bug if this assumption not correct and trying to record
		MPI_Comm_size(NETWORK_COMM_WORLD, &m_mpiSizeLocal);
		MPI_Comm_rank(NETWORK_COMM_WORLD, &m_mpiRankLocal);
	}
	else
	{
		m_mpiCommLayer = new_comm;
		MPI_Group_incl(orig_group, mpiProcsUsed.size(), &(mpiProcsUsed[0]), &new_group);
		MPI_Comm_create(NETWORK_COMM_WORLD, new_group, new_comm);
		if(binary_search(mpiProcsUsed.begin(),mpiProcsUsed.end(),m_population->network()->MPIGetNodeId()))
		{
			MPI_Comm_size(*m_mpiCommLayer, &m_mpiSizeLocal);
			MPI_Comm_rank(*m_mpiCommLayer, &m_mpiRankLocal);
		}
	}

	if(m_population->network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"done.\n";cout.flush();
	}

	m_commsLayersCreated = true;
}

vector<float> MPIDistribution::MPILayerReduceVariable(vector<float> data)
{
	MPI_Comm* comm = MPIGetCommLayer();

	vector<float> reduced(data.size());
	MPI_Allreduce(&data[0],&reduced[0],data.size(),MPI_FLOAT,MPI_SUM,*comm);

	return reduced;
}

// pop data after each transfer to reduce memory usage
vector<vector<float> > MPIDistribution::MPILayerReduceVariable(vector<vector<float> > data)
{
	MPI_Comm* comm = MPIGetCommLayer();

	vector<vector<float> > reduced(data.size());//(data.size(),vector<float>(data[0].size()));
	
	int size;
	//cout<<data.size();cout.flush();
	for(int i=0;i<data.size();i++)
	{
		vector<float> r(data[i].size());
		size = data[i].size();

		// Fix bug - should not need to use NETWORK_COMM_WORLD
		MPI_Allreduce(&data[i][0],&r[0],size,MPI_FLOAT,MPI_SUM,NETWORK_COMM_WORLD);//*comm);
		reduced[i] = r;
	}

	return reduced;
}

vector<vector<int> > MPIDistribution::MPILayerGatherVariables(vector<vector<int> > data)
{
	MPI_Comm* comm = MPIGetCommLayer(); // communicator for nodes in layer
	vector<vector<int> > gatheredData(data.size());
	
	for(int i=0;i<data.size();i++)
	{
		data[i].push_back(0); // will be removed later, allgatherv may not like zero-sized vectors
		int sizeLocal = data[i].size();
		int totalSize;
		MPI_Allreduce(&sizeLocal,&totalSize, 1,MPI_INT,MPI_SUM,*comm);

		vector<int> totalIndexes(totalSize);

		vector<int> rcounts(m_mpiSizeLocal);
		vector<int> displs(m_mpiSizeLocal);

		// counts from each node
		MPI_Allgather( &sizeLocal, 1, MPI_INT, &rcounts[0], 1, MPI_INT, *comm);

		for(int j=1;j<rcounts.size();j++)
		{
			displs[j] = rcounts[j-1]+displs[j-1];
		}

		MPI_Allgatherv( &data[i][0],data[i].size(),MPI_INT,&totalIndexes[0],&rcounts[0],&displs[0],MPI_INT,*comm);

		// remove all artificial values
		for(int j=rcounts.size()-1;j>-1;j--)
		{
			totalIndexes.erase(totalIndexes.begin() + displs[j] + rcounts[j] - 1);
		}

		gatheredData[i] = totalIndexes;
	}

	return gatheredData;
}

void Simulate(); // overridden from PopulationModifier and called every simulation time step
	vector<float> wta(vector<float> data); // the actual winner-take-all function,
												// which takes the input vector data and
												//  returns a vector with a 1 on the position 
												//  of the highest value in data and 0 on 
												//  all other positions

/// <summary>	 Makes all activity values in population (fullLayer = true) or columns (fullLayer = false)
///				available to the other processes also simulating parts of this population/column. </summary>
///
/// <param name="fullLayer">	true to make all activity values of the population local. </param>

void MPIDistribution::MakeActivityValuesLocal(bool fullLayer)
{
	if(m_commsHCsCreated == false && fullLayer == false)
	{
		this->MPICreateCommsHypercolumns();
	}
	else if(m_commsLayersCreated == false && fullLayer == true)
	{
		this->MPICreateCommLayer();	
	}

	vector<int> localHypercolumns = ((PopulationColumns*)m_population)->GetLocalHypercolumnIndexes();
	vector<Unit*> units = ((PopulationColumns*)m_population)->GetUnits();

	int totMcs =  ((PopulationColumns*)m_population)->GetNrUnitsTotal();

	vector<float> data;
	MPI_Comm* comm;

	bool allLocal = true;

	if(fullLayer == true && localHypercolumns.size()>0)
	{
		data = vector<float>(totMcs);
	
		//
		for(int i=0;i<units.size();i++)
			data[units[i]->GetUnitIdLocal()] = units[i]->GetValue();

		vector<float> recData(data.size());
		
		if(recData.size()>0)
		{
			comm = m_mpiCommLayer;
			MPI_Allreduce(&data[0],&recData[0],data.size(),MPI_FLOAT,MPI_SUM,*comm);

			m_population->SetValuesBuffer(recData);			
		}
	}
	else
	{
		for(int i=0;i<localHypercolumns.size();i++)
		{
			Hypercolumn* h = (Hypercolumn*)((PopulationColumns*)m_population)->GetHypercolumn(i);
			vector<RateUnit*> mcs = h->GetRateUnits();

			if(fullLayer == false)
			{
				vector<MPI_Comm*> mpiComms = MPIGetCommHCs();
				if(mpiComms.size()>0)
				{
					if(mpiComms.size()!=1)
					{
						//cerr<<"Assumption about locality of hypercolumns not correct!";
					}

					//assuming only one non-local hypercolumn on process
					comm = mpiComms[localHypercolumns[i]];
					data = vector<float>(h->GetTotalNrRateUnits());
					allLocal = true;
				}
			}
			else
			{
				data = vector<float>(totMcs);//mcsIndexes.size());
			}

			for(int j=0;j<mcs.size();j++)
			{
				data[mcs[j]->GetUnitIdLocalInHypercolumn()] = mcs[j]->GetValue();
			}

			if(units.size() < h->GetTotalNrRateUnits())
				allLocal = false;

			if(allLocal == false && fullLayer == false) // all local, no need to update
			{
				vector<float> recData(data.size());
				// using allreduce
				MPI_Allreduce(&data[0],&recData[0],data.size(),MPI_FLOAT,MPI_SUM,*comm);
				
				// using allgatherv
				/*for(int j=1;j<rcounts.size();j++)
				{
					displs[j] = rcounts[j-1]+displs[j-1];
				}

				MPI_Allgatherv( &data[i][0],data[i].size(),MPI_INT,&totalIndexes[0],&rcounts[0],&displs[0],MPI_INT,*comm);*/
				
				h->SetValues(recData);
			}
			else if(allLocal == true)
				h->SetValues(data);
		}
	}
}

// All nodes part of the layer gather the unit values
void MPIDistribution::MPIMakeLayerValuesLocal()
{
	this->MakeActivityValuesLocal(true);	
}

// All nodes part of the hypercolumns gather the unit values
void MPIDistribution::MPIMakeHypercolumnsValuesLocal()
{
	this->MakeActivityValuesLocal(false);	
}

void MPIDistribution::MPICreateCommsHypercolumns()
{
	// Create MPI handles for the hypercolumns that are non-locally distributed (used in case of send/receive values within a hypercolumn if minicolumns are on different processes)
	if(m_population->network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"Creating MPI Communicator for hypercolumns...";cout.flush();
	}

	MPI_Group orig_group, new_group;
	MPI_Comm_group(NETWORK_COMM_WORLD, &orig_group);

	vector<vector<int> > nodeHcIndexes = ((PopulationColumns*)m_population)->GetAllNodeIndexes();

	// TODO: check we do not have a memory leak for multiple runs and ok to do like this
	for(int i=0;i<m_mpiCommHCs.size();i++)
		delete m_mpiCommHCs[i];

	m_mpiCommHCs.clear();

	// this will currently create the extra communicators needed for shared hypercolumns on all processes
	// - only if more than one process uses them it will be created

	for(int i=0;i<nodeHcIndexes.size();i++)
	{
		if(nodeHcIndexes[i].size()>1)
		{
			MPI_Comm* new_comm = new MPI_Comm();
			MPI_Group_incl(orig_group, nodeHcIndexes[i].size(), &(nodeHcIndexes[i][0]), &new_group);
			MPI_Comm_create(NETWORK_COMM_WORLD, new_group, new_comm);

			m_mpiCommHCs.push_back(new_comm);
		}
		else
			m_mpiCommHCs.push_back(NULL);
	}

	if(m_population->network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"done. (commHCs=="<<m_mpiCommHCs.size()<<", sizeof=="<<sizeof(MPI_Comm)*m_mpiCommHCs.size()<<").\n"; cout.flush();
	}

	m_commsHCsCreated = true;
}

/// <summary>	Generates a division of a columnar populatin based on which parallelization scheme has been selected.
/// 			TODO: Simplify the way to extend to new parallelization schemes. </summary>
///
/// <param name="nrRateUnits">				Nr of rate units in each column. </param>
/// <param name="parallelizationScheme">	The parallelization scheme. </param>
/// <param name="mpiRank">					This mpi process. </param>
/// <param name="mpiSize">					Nr total processes. </param>
///
/// <returns>	Distribution. </returns>

vector<vector<vector<long> > > MPIDistribution::DivideEqualByHypercolumn(vector<int> nrRateUnits, ParallelizationSchemeLayer parallelizationScheme, int mpiRank,int mpiSize)
{
	vector<vector<vector<long> > > out;
	long nrUnits = nrRateUnits.size();

	float part = (float)nrUnits/(float)mpiSize;
	vector<vector<long> > hcNodeIndexes(nrRateUnits.size());

	if(part<1.0 || parallelizationScheme == this->ParallelizationRateUnits) // if we cannot divide by hypercolumn, divide by minicolumns
	{
		nrUnits = 0;
		for(unsigned int i=0;i<nrRateUnits.size();i++)
			nrUnits+=nrRateUnits[i];

		vector<vector<long> > div = DivideEqualByUnitsSingleHC(nrRateUnits,mpiRank,mpiSize);//DivideEqualByUnits(nrUnits,mpiRank,mpiSize);

		int startPos = 0;
		int currentUnitPos = 0;

		for(int i=0;i<nrRateUnits.size();i++) // make local?
		{
			vector<long> nIndexes; // int

			for(int j=startPos;j<div.size();j++)
			{
				if(div[j][0] != div[j][1]) // will be if more processes than units
				{
					if(div[j][0]<currentUnitPos+nrRateUnits[i] && div[j][1] > currentUnitPos)
						nIndexes.push_back(j);
					else if(div[j][1] == currentUnitPos)
					{
						// not divided hypercolumn, take one step forward
					}
					else
					{
						startPos = j-1;
						break;
					}
				}
			}

			currentUnitPos+=nrRateUnits[i];
			
			/*bool onThisProc = false;
			for(int j=0;j<nIndexes.size();j++)
			{
				if(nIndexes[j] == mpiRank)
				{
					onThisProc = true;
					break;
				}
			}

			if(onThisProc==true)*/
				hcNodeIndexes[i] = nIndexes;
		}

		out.push_back(hcNodeIndexes);
		div[0] = div[mpiRank]; if(div.size()>1) div.erase(div.begin()+1,div.end());

		out.push_back(div);
		return out;
	}

	vector<vector<long> > indexesInterval;

	// change to only do locally
	int i=mpiRank;
	//for(int i=0;i<mpiSize;i++)
	//{
		long lb = (long)(part * (float)i);
		long ub = (long)(part * (float)(i+1));

		if(i == mpiSize-1)
			ub = nrUnits;

		for(long m=lb;m<ub;m++)
			hcNodeIndexes[m].push_back(i);

		long tlb = 0;
		long tub = 0;

		for(int i=0;i<lb;i++)
			tlb+=nrRateUnits[i];

		for(int i=0;i<ub;i++)
			tub+=nrRateUnits[i];

		vector<long> v;
		v.push_back(tlb);
		v.push_back(tub);

		indexesInterval.push_back(v);
	//}

	out.push_back(hcNodeIndexes);
	out.push_back(indexesInterval);

	return out;//indexesInterval;
}


vector<vector<vector<long> > > MPIDistribution::DivideEqualByTotalUnits(vector<int> nrRateUnits, long currentUnit, long totalNrUnits, int mpiRank,int mpiSize)
{
	long nrUnits = 0;
	for(int i=0;i<nrRateUnits.size();i++)
	{
		nrUnits+=nrRateUnits[i];
	}

	vector<vector<vector<long> > > out;

	float part = (float)totalNrUnits/(float)mpiSize;
	vector<vector<long> > nodeIndexes(nrRateUnits.size());

	// does not happen
	
	vector<vector<long> > indexesInterval;

	// change to only do locally
	int i=mpiRank;
	//for(int i=0;i<mpiSize;i++)
	//{
		long lb = (long)(part * (float)i);
		long ub = (long)(part * (float)(i+1));
		if(i == mpiSize-1)
			ub = totalNrUnits;

		bool cont = false;
		/*if(lb<currentUnit+nrUnits && lb>=currentUnit)
			cont = true;
		if(ub<currentUnit+nrUnits && ub>=currentUnit)
			cont = true;
			*/
		if(lb<currentUnit+nrUnits && ub>=currentUnit)
			cont = true;

		if(cont == true) // this process should take part in this population
		{
			long currentPartUnit = 0;

			for(int m=0;m<nrRateUnits.size();m++)
			{
				for(int n=0;n<nrRateUnits[m];n++)
				{
					long unitNr = currentUnit+currentPartUnit;
					if(lb>=unitNr && ub>unitNr)
					{
						nodeIndexes[m].push_back(i);
						break;
					}
						currentPartUnit++;
				}
			}

			//for(long m=lb;m<ub;m++)
			//	nodeIndexes[m].push_back(i);

			// get local locations
			long tlb = 0;
			long tub = 0;

			if(lb<=currentUnit)
				tlb = 0;
			else
				tlb = lb-currentUnit;
				
			if(ub>currentUnit+nrUnits)
				tub = currentUnit+nrUnits;
			else
				tub = ub-currentUnit;

			vector<long> v;
			v.push_back(tlb);
			v.push_back(tub);

			indexesInterval.push_back(v);
		}
		else
		{
			vector<long> v;
			v.push_back(0);
			v.push_back(0);

			indexesInterval.push_back(v);
		}
	//}

	out.push_back(nodeIndexes);
	out.push_back(indexesInterval);

	return out;//indexesInterval;
}


// As divide equal by total units but does not allow one process to be in several populations if possible
vector<vector<vector<long> > > MPIDistribution::DivideEqualByPopulations(vector<int> nrRateUnits, long currentUnit, long totalNrUnits, int mpiRank,int mpiSize)
{
	long nrUnits = 0;
	for(int i=0;i<nrRateUnits.size();i++)
	{
		nrUnits+=nrRateUnits[i];
	}

	vector<vector<vector<long> > > out;

	float part = (float)totalNrUnits/(float)mpiSize;
	vector<vector<long> > nodeIndexes(nrRateUnits.size());

	// does not happen

	vector<vector<long> > indexesInterval;

	int i=mpiRank;
	long lb = (long)(part * (float)i);
	long ub = (long)(part * (float)(i+1));
	if(i == mpiSize-1)
		ub = totalNrUnits;

	bool cont = false;

	if(lb<currentUnit+nrUnits && ub>=currentUnit)
		cont = true;

	if(cont == true) // this process should take part in this population
	{
		long currentPartUnit = 0;

		for(int m=0;m<nrRateUnits.size();m++)
		{
			for(int n=0;n<nrRateUnits[m];n++)
			{
				long unitNr = currentUnit+currentPartUnit;
				if(lb>=unitNr && ub>unitNr)
				{
					nodeIndexes[m].push_back(i);
					break;
				}
				currentPartUnit++;
			}
		}

		// get local locations
		long tlb = 0;
		long tub = 0;

		if(lb<=currentUnit)
			tlb = 0;
		else
			tlb = lb-currentUnit;

		if(ub>currentUnit+nrUnits)
			tub = currentUnit+nrUnits;
		else
			tub = ub-currentUnit;

		vector<long> v;
		v.push_back(tlb);
		v.push_back(tub);

		indexesInterval.push_back(v);
	}
	else
	{
		vector<long> v;
		v.push_back(0);
		v.push_back(0);

		indexesInterval.push_back(v);
	}

	out.push_back(nodeIndexes);
	out.push_back(indexesInterval);

	return out;
}