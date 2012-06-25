#include <mpi.h>
#include "MPIDistribution.h"
#include "NetworkPopulation.h"

void MPIDistribution::MPICreateCommLayer()
{
	// Create MPI handle for entire network layer (used to gather data from all units in the layer)
	//MPI_Barrier(NETWORK_COMM_WORLD); //debug

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

	// currently debug off
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
	//MPI_Barrier(NETWORK_COMM_WORLD);
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

void MPIDistribution::MPIMakeValuesLocal(bool fullLayer)
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

	int totMcs =  ((PopulationColumns*)m_population)->GetNrUnitsTotal();//((PopulationColumns*)m_population)->GetRateUnits().size(); // slow

	vector<float> data;//mcsIndexes.size());
	MPI_Comm* comm;

	bool allLocal = true;

	if(fullLayer == true && localHypercolumns.size()>0)
	{
		data = vector<float>(totMcs);//mcsIndexes.size());
	
		//
		for(int i=0;i<units.size();i++)
			data[units[i]->GetUnitIdLocal()] = units[i]->GetValue();//((RateUnit*)units[i])->GetSubThresholdValue();//

		vector<float> recData(data.size());
		
		if(recData.size()>0)
		{
			comm = m_mpiCommLayer;
			MPI_Allreduce(&data[0],&recData[0],data.size(),MPI_FLOAT,MPI_SUM,*comm);

			m_population->SetValuesBuffer(recData);
			/*for(int j=0;j<units.size();j++)//data.size();j++)
			{
				((RateUnit*)units[j])->SetValue(recData[units[j]->GetUnitIdLocal()]);//recData[j]);//mcsIndexes[j]])->SetValue(recData[j]);
			}*/
		}
	}
	else
	{
		for(int i=0;i<localHypercolumns.size();i++)
		{
			//vector<int> mcsIndexes = ((PopulationColumns*)m_population)->GetRateUnitsIndexes(localHypercolumns[i]);
			Hypercolumn* h = (Hypercolumn*)((PopulationColumns*)m_population)->GetHypercolumn(i);
			vector<RateUnit*> mcs = h->GetRateUnits();

			if(fullLayer == false)
			{
				vector<MPI_Comm*> mpiComms = MPIGetCommHCs();
				if(mpiComms.size()>0)
				{
					if(mpiComms.size()!=1)
					{
						// (!)
						//cerr<<"W";//cerr<<"Assumption about locality of hypercolumns not correct!";
					}

					//assuming only one non-local hypercolumn on process
					comm = mpiComms[localHypercolumns[i]];//mpiComms[i];//mpiComms[0];//(MPIGetCommHCs())[localHypercolumns[i]];
					data = vector<float>(h->GetTotalNrRateUnits());//mcs.size());//mcsIndexes.size());//totMcs);
					allLocal = true;
				}
			}
			else
			{
				data = vector<float>(totMcs);//mcsIndexes.size());
			}
			

			//Hypercolumn* h = ((PopulationColumns*)m_population)->GetHypercolumn(i);
			//vector<RateUnit*> mcs = h->GetRateUnits();

			for(int j=0;j<mcs.size();j++)//mcsIndexes.size();j++)
			{
				data[mcs[j]->GetUnitIdLocalInHypercolumn()] = mcs[j]->GetValue();

				//if(units[j]->IsLocal() == true)//mcsIndexes[j]]->IsLocal() == true)
				//	data[((RateUnit*)units[j])->GetUnitIdLocalInHypercolumn()] = ((Unit*)units[j])->GetValue();//data[((RateUnit*)units[j])->GetUnitIdLocal()] = ((Unit*)units[j])->GetValue();//mcsIndexes[j]])->GetValue();
				//else
				//{
					//units[mcsIndexes[j]]->SetLocal(true); // will retrieve data
//					allLocal = false;
//				}
			}

			if(units.size() < h->GetTotalNrRateUnits())//mcsIndexes.size())
				allLocal = false;

			if(allLocal == false && fullLayer == false) // all local, no need to update
			{
				vector<float> recData(data.size());
				MPI_Allreduce(&data[0],&recData[0],data.size(),MPI_FLOAT,MPI_SUM,*comm);

				//Hypercolumn* h = ((PopulationColumns*)m_population)->GetHypercolumn(i);
				h->SetValues(recData);//values);
				/*for(int j=0;j<mcsIndexes.size();j++)
				{
				((RateUnit*)units[mcsIndexes[j]])->SetValue(recData[mcsIndexes[j]]);//mcsIndexes[j]])->SetValue(recData[j]);
				}*/
			}
			else if(allLocal == true)
				h->SetValues(data);

			//vector<int> nodeIndexes = GetNodeIndexes(localHypercolumns[i]); // other nodes utilizing this hypercolumn
		}
	}
}

// All nodes part of the layer gather the unit values
void MPIDistribution::MPIMakeLayerValuesLocal()
{
	this->MPIMakeValuesLocal(true);	
}

// not distr atm. - create fcn using mpi min instead
// assuming one hypercolumn in layer (or first hypercolumn in layer)
/*vector<long> MPIDistribution::MPIGetMaxIdInHypercolumns()
{	
	MPIMakeHypercolumnsValuesLocal();

	vector<int> localHypercolumns = ((PopulationColumns*)m_population)->GetLocalHypercolumnIndexes();
	vector<long> ids(localHypercolumns.size());

	for(int k=0;k<localHypercolumns.size();k++)
	{
//		vector<int> mcsIndexes = ((PopulationColumns*)m_population)->GetRateUnitsIndexes(localHypercolumns[k]);

		//vector<int> mcsIndexes = ((PopulationColumns*)m_population)->GetRateUnitsIndexes(0); // first hc
		vector<double> data(mcsIndexes.size());
		vector<Unit*> units = ((PopulationColumns*)m_population)->GetUnits(); // make pointer

		for(int j=mcsIndexes.size()-1;j>-1;j--)
		{
			data[j] = ((RateUnit*)units[mcsIndexes[j]])->GetValue();
		}

		float maxVal = -1e8;
		int maxIndex = -1;
		for(int i=0;i<data.size();i++)
		{
			if(data[i]>maxVal)
			{
				maxVal = data[i];
				maxIndex = i;
			}
		}

		ids[k] = units[mcsIndexes[maxIndex]]->GetUnitId();
	}

	return ids;
}*/


// All nodes part of the hypercolumns gather the unit values
void MPIDistribution::MPIMakeHypercolumnsValuesLocal()
{
	this->MPIMakeValuesLocal(false);	
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

	//for(int i=0;i<m_mpiCommHCs.size();i++)
//		delete m_mpiCommHCs[i];

	// watch memory leak (!)
	for(int i=0;i<m_mpiCommHCs.size();i++)
		delete m_mpiCommHCs[i];

	m_mpiCommHCs.clear();

	// this will currently create the extra communicators needed for shared hypercolumns on all processes
	// - change so that they will only be created on the processes that will use them (!)

	for(int i=0;i<nodeHcIndexes.size();i++)
	{
		if(nodeHcIndexes[i].size()>1)
		{
			// all processes currently need to go in here and all need to know the nodeHcIndexes (global) list - could be changed (!)
//			for(int j=0;j<nodeHcIndexes[i].size();j++)
//			{
				

				//if(nodeHcIndexes[i][j] == this->m_population->network()->MPIGetNodeId())
				//{
					MPI_Comm* new_comm = new MPI_Comm();
					MPI_Group_incl(orig_group, nodeHcIndexes[i].size(), &(nodeHcIndexes[i][0]), &new_group);
					MPI_Comm_create(NETWORK_COMM_WORLD, new_group, new_comm);

					m_mpiCommHCs.push_back(new_comm);
				//}
//			}
		}
		else
			m_mpiCommHCs.push_back(NULL);
	}

	/*for(int i=0;i<nodeHcIndexes.size();i++)
	{
		MPI_Comm* new_comm = new MPI_Comm();
		MPI_Group_incl(orig_group, nodeHcIndexes[i].size(), &(nodeHcIndexes[i][0]), &new_group);
		MPI_Comm_create(NETWORK_COMM_WORLD, new_group, new_comm);

		m_mpiCommHCs.push_back(new_comm);
	}*/

	if(m_population->network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"done. (commHCs=="<<m_mpiCommHCs.size()<<", sizeof=="<<sizeof(MPI_Comm)*m_mpiCommHCs.size()<<").\n"; cout.flush();
	}

	m_commsHCsCreated = true;
}

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