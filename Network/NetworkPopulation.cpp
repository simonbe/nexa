#include <math.h>
#include "Network.h"
#include "NetworkPopulation.h"
#include "MPIDistribution.h"
#include "NetworkUnitsIF.h"


Population::Population() 
{
	m_mpiDistribution = new MPIDistribution(this);
	m_keepValues = false;

	/*stringstream ss;
	ss<<"Layer"<<this->GetLayerId();
	m_name = ss.str();*/
	m_parallelizationScheme = MPIDistribution::ParallelizationDefault;
}

void Population::ResetLocalities()
{
	TimingStart("ResetLocalities");
	/*for(int i=0;i<m_units.size();i++)
	{
		if(m_units[i]->GetUnitIdLocal() >= m_localUnitIndexesInterval[0] &&
			m_units[i]->GetUnitIdLocal() < m_localUnitIndexesInterval[1])
			m_units[i]->SetLocal(true);
		else
			m_units[i]->SetLocal(false);

		m_units[i]->IsNewEvent(true); // ok?
	}*/
	for(int i=0;i<m_units.size();i++)
	{
		m_units[i]->SetLocal(true);	
		m_units[i]->IsNewEvent(true); // ok?
	}

	TimingStop("ResetLocalities");
}

Population::~Population()
{
	for(int i=0;i<m_units.size();i++)
		delete m_units[i];

	for(int i=0;i<m_connectionIncoming.size();i++)
	{
		if(m_connectionIncoming[i] != 0)
		{
			delete m_connectionIncoming[i];
			m_connectionIncoming[i] = NULL;
		}
	}

	// network now responsible for this deletion
	/*for(int i=0;i<m_eventLayers.size();i++)
	{
		if(m_eventLayers[i]!=0)
		{
			delete m_eventLayers[i];
			m_eventLayers[i] = 0;
		}
	}*/ 

	for(int i=0;i<m_preConnectivitys.size();i++)
	{
		if(m_preConnectivitys[i] != NULL)
		{
			delete m_preConnectivitys[i];
			m_preConnectivitys[i] = NULL;
		}
	}
}

void Population::Modify()
{
//	TimingStart(m_name);
//	TimingStart("Modify");

	// bcpnn etc
	//for(int i=0;i<m_connectionOutgoing.size();i++)
	for(int i=0;i<m_connectionIncoming.size();i++)
	{
		m_connectionIncoming[i]->ModifyConnection();
	//	m_connectionOutgoing[i]->ModifyConnection();
	}

	for(int i=0;i<m_eventLayers.size();i++)
		m_eventLayers[i]->Modify(); // typically seldom utilized

//	TimingStop("Modify");
//	TimingStop(m_name);
}

void PopulationColumns::ResetLocalities()
{
	TimingStart("ResetLocalities");
	// minicolumns

	for(int i=0;i<m_minicolumns.size();i++)
	{
/*		if(m_minicolumns[i]->GetUnitIdLocal() >= m_localUnitIndexesInterval[0] &&
			m_minicolumns[i]->GetUnitIdLocal() < m_localUnitIndexesInterval[1])
			m_minicolumns[i]->SetLocal(true);
		else
			m_minicolumns[i]->SetLocal(false);
*/
		m_minicolumns[i]->SetLocal(true);
		m_minicolumns[i]->IsNewEvent(true); // ok?
	}

	// hypercolumns - all by default local now
	for(int i=0;i<m_hypercolumns.size();i++)
	{
		m_hypercolumns[i]->SetLocal(true);//false);
		m_hypercolumns[i]->IsNewEvent(true); // ok?
	}

	/*for(int i=0;i<m_localHypercolumnIndexes.size();i++)
	{
		m_hypercolumns[m_localHypercolumnIndexes[i]]->SetLocal(true);
	}*/

	TimingStop("ResetLocalities");
}

// not used in this way atm.
Population::Population(Network* net, unsigned long nrUnits)
{
	/*m_networkType = Population::FixedWeightsAndRateUnits;
	m_nrUnits = nrUnits;
	m_network = net;

	m_unitIndexesInterval = m_mpiDistribution->DivideEqual(nrUnits,m_network->MPIGetNodeId(),m_network->MPIGetNrProcs());
	m_localUnitIndexesInterval = m_unitIndexesInterval[m_network->MPIGetNodeId()];
	*/
}

// Will add a units property for all incoming connections for layer
void Population::AddUnitsPropertyToInitialize(UnitModifier* p)
{
	if(m_units.size()>0)	// called after initialize
	{
		for(int i=0;i<m_units.size();i++)
		{
			m_units[i]->AddUnitModifier(p);
		}
	}
	else					// called before initialize
	{
		m_unitPropertiesToInitialize.push_back(p);
	}
}


/*void Population::AddUnitsPropertyIndependent(UnitModifier* p)
{
	for(int i=0;i<m_units.size();i++)
	{
		m_units[i]->AddUnitModifier(p);
	}
}*/

PopulationColumns::PopulationColumns(Network* net, unsigned long nrHypercolumns, unsigned long nrRateUnits, UnitType unitType, MPIDistribution::ParallelizationSchemeLayer parallelizationScheme, bool useSilentHypercolumns, float silentHypercolumnsThreshold)
{
	m_network = net;
	m_keepValues = false;

	m_unitType = unitType;
	m_parallelizationScheme = parallelizationScheme;
	m_useSilentHypercolumns = useSilentHypercolumns;
	m_silentHypercolumnsThreshold = silentHypercolumnsThreshold;
	m_firstRun = true;
	m_initialized = false;
	m_nrUnits = 0;
	m_nrRateUnits.clear();
	m_localHypercolumnIndexes.clear();
	m_localRateUnits.clear();
	m_mcHcIndexes.clear();
	m_nodeHcIndexes.clear();
	m_nodeHcIds.clear();

	if(parallelizationScheme == MPIDistribution::ParallelizationDefault || parallelizationScheme == MPIDistribution::ParallelizationSplitPopulations)
	{
		m_network->MPIAddNrUnitsParallelizationDefault(nrHypercolumns*nrRateUnits);
	}

	m_nrHypercolumns = nrHypercolumns;
	m_nrRateUnitsInHypercolumn = nrRateUnits;	
}

// not guaranteed to scale properly - change loop around m_nrRateUnits.size()...
void PopulationColumns::Initialize()
{
	///////// Prevously put in constructor

	double timeStart;
	if(this->network()->MPIGetNodeId() == 0)
	{
		timeStart = MPI_Wtime();
		cout<<"Setting up object PopulationColumns...";
		cout.flush();
	}

	vector<int> mcsInHc;

	for(int i=0;i<m_nrHypercolumns;i++)
	{
		vector<int> mcIndexes;
		/*for(int j=m_nrUnits;j<m_nrUnits+nrRateUnits;j++)
		{
			mcIndexes.push_back(j);
		}*/

		//m_mcHcIndexes.push_back(mcIndexes);

		m_nrRateUnits.push_back(m_nrRateUnitsInHypercolumn);
		m_nrUnits+=m_nrRateUnitsInHypercolumn;
		mcsInHc.push_back(m_nrRateUnitsInHypercolumn);
	}

	m_populationType = Population::FixedWeightsAndRateUnits;
	TimingStart("LayerInitDiv");
	vector<vector<vector<long> > > div;

	if(this->m_parallelizationScheme == MPIDistribution::ParallelizationDefault || this->m_parallelizationScheme == MPIDistribution::ParallelizationSplitPopulations)
	{
		div = m_mpiDistribution->DivideEqualByTotalUnits(mcsInHc,network()->MPIGetCurrentUnitParallelizationDefault(),network()->MPIGetNrUnitsParallelizationDefault(),network()->MPIGetNodeId(),network()->MPIGetNrProcs());
		network()->MPIAddCurrentUnitParallelizationDefault(m_nrUnits);

		vector<int> allUseProcess(network()->MPIGetNrProcs());
		int useProcess = 1;
		if( div[1][0][0] == 0 &&  div[1][0][1] == 0)
			useProcess = 0;

		MPI_Allgather(&useProcess,1,MPI_INT,&allUseProcess[0],1,MPI_INT,NETWORK_COMM_WORLD);

		for(int i=0;i<allUseProcess.size();i++)
		{
			if(allUseProcess[i] == 1)
				m_mpiProcessesUsed.push_back(i);
		}
	}
	else // DivideByHypercolumns
		div = m_mpiDistribution->DivideEqualByHypercolumn(mcsInHc,m_parallelizationScheme,network()->MPIGetNodeId(),network()->MPIGetNrProcs());//m_nrUnits,network()->MPIGetNodeId(),network()->MPIGetNrProcs());//(m_nrRateUnits,m_network->MPIGetNodeId(),m_network->MPIGetNrProcs());
	

	TimingStop("LayerInitDiv");
	
	m_nodeHcIndexes = div[0];
	//m_unitIndexesInterval = div[1];
	m_localUnitIndexesInterval = div[1][0];
	div.clear();
	mcsInHc.clear();
	//m_localUnitIndexesInterval = m_unitIndexesInterval[network()->MPIGetNodeId()];

	//m_hashMcHcIndexes

	int currentTotUnit = 0;
	vector<int> currentNodeIndexesLayer;
	int lastMaxUsedInterval = 0;

	currentTotUnit = m_localUnitIndexesInterval[0];//m_unitIndexesInterval[m_network->MPIGetNodeId()][0];

	// assuming all processes are involved in this layer
	//for(int i=0;i<m_network->MPIGetNrProcs();i++)
	//	m_nodeLayerIndexes.push_back(i);
	TimingStart("LayerInitHcIndexes");
	for(int i=0;i<m_nodeHcIndexes.size();i++)
	{
		for(int j=0;j<m_nodeHcIndexes[i].size();j++)
		{
			bool exist = false;
			for(int k=0;k<m_nodeLayerIndexes.size();k++)
			{
				if(m_nodeLayerIndexes[k] == m_nodeHcIndexes[i][j])
				{
					exist = true;
					break;
				}
			}

			if(exist == false)
				m_nodeLayerIndexes.push_back(m_nodeHcIndexes[i][j]);
		}
	}

	TimingStop("LayerInitHcIndexes");

	//for(int i=m_unitIndexesInterval[m_network->MPIGetNodeId()][0];i<m_unitIndexesInterval[m_network->MPIGetNodeId()][1];i++)
	TimingStart("LayerInitLocalUnitIndexes");
	for(int i=m_localUnitIndexesInterval[0];i<m_localUnitIndexesInterval[1];i++)
	{
		int cHc = i/m_nrRateUnitsInHypercolumn;

		if(find(m_localHypercolumnIndexes.begin(),m_localHypercolumnIndexes.end(),cHc) == m_localHypercolumnIndexes.end())
			m_localHypercolumnIndexes.push_back(cHc);
	}
	
	sort(m_localHypercolumnIndexes.begin(), m_localHypercolumnIndexes.end()); // not needed with current division
	TimingStop("LayerInitLocalUnitIndexes");
	
	//m_nodeHcIndexes.clear(); // fix, change logic
	m_unitIndexesInterval.clear(); // fix, change logic

	if(m_network->MPIGetNodeId() == 0)
	{
		cout<<"\n(process 0) LocalHcIndex: "<<m_localHypercolumnIndexes.size()<<", localHcIndexes: "<<m_localHypercolumnIndexes.size() <<", nodeLayerIndexes: "<<m_nodeLayerIndexes.size()<<".\n";
		cout.flush();
	}
	/*
	for(int i=0;i<nrHypercolumns;i++)
	//for(int i=m_unitIndexesInterval[m_network->MPIGetNodeId()][0];i<m_unitIndexesInterval[m_network->MPIGetNodeId()][1];i++)
	{
		bool contains = false;
		int currentUnit = currentTotUnit;
		int nrRateUnitsLocal = 0;
		
		int itemsErased = 0;

		vector<int> currentNodeIndexes;
		vector<vector<long> > unitIntervalsTemp = m_unitIndexesInterval;
		int usedMaxInterval = 0;
		
		///// speedup
		for(int k=lastMaxUsedInterval;k<unitIntervalsTemp.size();k++)
		{
			if(currentUnit>= unitIntervalsTemp[k][0] &&
				currentUnit< unitIntervalsTemp[k][1])
			{
				usedMaxInterval = k;
				k=unitIntervalsTemp.size();
			}
		}

		if(usedMaxInterval>0)
		{
			unitIntervalsTemp.erase(unitIntervalsTemp.begin(),unitIntervalsTemp.begin()+usedMaxInterval);
			lastMaxUsedInterval += usedMaxInterval;
			itemsErased+=usedMaxInterval;
		}

		///// end 
		usedMaxInterval = 0;

		for(int j=0;j<nrRateUnits;j++)
		{
			//for(int k=0;k<m_unitIndexesInterval.size();k++) // does not scale (!)
			for(int k=0;k<unitIntervalsTemp.size();k++)
			{
				if(currentUnit>= unitIntervalsTemp[k][0] &&
					currentUnit< unitIntervalsTemp[k][1])
				{
					int kpe = k + itemsErased;
					if(kpe==network()->MPIGetNodeId())
						contains = true;

					bool alreadyContains = false;
					bool alreadyContainsLayer = false;

					for(int m=0;m<currentNodeIndexes.size();m++)
					{
						if(currentNodeIndexes[m] == kpe)
							alreadyContains = true;
					}
					for(int m=0;m<currentNodeIndexesLayer.size();m++)
					{
						if(currentNodeIndexesLayer[m] == kpe)
							alreadyContainsLayer = true;
					}

					if(alreadyContains == false)
						currentNodeIndexes.push_back(kpe);
					if(alreadyContainsLayer == false)
						currentNodeIndexesLayer.push_back(kpe);

					nrRateUnitsLocal++;
					usedMaxInterval = k;
					k = unitIntervalsTemp.size();
					break;
				}
			}

			if(usedMaxInterval>0)
			{
				unitIntervalsTemp.erase(unitIntervalsTemp.begin(),unitIntervalsTemp.begin()+usedMaxInterval);
				itemsErased+=usedMaxInterval;
			}

			currentUnit++;
		}

		if(nrRateUnitsLocal>0 && nrRateUnitsLocal != nrRateUnits) // this hypercolumn is split on multiple processes (not default)
		{
			m_nodeHcIndexes.push_back(currentNodeIndexes);
			//m_nodeHcIds.push_back(
		}

		if(contains == true)
			m_localHypercolumnIndexes.push_back(i); // used ? yes

		currentTotUnit+=nrRateUnits;
	}

	m_nodeLayerIndexes = currentNodeIndexesLayer;*/

	if(this->network()->MPIGetNodeId() == 0)
	{
		cout<<"ok. (Time (process 0): "<<MPI_Wtime()-timeStart<<")\n";
		cout.flush();
	}



	/// end previously in constructor
	////////////////////////////////////////////////////




	//// Initialize the individual units



	if(this->network()->MPIGetNodeId() == 0)
	{
		cout<<"Initializing local hypercolumns, ";
		cout.flush();
	}

	if(m_initialized == false || m_initialized == true) // re-initializing (as a result of added populations etc) should only initialize certain parts, so may need to change
	{
		long localId = 0;

		TimingStart("HypercolumnsInit");

		for(long l=0;l<m_nrRateUnits.size();l++)
		{
			bool makeLocal = false;

			for(int i=0;i<m_nrRateUnits[l];i++)
			{
				int lId = localId+i;
				if(lId >= m_localUnitIndexesInterval[0] && lId < m_localUnitIndexesInterval[1])
				{
					makeLocal = true;
					break;
				}
			}

			long unitId = network()->GetNextUnitId();
			if(l==0)
				m_unitsStartId = unitId + m_nrRateUnits.size(); // assumes the startid starts first with the minicolumns and there is a strict ordering

			int hcId = network()->GetNextHypercolumnId(); // not used anymore?

			if(makeLocal == true)
			{
				Hypercolumn* h;
				
				//m_localHypercolumnIds.push_back(unitId);
				if(this->m_useSilentHypercolumns == true)
					h = new Hypercolumn(this->m_silentHypercolumnsThreshold);
				else
					h = new Hypercolumn();

				
				h->SetUnitId(unitId);//hcId);
				h->SetTotalNrRateUnits(m_nrRateUnits[l]);

				//h->SetNodeId(this->MPIBelongsToNode(localId)); // will be set to same node as the first minicolumn belongs to
				h->SetUnitIdLocal(l); // will be ordered first - may be problematic in vector (local) index based lookups
				h->network(this->network());
				
				this->AddUnit(h);

				localId+=m_nrRateUnits[l];
			}
			else
			{	
				localId+=m_nrRateUnits[l];
			}
		}

		TimingStop("HypercolumnsInit");
		
		localId = 0;

		if(this->network()->MPIGetNodeId() == 0)
		{
			cout<<"minicolumns...";
			cout.flush();
		}

		int currentHc = 0;
			
		TimingStart("RateUnitsInit");

		if(this->network()->MPIGetNodeId() == 0)
		{
			cout<<"Interval process 0 = ["<<m_localUnitIndexesInterval[0]<<", "<<m_localUnitIndexesInterval[1]<<"] ";
			cout.flush();
		}
		
		for(long l=0;l<m_nrRateUnits.size();l++)
		{
			//h->SetUnitIdLocal(localId);
			//localId++;

			bool hcUsed = false;

			for(int i=0;i<m_nrRateUnits[l];i++)
			{
				int lId = localId+i;
				if(lId >= m_localUnitIndexesInterval[0] && lId < m_localUnitIndexesInterval[1])
				{
					hcUsed = true;
					break;
				}
			}

			for(int i=0;i<m_nrRateUnits[l];i++)
			{
				bool makeLocal = false;

				if(localId >= m_localUnitIndexesInterval[0] && localId < m_localUnitIndexesInterval[1])
					makeLocal = true;

				long id = network()->GetNextUnitId();

				if(l==0 && i==0)
					m_startRateUnitId = id;

				if(makeLocal == true)
				{
					Hypercolumn* h = this->GetHypercolumns()[currentHc];

					if(l==0 && i==0)
						m_startRateUnitId = id;

					// reduce to T_ etc

					if(m_unitType == PopulationColumns::Graded || m_unitType == PopulationColumns::GradedThresholded)
					{
						RateUnit* m;
						if(m_unitType == PopulationColumns::Graded)
							m = new RateUnit(false);
						else if(m_unitType == PopulationColumns::GradedThresholded)
							m = new RateUnit(true);
						//TransferBcpnnOnline* t = new TransferBcpnnOnline(); // transfer function

						m->SetUnitId(id);
						m->SetUnitIdLocal(localId);
						m->SetUnitIdLocalHypercolumn(i);
						m->SetHypercolumnId(h->GetUnitId());//hcId);
						this->SetUnitIdLocalId(id,localId);
						//m->AddUnitModifier(t);

						//m->SetNodeId(this->MPIBelongsToNode(localId));
						m->SetPopulation(this);
						m->network(this->network());

						if(m_useTrace == true)
						{
							m->SetTrace(true,m_traceStrength);
						}

						this->AddUnit(m);
						h->AddRateUnit(m);
						m->AddHypercolumn(h);
					}
					else if(m_unitType == PopulationColumns::adEIF)
					{
#if GSL_AVAILABLE==1
						UnitadEIF* m = new UnitadEIF();
						m->SetUnitId(id);//m_network->GetNextUnitId());
						m->SetUnitIdLocal(localId);
						m->SetUnitIdLocalHypercolumn(i);
						m->SetHypercolumnId(h->GetUnitId());
						//m->AddUnitModifier(t);
						this->SetUnitIdLocalId(id,localId);

						m->SetNodeId(this->MPIBelongsToNode(localId));
						m->Population(this);
						m->network(this->network());

						this->AddUnit(m);
						h->AddRateUnit(m);
						m->AddHypercolumn(h);
#endif
					}
					else if(m_unitType == PopulationColumns::IF)
					{
						UnitIF* m = new UnitIF();
						m->SetUnitId(id);//network()->GetNextUnitId());
						m->SetUnitIdLocal(localId);
						m->SetUnitIdLocalHypercolumn(i);
						m->SetHypercolumnId(h->GetUnitId());//hcId);
						this->SetUnitIdLocalId(id,localId);
						//m->AddUnitModifier(t);

						m->SetNodeId(this->MPIBelongsToNode(localId));
						m->SetPopulation(this);
						m->network(this->network());

						this->AddUnit(m);
						h->AddRateUnit(m);
						m->AddHypercolumn(h);
					}
				}

				localId++;
			}

			if(hcUsed == true)
				currentHc++;
		}

		TimingStop("RateUnitsInit");
		
		// could put in unit properties here

		ResetLocalities();

		for(int i=0;i<m_minicolumns.size();i++)
		{
			// global/non-connection-specific unit properties/transfer functions
			for(int j=0;j<m_unitPropertiesToInitialize.size();j++)
				m_minicolumns[i]->AddUnitModifier(m_unitPropertiesToInitialize[j]);
			
			if(m_minicolumns[i]->IsLocal())
				m_localRateUnits.push_back(m_minicolumns[i]);
		}
	}

	if(this->network()->MPIGetNodeId() == 0)
	{
		cout<<"ok ("<<this->GetRateUnits().size()<<" nr local mcs, "<<this->GetHypercolumns().size()<<" nr local hcs, "<<this->GetUnitsAll()->size()<<" tot units). ";
		cout.flush();
	}

	MPI_Barrier(NETWORK_COMM_WORLD); // for debug, remove (!)
}

void Population::InitializeParallelization()
{
/*	if(m_parallelizationScheme == MPIDistribution::ParallelizationDefault || m_parallelizationScheme == MPIDistribution::ParallelizationSplitPopulations)
	{
		vector<int> allUseProcess(network()->MPIGetNrProcs());
		int useProcess = 1;
		if(m_localUnitIndexesInterval[0] == 0 && m_localUnitIndexesInterval[1] == 0)
			useProcess = 0;

		MPI_Allgather(&useProcess,1,MPI_INT,&allUseProcess[0],1,MPI_INT,NETWORK_COMM_WORLD);

		for(int i=0;i<allUseProcess.size();i++)
		{
			if(allUseProcess[i] == 1)
				m_mpiProcessesUsed.push_back(i);
		}
	}*/
}

void PopulationColumns::InitializeConnectionsEventsAndParameters()
{
	
	for(int i=0;i<m_pres.size();i++)
	{
		m_preConnectivitys[i]->network(this->network());
		m_preConnectivitys[i]->Initialize();

		// Multiple parameters implementation; adds vector and functions to parameter class
		//m_preConnectivitys[i]->AddParameters();
	}

	for(int i=0;i<m_eventLayers.size();i++)
	{
		m_eventLayers[i]->network(this->network());
		m_eventLayers[i]->Initialize(this);
	}

	//GenerateHashTables();
	//if(m_initialized == false)
	//	MPI()->MPICreateCommLayer();

	m_initialized = true;
}

void Population::GenerateHashTables()
{
	if(network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"Generating hash tables (network layer)...";
	}

/*	for(int i=0;i<m_units.size();i++)
	{
		if(m_units[i]->IsLocal() == true)
		{
			vector<Connection*> preConns = this->GetIncomingConnections();//m_units[i]->GetPreConnections();

			if(preConns.size()>0)
			{
				for(int k=0;k<preConns.size();k++)
				{
					vector<long> preIds = preConns[k]->GetPreIds(m_units[i]->GetUnitId());

					//vector<Unit*> pres =	m_units[i]->GetPres();
					for(int j=0;j<preIds.size();j++)
					{
						long unitId = preIds[j];//pres[j]->GetUnitId();
						map<long,vector<Unit*> >::iterator itr;

						vector<Unit*> conUnits;
						if ( (itr = network()->HashIdPreUnit()->find(unitId)) != network()->HashIdPreUnit()->end())
							conUnits = (*network()->HashIdPreUnit())[unitId];//m_hashIdPreUnit[unitId];

						conUnits.push_back(m_units[i]);

						(*network()->HashIdPreUnit())[unitId] = conUnits;//m_hashIdPreUnit[unitId] = conUnits;
					}
				}
			}
		}
	}

	if(network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"done.\n";
	}

	*/
}

/*void PopulationColumns::SendAndReceiveVersionISend()
{
	// NOTE: The MPI parts will be moved into general functions in parent class, but calling the functions will be the same

	vector<int> postNodes = this->GetPostNodes();
	vector<int> preNodes = this->GetPreNodes();

	int maxIter = max(postNodes.size(),preNodes.size());

	for(int i=0;i<maxIter;i++)
	{
		// irecv 1 - 
		// a) could use MPI_ANY_SOURCE
		// b) could change communicator to specific group

		int tagAct = 1;
		int nrRecv = 0;

		MPI_Request mpiRequest1,mpiRequest2,mpiRequest3,mpiRequest1s,mpiRequest2s,mpiRequest3s;
		vector<long> recvIds;
		vector<float> recvValues;

		if(i<=preNodes.size())
		{
			if(preNodes[i]!=this->network()->MPIGetNodeId())
				MPI_Irecv(&nrRecv,1,MPI_INT,preNodes[i],tagAct,NETWORK_COMM_WORLD,&mpiRequest1);
		}

		// isend
		vector<UnitModifier*> sendEvents;
		int nrEvents;

		if(i<=postNodes.size() && postNodes.size()>0)
		{
			if(postNodes[i]!=this->network()->MPIGetNodeId())
			{
				for(int j=0;j<m_minicolumns.size();j++)
				{
					//if this is post node
					if(m_minicolumns[j]->IsPostNode(postNodes[i])) // not most efficient
					{
						sendEvents.push_back(m_minicolumns[j]->CreateEvent());
					}
				}

				tagAct = 1;
				nrEvents = sendEvents.size();
				MPI_Isend(&nrEvents,1,MPI_INT,postNodes[i],tagAct,NETWORK_COMM_WORLD,&mpiRequest1s);
			}
			else
			{
				// add local events
				for(int j=0;j<m_minicolumns.size();j++)
				{
					vector<Unit*> recvUnits = (*m_network->HashIdPreUnit())[m_minicolumns[j]->GetUnitId()];//m_hashIdPreUnit[m_minicolumns[j]->GetUnitId()];
					UnitModifier* e = m_minicolumns[j]->CreateEvent();
					m_allEventsIncoming.push_back(e);

					for(int k=0;k<recvUnits.size();k++)
					{
						recvUnits[k]->AddEventIncoming(e);	
					}
				}
			}
		}

		// irecv 2
		if(i<=preNodes.size() && preNodes.size()>0)
		{
			if(preNodes[i]!=this->network()->MPIGetNodeId())
			{
				MPI_Wait(&mpiRequest1,MPI_STATUS_IGNORE);

				recvIds = vector<long>(nrRecv);
				recvValues = vector<float>(nrRecv);
				tagAct = 2;
				MPI_Irecv(&recvIds[0],nrRecv,MPI_LONG,preNodes[i],tagAct,NETWORK_COMM_WORLD,&mpiRequest2);
				tagAct = 3;
				MPI_Irecv(&recvValues[0],nrRecv,MPI_FLOAT,preNodes[i],tagAct,NETWORK_COMM_WORLD,&mpiRequest3);
			}
		}
		
		// isend 2
		if(i<=postNodes.size() && postNodes.size()>0)
		{
			if(postNodes[i]!=this->network()->MPIGetNodeId())
			{
				tagAct = 2;
				vector<long> ids;
				vector<float> values;
				for(int j=0;j<sendEvents.size();j++)
				{
					ids.push_back(sendEvents[j]->GetFromUnitId());
					values.push_back(sendEvents[j]->GetValue());
					delete sendEvents[j];
				}

				sendEvents.clear();

				MPI_Send(&ids[0],nrEvents,MPI_LONG,postNodes[i],tagAct,NETWORK_COMM_WORLD);
				tagAct = 3;
				MPI_Send(&values[0],nrEvents,MPI_FLOAT,postNodes[i],tagAct,NETWORK_COMM_WORLD);
			}
		}

		// irecv 3
		if(i<=preNodes.size() && preNodes.size()>0)
		{
			if(preNodes[i]!=this->network()->MPIGetNodeId())
			{
				MPI_Wait(&mpiRequest2,MPI_STATUS_IGNORE);
				MPI_Wait(&mpiRequest3,MPI_STATUS_IGNORE);

				for(int j=0;j<nrRecv;j++)
				{
					Unit* unit = m_network->GetUnitFromId(recvIds[j]);//m_hashIdUnit[recvIds[j]];
					UnitModifier* e = unit->CreateEvent();
					e->SetValue(recvValues[j]);
					vector<Unit*> recvUnits = *m_network->GetPreUnitsFromId(recvIds[j]);//m_hashIdPreUnit[recvIds[j]]; // hash for the units on this node that has incoming unit presynaptically
					m_allEventsIncoming.push_back(e);

					for(int k=0;k<recvUnits.size();k++)
						recvUnits[k]->AddEventIncoming(e);
					//unit->SetEventOutgoing(e);
					unit->SetLocal(true);// fetch locally from now on
				}
			}
		}
	}
}
*/

void PopulationColumns::Simulate()
{
	TimingStart(m_name);

	/*for(int i=0;i<m_localRateUnits.size();i++)//m_minicolumns.size();i++)
	{
		//if(m_minicolumns[i]->IsLocal())
		//m_mi
		nicolumns[i]->SendReceive();
		m_localRateUnits[i]->SendReceiveNext();
	}*/

	/*if(this->network()->MPIGetNodeId() == 0)
	{
		cout<<"-";
		cout.flush();
	}*/

	//SendAndReceiveVersionISend();
	//SendAndReceiveVersionAllgather();

	if(IsOn())//m_unitsFixed == false)
	{
		/*if(this->network()->MPIGetNodeId() == 0)
		{
			cout<<"{";
			cout.flush();
		}*/

		//ResetLocalities(); // not necessary

		TimingStart("SimulateUnits");

		int count = 0;

		/*for(int i=0;i<m_hypercolumns.size();i++)
		{
			if(m_hypercolumns[i]->IsLocal() == true)
				m_hypercolumns[i]->Simulate();

			vector<RateUnit*> mcs = ((Hypercolumn*)m_hypercolumns[i])->GetRateUnits();

			for(int j=0;j<mcs.size();j++)
			{
				if(mcs[j]->IsLocal() == true)
				{
					mcs[j]->Simulate();
					count++;
				}	
			}
		}*/

		//for(int i=m_localUnitIndexesInterval[0];i< m_localUnitIndexesInterval[1]; i++)
		for(int i=0;i<m_minicolumns.size();i++)
		{
			m_minicolumns[i]->Simulate();
			count++;
		}

		for(int i=0;i<m_hypercolumns.size();i++)
		{
			m_hypercolumns[i]->Simulate();
			count++;
		}

		/*if(network()->MPIGetNodeId() == 0)
			cout<<count;*/

		TimingStop("SimulateUnits");

		/*
		for(int i=0;i<m_minicolumns.size();i++)
		{
			this->network()->SendAndReceiveVersionAllgather(); // slow to have it here

			if(m_minicolumns[i]->IsLocal() == true)
				m_minicolumns[i]->Simulate();
		}
		*/

/*		for(int i=0;i<m_minicolumns.size();i++)
			if(m_minicolumns[i]->IsLocal() == true)
				m_minicolumns[i]->ClearEventsIncoming(); // must be done after all simulate since units share incoming events
*/


/*		for(int i=0;i<m_hypercolumns.size();i++)
			if(m_hypercolumns[i]->IsLocal() == true)
				m_hypercolumns[i]->ClearEventsIncoming();
*/
		// spiking units

		/*if(this->network()->MPIGetNodeId() == 0)
		{
			cout<<"'";
			cout.flush();
		}*/

		/*for(int i=0;i<m_allEventsIncoming.size();i++)
		{
		delete m_allEventsIncoming[i];
		}

		m_allEventsIncoming.clear();*/

		//ResetLocalities(); // ok?

		/*	if(this->network()->MPIGetNodeId() == 0)
		{
		cout<<"/";
		cout.flush();
		}*/

		/*	if(this->network()->MPIGetNodeId() == 0)
		{
		cout<<"\\";
		cout.flush();
		}*/

		//ResetLocalities();
	}
	else
	{
		//TimingStart("IsNotOnClearEvents");
		
		// necessary to do? also check memory
		// also currently being done in Unit::Simulate
		
		/*for(int i=m_localUnitIndexesInterval[0];i< m_localUnitIndexesInterval[1]; i++)
		{
			m_minicolumns[i]->ClearEventsIncoming();
		}*/

		for(int i=0;i<m_minicolumns.size();i++)
			m_minicolumns[i]->ClearEventsIncoming();

		//TimingStop("IsNotOnClearEvents");
		/*for(int i=0;i<m_minicolumns.size();i++)
			if(m_minicolumns[i]->IsLocal() == true)
				m_minicolumns[i]->ClearEventsIncoming();
		for(int i=0;i<m_hypercolumns.size();i++)
			if(m_hypercolumns[i]->IsLocal() == true)
				m_hypercolumns[i]->ClearEventsIncoming();
				*/
	}

	// Note: PopulationModifiers are run even if layer is switched off (is overlay on layer, not "inside")
	TimingStart("PopulationModifiers");
	for(int i=0;i<m_eventLayers.size();i++)
	{
		if(m_eventLayers[i]->IsOn())
			m_eventLayers[i]->Simulate();
	}
	TimingStop("PopulationModifiers");

	//TimingStart("IsNewEvent");
	// correct to have this here for all types of units?
	/*for(int i=0;i<m_hypercolumns.size();i++)
	{
		vector<RateUnit*> mcs = ((Hypercolumn*)m_hypercolumns[i])->GetRateUnits();

		for(int j=0;j<mcs.size();j++)
		{
			if(mcs[j]->IsLocal() == true)
			{
				if(mcs[j]->GetValue() != 0.0)
				{
					mcs[j]->IsNewEvent(true); // currently removes prev existing events (ok?)
					mcs[j]->AddEventOutgoing();
				}
				else
					mcs[j]->IsNewEvent(false);
			}
		}
	}*/

	for(int j=0;j<m_minicolumns.size();j++)//for(int j=m_localUnitIndexesInterval[0];j< m_localUnitIndexesInterval[1]; j++)
	{
		if(m_minicolumns[j]->GetValue() != 0.0)
		{
			m_minicolumns[j]->IsNewEvent(true); // currently removes prev existing events (ok?)
			m_minicolumns[j]->AddEventOutgoing();
		}
		else
			m_minicolumns[j]->IsNewEvent(false);
	}

	//TimingStop("IsNewEvent");

	if(IsRecording())
	{
		TimingStart("RecordingLayer");

		if(m_meterLayer->GetSaveType() == Storage::MPI_Binary)//m_network->GetFilePreference() == Network::MPI_Binary)
		{
			vector<float> values = vector<float>(m_minicolumns.size());
			for(int j=0;j<m_minicolumns.size();j++)
			{
				values[j] = m_minicolumns[j]->GetValue();
			}

			m_recordedValues.push_back(values);
		}
		else
		{
			vector<float> values = vector<float>(m_minicolumns.size());
			for(int j=0;j<m_minicolumns.size();j++)
			{
				values[j] = m_minicolumns[j]->GetValue();
			}

			m_recordedValues.push_back(values);

			/*
			this->MPI()->MPIMakeLayerValuesLocal();//MPIMakeHypercolumnsValuesLocal(); // fix for recording

			if(this->network()->MPIGetNodeId() == 0)
			{
				//if(empty == false)
				vector<float> values = this->GetValuesBuffer();
				//for(int i=0;i<m_minicolumns.size();i++)
				//	values[i] = ((RateUnit*)m_minicolumns[i])->GetSubThresholdValue();

				values.insert(values.begin(),this->network()->GetCurrentTimeStep());
				m_recordedValues.push_back(values);//values);

			}
			*/
		}

		TimingStop("RecordingLayer");
	}

	// for graded, go through all incoming connections and retrieve events
	// for spiking, could allow any event to come during this timestep

	TimingStop(m_name);
}

void Population::ClearEventsIncoming()
{
	// Unit properties reset history
	/*vector<UnitModifier*> layerUnitProperties = this->GetUnitPropertiesLayer();

	for(int i=0;i<layerUnitProperties.size();i++)
	{
		layerUnitProperties[i]->Clear();
	}*/
	

	// assumes even distribution
	for(int i=m_localUnitIndexesInterval[0];i< m_localUnitIndexesInterval[1]; i++)
	{
		m_units[i]->ClearEventsIncoming();
	}
}

void PopulationColumns::ClearEventsIncoming()
{
	// assumes even distribution
	//for(int i=m_localUnitIndexesInterval[0];i< m_localUnitIndexesInterval[1]; i++)
	for(int i=0;i<m_minicolumns.size();i++)
	{
		m_minicolumns[i]->ClearEventsIncoming();
	}
}

void PopulationColumns::SetValue(int localIndex, float value)
{
	for(int i=0;i<m_minicolumns.size();i++)
	{
		if(m_minicolumns[i]->GetUnitIdLocal() == localIndex)
		{
			m_minicolumns[i]->SetValue(value);
			break;
		}
	}
}

void PopulationColumns::SetValuesAll(vector<float> values, bool alsoSimulate)
{
	if(m_minicolumns.size() == values.size())
	{
		for(int i=0;i<m_minicolumns.size();i++)
		{
			m_minicolumns[i]->SetValue(values[i]);
			if(alsoSimulate == false)
				((RateUnit*)m_minicolumns[i])->SetNoUpdatingCurrentTimeStep(true); // test this
		}
	}
	else
	{
		int index = 0;
		for(int i=0;i<m_minicolumns.size();i++)//values.size();i++)
		{
			//		if(index==values.size())
			//			index-=values.size(); // repeat pattern

			if(m_minicolumns[index]->GetUnitIdLocal()<values.size())
			{
				((RateUnit*)m_minicolumns[i])->SetValue(values[m_minicolumns[index]->GetUnitIdLocal()]);//index]);
			}
			
			if(alsoSimulate == false)
				((RateUnit*)m_minicolumns[i])->SetNoUpdatingCurrentTimeStep(true); // test this

			index++;
		}
	}

	if(alsoSimulate == true)
		this->m_keepValues = true;
	else
		this->m_keepValues = false;
	// clear also all incoming values
}

vector<float> PopulationColumns::GetValuesLocal()
{
	vector<float> f(m_minicolumns.size());

	for(int i=0;i<m_minicolumns.size();i++)
	{
		f[i] = ((RateUnit*)m_minicolumns[i])->GetValue();
	}

	return f;
}

/*void PopulationColumns::SetValuesAll(vector<short> values)
{
	for(int i=0;i<m_minicolumns.size();i++)
	{
		((RateUnit*)m_minicolumns[i])->SetValue((float)values[i]);
	}
}*/

void PopulationColumns::Reset()
{
	for(int i=0;i<m_minicolumns.size();i++)
	{
		((RateUnit*)m_minicolumns[i])->SetValue(0.0);
	}
}

void PopulationColumns::SetValuesLocal(vector<float> values)
{
	//	if(values.size()!=m_units.size()) // should always be the same size as values
	//		Logger::print("Warning: values not same size as units.");

	long start = m_localUnitIndexesInterval[0]; // relative
	long end = m_localUnitIndexesInterval[1];

	if(values.size()>end-start) // absolute
	{
		start = 0;
		end = values.size();
	}

	for(long i=start;i<end;i++) 
	{
		//((RateUnit*)m_units[i])->SetValue(values[i-start]);//SetEventOutgoing(m_units[i]->CreateEvent(values[i-start]));
		((RateUnit*)m_minicolumns[i])->SetValue(values[i-start]);
	}
}


/*vector<vector<float> > Population::GetLocalWeights()
{
	vector<vector<float> > weights;
	vector<Unit*> localUnits = GetLocalUnits();

	for(int i=0;i<localUnits.size();i++)
	{
		vector<Connection*> pre = localUnits[i]->GetPreConnections();
		vector<float> w;

		for(int j=0;j<pre.size();j++)
		{
			ConnectionFixed* fixed = (ConnectionFixed*)pre[j];
			w.push_back(fixed->GetWeight(j));
		}

		weights.push_back(w);
	}

	return weights;
}*/

/*vector<vector<float> > Population::GetWeights()
{
	vector<vector<float> > weights;

	for(int i=0;i<m_units.size();i++)
	{
		vector<Connection*> pre = m_units[i]->GetPreConnections();
		vector<float> w;

		for(int j=0;j<pre.size();j++)
		{
			ConnectionFixed* fixed = (ConnectionFixed*)pre[j];
			w.push_back(fixed->GetWeight(j));
		}

		weights.push_back(w);
	}

	return weights;
}*/

void Population::AddPre(Population* layer, Connectivity* connectionType)
{
	m_pres.push_back(layer);
	m_preConnectivitys.push_back(connectionType);
	connectionType->SetPreAndPost(layer,this);
	connectionType->network(network());

	layer->AddPost(this, connectionType);
}

void Population::AddPost(Population* layer, Connectivity* connectionType)
{
	m_posts.push_back(layer);
	m_postConnectivitys.push_back(connectionType);
	connectionType->SetPreAndPost(this,layer);
	connectionType->network(network());

	//Connectivity* c2 = new Connectivity(*connectionType); // copy
	//FullConnectivity* fullConnectivity = new FullConnectivity();
	//layer->AddPre(this,(Connectivity*)fullConnectivity);
}

vector<Unit*> Population::GetLocalUnits()
{
	vector<Unit*> localUnits;

	for(int i=0;i<m_units.size();i++)
		if(m_units[i]->IsLocal())
			localUnits.push_back(m_units[i]);
	return localUnits;
}

vector<int> Population::GetMPIDistribution(int nodeId)
{
	vector<int> unitIndexes;

	for(int i=0;i<m_units.size();i++)
	{
		if(m_units[i]->GetNodeId() == nodeId)
			unitIndexes.push_back(i);
	}

	return unitIndexes;
}

vector<int> PopulationColumns::GetMPIDistribution()//int nodeId)
{
	vector<int> unitIndexes(m_minicolumns.size());

	for(int i=0;i<m_minicolumns.size();i++)
	{
		unitIndexes[i] = m_minicolumns[i]->GetUnitIdLocal();
		//if(m_minicolumns[i]->GetNodeId() == nodeId)
		//	unitIndexes.push_back(i);
	}

	return unitIndexes;
}

vector<int> PopulationColumns::GetMPIDistributionHypercolumns(int nodeId)
{
	vector<int> unitIndexes;

	for(int i=0;i<m_hypercolumns.size();i++)
	{
		if(m_hypercolumns[i]->GetNodeId() == nodeId)
			unitIndexes.push_back(i);
	}

	return unitIndexes;
}

void PopulationColumns::AddUnit(Unit* unit)
{
	m_units.push_back(unit);
	network()->AddUnit(unit);//network()->AddUnitToHash(unit);// m_hashIdUnit[unit->GetUnitId()] = unit;

	if(unit->GetType().compare("minicolumn")==0)
	{
		m_minicolumns.push_back((RateUnit*)unit);
	}
	else if(unit->GetType().compare("hypercolumn")==0)
	{
		m_hypercolumns.push_back((Hypercolumn*)unit);
	}
}

void Population::AddLayerEvent(PopulationModifier* eventLayer)
{
	eventLayer->network(this->network()); // move?
	eventLayer->Initialize(this);
	m_eventLayers.push_back(eventLayer);
	m_network->AddNetworkObjectToDelete(eventLayer); // for avoiding duplicate deletion - network responsible for this object
}

void Population::Dispose()
{
	DisposeParent();
}

void Population::DisposeParent()
{
	for(int i=0;i<m_units.size();i++)
	{
		m_units[i]->Dispose(); // can put in destructor
		delete m_units[i];
	}

	for(int i=0;i<m_connectionOutgoing.size();i++) // contains hash weights, will be rebuilt
	{
		vector<ConnectionModifier*> eventConnections = m_connectionOutgoing[i]->GetConnectionModifiers();
		for(int j=0;j<eventConnections.size();j++)
		{
//			eventConnections[j]->Reset();			
		}

		m_connectionOutgoing[i]->Dispose();
		delete m_connectionOutgoing[i];
	}

	// needs to be rebuilt in case of reset
	for(int i=0;i<m_eventLayers.size();i++)
	{
	/*	if(m_eventLayers[i] != NULL)
			m_eventLayers[i]->Reset(); // will work as a reset for e.g. the mds case (!) ok?
			*/
		//delete m_eventLayers[i];
	}

	//m_eventLayers.clear();
	
	m_connectionOutgoing.clear();
	m_connectionIncoming.clear();
	m_units.clear();

	for(int i=0;i<m_preConnectivitys.size();i++)
	{
		m_preConnectivitys[i]->SetInitialized(false);
	}

	m_nrUnits = 0;
}

void PopulationColumns::Dispose()
{
	DisposeParent();
	m_minicolumns.clear();
	m_hypercolumns.clear();
	m_localRateUnits.clear();
}

vector<int> PopulationColumns::GetUnitIdLocals()
{
	vector<int> ids(m_minicolumns.size());
	for(int i=0;i<m_minicolumns.size();i++)
	{
		ids[i] = ((RateUnit*)m_minicolumns[i])->GetUnitIdLocal();
	}

	return ids;
}


void PopulationColumns::SetSilentHypercolumns(bool useSilent, float threshold)
{
	m_useSilentHypercolumns = useSilent;
	m_silentHypercolumnsThreshold = threshold;
	for(int i=0;i<m_hypercolumns.size();i++)
	{
		((Hypercolumn*)m_hypercolumns[i])->SetUseSilent(useSilent);
		((Hypercolumn*)m_hypercolumns[i])->SetSilentThreshold(threshold);
	}
}