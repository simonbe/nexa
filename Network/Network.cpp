#include <mpi.h>
#include <math.h>
#include <ctime>
#include <cstring>
#include <fstream>
#include <deque>
#include "Network.h"
#include "Meter.h"

#if MUSIC_AVAILABLE == 1
MPI_Comm Network::MPIComm;
#endif

// Standard initialization
Network::Network()
{
	// MPI_Init needs to be run before Network object created

	MPI_Comm_size(NETWORK_COMM_WORLD, &m_mpiNrProcs);
	MPI_Comm_rank(NETWORK_COMM_WORLD, &m_mpiNodeId);

	m_currentUnitId = 0;
	m_currentHypercolumnId = 0;
	m_currentLayerId = 0;
	m_currentTimeStep = 0;
	m_simulationResolution = 0.1; //[ms]
	m_firstRun = true;
	m_networkDetailsStored = false;
	m_useTiming = true;
	m_isUsingDelays = false; // gets turned on if IF units are used (or change manually with SetUsingDelays)
	m_keepCommunicationBuffer = false; // can be turned on to reduce communication by saving a part of the communication buffer (manually specfied) from prev time step
	m_isTrackingHypercolumnIds = false; // gets turned on if bcpnn is added

	m_filenameNetworkDetails = new char[50];
	m_filenameTimings = new char[50];
	m_filenameAnalysis = new char[50];
	sprintf(m_filenameNetworkDetails,"NetworkDetails%d.txt",m_mpiNrProcs);
	sprintf(m_filenameTimings,"NetworkTimings%d.txt",m_mpiNrProcs);
	sprintf(m_filenameAnalysis,"NetworkAnalysis%d.txt",m_mpiNrProcs);
	m_savePreference = Storage::MPI_Binary;

#if COMMUNICATION_BUFFER_HASH==1 // compiler specific
	m_communicationBufferType = this->CommunicationBufferHash;
#else
	m_communicationBufferType = this->CommunicationBufferVector;
#endif

	m_synapseModel = this->SynapsesStandard;

	// srand(this->MPIGetNodeId()); // default set in SetSeed
	m_name = "Network";
	m_useExtraFilenameString = false;
	m_networkParameters = new NetworkParameters();
	
	////// Parallelization scheme trackers
	m_mpiNrUnitsParallelizationDefault = 0; // keep track of division for default parallelization scheme
	m_mpiCurrentUnitParallelizationDefault = 0;

	this->SetSeed(false);//false); // default seeds (different important for e.g. random connectivities)
}

// Additional standard filenames for recording
void Network::SetExtraFilenameString(char* extraString)
{
	m_useExtraFilenameString = true;
	m_extraFilenameString = extraString;

	delete[] m_filenameNetworkDetails;
	delete[] m_filenameTimings;
	delete[] m_filenameAnalysis;
	
	m_filenameNetworkDetails = new char[100];
	m_filenameTimings = new char[100];
	m_filenameAnalysis = new char[100];
	
	sprintf(m_filenameNetworkDetails,"NetworkDetails%d_%s.txt",m_mpiNrProcs,extraString);
	sprintf(m_filenameTimings,"NetworkTimings%d_%s.txt",m_mpiNrProcs,extraString);
	sprintf(m_filenameAnalysis,"NetworkAnalysis%d_%s.txt",m_mpiNrProcs,extraString);
}

// Depending on connectivity pattern, may want same seed on all processes
void Network::SetSeed(bool allSame)//, float x);
{
	if(allSame == false)
		srand(this->MPIGetNodeId() + m_runId);
	else
		srand(8);//10);//x);
}

Network::~Network()
{
	this->ClearEventsIncoming();

	//delete m_filenameNetworkDetails;
	//delete m_filenameTimings;
	//delete m_filenameAnalysis;

	m_hashSynapses.clear();

	for(int i=0;i<m_populations.size();i++)
	{
		delete m_populations[i];
	}

	for(int i=0;i<m_meters.size();i++)
		delete m_meters[i];

	for(int i=0;i<m_analysis.size();i++)
		delete m_analysis[i];

	for(int i=0;i<m_connectivityTypes.size();i++)
		delete m_connectivityTypes[i];

	for(int i=0;i<m_networkObjectsToDelete.size();i++)
	{
		delete m_networkObjectsToDelete[i];
		m_networkObjectsToDelete[i] = NULL;
	}
}

// Decide if record timing from the main class
// includes timing of
// - Communication
// - Initializations
void Network::SetUseTiming(bool useTiming)
{
	m_useTiming = useTiming;
}

void Network::SetFilenamesAdditional(string additional)
{
	m_filenameNetworkDetails = new char[100];
	m_filenameTimings = new char[100];
	m_filenameAnalysis = new char[100];
	sprintf(m_filenameNetworkDetails,"NetworkDetails%d_%s.txt",m_mpiNrProcs, additional.c_str());
	sprintf(m_filenameTimings,"NetworkTimings%d_%s.txt",m_mpiNrProcs,additional.c_str());
	sprintf(m_filenameAnalysis,"NetworkAnalysis%d_%s.txt",m_mpiNrProcs,additional.c_str());
}

// Builds network structures and initializes all units
void Network::Initialize()
{
	MPI_Barrier(NETWORK_COMM_WORLD);
	TimingStart("InitializeNetwork");
	int nodeId = this->MPIGetNodeId();

	for(int i=0;i<m_populations.size();i++)
	{
	//	if(m_populations[i]->IsInitialized() == false) // now taken care of in each network object
			m_populations[i]->Initialize();
	}

	for(int i=0;i<m_populations.size();i++)
	{
		m_populations[i]->InitializeParallelization();
	}

	// create communicators for populations
	for(int i=0;i<m_populations.size();i++)
	{
		m_populations[i]->MPI()->MPICreateCommLayer();
	}

	for(int i=0;i<m_populations.size();i++)
	{
		m_populations[i]->InitializeConnectionsEventsAndParameters();
	}

	// optimization cache of all possible input ids (derived from hashed synapses)
	if(this->MPIGetNodeId() == 0)
	{
		cout<<"Create tables all pre ids union...";
		cout.flush();
	}

	this->CreateAllPreIdsUnion();

	if(this->MPIGetNodeId() == 0)
	{
		cout<<"ok.\n";
		cout.flush();
	}

	for(int i=0;i<m_populations.size();i++)
	{
		if(m_populations[i]->GetUnitsAll()->size() != 0)//m_populations[i]->GetNrUnits() != 0)
			m_populationIndexesThisProc.push_back(i);
	}
	
	TimingStart("CreateAllPostProcs");	

	if(this->MPIGetNodeId() == 0)
	{
		cout<<"Create tables all post processes...";
		cout.flush();
	}

	// create which post processes to send each event output from a unit
	// not needed for all setups (change to check if needed) (not correctly scalable)
	this->CreateAllPostProcs();

	TimingStop("CreateAllPostProcs");

	if(this->MPIGetNodeId() == 0)
	{
		cout<<"ok.\n";
		cout.flush();
	}

	if(this->MPIGetNodeId() == 0)
	{
		cout<<"(local sizeof hashSynapses == "<<sizeof((*network()->GetSynapses()))<< ")";
		cout.flush();
	}

	TimingStop("InitializeNetwork");
	MPI_Barrier(NETWORK_COMM_WORLD);
}



void Network::Simulate(int nrTimesteps)
{
	for(int i=0;i<nrTimesteps;i++)
		Simulate();
}

// Main simulation routine for advancing one time step
void Network::Simulate()
{
	if(m_firstRun == true) // can move
	{
		// store network details to disk
		if(this->MPIGetNodeId() == 0)
		{
			StoreNetworkDetails();
		}

		MPI_Barrier(NETWORK_COMM_WORLD); // to get start of timing right
	}

	MPI_Barrier(NETWORK_COMM_WORLD);

	TimingStart(m_name);

	// Action

	TimingStart("Simulate");

	for(int i=0;i<m_populationIndexesThisProc.size();i++)//m_populations.size();i++)
	{
		// distribute all events, allgather-version
		//		SendAndReceiveVersionAllgather();
		//m_populations[i]->ResetLocalities(); // needed?
		//	if(m_populations[i]->GetNrUnits() != 0)
		m_populations[m_populationIndexesThisProc[i]]->Simulate();//m_populations[i]->Simulate();
	}
	TimingStop("Simulate");

//	for(int i=0;i<m_minicolumns.size();i++)
	//((RateUnit*)m_minicolumns[i])->SetBeta(0.0); // special case reset - if >1 connections, reset needs to be outside 

	
	// put in correct place ?
	MPI_Barrier(NETWORK_COMM_WORLD);

	// additional timing for comm
//	TimingStart("SendAndReceiveVersionAllgather");
	
#if USE_COMMUNICATION_ALLTOALL == 1 // not needed to be on compile level
	CommunicationVersionAlltoall();
#else
	CommunicationVersionAllgather(); // Communicate all messages
#endif

	
//	TimingStop("SendAndReceiveVersionAllgather");

	MPI_Barrier(NETWORK_COMM_WORLD);

	// Modify
	TimingStart("Modify");
	for(int i=0;i<m_populationIndexesThisProc.size();i++)//m_populations.size();i++)
	{
		//if(m_populations[i]->GetNrUnits() != 0)
		//	m_populations[i]->Modify();
		m_populations[m_populationIndexesThisProc[i]]->Modify();
//		m_populations[i]->ResetLocalities(); // needed?
	}
	TimingStop("Modify");

	MPI_Barrier(NETWORK_COMM_WORLD);

	// Analysis
	TimingStart("Analysis");
	for(int i=0;i<this->m_analysis.size();i++)
	{
		if(m_analysis[i]->IsOn())
			m_analysis[i]->Simulate();
	}
	TimingStop("Analysis");

	//this->ClearMemory();
	m_currentTimeStep+=this->GetTimeResolution();


	TimingStop(m_name);

	if(m_firstRun)
		m_firstRun = false;
}


// Communication routine
// - uses Allgather as main mpi routine
// - will transmit outgoing messages to all processes.
void Network::CommunicationVersionAllgather()
{
	TimingStart("SendAndReceiveVersionAllgather");
	
	// go through all layers and transfer events data
	vector<short> eventTypes;
	vector<long> eventIds;
	vector<vector<long> > totalEventIds(m_mpiNrProcs);
	vector<long> eventHypercolumnIds;
	vector<vector<long> > totalEventHypercolumnIds(m_mpiNrProcs);
	vector<float> eventData; // currently only floats (will get changed, MPI supports sending varying variables by the w versions, e.g. MPI_Alltoallw)
	vector<vector<float> > totalEventData(m_mpiNrProcs);

	vector<vector<int> > eventDestinationsProcs(m_mpiNrProcs);

	vector<float> localEventData;
	vector<long> localEventIds;
	vector<long> localEventHypercolumnIds;

	TimingStart("CommunicationGetEvents");

	/*//	- (currently) getting all possible events
	for(int i=0;i<m_populations.size();i++)
	{
		vector<Unit*> localUnits = ((PopulationColumns*)m_populations[i])->GetLocalRateUnits(); // "minicolumn"...

		for(int j=0;j<localUnits.size();j++)
		{
			if(localUnits[j]->IsNewEvent() == true)
			{
				// also sending type of event info atm, not necessary in case of one-type network (+1 short)
				vector<UnitModifier*> events = localUnits[j]->GetEvents();

				for(int k=0;k<events.size();k++)
				{
					eventTypes.push_back(events[k]->GetEventTypeId());
					eventIds.push_back(events[k]->GetFromUnitId());
					eventHypercolumnIds.push_back(events[k]->GetFromHypercolumnId());
					vector<float> eData = events[k]->GetEventData();
					for(int m=0;m<eData.size();m++)
						eventData.push_back(eData[m]);
				}

				localUnits[j]->IsNewEvent(false); // will be set true again by by unit's simulateeventque fcn or CreateEvent fcn
			}
		}
	}*/

	//	get all events and where to send them
	int currentIndex = 0;
	int nrEvents = 0;
	//for(int i=0;i<m_populations.size();i++)
	for(int i=0;i<m_populationIndexesThisProc.size();i++)
	{
		vector<Unit*> localUnits = ((PopulationColumns*)m_populations[m_populationIndexesThisProc[i]])->GetLocalRateUnits(); // "minicolumn"...

		vector<int> destinations;
		bool firstDest = true;

		if(m_firstRun == true) // build which processes to send to only first simulation step
		{
			vector<int> toProcesses;

			for(int m=0;m<m_populations[m_populationIndexesThisProc[i]]->GetOutgoingConnections().size();m++)
			{
				vector<int> prcs = m_populations[m_populationIndexesThisProc[i]]->GetOutgoingConnections()[m]->PostLayer()->MPIGetProcessesUsed();

				if(prcs.size()==0) // all procs take part in population so send to all
				{
					for(int n=0;n<this->MPIGetNrProcs();n++)
					{
						toProcesses.push_back(n);
					}
				}
				else // only to subset
				{
					for(int n=0;n<prcs.size();n++)
					{
						toProcesses.push_back(prcs[n]);
					}
				}
			}

			sort(toProcesses.begin(),toProcesses.end());
			toProcesses.erase(unique(toProcesses.begin(),toProcesses.end()),toProcesses.end()); // remove duplicates

			m_communicationToProcs.push_back(toProcesses);
		}

		for(int j=0;j<localUnits.size();j++)
		{
			if(localUnits[j]->IsNewEvent() == true)
			{
				// also sending type of event info atm, not necessary in case of one-type network (+1 short)
				vector<UnitModifier*> events = localUnits[j]->GetEvents();
				vector<int>* postProcesses = localUnits[j]->GetPostProcesses();

				long eventId;
				short eventType;
				long eventHypercolumnId;

				for(int k=0;k<events.size();k++)
				{
					eventType = events[k]->GetEventTypeId();
					eventId = events[k]->GetFromUnitId();
					eventHypercolumnId = events[k]->GetFromHypercolumnId();

					//eventTypes.push_back(eventType);
					//eventIds.push_back(eventId);
					//eventHypercolumnIds.push_back(eventHypercolumnId);

					vector<float> eData = events[k]->GetEventData();
					//for(int m=0;m<eData.size();m++)
					//	eventData.push_back(eData[m]);

					bool onlyLocal = true;

					if( postProcesses->size()>0)
						onlyLocal = false;

					/*for(int n=0;n<m_communicationToProcs[i].size();n++)
					{
						if(m_communicationToProcs[i][n] == m_mpiNodeId) // local data, no need to send
						{
							//localEventIds.push_back(eventId);
							//localEventData.push_back(eData[0]);
							//if(m_isTrackingHypercolumnIds)
							//	localEventHypercolumnIds.push_back(eventHypercolumnId);
						}
						else // non-local data
						{
							onlyLocal = false;
							//eventDestinationsProcs[m_communicationToProcs[i][n]].push_back(currentIndex);
							//totalEventIds[m_communicationToProcs[i][n]].push_back(eventId);
							//if(m_isTrackingHypercolumnIds)
							//	totalEventHypercolumnIds[m_communicationToProcs[i][n]].push_back(eventHypercolumnId);

							// unnecessary atm? all models currently only use one data point
							//for(int p=0;p<eData.size();p++)
							//	totalEventData[m_communicationToProcs[i][n]].push_back(eData[p]);
						}
					}*/

					if(onlyLocal == false)
					{
						eventIds.push_back(eventId);
						eventData.push_back(eData[0]);
						if(m_isTrackingHypercolumnIds)
							eventHypercolumnIds.push_back(eventHypercolumnId);

						nrEvents++;
					}
					else // do not communicate, only local process needs it
					{
						localEventIds.push_back(eventId);
						localEventData.push_back(eData[0]);
						if(m_isTrackingHypercolumnIds)
							localEventHypercolumnIds.push_back(eventHypercolumnId);
					}

					currentIndex++;
				}

				localUnits[j]->IsNewEvent(false); // will be set true again by by unit's simulateeventque fcn or CreateEvent fcn
			}
		}
	}

	TimingStop("CommunicationGetEvents");

	// redistribute nr events
//	int nrEvents = eventIds.size();
	int totalNrEvents = 0;
	//TimingStart("CommunicationBarrierTest");
	//MPI_Barrier(NETWORK_COMM_WORLD);
	//TimingStop("CommunicationBarrierTest");
	
	//TimingStart("CommunicationReduceEvents");
	TimingStart("Allreduce");

	MPI_Allreduce(&nrEvents,&totalNrEvents, 1,MPI_INT,MPI_SUM,NETWORK_COMM_WORLD); // could switch to multiple MPI communicators as network becomes large (and is modular enough to benefit).

	TimingStop("Allreduce");

	// debug - will be removed
	if(this->MPIGetNodeId() == 0)
		cout<<"{"<<localEventIds.size()<<"/"<<totalNrEvents<<"} ";
	
	//TimingStop("CommunicationReduceEvents");
	vector<long> allIds(totalNrEvents);
	vector<long> allHypercolumnIds(totalNrEvents);
	vector<short> allTypes(totalNrEvents);
	vector<float> allData(totalNrEvents);

	if(totalNrEvents>0)
	{
		// graded minicolumn distribution
		//	- (currently) reconstructing the events
		vector<float>	values(eventIds.size());
		vector<long>	ids(eventIds.size());
		vector<short>	types(eventIds.size());
		vector<short>	nrValues(eventIds.size());
		vector<vector<float> > data(eventIds.size());

		vector<int> rcounts(m_mpiNrProcs);// = new int[m_mpiNrProcs];
		vector<int> displs(m_mpiNrProcs);

		// counts from each node
		TimingStart("Allgather");

		MPI_Allgather( &nrEvents, 1, MPI_INT, &rcounts[0], 1, MPI_INT, NETWORK_COMM_WORLD);

		for(int i=1;i<rcounts.size();i++)
		{
			displs[i] = rcounts[i-1]+displs[i-1];
		}

		if(eventIds.size()>0)
		{
			MPI_Allgatherv( &eventIds[0],nrEvents,MPI_LONG,&allIds[0],&rcounts[0],&displs[0],MPI_LONG,NETWORK_COMM_WORLD);
			if(m_isTrackingHypercolumnIds)
				MPI_Allgatherv( &eventHypercolumnIds[0],eventIds.size(),MPI_LONG,&allHypercolumnIds[0],&rcounts[0],&displs[0],MPI_LONG,NETWORK_COMM_WORLD); // not necessary always
			//MPI_Allgatherv( &eventTypes[0],eventTypes.size(),MPI_SHORT,&allTypes[0],&rcounts[0],&displs[0],MPI_SHORT,NETWORK_COMM_WORLD);
			
			// assuming only type 1 or 2, otherwise change totDataSize, rcounts and displs
			int totDataSize = eventData.size()*1;
			MPI_Allgatherv( &eventData[0],totDataSize,MPI_FLOAT,&allData[0],&rcounts[0],&displs[0],MPI_FLOAT,NETWORK_COMM_WORLD);
		}
		else
		{
			MPI_Allgatherv( NULL,0,MPI_INT,&allIds[0],&rcounts[0],&displs[0],MPI_LONG,NETWORK_COMM_WORLD);

			if(m_isTrackingHypercolumnIds)
				MPI_Allgatherv( NULL,0,MPI_INT,&allHypercolumnIds[0],&rcounts[0],&displs[0],MPI_LONG,NETWORK_COMM_WORLD);
			//MPI_Allgatherv( NULL,0,MPI_INT,&allTypes[0],&rcounts[0],&displs[0],MPI_SHORT,NETWORK_COMM_WORLD);
			MPI_Allgatherv( NULL,0,MPI_FLOAT,&allData[0],&rcounts[0],&displs[0],MPI_FLOAT,NETWORK_COMM_WORLD);
		}

		TimingStop("Allgather");
		//TimingStart("CommunicationAddingEvents");

		//clear buffers
	}

#if COMMUNICATION_BUFFER_HASH == 0
		//if(m_communicationBufferType == this->CommunicationBufferVector)
		//{
			// IF delays and delay lines
			if(this->m_isUsingDelays == true)
			{
				for(int j=0;j<m_incomingBufferData.size();j++)
				{
					m_incomingBufferDataDelays[j].clear();
					//m_incomingBufferHypercolumnIds.clear();
				}

				for(int j=0;j<allIds.size();j++)
				{
					while(m_incomingBufferData.size()<allIds[j]+1)
					{
						m_incomingBufferDataDelays.push_back(vector<float>(0));
						m_incomingBufferHypercolumnIds.push_back(-1);
					}

					m_incomingBufferDataDelays[allIds[j]].push_back(allData[j]);
					m_incomingBufferHypercolumnIds[allIds[j]] = allHypercolumnIds[j];
				}
			}
			else // default, faster access
			{
				for(int j=0;j<m_allPreIds.size();j++)
					m_incomingBufferData[m_allPreIds[j]] = 0; // allPreIds[j] - allPreIds[0]

				for(int j=0;j<allIds.size();j++)
				{
					/*				while(m_incomingBufferData.size()<allIds[j]+1)
					{
					m_incomingBufferData.push_back(vector<float>(0));
					m_incomingBufferHypercolumnIds.push_back(-1);
					}*/

					if(allIds[j] < m_incomingBufferData.size()) // only buffer those that this process may need to access
					{
						m_incomingBufferData[allIds[j]] = allData[j];
						m_incomingBufferHypercolumnIds[allIds[j]] = allHypercolumnIds[j];
					}
				}
			}
#else		//}
		//else if(m_communicationBufferType == this->CommunicationBufferHash)
		//{
			// atm using map instead of hash impl
			//this->m_incomingBufferDataHash = map<long,float>(allIds.size());
	//m_incomingBufferDataHash = unordered_map<long,float>(allIds.size());


	// remove unused 
	// (not optimized, slow - remove lowest/highest mixed, switch search)
	TimingStart("HashBuffersPrepare");

/*	TimingStart("HashBuffersPrepare00");
	sort(allIds.begin(),allIds.end());
	TimingStop("HashBuffersPrepare00");
	TimingStart("HashBuffersPrepare01");
	//sort(allData.begin(),allData.end());
	TimingStop("HashBuffersPrepare01");
	if(m_isTrackingHypercolumnIds)
		sort(allHypercolumnIds.begin(),allHypercolumnIds.end());
		*/
	
	vector<long> tempIds = allIds;
	vector<int> useIds;

	long offset=0;
	if(m_allPreIds.size()>0 && tempIds.size()>0)
	{
		vector<long>::iterator lb = m_allPreIds.begin();
		if(m_allPreIds[0]<tempIds[0])
			lb = lower_bound(m_allPreIds.begin(),m_allPreIds.end(),tempIds[0]);//m_allPreIds.begin();

		vector<long>::iterator ub = m_allPreIds.end();
		if(m_allPreIds[m_allPreIds.size()-1]>tempIds[tempIds.size()-1])
			ub = upper_bound(m_allPreIds.begin(),m_allPreIds.end(),tempIds[tempIds.size()-1]);//m_allPreIds.end();

		int ubInt,lbInt;
		ubInt = int(ub-m_allPreIds.begin());
		lbInt = int(lb-m_allPreIds.begin());
		//useIds.reserve(m_allPreIds.size());//1.5*m_nrLocalMessagesLastTimestep); // optimization

		//TimingStart("HashBuffersPrepare2");

		for(int i=tempIds.size()-1;i>-1;i--)//0;i<tempIds.size();i++)//
		{
			if(tempIds[i] <= m_allPreIds[ubInt-1] && tempIds[i]>= m_allPreIds[lbInt])//minId) // often enough
			{
				vector<long>::iterator it;//,it2;
				//it2 = upper_bound(lb,ub,tempIds[i]);
				//ub = it2;
				it = lower_bound(lb,ub,tempIds[i]);//upper_bound(lb,ub,tempIds[i]);//lower_bound(lb,ub, tempIds[i]);
				
				//ub = it;
				//if(it!=ub)
				//{
				
				//lb = it;
				//ubInt = int(ub-m_allPreIds.begin());
				//lbInt = int(lb-m_allPreIds.begin());

				if(tempIds[i]==*it)//tempIds[i]//m_allPreIds[int(it-m_allPreIds.begin()-1)])// && !(tempIds[i]>*it))
				{
					//	if(m_allPreIds[int(lb-m_allPreIds.begin())] == tempIds[i])	// change check
					useIds.push_back(i);
				}
				//}
			}
		}
		
		//TimingStop("HashBuffersPrepare2");

		/*for(int i=0;i<m_populations.size();i++)
		{
			long startId,endId;
			if(m_populations[i]->GetUnitsAll()->size()>0 && m_populations[i]->IsOn()) // only continue if this process is part of population
			{
				startId = m_populations[i]->GetUnitsAll()->at(0)->GetUnitId();
				endId = m_populations[i]->GetUnitsAll()->at(m_populations[i]->GetUnitsAll()->size()-1)->GetUnitId();

				vector<Connection*> conns = m_populations[i]->GetIncomingConnections();
				for(int j=0;j<conns.size();j++)
				{
					if(conns[j]->KeepActiveBuffer() == false)
					{
						conns[j]->ClearActiveBuffer();

						//for(int k=0;k<allIds.size();k++)
						//	m_populations[i]->GetIncomingConnections()[j]->AddActiveEvent(allIds[k],allData[k]);//,allData[k]);
						if(allIds.size()>0)
							conns[j]->AddActiveEvents(allIds,allData);
					}
				}
			}
		}*/
	}

	//TimingStart("HashBuffersPrepare3");
	vector<long> tempIds2(useIds.size());
	vector<long> tempHypercolumnIds;
	vector<float> tempData(useIds.size());
	if(m_isTrackingHypercolumnIds)
		tempHypercolumnIds = vector<long>(useIds.size());

	reverse(useIds.begin(),useIds.end());
	//sort(useIds.begin(),useIds.end());

	for(int i=0;i<useIds.size();i++)
	{
		tempIds2[i] = allIds[useIds[i]];
		tempData[i] = allData[useIds[i]];

		if(m_isTrackingHypercolumnIds)
			tempHypercolumnIds[i] = allHypercolumnIds[useIds[i]];
	}

	allIds = tempIds2;
	allData = tempData;
	if(m_isTrackingHypercolumnIds)
		allHypercolumnIds = tempHypercolumnIds;
	//TimingStop("HashBuffersPrepare3");

	// also add the locally stored ones
	for(int i=0;i<localEventIds.size();i++)
	{
		allIds.push_back(localEventIds[i]);
		allData.push_back(localEventData[i]);
	}

	if(m_isTrackingHypercolumnIds)
	{
		for(int i=0;i<localEventHypercolumnIds.size();i++)
		{
			allHypercolumnIds.push_back(localEventHypercolumnIds[i]);
		}
	}


	TimingStop("HashBuffersPrepare");

	TimingStart("HashBuffersFilling");

	m_incomingBufferDataHash.clear();
	m_incomingBufferHypercolumnIdsHash.clear();

	/*if(this->MPIGetNodeId() == 0)
	{
		cout<<"["<<allIds.size()<<"]";
	}
	else if(this->MPIGetNodeId() == 1)
	{
		cout<<"/"<<allIds.size()<<"/";
	}
	else if(this->MPIGetNodeId() == this->MPIGetNrProcs()-1)
	{
		cout<<"|"<<allIds.size()<<"|";
	}
	else if(this->MPIGetNodeId() == this->MPIGetNrProcs()-2)
	{
		cout<<"*"<<allIds.size()<<"*";
	}

	cout.flush();*/

	if(this->m_keepCommunicationBuffer == true) // used to lower communication in case of static input, e.g. input data from file over some time steps
	{
		m_incomingBufferDataHash = this->m_bufferToKeepData;
		m_incomingBufferHypercolumnIdsHash = this->m_bufferToKeepHypercolumnIds; 
	}

	if(m_isTrackingHypercolumnIds) // set on by bcpnn
	{
#if USE_UNORDERED_MAP == 1
		m_incomingBufferDataHash.rehash(allIds.size() + m_incomingBufferDataHash.size());
		m_incomingBufferHypercolumnIdsHash.rehash(allData.size() + m_incomingBufferHypercolumnIdsHash.size());
#endif

		for(int i=0;i<allIds.size();i++) // check how slow (?)
		{
			m_incomingBufferDataHash[allIds[i]] = allData[i];
			m_incomingBufferHypercolumnIdsHash[allIds[i]] = allHypercolumnIds[i]; // when used?
		}
	}
	else
	{
#if USE_UNORDERED_MAP == 1
		m_incomingBufferDataHash.rehash(allData.size() + m_incomingBufferDataHash.size());
#endif

		for(int i=0;i<allIds.size();i++) // check how slow (?)
		{
			m_incomingBufferDataHash[allIds[i]] = allData[i];
		}
	}

	TimingStop("HashBuffersFilling");

	TimingStart("HashBuffersAddEvent");

#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	for(int i=0;i<m_populations.size();i++)
	{
		if(m_populations[i]->GetUnitsAll()->size()>0 && m_populations[i]->IsOn()) // only continue if this process is part of population
		{
			vector<Connection*> conns = m_populations[i]->GetIncomingConnections();
			for(int j=0;j<conns.size();j++)
			{
				if(conns[j]->KeepActiveBuffer() == false)
				{
					conns[j]->ClearActiveBuffer();

					//for(int k=0;k<allIds.size();k++)
					//	m_populations[i]->GetIncomingConnections()[j]->AddActiveEvent(allIds[k],allData[k]);//,allData[k]);
					if(allIds.size()>0)
						conns[j]->AddActiveEvents(allIds,allData);
				}
			}
		}
	}
#endif
	TimingStop("HashBuffersAddEvent");

#endif
		//TimingStop("CommunicationAddingEvents");
	
	TimingStop("SendAndReceiveVersionAllgather");
}


// Communication routine
// - uses MPI_Alltoall as main mpi routine
// - reduces total nr messages compared to allgather version. Will not transmit messages from a population to another which it is not connected to.
void Network::CommunicationVersionAlltoall()
{
	TimingStart("SendAndReceiveVersionAlltoall");
	
	// go through all layers and transfer events data
	vector<short> eventTypes;
	vector<long> eventIds;
	vector<vector<long> > totalEventIds(m_mpiNrProcs);
	vector<long> eventHypercolumnIds;
	vector<vector<long> > totalEventHypercolumnIds(m_mpiNrProcs);
	vector<float> eventData; // currently only floats (will get changed)
	vector<vector<float> > totalEventData(m_mpiNrProcs);
	vector<vector<int> > eventDestinationsProcs(m_mpiNrProcs);

	vector<float> localEventData;
	vector<long> localEventIds;
	vector<long> localEventHypercolumnIds;

	TimingStart("CommunicationGetEvents");

	//	get all events and where to send them
	int currentIndex = 0;
	int nrEvents = 0;
	//for(int i=0;i<m_populations.size();i++)
	for(int i=0;i<m_populationIndexesThisProc.size();i++)
	{
		vector<Unit*> localUnits = ((PopulationColumns*)m_populations[m_populationIndexesThisProc[i]])->GetLocalRateUnits(); // "minicolumn"...

		vector<int> destinations;
		bool firstDest = true;

		if(m_firstRun == true) // build which processes to send to only first time
		{
			vector<int> toProcesses;

			for(int m=0;m<m_populations[m_populationIndexesThisProc[i]]->GetOutgoingConnections().size();m++)
			{
				vector<int> prcs = m_populations[m_populationIndexesThisProc[i]]->GetOutgoingConnections()[m]->PostLayer()->MPIGetProcessesUsed();

				if(prcs.size()==0) // all procs take part in population so send to all
				{
					for(int n=0;n<this->MPIGetNrProcs();n++)
					{
						toProcesses.push_back(n);
					}
				}
				else // only to subset
				{
					for(int n=0;n<prcs.size();n++)
					{
						toProcesses.push_back(prcs[n]);
					}
				}
			}

			sort(toProcesses.begin(),toProcesses.end());
			toProcesses.erase(unique(toProcesses.begin(),toProcesses.end()),toProcesses.end()); // remove duplicates

			m_communicationToProcs.push_back(toProcesses);
		}

		for(int j=0;j<localUnits.size();j++)
		{
			if(localUnits[j]->IsNewEvent() == true)
			{
				// also sending type of event info atm, not necessary in case of one-type network (+1 short)
				vector<UnitModifier*> events = localUnits[j]->GetEvents();
				
				vector<int>* postProcesses = localUnits[j]->GetPostProcesses();
				
				long eventId;
				short eventType;
				long eventHypercolumnId;

				for(int k=0;k<events.size();k++)
				{
					eventType = events[k]->GetEventTypeId();
					eventId = events[k]->GetFromUnitId();
					eventHypercolumnId = events[k]->GetFromHypercolumnId();

					//eventTypes.push_back(eventType);
					//eventIds.push_back(eventId);
					//eventHypercolumnIds.push_back(eventHypercolumnId);

					vector<float> eData = events[k]->GetEventData();
					//for(int m=0;m<eData.size();m++)
					//	eventData.push_back(eData[m]);

					//for(int n=0;n<m_communicationToProcs[i].size();n++)
					//{
					//if(postProcesses->size()==0)//m_communicationToProcs[i][n] == m_mpiNodeId) // local data, no need to send
					//{
						localEventIds.push_back(eventId);
						localEventData.push_back(eData[0]);
						if(m_isTrackingHypercolumnIds)
							localEventHypercolumnIds.push_back(eventHypercolumnId);
					//}
					if(postProcesses->size() != 0)//else // non-local data
					{
						for(int m = 0;m<postProcesses->size();m++)
						{
							eventDestinationsProcs[postProcesses->at(m)].push_back(currentIndex);
							totalEventIds[postProcesses->at(m)].push_back(eventId);
							if(m_isTrackingHypercolumnIds)
								totalEventHypercolumnIds[postProcesses->at(m)].push_back(eventHypercolumnId);

							// unnecessary atm? all models currently only use one data point
							for(int p=0;p<eData.size();p++)
								totalEventData[postProcesses->at(m)].push_back(eData[p]);
							/*
							eventDestinationsProcs[m_communicationToProcs[i][n]].push_back(currentIndex);
							totalEventIds[m_communicationToProcs[i][n]].push_back(eventId);
							if(m_isTrackingHypercolumnIds)
								totalEventHypercolumnIds[m_communicationToProcs[i][n]].push_back(eventHypercolumnId);

							// unnecessary atm? all models currently only use one data point
							for(int p=0;p<eData.size();p++)
								totalEventData[m_communicationToProcs[i][n]].push_back(eData[p]);*/
						}
						nrEvents++;
					}
					//}

					currentIndex++;
				}

				localUnits[j]->IsNewEvent(false); // will be set true again by by unit's simulateeventque fcn or CreateEvent fcn
			}
		}
	}

	TimingStop("CommunicationGetEvents");

	// change to not send the local events, only external

	// redistribute nr events
	//int nrEvents = eventIds.size();
	int totalNrEvents = 0;
	
	TimingStart("Allreduce");

	MPI_Allreduce(&nrEvents,&totalNrEvents, 1,MPI_INT,MPI_SUM,NETWORK_COMM_WORLD); // could switch to multiple MPI communicators as network becomes large (and is modular enough to benefit).

	TimingStop("Allreduce");

	// will be removed
	if(this->MPIGetNodeId() == 0)
		cout<<"{"<<localEventIds.size()<<"/"<<totalNrEvents<<"} ";

	vector<long> allIds;
	vector<float> allData;
	vector<long> allHypercolumnIds;

	if(totalNrEvents>0)
	{
		//TimingStart("CommunicationBarrierTest");
		//MPI_Barrier(NETWORK_COMM_WORLD);
		//TimingStop("CommunicationBarrierTest");

		//TimingStart("CommunicationReduceEvents");

		// send nr this process wants to send to everyone
		vector<int> nrMessagesSend(m_mpiNrProcs);
		vector<int> nrMessagesRecv(m_mpiNrProcs);
		vector<int> strideSend(m_mpiNrProcs);
		vector<int> strideRecv(nrMessagesRecv.size());
		//int nrMessages = eventIds.size();

		int totMessToSend = 0;
		for(int i=0;i<eventDestinationsProcs.size();i++)
		{
			nrMessagesSend[i] = eventDestinationsProcs[i].size();
			totMessToSend += nrMessagesSend[i];

			if(i>0)
				strideSend[i] = eventDestinationsProcs[i-1].size() + strideSend[i-1];
			else
				strideSend[0] = 0;
		}
		TimingStart("Alltoall1");
		// retrieve nr all processes want to receive
		MPI_Alltoall(&nrMessagesSend[0],1,MPI_INT,&nrMessagesRecv[0],1,MPI_INT,NETWORK_COMM_WORLD);//m_mpiNrProcs
		TimingStop("Alltoall1");

		// use to build up strides etc
		int totMessToRecv = 0;

		for(int i=0;i<nrMessagesRecv.size();i++)
		{
			totMessToRecv+=nrMessagesRecv[i];
			if(i>0)
				strideRecv[i] = nrMessagesRecv[i-1] + strideRecv[i-1];
			else
				strideRecv[i] = 0;
		}

		vector<long> sendIds(totMessToSend);
		vector<float> sendData(totMessToSend); // assumes only one message sent
		vector<long> sendHypercolumnIds(totMessToSend);

		currentIndex = 0;
		for(int i=0;i<totalEventIds.size();i++)
		{
			for(int j=0;j<totalEventIds[i].size();j++)
			{
				sendIds[currentIndex] = totalEventIds[i][j];
				sendData[currentIndex] = totalEventData[i][j];
				if(m_isTrackingHypercolumnIds)
					sendHypercolumnIds[currentIndex] = totalEventHypercolumnIds[i][j];
				currentIndex++;
			}
		}

		allIds = vector<long>(totMessToRecv);
		allData = vector<float>(totMessToRecv);
		
		TimingStart("Alltoall2");
		// all-scatter/-gather ids
		if(sendIds.size()==0)
			MPI_Alltoallv(NULL, &nrMessagesSend[0], &strideSend[0], MPI_LONG, &allIds[0], &nrMessagesRecv[0], &strideRecv[0], MPI_LONG, NETWORK_COMM_WORLD);
		else
			MPI_Alltoallv(&sendIds[0], &nrMessagesSend[0], &strideSend[0], MPI_LONG, &allIds[0], &nrMessagesRecv[0], &strideRecv[0], MPI_LONG, NETWORK_COMM_WORLD);

		// all-scatter/-gather data
		if(sendData.size()==0)
			MPI_Alltoallv(NULL, &nrMessagesSend[0], &strideSend[0], MPI_FLOAT, &allData[0], &nrMessagesRecv[0], &strideRecv[0], MPI_FLOAT, NETWORK_COMM_WORLD);
		else
			MPI_Alltoallv(&sendData[0], &nrMessagesSend[0], &strideSend[0], MPI_FLOAT, &allData[0], &nrMessagesRecv[0], &strideRecv[0], MPI_FLOAT, NETWORK_COMM_WORLD);

		// all-scatter/-gather hypercolumn ids if applicable (bcpnn uses)
		if(m_isTrackingHypercolumnIds)
		{
			allHypercolumnIds = vector<long>(totMessToRecv);
			if(sendHypercolumnIds.size()==0)
				MPI_Alltoallv(NULL, &nrMessagesSend[0], &strideSend[0], MPI_LONG, &allHypercolumnIds[0], &nrMessagesRecv[0], &strideRecv[0], MPI_LONG, NETWORK_COMM_WORLD);
			else
				MPI_Alltoallv(&sendHypercolumnIds[0], &nrMessagesSend[0], &strideSend[0], MPI_LONG, &allHypercolumnIds[0], &nrMessagesRecv[0], &strideRecv[0], MPI_LONG, NETWORK_COMM_WORLD);
		}

		TimingStop("Alltoall2");

	}

	// add local data to buffers
	for(int i=0;i<localEventIds.size();i++)
	{
		allIds.push_back(localEventIds[i]);
		allData.push_back(localEventData[i]);
		if(m_isTrackingHypercolumnIds)
			allHypercolumnIds.push_back(localEventHypercolumnIds[i]);
	}


	TimingStart("HashBuffersFilling");

	m_incomingBufferDataHash.clear();
	m_incomingBufferHypercolumnIdsHash.clear();

	if(this->m_keepCommunicationBuffer == true) // used to lower communication in case of static input, e.g. input data from file over some time steps
	{
		m_incomingBufferDataHash = this->m_bufferToKeepData;
		m_incomingBufferHypercolumnIdsHash = this->m_bufferToKeepHypercolumnIds; 
	}

	if(m_isTrackingHypercolumnIds) // set on by bcpnn
	{
#if USE_UNORDERED_MAP == 1
		m_incomingBufferDataHash.rehash(allIds.size() + m_incomingBufferDataHash.size());
		m_incomingBufferHypercolumnIdsHash.rehash(allData.size() + m_incomingBufferHypercolumnIdsHash.size());
#endif

		for(int i=0;i<allIds.size();i++) // check how slow (?)
		{
			m_incomingBufferDataHash[allIds[i]] = allData[i];
			m_incomingBufferHypercolumnIdsHash[allIds[i]] = allHypercolumnIds[i]; // when used?
		}
	}
	else
	{
#if USE_UNORDERED_MAP == 1
		m_incomingBufferDataHash.rehash(allData.size() + m_incomingBufferDataHash.size());
#endif

		for(int i=0;i<allIds.size();i++) // check how slow (?)
		{
			m_incomingBufferDataHash[allIds[i]] = allData[i];
		}
	}

	TimingStop("HashBuffersFilling");

	TimingStart("HashBuffersAddEvent");

#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	for(int i=0;i<m_populations.size();i++)
	{
		if(m_populations[i]->GetUnitsAll()->size()>0 && m_populations[i]->IsOn()) // only continue if this process is part of population
		{
			vector<Connection*> conns = m_populations[i]->GetIncomingConnections();

			for(int j=0;j<conns.size();j++)
			{
				if(conns[j]->KeepActiveBuffer() == false)
				{
					conns[j]->ClearActiveBuffer();

					//for(int k=0;k<allIds.size();k++)
					//	m_populations[i]->GetIncomingConnections()[j]->AddActiveEvent(allIds[k],allData[k]);//,allData[k]);
					if(allIds.size()>0)
						conns[j]->AddActiveEvents(allIds,allData);
				}
			}
		}
	}
#endif
	TimingStop("HashBuffersAddEvent");

	TimingStop("SendAndReceiveVersionAlltoall");
}



long Network::AddUnitModifierIncoming(UnitModifier* e)
{
	long id = e->GetFromUnitId();
	delete m_eventsUnitIncoming[id];
	m_eventsUnitIncoming[id] = e;
	return id;
	//m_eventsUnitIncoming.push_back(e);
	//return m_eventsUnitIncoming.size() - 1;
}

// Recording routine
void Network::AddMeter(Meter* meter)
{
	m_meters.push_back(meter);
	meter->network(this);
}

UnitModifier* Network::GetUnitModifierIncoming(int id)
{
	return m_eventsUnitIncoming[id];
	//return m_eventsUnitIncoming[index];
}

// not used atm.
void Network::ClearMemory()
{
	// moved to dispose

/*	int nr = m_eventsUnitIncoming.size();
	for(int i=0;i<nr;i++)
	{
		delete m_eventsUnitIncoming[i];
		m_eventsUnitIncoming[i] = NULL;
	}

	m_eventsUnitIncoming.clear();
	*/
}

void Network::Dispose() // also used for reset
{
	this->ClearEventsIncoming();

	//delete m_filenameNetworkDetails;
	//delete m_filenameTimings;
	//delete m_filenameAnalysis;

	m_hashSynapses.clear();

	for(int i=0;i<m_populations.size();i++)
	{
		delete m_populations[i];
	}

	for(int i=0;i<m_meters.size();i++)
		delete m_meters[i];

	for(int i=0;i<m_analysis.size();i++)
		delete m_analysis[i];

	for(int i=0;i<m_connectivityTypes.size();i++)
		delete m_connectivityTypes[i];

	for(int i=0;i<m_networkObjectsToDelete.size();i++)
	{
		delete m_networkObjectsToDelete[i];
		m_networkObjectsToDelete[i] = NULL;
	}

	m_populations.clear();
	m_meters.clear();
	m_analysis.clear();
	m_connectivityTypes.clear();
	m_networkObjectsToDelete.clear();

	/*for(int i=0;i<m_populations.size();i++)
	{
		m_populations[i]->Dispose();
	}

	m_firstRun = true;
	m_currentUnitId = 0;
	m_currentHypercolumnId = 0;
	m_currentLayerId = 0;
	m_currentTimeStep = 0;

	// meters - should they be disposed here?

	for(int i=0;i<m_meters.size();i++)
	{
		delete m_meters[i];
	}
	
	m_meters.clear();*/

}

void Network::Reset()
{
	this->ClearEventsIncoming();

	//delete m_filenameNetworkDetails;
	//delete m_filenameTimings;
	//delete m_filenameAnalysis;

	m_hashSynapses.clear();

	for(int i=0;i<m_populations.size();i++)
	{
		delete m_populations[i];
	}

	for(int i=0;i<m_meters.size();i++)
		delete m_meters[i];

//	for(int i=0;i<m_analysis.size();i++)
//		delete m_analysis[i];
	m_analysis.clear(); // Error in timing if deleting, change

	for(int i=0;i<m_connectivityTypes.size();i++)
		delete m_connectivityTypes[i];

	for(int i=0;i<m_networkObjectsToDelete.size();i++)
	{
		delete m_networkObjectsToDelete[i];
		m_networkObjectsToDelete[i] = NULL;
	}

	m_populations.clear();
	m_meters.clear();
	m_analysis.clear();
	m_connectivityTypes.clear();
	m_networkObjectsToDelete.clear();

	this->Parameters()->Reset();
}

// Remove all the current events/messages in this time step which are about to arrive
// - use to reset the communication in a network, e.g. if stimuli is changed or something
void Network::ClearEventsIncoming()
{
	for(int i=0;i<m_populations.size();i++)
		m_populations[i]->ClearEventsIncoming();

	// necessary?
	map<long,UnitModifier*>::iterator p;

	for(p = m_eventsUnitIncoming.begin(); p != m_eventsUnitIncoming.end(); p++)
	{
		delete m_eventsUnitIncoming[p->first];
	}

	m_eventsUnitIncoming.clear();
#if COMMUNICATION_BUFFER_HASH == 1
	this->m_incomingBufferDataHash.clear();
	this->m_incomingBufferDataDelays.clear();
	this->m_incomingBufferHypercolumnIdsHash.clear();
#else
	this->m_incomingBufferData = vector<float>(m_incomingBufferData.size(),0.0);
	this->m_incomingBufferDataDelays = vector<vector<float> >(m_incomingBufferDataDelays.size(),vector<float>(0));
#endif
}

void Network::ClearEventsIncoming(vector<long> fromPreIds)
{
	for(int i=0;i<fromPreIds.size();i++)
	{
#if COMMUNICATION_BUFFER_HASH == 1
		m_incomingBufferDataHash[fromPreIds[i]] = 0.0;
#else
		m_incomingBufferData[fromPreIds[i]] = 0.0;
#endif
	}
}

void ConnectionModifier::SetConnection(Connection* c)
{
	m_connection = c;
}

// Adds population of units to this network
void Network::AddLayer(Population* layer)
{
	layer->network(this);

	layer->SetLayerId(m_currentLayerId);
	stringstream ss;
	ss<<"Layer"<<m_currentLayerId;
	if(layer->GetName().size() == 0)
		layer->SetName(ss.str());

	m_populations.push_back(layer);
	m_currentLayerId++;
}

/*void Network::AddUnitToHash(Unit* unit)
{
	m_hashIdUnit[unit->GetUnitId()] = unit;
}*/

// Hashes a unit
void Network::AddUnit(Unit* unit)
{
	m_hashIdUnit[unit->GetUnitId()] = unit;
	/*if(m_listIdUnit.size()<unit->GetUnitId()+1)
	{
		for(int i=m_listIdUnit.size();i<unit->GetUnitId()+1;i++)
			m_listIdUnit.push_back(NULL);
	}

	m_listIdUnit[unit->GetUnitId()] = unit;*/
}

int Network::GetSeed()
{
	int seed;

	if(MPIGetNodeId() == 0)
	{
		seed = (unsigned)time(NULL);
		MPI_Bcast(&seed,1,MPI_UNSIGNED,0,NETWORK_COMM_WORLD);
	}
	else
	{
		MPI_Bcast(&seed,1,MPI_UNSIGNED,0,NETWORK_COMM_WORLD);
	}

	return seed;
}

// In-between simulation steps or post-run: Save all output data to new files or add to existing files (depending on how meter was set up)
void Network::RecordAll()
{
	TimingStart("RecordAll");

	// resaving network details (filenames, network structure may have changed)
	if(this->MPIGetNodeId() == 0)
		StoreNetworkDetails();

	for(int i=0;i<m_meters.size();i++)
	{
		if(this->m_useExtraFilenameString == true)
			m_meters[i]->SetExtraFilenameString(m_extraFilenameString);
		m_meters[i]->RecordAll(this->GetCurrentTimeStep());
	}

	TimingStop("RecordAll");
}

// Stores the details about the network setup to file, resulting file currently not used to load in details at a later time.
void Network::StoreNetworkDetails()
{
	if(true)//m_networkDetailsStored == false)
	{
		ofstream myfile;
		myfile.open (m_filenameNetworkDetails);
		myfile << "Network Details - do not change manually\n";

		// computer info
		myfile<<"Nodes: "<< this->MPIGetNrProcs()<<"\n"; // imp for nr connections files

		// Layers
		myfile << "Layers: "<< m_populations.size() << "\n";
		//	myfile << "Connections: not ava\n";

		ofstream geomfile;
		bool isGeomOpen = false;

		for(int i=0;i<m_populations.size();i++)
		{
			myfile<<"Id: "<<m_populations[i]->GetLayerId();
			myfile<<" Type: "<<m_populations[i]->GetPopulationType();
			myfile<<" NrUnits: "<<m_populations[i]->GetNrUnitsTotal();

			// should be layer type specific
			int nrHcs = ((PopulationColumns*)m_populations[i])->GetHypercolumns().size();
			myfile<<" NrHypercolumns: "<<nrHcs;
			myfile<<" NrRateUnits:";
			for(int j=0;j<nrHcs;j++)
			{
				int nrMcs = ((PopulationColumns*)m_populations[i])->GetHypercolumns()[j]->GetRateUnits().size();
				myfile<<" "<<nrMcs;
			}

			// replace to get from connection
			if(m_populations[i]->GetUnits().size()>0 && m_populations[i]->GetIncomingConnections().size()>0)
				if(m_populations[i]->GetIncomingConnections()[0]->GetUnitModifier("geometry") != NULL)//if(m_populations[i]->GetUnits()[0]->GetUnitModifier("geometry") != NULL)
				{
					if(isGeomOpen == false)
					{
						geomfile.open ("geometry.txt");
						isGeomOpen = true;
					}

					//	myfile<<" Geometries: ";
					vector<RateUnit*> mcs = ((PopulationColumns*)m_populations[i])->GetRateUnits();
					for(int j=0;j<mcs.size();j++)
					{
						geomfile<<((GeometryUnit*)m_populations[i]->GetIncomingConnections()[0]->GetUnitModifier(7))->GetValuesAsString()<<"\n";//((GeometryUnit*)mcs[j]->GetUnitModifier(7))->GetValuesAsString()<<"\n";
						//		myfile<<((GeometryUnit*)m_populations[i]->GetUnits()[j]->GetUnitModifier(7))->GetValuesAsString()<<" ";
					}
				}

				myfile<<"\n";
		}

		if(isGeomOpen)
			geomfile.close();

		// Meters info
		myfile << "Meters: "<<m_meters.size()<<"\n";

		for(int i=0;i<m_meters.size();i++)
		{
			vector<Meter::MeterType> meterTypes = m_meters[i]->GetMeterTypes();
			myfile << "Filename: "<<m_meters[i]->GetFilename()<<" Type: ";
			for(int j=0;j<meterTypes.size();j++)
			{
				myfile<<meterTypes[j];
				if(j!=meterTypes.size()-1)
					myfile<<" ";
			}

			vector<int> additionalInfo = m_meters[i]->GetAdditionalInfo();
			if(additionalInfo.size()>0)
			{
				myfile<<" Info: ";

				for(int j=0;j<additionalInfo.size();j++)
				{
					myfile<<additionalInfo[j];
					if(j!=additionalInfo.size()-1)
						myfile<<" ";
				}
			}

			myfile<<"\n";
		}

		myfile.close();

		m_networkDetailsStored = true;
	}
}

vector<float> Network::PrintNetworkDetails()
{
	vector<float> out;

	if(this->MPIGetNodeId() == 0)
	{
		cout<<"\nNetwork Details\n";
		cout<<"----------------\n";
		cout<<"Total nr nodes = "<< this->MPIGetNrProcs()<<"\n";
		cout<<"Total non-virtual units = "<< m_hashIdUnit.size()<<"\n";//m_listIdUnit.size()<<"\n"; // will change to only local
		out.push_back(this->MPIGetNrProcs());
		out.push_back(m_hashIdUnit.size());//m_listIdUnit.size());
	}	

	long sum = 0;
	long nrUnits = 0;

	for(int i=0;i<this->m_populations.size();i++)
	{
		vector<Connection*> inconns = m_populations[i]->GetIncomingConnections();
		vector<Unit*> units = m_populations[i]->GetLocalUnits();
		nrUnits+=units.size();

		/*for(int j=0;j<inconns.size();j++)
		{
			vector<vector<long>* > ids = inconns[j]->GetPreIds();
			for(int k=0;k<ids.size();k++)
				sum+= (*ids[k]).size();
		}*/
	}

	long totSum = 0;
	long totUnits = 0;
	MPI_Reduce(&sum,&totSum,1,MPI_LONG,MPI_SUM,0,NETWORK_COMM_WORLD);
	MPI_Reduce(&nrUnits,&totUnits,1,MPI_LONG,MPI_SUM,0,NETWORK_COMM_WORLD);

	if(this->MPIGetNodeId() == 0)
	{
		cout<<"Total nr units = "<<totUnits<<"\n";
		cout<<"Av units/node = "<<(float)totUnits/(float)this->MPIGetNrProcs()<<"\n";
		cout<<"Connections = "<<totSum<<"\n";
		cout<<"Av connections/node = "<<(float)totSum/(float)this->MPIGetNrProcs()<<"\n";
		cout<<"Av connections/unit = "<<(float)totSum/(float)totUnits<<"\n";

		out.push_back(totUnits);
		out.push_back(totSum);
		cout.flush();
	}

	MPI_Barrier(NETWORK_COMM_WORLD);
	
	cout<<this->MPIGetNodeId()<<":"<<nrUnits<<" ";
	cout.flush();
	
	MPI_Barrier(NETWORK_COMM_WORLD);

	if(this->MPIGetNodeId() == 0)
	{
		cout<<"\n\n";
		cout.flush();
	}

	return out;
}

void Network::StoreAnalysis()
{
	ofstream myfile;

	if(this->MPIGetNodeId() == 0)
	{
		cout<<"\nStoring analysis data...";
		cout.flush();

		stringstream ss;
		myfile.open (m_filenameAnalysis);
	}

	vector<vector<float> > totData;
	vector<string> names;

	for(int i=0;i<m_analysis.size();i++)
	{
		m_analysis[i]->Finalize();	// could move to global finalizer as it may be used not only in storing results.
		// either finalizes or finalizes + writes data
		string name = m_analysis[i]->GetName();

		if(strcmp(name.c_str(),"Timing")!=0) // special case
		{
			names.push_back(name);
			if(this->MPIGetNodeId() == 0)
			{
				cout<<name<<"..";cout.flush();
			}
		}

		if(this->MPIGetNodeId() == 0)
		{
			vector<vector<float> > results = m_analysis[i]->GetResults();
			cout<<"#="<<results.size()<<"..";cout.flush();
			// first column is typically time ?
			if(results.size()>0)
			{
				for(int m=0;m<results.size();m++)
				{
					if(totData.size()<=m)
						totData.push_back(vector<float>());

					for(int n=0;n<results[m].size();n++)
					{
						totData[m].push_back(results[m][n]);
					}
				}
			}
		}
	}

	if(this->MPIGetNodeId() == 0)
	{
		if(totData.size()>0)
		{
			for(int i=0;i<names.size();i++)
			{
				myfile<<names[i];
				if(i!=names.size()-1)
					myfile<<",";
				else
					myfile<<"\n";
			}

			for(int i=0;i<totData.size();i++)
			{
				for(int j=0;j<totData[i].size();j++)
				{
					myfile<<totData[i][j];
					if(j!=totData[i].size()-1)
						myfile<<",";
				}
				if(i!=totData.size()-1)
					myfile<<"\n";
			}

			cout<<"ok.";cout.flush();
		}
		else cout<<"no data.";

		
		myfile.close();
	}

	// additional analysis classes here

	StoreTimings();
}

void Network::StoreTimings()
{
	if(this->MPIGetNodeId() == 0)
	{
		cout<<"\nPrinting network timings...";

		ofstream myfile;
		stringstream ss;
		myfile.open (m_filenameTimings);

		if(this->IsTiming())
		{
			myfile<<this->Timing()->GetTimingString();/*"Network, ";
											for(int i=0;i<this->GetTimingTotal().size();i++)
											{
											myfile<<this->GetTimingTotal()[i];
											if(i!=this->GetTimingTotal().size()-1)
											myfile<<", ";
											}*/
		}

		for(int i=0;i<m_populations.size();i++)
		{
			if(m_populations[i]->IsTiming())
			{
				string s = m_populations[i]->Timing()->GetTimingString();
				myfile<<i+1<<". "<<s;
			}

			for(int j=0;j<m_populations[i]->GetPopulationModifiers().size();j++)
			{
				if(m_populations[i]->GetPopulationModifiers()[j]->IsTiming())
				{
					string s = m_populations[i]->GetPopulationModifiers()[j]->Timing()->GetTimingString();
					myfile<<i+1<<".A"<<j+1<<" "<<s;
				}
			}

			for(int j=0;j<m_populations[i]->GetPreConnectivitys().size();j++)
			{
				if(m_populations[i]->GetPreConnectivitys()[j]->IsTiming())
				{
					string s = m_populations[i]->GetPreConnectivitys()[j]->Timing()->GetTimingString();
					myfile<<i+1<<".B"<<j+1<<" "<<s;
				}
			}

			for(int j=0;j<m_populations[i]->GetIncomingConnections().size();j++)
			{
				for(int k=0;k<m_populations[i]->GetIncomingConnections()[j]->GetConnectionModifiers().size();k++)
				{
					if(m_populations[i]->GetIncomingConnections()[j]->GetConnectionModifiers()[k]->IsTiming())
					{
						string s = m_populations[i]->GetIncomingConnections()[j]->GetConnectionModifiers()[k]->Timing()->GetTimingString();
						myfile<<i+1<<".C"<<j+1<<" "<<s;
					}
				}
			}
		}

		myfile.close();

		cout<<"ok.\n";
	}
}

// Hashed weights and delays put in main network
void Network::SetWeight(float weight, long preId, long postId)
{
	/*if(m_hashSynapses.size()<postId+1)
	{
		while(m_hashSynapses.size()<postId+1)
			m_hashSynapses.push_back(map<long,SynapseStandard>());
	}*/

	m_hashSynapses[postId][preId].weight = weight;
}

// may get replaced
void Network::SetDelay(float delay, long preId, long postId)
{
	/*if(m_hashSynapses.size()<postId+1)
	{
		while(m_hashSynapses.size()<postId+1)
			m_hashSynapses.push_back(map<long,SynapseStandard>());
	}*/
#if USE_DELAYS==1
	m_hashSynapses[postId][preId].delay = delay;
#endif
}

float Network::GetWeight(long preId, long postId)
{
	return m_hashSynapses[postId][preId].weight;
}

float Network::GetDelay(long preId, long postId)
{
#if USE_DELAYS==1
	return m_hashSynapses[postId][preId].delay;
#else
	return 0;
#endif
}

// not used anymore - remove
bool Network::ConnectionExists(long preId, long postId)
{
	/*if(m_hashSynapses.count(postId)>0)
		if(m_hashSynapses[postId].count(preId)>0)
			return true;
	*/
	return false;
}


/*void Network::AddPostIds(vector<long> preIds, long postId)
{
	for(int i=0;i<preIds.size();i++)
	{
		if(m_listPosts.size()<preIds[i]+1) // change
		{
			for(int j=m_listPosts.size();j<preIds[i]+1;j++)
				m_listPosts.push_back(vector<long>());
		}

		//used? remove
		m_listPosts[preIds[i]].push_back(postId);
	}
}*/	

// Run in cases we need to have easy access to all pre unit-ids
// can optimize if needed
// need to call every rebuild
// (or add a function that can add the new connection pre id if needed)
void Network::CreateAllPreIdsUnion()
{
	m_allPreIds.clear();

#if USE_UNORDERED_MAP == 1
	unordered_map<long,int> tempMap;
	unordered_map<long, unordered_map<long, SynapseStandard> >::iterator it1; // post
	unordered_map<long, SynapseStandard>::iterator it2; //pre
#else
	map<long,int> tempMap;
	map<long, map<long, SynapseStandard> >::iterator it1; // post
	map<long, SynapseStandard>::iterator it2; //pre
#endif
	/*for(it1 = m_hashSynapses.begin();it1!=m_hashSynapses.end();it1++)
	{
		for(it2 = it1->second.begin();it2!=it1->second.end();it2++)
		{*/

	// high memory usage (!)
	long totalLocalSynapses = 0;
//	for(int i=0;i<m_hashSynapses.size();i++)
//	{
	for(it1 = m_hashSynapses.begin();it1!=m_hashSynapses.end();it1++)
	{
		for(it2= it1->second.begin();it2!=it1->second.end();it2++)//it2 = m_hashSynapses[i].begin();it2!=m_hashSynapses[i].end();it2++)
		{
			// no duplicates (not using nr times)
			//if(find(m_allPreIds.begin(),m_allPreIds.end(),it2->first) == m_allPreIds.end())
//				m_allPreIds.push_back(it2->first);
			tempMap[it2->first] = tempMap[it2->first]+1;
			totalLocalSynapses++;
		}
	}
#if USE_UNORDERED_MAP == 1	
	unordered_map<long,int>::iterator it3;
#else
	map<long,int>::iterator it3;
#endif

	for(it3 = tempMap.begin();it3!=tempMap.end();it3++)
		m_allPreIds.push_back(it3->first);

	// sort
	sort(m_allPreIds.begin(),m_allPreIds.end());

	if(m_allPreIds.size()>0)
	{
#if COMMUNICATION_BUFFER_HASH == 0
		m_incomingBufferData.resize(m_allPreIds[m_allPreIds.size()-1]+1);// = vector<float>(m_allPreIds[m_allPreIds.end()-1],0.0); // change to ... -m_allPreIds[0] to only get necessary nr
		m_incomingBufferHypercolumnIds.resize(m_allPreIds[m_allPreIds.size()-1]+1);
#endif
	}

	if(this->MPIGetNodeId()==0)
	{
		cout<<"Total nr local synapses (process 0) == "<<totalLocalSynapses;
	}

	if(this->MPIGetNodeId()==1)
	{
		cout<<"Total nr local synapses (process 1) == "<<totalLocalSynapses;
	}

	if(this->MPIGetNodeId()==this->MPIGetNrProcs()-1)
	{
		cout<<"Total nr local synapses (process "<<this->MPIGetNrProcs()-1<<") == "<<totalLocalSynapses;
	}
}

// optimize (!)
void Network::CreateAllPostProcs()
{
	for(int i=0;i<this->MPIGetNrProcs();i++)
	{
		if(i == this->MPIGetNodeId())
		{
			int nrPreIds = m_allPreIds.size();
			MPI_Bcast(&nrPreIds,1,MPI_INT,i,NETWORK_COMM_WORLD);

			MPI_Bcast(&m_allPreIds[0],nrPreIds,MPI_LONG,i,NETWORK_COMM_WORLD);
		}
		else
		{
			int nrPreIds;
			MPI_Bcast(&nrPreIds,1,MPI_INT,i,NETWORK_COMM_WORLD);

			vector<long> preIds(nrPreIds);
			MPI_Bcast(&preIds[0],nrPreIds,MPI_LONG,i,NETWORK_COMM_WORLD);

			//for(int j=0;j<preIds.size();j++)
			//{
				for(int k=0;k<m_populationIndexesThisProc.size();k++)
				{
					// if no connectivity in between populations, skip checking
					//bool isConnectivityInBetween = true;
					//for(int n=0;n<m_populations[m_populationIndexesThisProc[k]]->Float

					//if(isConnectivityInBetween)
					//{
						vector<Unit*>* units = m_populations[m_populationIndexesThisProc[k]]->GetUnitsAll();
						long unitId;
						for(int m=0;m<units->size();m++)
						{
							unitId = units->at(m)->GetUnitId();
							if(binary_search(preIds.begin(),preIds.end(),unitId))//find(preIds.begin(),preIds.end(),units->at(m)->GetUnitId()) != preIds.end())
							//if(preIds[j] == units->at(m)->GetUnitId())
							{
								units->at(m)->AddPostProcess(i);
							}
						}
					//}
				}
			//}
		}
	}
}
/*
void Network::AddMeter(char* filename, NetworkObject* object, Meter::MeterType type)
{
	Meter* meter = new Meter(filename, Storage::CSV);
	meter->AttachObject(object,type);
	m_meters.push_back(meter);
}
*/

// Copies post synapses, same pre
// - useful to copy a part of the network to another part (learnt set of connections to an empty set of connections at another place)
void Network::CopyConnectionsPost(Connection* from, Connection* to, bool copyValues)
{
	// assumes linearly distributed / sorted unit ids
	map<long, SynapseStandard>::iterator it;

	vector<long> postIds = from->GetPostIds();
	vector<long> newPostIds = to->GetPostIds();
	vector<long> newPostLocalIds = to->GetPostLocalIds();
	to->Clear(); // erases synapses

	for(int i=0;i<postIds.size();i++)
	{
		long postId = from->GetPostIds()[i];
		vector<long> preIds = from->GetPreIds(postId);//this->GetPreSynapses(postId);
		to->AddConnections(preIds,newPostIds[i],newPostLocalIds[i]);

		for(int j=0;j<preIds.size();j++)
		{
			if(copyValues)
			{
				float weight = GetWeight(preIds[j],postIds[i]);
				this->SetWeight(weight,preIds[j],newPostIds[i]);
			}
			else
			{
				this->SetWeight(0.0,preIds[j],newPostIds[i]);
			}
		}
	}

	to->CreatePreIdsUnion();
}


/// To reduce communication, keep buffer (only!) next time step
void Network::KeepCommunicationBuffer(long startId, long endId)
{
	this->m_keepCommunicationBuffer = true;

#if USE_UNORDERED_MAP == 1
	unordered_map<long,float>::iterator it1;
	unordered_map<long,long>::iterator it2;
#else
	map<long,float>::iterator it1;
	map<long,long>::iterator it2;
#endif

	m_bufferToKeepData.clear();
	m_bufferToKeepHypercolumnIds.clear();

	for(it1 = m_incomingBufferDataHash.begin();it1!=m_incomingBufferDataHash.end();it1++)
	{
		if(it1->first>=startId && it1->first <=endId)
		{
			m_bufferToKeepData[it1->first] = it1->second;
		}
	}

	for(it2 = m_incomingBufferHypercolumnIdsHash.begin();it2!=m_incomingBufferHypercolumnIdsHash.end();it2++)
	{
		if(it2->first>=startId && it2->first <=endId)
		{
			m_bufferToKeepHypercolumnIds[it2->first] = it2->second;
		}
	}
}

// Main run file, overridable
void Network::Run()
{
	this->NetworkSetupStructure();

	// retrieve extra parameters here
	this->NetworkSetupParameters();

	m_runId =0;
	char* extraString;
	extraString = new char[50];
	

	if(this->Parameters()->ParametersLeft())
	{	
		sprintf(extraString,"%d",m_runId);
		this->SetExtraFilenameString(extraString);
		this->Parameters()->SetNextParameters();
		m_runId++;
	}

	this->Initialize();
	this->NetworkSetupMeters();

	this->NetworkRun();

	
	// Rerun network if multiple parameters specified in setup
	while(this->Parameters()->ParametersLeft() == true)
	{
		this->Reset();
		this->NetworkSetupStructure();

		// retrieve extra parameters here
		this->NetworkSetupParameters();

		sprintf(extraString,"%d",m_runId);
		this->SetExtraFilenameString(extraString);
		this->Parameters()->SetNextParameters();

		this->Initialize();
		this->NetworkSetupMeters();

		this->NetworkRun();
		m_runId++;
	}

/*	if(true)//m_parameterValues.size() == 0) // only run for original parameters setup in NetworkSetupStructure
	{
		this->Initialize();
		this->NetworkSetupMeters();
		this->NetworkRun();

		//~this();
	}
	else
	{
		// run as defined in SetupParameters
		for(int i=0;i<1;i++)//m_parameterValues.size();i++)
		{
			// set current file names as defined from parameters
			//this->SetExtraFilenameString(...);

			this->Initialize();
			this->NetworkSetupMeters();
			this->NetworkRun();

			//~this();
		}
	}*/
}