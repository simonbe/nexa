#include <mpi.h>
#include <math.h>
#include <ctime>
#include <cstring>
#include <fstream>
#include <deque>
#include "Network.h"
#include "Meter.h"

// MUSIC specifics
#if MUSIC_AVAILABLE == 1
MPI_Comm Network::MPIComm;
#endif

/// <summary>	Network constructor. Sets start values of ids, flags, default filenames etc. </summary>

Network::Network()
{
	// MPI_Init needs to be run before Network object created

	MPI_Comm_size(NETWORK_COMM_WORLD, &m_mpiNrProcs);
	MPI_Comm_rank(NETWORK_COMM_WORLD, &m_mpiNodeId);

	m_currentUnitId = 0;
	m_currentHypercolumnId = 0;
	m_currentLayerId = 0;
	m_currentTimeStep = 0;
	m_simulationResolution = 0.1;
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

	m_name = "Network";
	m_useExtraFilenameString = false;
	m_networkParameters = new NetworkParameters();
	
	////// Parallelization scheme trackers
	m_mpiNrUnitsParallelizationDefault = 0; // keep track of division for default parallelization scheme
	m_mpiCurrentUnitParallelizationDefault = 0;

	this->SetSeed(false, GetSeed()); // default seeds (different important for e.g. random connectivities),
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

/// <summary>	Sets random number generator seed. 
/// 			- Can use srand anywhere in code instead.
/// 			- Some connectivitys use it to set same seed on all processes during network initialization. </summary>
///
/// <param name="allSame">	true to set same seed on all processes. </param>
/// <param name="x">	  	Seed value. </param>

void Network::SetSeed(bool allSame, float x)
{
	if(allSame == false)
		srand(this->MPIGetNodeId() + m_runId + x);
	else
		srand(x);
}

/// <summary>	Destructor. Called during simulation for multiple parameters runs.
/// 			TODO: Check difference to Dispose. </summary>

Network::~Network()
{
	this->ClearEventsIncoming();

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

/// <summary>	 Adds extra string to files that are outputted by default (such as NetworkTimings.txt if timing is used).
/// 			 For multiple parameters runs, this is called for each new independent run to get different filenames. 
/// 			 TODO: May get moved to meter or some utils class instead.</summary>
///
/// <param name="additional">	Additional filename string. </param>

void Network::SetFilenamesAdditional(string additional)
{
	m_filenameNetworkDetails = new char[100];
	m_filenameTimings = new char[100];
	m_filenameAnalysis = new char[100];
	sprintf(m_filenameNetworkDetails,"NetworkDetails%d_%s.txt",m_mpiNrProcs, additional.c_str());
	sprintf(m_filenameTimings,"NetworkTimings%d_%s.txt",m_mpiNrProcs,additional.c_str());
	sprintf(m_filenameAnalysis,"NetworkAnalysis%d_%s.txt",m_mpiNrProcs,additional.c_str());
}

/// <summary>	Builds network structures and initializes all units. </summary>

void Network::Initialize()
{
	MPI_Barrier(NETWORK_COMM_WORLD);
	TimingStart("InitializeNetwork");
	int processId = this->MPIGetNodeId();

	for(int i=0;i<m_populations.size();i++)
	{
		m_populations[i]->Initialize();
	}

	// create communicators for populations 
	// TODO: Check if this has any effect any longer (called anyway when asked to be used in a simulation)
	for(int i=0;i<m_populations.size();i++)
	{
		m_populations[i]->MPI()->MPICreateCommLayer();
	}

	for(int i=0;i<m_populations.size();i++)
	{
		m_populations[i]->InitializeProjectionsEventsAndParameters();
	}

	// optimization cache of all possible input ids (derived from hashed synapses)
	// TODO: not necessary to have in most cases (and uses extra memory) so could be made optional
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
		if(m_populations[i]->GetUnitsAll()->size() != 0)
			m_populationIndexesThisProc.push_back(i);
	}
	
	TimingStart("CreateAllPostProcs");	

	if(this->MPIGetNodeId() == 0)
	{
		cout<<"Create tables all post processes...";
		cout.flush();
	}

	// create which post processes to send each event output from a unit
	// TODO: could be made optional, not needed for all setups
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

/// <summary>	Simulates a number of timesteps after each other. </summary>
///
/// <param name="nrTimesteps">	Number of timesteps to simulate. </param>

void Network::Simulate(int nrTimesteps)
{
	for(int i=0;i<nrTimesteps;i++)
		Simulate();
}

/// <summary>	Advances simulation one time step </summary>

void Network::Simulate()
{
	if(m_firstRun == true)
	{
		// run first time step in simulation

		// store network details to disk
		if(this->MPIGetNodeId() == 0)
		{
			StoreNetworkDetails();
		}

		MPI_Barrier(NETWORK_COMM_WORLD); // to get start of timing right, may be removed
	}

	TimingStart(m_name);

	// Action

	TimingStart("Simulate");

	for(int i=0;i<m_populationIndexesThisProc.size();i++)
	{
		m_populations[m_populationIndexesThisProc[i]]->Simulate();
	}
	TimingStop("Simulate");

	// additional timing for comm
	
#if USE_COMMUNICATION_ALLTOALL == 1 // TODO: not needed to be on compile level
	CommunicationVersionAlltoall();  // Communicates outgoing activity only to receiving populations
#else
	CommunicationVersionAllgather(); // Communicates all outgoing activity to all processes
#endif

	// Call population and Projection modifiers, e.g. synaptic plasticity
	TimingStart("Modify");
	for(int i=0;i<m_populationIndexesThisProc.size();i++)//m_populations.size();i++)
	{
		m_populations[m_populationIndexesThisProc[i]]->Modify();
	}
	TimingStop("Modify");

	// Call analysis classes
	TimingStart("Analysis");
	for(int i=0;i<this->m_analysis.size();i++)
	{
		if(m_analysis[i]->IsOn())
			m_analysis[i]->Simulate();
	}
	TimingStop("Analysis");

	m_currentTimeStep+=this->GetTimeResolution();

	TimingStop(m_name);

	if(m_firstRun)
		m_firstRun = false;
}

/// <summary>	Communication routine which uses Allgather as main mpi routine
/// 			 - will transmit outgoing messages to all processes. </summary>

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

			for(int m=0;m<m_populations[m_populationIndexesThisProc[i]]->GetOutgoingProjections().size();m++)
			{
				vector<int> prcs = m_populations[m_populationIndexesThisProc[i]]->GetOutgoingProjections()[m]->PostLayer()->MPIGetProcessesUsed();

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

					vector<float> eData = events[k]->GetEventData();

					bool onlyLocal = true;

					if( postProcesses->size()>0)
						onlyLocal = false;
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
	int totalNrEvents = 0;
	TimingStart("Allreduce");

	MPI_Allreduce(&nrEvents,&totalNrEvents, 1,MPI_INT,MPI_SUM,NETWORK_COMM_WORLD); // could switch to multiple MPI communicators as network becomes large (and is modular enough to benefit).

	TimingStop("Allreduce");

	// TODO: may get removed
	if(this->MPIGetNodeId() == 0)
		cout<<"{"<<localEventIds.size()<<"/"<<totalNrEvents<<"} ";
	
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

	// TODO: (not optimized - could remove lowest/highest mixed, switch search method)
	TimingStart("HashBuffersPrepare");
	
	vector<long> tempIds = allIds;
	vector<int> useIds;

	long offset=0;
	if(m_allPreIds.size()>0 && tempIds.size()>0)
	{
		vector<long>::iterator lb = m_allPreIds.begin();
		if(m_allPreIds[0]<tempIds[0])
			lb = lower_bound(m_allPreIds.begin(),m_allPreIds.end(),tempIds[0]);

		vector<long>::iterator ub = m_allPreIds.end();
		if(m_allPreIds[m_allPreIds.size()-1]>tempIds[tempIds.size()-1])
			ub = upper_bound(m_allPreIds.begin(),m_allPreIds.end(),tempIds[tempIds.size()-1]);

		int ubInt,lbInt;
		ubInt = int(ub-m_allPreIds.begin());
		lbInt = int(lb-m_allPreIds.begin());
		//useIds.reserve(m_allPreIds.size());//1.5*m_nrLocalMessagesLastTimestep); // optimization

		//TimingStart("HashBuffersPrepare2");

		for(int i=tempIds.size()-1;i>-1;i--)
		{
			if(tempIds[i] <= m_allPreIds[ubInt-1] && tempIds[i]>= m_allPreIds[lbInt])// often enough
			{
				vector<long>::iterator it;
				it = lower_bound(lb,ub,tempIds[i]);
				
				if(tempIds[i]==*it)
				{
					//	if(m_allPreIds[int(lb-m_allPreIds.begin())] == tempIds[i])	// change check
					useIds.push_back(i);
				}
			}
		}
		
		//TimingStop("HashBuffersPrepare2");

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

	// Add incoming events to receiving populations
	#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	for(int i=0;i<m_populations.size();i++)
	{
		if(m_populations[i]->GetUnitsAll()->size()>0 && m_populations[i]->IsOn()) // only continue if this process is part of population
		{
			vector<Projection*> conns = m_populations[i]->GetIncomingProjections();
			for(int j=0;j<conns.size();j++)
			{
				if(conns[j]->KeepActiveBuffer() == false)
				{
					conns[j]->ClearActiveBuffer();

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

/// <summary>	Communication routine which uses MPI_Alltoall as main mpi routine
/// // - reduces total nr messages compared to allgather version. Will not transmit messages from a population to another which it is not connected to.
/// </summary>

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

			for(int m=0;m<m_populations[m_populationIndexesThisProc[i]]->GetOutgoingProjections().size();m++)
			{
				vector<int> prcs = m_populations[m_populationIndexesThisProc[i]]->GetOutgoingProjections()[m]->PostLayer()->MPIGetProcessesUsed();

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

					vector<float> eData = events[k]->GetEventData();
					localEventIds.push_back(eventId);
					localEventData.push_back(eData[0]);
					if(m_isTrackingHypercolumnIds)
						localEventHypercolumnIds.push_back(eventHypercolumnId);

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
						}
						nrEvents++;
					}

					currentIndex++;
				}

				localUnits[j]->IsNewEvent(false); // will be set true again by by unit's simulateeventque fcn or CreateEvent fcn
			}
		}
	}

	TimingStop("CommunicationGetEvents");

	// TODO: change to not send the local events, only external

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
		// send nr this process wants to send to everyone
		vector<int> nrMessagesSend(m_mpiNrProcs);
		vector<int> nrMessagesRecv(m_mpiNrProcs);
		vector<int> strideSend(m_mpiNrProcs);
		vector<int> strideRecv(nrMessagesRecv.size());
		
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
		
		// all-to-all ids
		// all-to-all data
		if(sendIds.size()==0)
		{	
			if(totMessToRecv == 0)
			{
				MPI_Alltoallv(NULL, &nrMessagesSend[0], &strideSend[0], MPI_LONG, NULL, &nrMessagesRecv[0], &strideRecv[0], MPI_LONG, NETWORK_COMM_WORLD);
				MPI_Alltoallv(NULL, &nrMessagesSend[0], &strideSend[0], MPI_FLOAT, NULL, &nrMessagesRecv[0], &strideRecv[0], MPI_FLOAT, NETWORK_COMM_WORLD);
			}
			else
			{
				MPI_Alltoallv(NULL, &nrMessagesSend[0], &strideSend[0], MPI_LONG, &allIds[0], &nrMessagesRecv[0], &strideRecv[0], MPI_LONG, NETWORK_COMM_WORLD);
				MPI_Alltoallv(NULL, &nrMessagesSend[0], &strideSend[0], MPI_FLOAT, &allData[0], &nrMessagesRecv[0], &strideRecv[0], MPI_FLOAT, NETWORK_COMM_WORLD);
			}
		}
		else
		{
			if(totMessToRecv == 0)
			{
				MPI_Alltoallv(&sendIds[0], &nrMessagesSend[0], &strideSend[0], MPI_LONG, NULL, &nrMessagesRecv[0], &strideRecv[0], MPI_LONG, NETWORK_COMM_WORLD);
				MPI_Alltoallv(&sendData[0], &nrMessagesSend[0], &strideSend[0], MPI_FLOAT, NULL, &nrMessagesRecv[0], &strideRecv[0], MPI_FLOAT, NETWORK_COMM_WORLD);
			}
			else
			{
				MPI_Alltoallv(&sendIds[0], &nrMessagesSend[0], &strideSend[0], MPI_LONG, &allIds[0], &nrMessagesRecv[0], &strideRecv[0], MPI_LONG, NETWORK_COMM_WORLD);
				MPI_Alltoallv(&sendData[0], &nrMessagesSend[0], &strideSend[0], MPI_FLOAT, &allData[0], &nrMessagesRecv[0], &strideRecv[0], MPI_FLOAT, NETWORK_COMM_WORLD);
			}
		}

		// all-to-all hypercolumn ids if applicable (bcpnn uses)
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
			vector<Projection*> conns = m_populations[i]->GetIncomingProjections();

			for(int j=0;j<conns.size();j++)
			{
				if(conns[j]->KeepActiveBuffer() == false)
				{
					conns[j]->ClearActiveBuffer();

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

/// <summary>	Adds a recording meter to network. </summary>
///
/// <param name="meter">	The meter. </param>

void Network::AddMeter(Meter* meter)
{
	m_meters.push_back(meter);
	meter->network(this);
}

/// <summary>
/// Performs application-defined tasks associated with freeing, releasing, or resetting unmanaged
/// resources.
/// (Resetting have been moved to Reset function)
/// </summary>

void Network::Dispose()
{
	this->ClearEventsIncoming();

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

}

/// <summary>	Resets network. Run in-between multiple indepent simulation runs (as decided from parameter class). </summary>

void Network::Reset()
{
	this->ClearEventsIncoming();

	m_hashSynapses.clear();

	for(int i=0;i<m_populations.size();i++)
	{
		delete m_populations[i];
	}

	for(int i=0;i<m_meters.size();i++)
		delete m_meters[i];

	m_analysis.clear(); // TODO: Check if could get error in timing if deleting and change

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

/// <summary>	 Removes all the current events/messages in this time step which are about to arrive
/// - use to reset the communication in a network, e.g. if an input stimuli is changed and we want a clean responding network.</summary>

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

void ProjectionModifier::SetProjection(Projection* c)
{
	m_projection = c;
}

/// <summary>	Adds population of units to this network. </summary>
///
/// <param name="layer">	The population. </param>

void Network::AddPopulation(Population* population)
{
	population->network(this);

	population->SetLayerId(m_currentLayerId);
	stringstream ss;
	ss<<"Layer"<<m_currentLayerId;
	if(population->GetName().size() == 0)
		population->SetName(ss.str());

	m_populations.push_back(population);
	m_currentLayerId++;
}

/// <summary>	Maps a unit id to unit. </summary>
///
/// <param name="unit">	The unit. </param>

void Network::AddUnit(Unit* unit)
{
	m_hashIdUnit[unit->GetUnitId()] = unit;
}

/// <summary>	Gets the same random seed on all processes. </summary>
///
/// <returns>	The random seed. </returns>

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

/// <summary>	Save all output data to new files or add to existing files (depending on how meter was set up). 
/// 			Can be called during a simulation or after all time steps have been completed. </summary>

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

/// <summary>	Stores the details about the network setup to default file.</summary>

void Network::StoreNetworkDetails()
{
	ofstream myfile;
	myfile.open (m_filenameNetworkDetails);
	myfile << "Network Details - do not change manually\n";

	// computer info
	myfile<<"NrProcesses: "<< this->MPIGetNrProcs()<<"\n"; // imp for nr Projections files

	// Populations
	myfile << "Populations: "<< m_populations.size() << "\n";

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

		// Save a geometry file if units have added geometry to have positions
		// TODO: could replace to get from Projection
		if(m_populations[i]->GetUnits().size()>0 && m_populations[i]->GetIncomingProjections().size()>0)
			if(m_populations[i]->GetIncomingProjections()[0]->GetUnitModifier("geometry") != NULL)
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
					geomfile<<((GeometryUnit*)m_populations[i]->GetIncomingProjections()[0]->GetUnitModifier(7))->GetValuesAsString()<<"\n";
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

/// <summary>	Print details about network to screen and put specific details in return vector. </summary>
///
/// <returns>	[nr processes, nr local units on process 0, total nr units in network, total nr connections in network]. </returns>

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
		vector<Projection*> inconns = m_populations[i]->GetIncomingProjections();
		vector<Unit*> units = m_populations[i]->GetLocalUnits();
		nrUnits+=units.size();
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
		cout<<"Av Connections/node = "<<(float)totSum/(float)this->MPIGetNrProcs()<<"\n";
		cout<<"Av Connections/unit = "<<(float)totSum/(float)totUnits<<"\n";

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

/// <summary>	Stores data that has been put in the analysis and
/// 			do not have their own file writing to a default csv file. </summary>

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

/// <summary>	Stores timing details about simulation to default file
/// 			(contains information for all network objects where timing has been turned on).  </summary>

void Network::StoreTimings()
{
	if(this->MPIGetNodeId() == 0)
	{
		cout<<"\nPrinting network timings...";

		ofstream myfile;
		stringstream ss;
		myfile.open (m_filenameTimings);

		myfile<<"Name / Average time / Std dev / Time max / Mean max time / Std dev max time /  Min time / Min process / Max time / Max process\n";
		myfile<<"------------------------------------------------------------------------------------------------------------------------------\n\n";

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

			for(int j=0;j<m_populations[i]->GetIncomingProjections().size();j++)
			{
				for(int k=0;k<m_populations[i]->GetIncomingProjections()[j]->GetProjectionModifiers().size();k++)
				{
					if(m_populations[i]->GetIncomingProjections()[j]->GetProjectionModifiers()[k]->IsTiming())
					{
						string s = m_populations[i]->GetIncomingProjections()[j]->GetProjectionModifiers()[k]->Timing()->GetTimingString();
						myfile<<i+1<<".C"<<j+1<<" "<<s;
					}
				}
			}
		}

		myfile.close();

		cout<<"ok.\n";
	}
}

/// <summary>	Sets a weight between two units. 
/// 			TODO: Move to be dependent on synapse class/structure. </summary>
///
/// <param name="weight">	The weight. </param>
/// <param name="preId"> 	Pre id of a unit. </param>
/// <param name="postId">	Post id of a unit. </param>

void Network::SetWeight(float weight, long preId, long postId)
{
	m_hashSynapses[postId][preId].weight = weight;
}

/// <summary>	Sets a delay between two units.
/// 			TODO: Move to be dependent on synapse class/structure. </summary>
///
/// <param name="delay"> 	The delay. </param>
/// <param name="preId"> 	Id of pre-unit. </param>
/// <param name="postId">	Id of post-unit. </param>

void Network::SetDelay(float delay, long preId, long postId)
{
#if USE_DELAYS==1
	m_hashSynapses[postId][preId].delay = delay;
#endif
}

/// <summary>	Gets a weight between two units. </summary>
///
/// <remarks>	Post Lazarus, 9/5/2012. </remarks>
///
/// <param name="preId"> 	Id of pre-unit. </param>
/// <param name="postId">	Id of post-unit. </param>
///
/// <returns>	The weight. </returns>

float Network::GetWeight(long preId, long postId)
{
	return m_hashSynapses[postId][preId].weight;
}

/// <summary>	Gets a delay between two units. </summary>
///
/// <remarks>	Post Lazarus, 9/5/2012. </remarks>
///
/// <param name="preId"> 	Id of pre-unit. </param>
/// <param name="postId">	Id of post-unit. </param>
///
/// <returns>	The delay. </returns>

float Network::GetDelay(long preId, long postId)
{
#if USE_DELAYS==1
	return m_hashSynapses[postId][preId].delay;
#else
	return 0;
#endif
}

/// <summary>	Creates all union of all pre ids on process.
/// 			Run in cases we need to have easy access to all pre unit-ids
///				Needs to be called every rebuild of connections.
///				TODO: optimize.
///				TODO: add a function that can add the new Projection pre id (if needed). </summary>

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

	// could lead to high memory usage
	// TODO: check memory and change logic if needed
	long totalLocalSynapses = 0;
	for(it1 = m_hashSynapses.begin();it1!=m_hashSynapses.end();it1++)
	{
		for(it2= it1->second.begin();it2!=it1->second.end();it2++)
		{
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

/// <summary>	Creates a vector for each unit to which processes to send an output
///				TODO: Make non-manually activated.
///				(depending on communication method, sparse send methods esp, this may or may not be needed) </summary>

void Network::CreateAllPostProcs()
{
	for(int i=0;i<this->MPIGetNrProcs();i++)
	{
		if(i == this->MPIGetNodeId())
		{
			int nrPreIds = m_allPreIds.size();
			MPI_Bcast(&nrPreIds,1,MPI_INT,i,NETWORK_COMM_WORLD);
			if(nrPreIds>0)
				MPI_Bcast(&m_allPreIds[0],nrPreIds,MPI_LONG,i,NETWORK_COMM_WORLD);
			else
				MPI_Bcast(NULL,nrPreIds,MPI_LONG,i,NETWORK_COMM_WORLD);
		}
		else
		{
			int nrPreIds;
			MPI_Bcast(&nrPreIds,1,MPI_INT,i,NETWORK_COMM_WORLD);

			vector<long> preIds(nrPreIds);
			if(nrPreIds>0)
				MPI_Bcast(&preIds[0],nrPreIds,MPI_LONG,i,NETWORK_COMM_WORLD);
			else
				MPI_Bcast(NULL,nrPreIds,MPI_LONG,i,NETWORK_COMM_WORLD);

				for(int k=0;k<m_populationIndexesThisProc.size();k++)
				{
						vector<Unit*>* units = m_populations[m_populationIndexesThisProc[k]]->GetUnitsAll();
						long unitId;
						for(int m=0;m<units->size();m++)
						{
							unitId = units->at(m)->GetUnitId();
							if(binary_search(preIds.begin(),preIds.end(),unitId))
							{
								units->at(m)->AddPostProcess(i);
							}
						}
				}
		}
	}
}

/// <summary>	Copies post synapses, same pre
///				- useful to copy a part of the network to another part   
///				(learnt set of Projections to an empty set of Projections at another place) 
///				TODO: make copy values synapse specific.</summary>
///
/// <param name="from">		 	Source projection. </param>
/// <param name="to">		 	Target projection. </param>
/// <param name="copyValues">	true if values in connections should be copied as well. </param>

void Network::CopyProjections(Projection* from, Projection* to, bool copyValues)
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
		to->AddProjections(preIds,newPostIds[i],newPostLocalIds[i]);

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

	// Re-map/-hash
	// not necessary if getpreidsall never used (TODO: check it gets recreated if it has not been run and remove here)
	to->CreatePreIdsUnion();
}

/// <summary>	Keep communication buffer. 
/// 			Can be used to reduce communication by keeping a buffer for the next time step (only one time step).</summary>
///				TODO: Change so not dependent on linear assumption of unit ids.
///
/// <param name="startId">	Start id of unit to keep incoming information for. </param>
/// <param name="endId">  	End id of unit to keep incoming information for. </param>

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

/// <summary>	Main run routine. Call to start initialization and run a network model. </summary>

void Network::Run()
{
	// retrieve description of network
	this->NetworkSetupStructure();

	// retrieve extra parameters here
	// TODO: check if this is natural to have at all
	this->NetworkSetupParameters();

	// Run for first parameter (special as most runs may only be done for one set of parameters)
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

	RecordAll(); // Write to disk anything that has not already been written at the end of a simulation
	StoreAnalysis(); // Write analysis at end by default
}