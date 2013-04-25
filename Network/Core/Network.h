// file:	Core\Network.h
//
// summary:	Declares Network class
// - Derive from to setup an own simulation (but can also call functions directly)
// - Has compiler flags for various libraries

#pragma once
#ifndef NETWORK_H
#define NETWORK_H


#define DEBUG_LEVEL 3 // not used, remove

// Libraries available
#define MUSIC_AVAILABLE 0
#define SPARSE_HASH_AVAILABLE 0
#define GSL_AVAILABLE 0
#define BOOST_AVAILABLE 0
#define VISIT_AVAILABLE 0
#define TR1_COMPILER 0
#define COMMUNICATION_BUFFER_HASH 1 // may be turned into normal flag as does not need to be on compile level
#define USE_HASHED_ACTIVE_COMMUNICATION 1	// faster simulation of units in most cases but uses more memory (extra hash scaling with nr of synapses)
											// does not need to be on compiler level
#define USE_DELAYS 0 // not needed to be put at compiler level, will change
#define USE_UNORDERED_MAP 1 // faster access to e.g. GetPreValue

#define USE_COMMUNICATION_ALLTOALL 1 // otherwise defaults to a communication pattern using allgatherv (all processes get all output messages), not needed to be put at compiler level

// Help pre-sets for different computer architectures
#define CRAY 3
#define JUGENE 2
#define BGL 1
#define PC 0

#define EPS 1e-9
#include <map>
#include <mpi.h>
#include <algorithm>
#include <vector>
#include <iostream>

// Depending on compiler and version, unordered_map is included differently
#if USE_UNORDERED_MAP == 1
#if TR1_COMPILER == 1
#include <tr1/unordered_map>
#else
#include <unordered_map>
#endif
#endif

// TODO: currently not used, check functionality
#if SPARSE_HASH_AVAILABLE == 1
	#include <google/type_traits.h>
	#include <google/sparsetable>

	using GOOGLE_NAMESPACE::sparsetable;
	using GOOGLE_NAMESPACE::sparse_hash_map;
	using GOOGLE_NAMESPACE::sparse_hash_set;
	using GOOGLE_NAMESPACE::dense_hash_map;
	using GOOGLE_NAMESPACE::dense_hash_set;
	using GOOGLE_NAMESPACE::HashtableInterface_SparseHashMap;
	using GOOGLE_NAMESPACE::HashtableInterface_SparseHashSet;
	using GOOGLE_NAMESPACE::HashtableInterface_SparseHashtable;
	using GOOGLE_NAMESPACE::HashtableInterface_DenseHashMap;
	using GOOGLE_NAMESPACE::HashtableInterface_DenseHashSet;
	using GOOGLE_NAMESPACE::HashtableInterface_DenseHashtable;
#endif

#include "NetworkObject.h"
#include "NetworkPopulation.h"
#include "Meter.h"
#include "DataSources.h"
#include "Storage.h"
#include "DataSources.h"
#include "NetworkParameters.h"


class ProjectionModifier;
class Population;
class Unit;
class Meter;
class NetworkParameters;

using namespace std;

/// <summary>	Network. Main simulation class. Recommended to inherit and extend for each new network model.</summary>

class Network : public NetworkObject
{
	#pragma message("Including Network")

public:

	Network();
	~Network();


// TODO: Synapse type atm decided on compile level, change
#if USE_DELAYS==1
	struct SynapseStandard // also allow alternate synapses
	{
		float weight;
		float delay;
	};
#else
	struct SynapseStandard
	{
		float weight;
	};
#endif

	enum SynapseModel 
	{
		SynapsesStandard
	};

	enum CommunicationBufferType
	{
		CommunicationBufferVector, // not good memory scaling, but fast (ok to use for small networks)
		CommunicationBufferHash // not memory bound, generally slower depending on hash/table implementation
	};

	// Random number generator seed, can also used srand straight away
	void SetSeed(bool allSame = false, float x = 8);

	void SetMPIParameters(int processId, int nrProcs)
	{
		m_mpiNodeId = processId;
		m_mpiNrProcs = nrProcs;
	}

	// TODO: Check why some classes use this
	void SetTimeResolution(float timeStep)
	{
		m_simulationResolution = timeStep;
	}

	float GetTimeResolution()
	{
		return m_simulationResolution;
	}

	// Used in network construction to get a unique global id of unit
	//  TODO: check why not move
	long GetNextUnitId()
	{
		long l = m_currentUnitId;
		m_currentUnitId++;
		return l;
	}

	// Used in network construction to get a unique global id of hypercolumn (in case of a columnar population)
	//  TODO: check why not move
	int GetNextHypercolumnId()
	{
		long l = m_currentHypercolumnId;
		m_currentHypercolumnId++;
		return l;
	}

	/// <summary>	Get the process id (0...N-1 for N processes). </summary>
	///
	/// <returns>	This MPI process id. </returns>

	int MPIGetNodeId()
	{
		return m_mpiNodeId;
	}

	/// <summary>	Gets nr of processes network is running on (communicator NETWORK_MPI_WORLD). </summary>
	///
	/// <returns>	Nr MPI processes. </returns>

	int MPIGetNrProcs()
	{
		return m_mpiNrProcs;
	}

	long MPIGetNrUnitsParallelizationDefault()
	{
		return m_mpiNrUnitsParallelizationDefault;
	}

	long MPIGetCurrentUnitParallelizationDefault()
	{
		return m_mpiCurrentUnitParallelizationDefault;
	}

	void MPIAddNrUnitsParallelizationDefault(long nrUnits)
	{
		m_mpiNrUnitsParallelizationDefault += nrUnits;
	}

	void MPIAddCurrentUnitParallelizationDefault(long nrUnits)
	{
		m_mpiCurrentUnitParallelizationDefault += nrUnits;
	}


	float GetCurrentTimeStep() { return m_currentTimeStep; }
	void SetCurrentTimeStep(float timeStep) { m_currentTimeStep = timeStep; }

	void AddPopulation(Population* population);
	Population* GetLayer(int index)
	{
		return m_populations[index];
	}

	void Initialize();
	void Simulate();
	void Simulate(int nrTimesteps);
	void CommunicationVersionAllgather();
	void CommunicationVersionAlltoall();

	void SetExtraFilenameString(char* extraString);

	void KeepCommunicationBuffer(long startId, long endId); // way to reduce communication
	void KeepCommunicationBuffer(bool keep) // way to reduce communication
	{
		this->m_keepCommunicationBuffer = keep;
	}

	void Reset();
	void AddUnit(Unit* unit);

	Unit* GetUnitFromId(long unitId) // only local atm?
	{
		return m_hashIdUnit[unitId];
	}

	float GetSimulationResolution()
	{
		return m_simulationResolution;
	}

	void SetSimulationResolution(float res)
	{
		m_simulationResolution = res;
	}

	void AddMeter(Meter* meter);

	Meter* GetMeter(int index)
	{
		return m_meters[index];
	}

	vector<Meter*> GetMeters()
	{
		return m_meters;
	}

	void ClearEventsIncoming();
	void ClearEventsIncoming(vector<long> fromPreIds);

	int GetSeed();

	void RecordAll();
	void StoreNetworkDetails();
	void StoreTimings();
	void StoreAnalysis();

	void Dispose();

	vector<float> PrintNetworkDetails();
	
	void SetWeight(float weight, long preId, long postId);
	void SetDelay(float delay, long preId, long postId);
	float GetWeight(long preId, long postId);
	float GetDelay(long preId, long postId);
	
	void SetUsingDelays(bool usingDelays) // note: not all neural units can use the delays without extending the classes
	{
		m_isUsingDelays = usingDelays;
	}
	
	bool IsUsingDelays()
	{
		return m_isUsingDelays;
	}

	void SetTrackingHypercolumnIds(bool trackingHypercolumnIds)
	{
		m_isTrackingHypercolumnIds = trackingHypercolumnIds;
	}

	bool IsTrackingHypercolumnIds()
	{
		return m_isTrackingHypercolumnIds;
	}

	int GetNrPopulations()
	{
		return m_populations.size();
	}

	void AddAnalysis(Analysis* analysis)
	{
		analysis->network(this);
		analysis->SwitchOnOff(true);
		m_analysis.push_back(analysis);
	}

	void AddTiming(NetworkObject* networkObject)
	{
		if(networkObject!=NULL)
		{
			networkObject->SetTiming(true, this);
			m_analysis.push_back(networkObject->GetAnalysisTiming());
		}
	}

	void SetFilenamesAdditional(string additional);


	float GetPreValue(long preId) // assumes no delays
	{
#if COMMUNICATION_BUFFER_HASH == 1
		return m_incomingBufferDataHash[preId];
#else
		return m_incomingBufferData[preId];
#endif
	}

#if USE_UNORDERED_MAP == 1
	unordered_map<long, SynapseStandard>* GetPreSynapses(long postId)
	{
		return &m_hashSynapses[postId];
	}
#else
	map<long, SynapseStandard>* GetPreSynapses(long postId)
	{
		return &m_hashSynapses[postId];
	}
#endif

	/// <summary>	Checks if a unit has any incoming connections. </summary>
	///
	/// <param name="postId">	id of post unit. </param>
	///
	/// <returns>	true if it has. </returns>

	bool PreSynapsesExist(long postId)
	{
		if(m_hashSynapses.find(postId) == m_hashSynapses.end())
			return false;
		else return true;
	}

#if USE_UNORDERED_MAP == 1
	unordered_map<long, unordered_map<long, SynapseStandard> >* GetSynapses()
	{
		return &m_hashSynapses;
	}
#else
	map<long, map<long, SynapseStandard> >* GetSynapses()
	{
		return &m_hashSynapses;
	}
#endif

	/// <summary>	Removes synapse values (if e.g. a connection is removed). </summary>
	///
	/// <param name="postId">	id of post unit. </param>
	/// <param name="preId"> 	id of pre unit. </param>

	void EraseSynapseValues(long postId, long preId)
	{
		m_hashSynapses[postId].erase(preId);
		if(m_hashSynapses[postId].size() == 0)
			m_hashSynapses.erase(postId);
	}

	long GetIncomingBufferHypercolumnIds(long preId)
	{
#if COMMUNICATION_BUFFER_HASH == 1
		return m_incomingBufferHypercolumnIdsHash[preId];
#else
		return m_incomingBufferHypercolumnIds[preId];
#endif
	}

	void CreateAllPreIdsUnion();
	void CreateAllPostProcs();

	vector<Analysis*> GetAnalysis()
	{
		return m_analysis;
	}

	void SetFilenameAnalysis(char* filename)
	{
		m_filenameAnalysis = filename;
	}

	void SetFilenameTimings(char* filename)
	{
		m_filenameTimings = filename;
	}

	void SetFilenameNetworkDetails(char* filename)
	{
		m_filenameNetworkDetails = filename;
	}

	// Copy Projections including weight values from one existing projection Projection to another (empty) projection
	void CopyProjections(Projection* from, Projection* to, bool copyValues);

	// Help function where objects the network should be responsible to delete can be put
	void AddNetworkObjectToDelete(void* obj) 
	{
		bool exists = false;
		for(int i=0;i<m_networkObjectsToDelete.size();i++)
		{
			if(m_networkObjectsToDelete[i] == obj)
			{
				exists = true;
				break;
			}
		}
		
		if(exists == false)
			m_networkObjectsToDelete.push_back(obj);
	}

	Storage::FilePreference GetFilePreference()
	{
		return m_savePreference;
	}

	void SetFilePreference(Storage::FilePreference savePreference)
	{
		m_savePreference = savePreference;
	}


	//////////// Interface for driving simulation from a network object

	void Run();

	//////////////////////////////////////////////////////////////////


	//////////// Interface for driving simulation from a network object
	//// See example class in NetworkTestsNetwork
	virtual void NetworkSetupStructure() { };
	virtual void NetworkSetupMeters() { };
	virtual void NetworkSetupParameters() { };
	virtual void NetworkRun() { };

	NetworkParameters* Parameters()
	{
		return m_networkParameters;
	}
	//////////////////////////////////////////////////////////////

	 

	/// <summary>	Sets run identifier. Used in multiple parameters runs to a different identifier for each independent run. </summary>
	///
	/// <param name="index">	Run id. </param>

	void SetRunId(int index)
	{
		m_runId = index;
	}

	/// <summary>	Gets run identifier. </summary>
	///
	/// <returns>	Run id. </returns>

	int GetRunId()
	{
		return m_runId;
	}

private:

	///////////////////////////////////
	// (currently) global settings

	CommunicationBufferType m_communicationBufferType; // what type of communication buffer incoming data should be kept in
	SynapseModel m_synapseModel; // TODO: extend

	///////////////////////////////////////////

	bool m_firstRun;					// keeps track of first simulation time step
	bool m_networkDetailsStored;
	bool m_useTiming;
	bool m_isUsingDelays;				// keeps track of if the network synapse model is using delays
	bool m_isTrackingHypercolumnIds;	// if pre hypercolumn ids are needed to be communicated etc (used in bcpnn)
	
	// extra filename strings to keep track of independent runs
	bool m_useExtraFilenameString;
	char* m_extraFilenameString;

	Storage::FilePreference m_savePreference;
	char* m_filenameNetworkDetails;
	char* m_filenameTimings;
	char* m_filenameAnalysis;

	float m_simulationResolution;
	
	int m_mpiNodeId; // mpi process id
	int m_mpiNrProcs; // nr of total processes available to network
	
	// used when building network to get unique ids
	long m_mpiCurrentUnitParallelizationDefault;
	long m_mpiNrUnitsParallelizationDefault;
	long m_currentUnitId;
	int m_currentHypercolumnId;
	int m_currentLayerId;

	// current time step used when writing data to disk
	float m_currentTimeStep;

	vector<Population*> m_populations;
	vector<int> m_populationIndexesThisProc;
	vector<vector<int> > m_communicationToProcs;
	vector<Connectivity*> m_connectivityTypes;
	vector<Meter*> m_meters;
	vector<Analysis*> m_analysis;

	// Note: for the large vectors, may want to sort the vector and then use the binary_search, lower_bound, or upper_bound algorithms

	map<long, Unit*> m_hashIdUnit;
	map<long,UnitModifier*> m_eventsUnitIncoming;
#if USE_UNORDERED_MAP == 1
	unordered_map<long, unordered_map<long, SynapseStandard> > m_hashSynapses;
#else
	map<long, map<long, SynapseStandard> > m_hashSynapses;
#endif

	vector<vector<float> > m_incomingBufferDataDelays;
	vector<float> m_incomingBufferData; // 

#if USE_UNORDERED_MAP == 1
	unordered_map<long,float> m_incomingBufferDataHash; // alternate version (slower, but not as memory bound - but could be postitive cache effects) - currently map and not spec hash impl
	unordered_map<long,long> m_incomingBufferHypercolumnIdsHash;
#else
	map<long,float> m_incomingBufferDataHash; // alternate version (slower, but not as memory bound - but could be postitive cache effects) - currently map and not spec hash impl
	map<long,long> m_incomingBufferHypercolumnIdsHash; 
#endif

	map<long, vector<float> > m_incomingBufferDataDelaysHash;
	vector<long> m_incomingBufferHypercolumnIds; // can store this locally to reduce communication and processing

	vector<long> m_allPreIds;
	vector<void*> m_networkObjectsToDelete; // network object the main network class is responsible to delete

	// turn on to keep communication buffer for next timestep to reduce communication
	bool m_keepCommunicationBuffer;

	// buffers keeping communicated (received + local) data
#if USE_UNORDERED_MAP == 1
	unordered_map<long,float> m_bufferToKeepData;
	unordered_map<long,long> m_bufferToKeepHypercolumnIds;
#else
	map<long,float> m_bufferToKeepData;
	map<long,long> m_bufferToKeepHypercolumnIds;
#endif
		
	// Handling of multiple parameters
	NetworkParameters* m_networkParameters;

	// Index for the current run if multiple parameters used (also used to set different random seeds for independent runs)
	int m_runId;

	// MUSIC specific
#if MUSIC_AVAILABLE == 1
	   public:
		   static MPI_Comm MPIComm;
#endif
};

	// MUSIC specific - need to use specific communicator if MUSIC is used as it takes up its own process/processes
#if MUSIC_AVAILABLE == 1
       //TODO: This could  be set to the setup->communicator() when music is used.
       #define NETWORK_COMM_WORLD Network::MPIComm
#else
       #define NETWORK_COMM_WORLD MPI_COMM_WORLD
#endif

#endif