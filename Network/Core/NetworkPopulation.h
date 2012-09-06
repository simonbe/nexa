#pragma once
#ifndef Population_H
#define Population_H

#include "Network.h"
#include "NetworkPopulationModifier.h"
#include "NetworkProjections.h"
#include "NetworkUnits.h"
#include "MPIDistribution.h"

class Network;
class Connectivity;
class Unit;
class Meter;
class MPIDistribution;

using namespace std;

/// <summary>	Population of neural units. </summary>

class Population : public NetworkObject
{
public:

	Population();

	virtual ~Population();

	void SetLayerId(int layerId)
	{
		m_layerId = layerId;
	}

	int GetLayerId()
	{
		return m_layerId;
	}

	enum PopulationType
	{
		FixedWeightsAndRateUnits,
		Data
	};

	// Initializations
	virtual void Initialize() = 0;
	virtual void InitializeProjectionsEventsAndParameters() = 0;
	
	// Run before global communication
	virtual void Simulate() = 0;
	
	// Run after global communication
	virtual void Modify();

	virtual void Reset() = 0;

	virtual void ResetLocalities();
	virtual void SetValuesLocal(vector<float> values) = 0;
	virtual void SetValuesAll(vector<float> values, bool alsoSimulate = false) = 0;

	void AddPre(Population* layer, Connectivity* ProjectionType, bool calledFromAddPost = false);
	void AddPost(Population* layer, Connectivity* ProjectionType, bool calledFromAddPre = false);

	virtual void AddUnit(Unit* unit)
	{
		m_units.push_back(unit);
	}

	virtual vector<Unit*> GetUnits()
	{
		return m_units;
	}

	vector<Unit*>* GetUnitsAll()
	{
		return &m_units;
	}

	void AddUnitsModifierToInitialize(UnitModifier* p);

	virtual vector<Unit*> GetUnits(string type)
	{
		return m_units;
	}

	vector<Unit*> GetLocalUnits();
	
	virtual long GetNrUnitsTotal()
	{
		return m_units.size();
	}

	void SetPopulationType(PopulationType type)
	{
		m_populationType = type;
	}

	PopulationType GetPopulationType()
	{
		return m_populationType;
	}

	int MPIBelongsToNode(unsigned long localUnitId) // default method, will be put in MPI class instead
	{
		for(int i=0;i<m_unitIndexesInterval.size();i++)
		{
			vector<long> localUnitIndexesInterval = m_unitIndexesInterval[i];

			if(localUnitId >= localUnitIndexesInterval[0] &&
			localUnitId < localUnitIndexesInterval[1])
			return i; // should return from hash m_hashIndexNode[i] - could have a number of nodes specialized in this layer then
		}

		return -1;
	}

	vector<vector<float> > GetWeights();
	vector<vector<float> > GetLocalWeights();

	void AddPopulationModifier(PopulationModifier* eventLayer);

	void AddPostNode(int node)
	{
		vector<int>::iterator result;
		result = find( m_postNodes.begin(), m_postNodes.end(), node );

		if( result == m_postNodes.end() )
			m_postNodes.push_back(node);
	}

	bool IsPostNode(int node)
	{
		vector<int>::iterator result;
		result = find( m_postNodes.begin(), m_postNodes.end(), node );

		if( result == m_postNodes.end() )
			return false;
		else return true;
	}

	void AddPreNode(int node)
	{
		vector<int>::iterator result;
		result = find( m_preNodes.begin(), m_preNodes.end(), node );

		if( result == m_preNodes.end() )
			m_preNodes.push_back(node);
	}

	bool IsPreNode(int node)
	{
		vector<int>::iterator result;
		result = find( m_preNodes.begin(), m_preNodes.end(), node );

		if( result == m_preNodes.end() )
			return false;
		else return true;
	}

	vector<int> GetPostNodes()
	{
		return m_postNodes;
	}

	vector<int> GetPreNodes()
	{
		return m_preNodes;
	}

	MPIDistribution* MPI()
	{
		return m_mpiDistribution;
	}

	vector<int> GetNodeLayerIndexes()
	{
		return m_nodeLayerIndexes;
	}

	vector<int> MPIGetProcessesUsed()
	{
		return m_mpiProcessesUsed;
	}

	/*Unit* GetUnitFromId(long id)
	{
		return m_hashIdUnit[id];
	}*/

	void AddOutgoingProjection(Projection* c)
	{
		m_projectionOutgoing.push_back(c);
	}

	void AddIncomingProjection(Projection* c)
	{
		m_projectionIncoming.push_back(c);
	}

	vector<Projection*> GetOutgoingProjections()
	{
		return m_projectionOutgoing;
	}

	vector<Projection*> GetIncomingProjections()
	{
		return m_projectionIncoming;
	}

	vector<PopulationModifier*> GetPopulationModifiers()
	{
		return m_populationModifiers;
	}

	bool IsFirstRun()
	{
		return m_firstRun;
	}

	void SetFirstRun(bool firstRun)
	{
		m_firstRun = firstRun;
	}

	vector<vector<float> > GetValuesToRecord()
	{
	//	if(IsRecording())
//		{
			return m_recordedValues;
//		}
//		else return vector<vector<float> >(0);
	}

	void ClearMemory()
	{
		for(int i=0;i<m_recordedValues.size();i++)
			m_recordedValues[i].clear();

		m_recordedValues.clear();
	}

	virtual void ClearEventsIncoming();

	void Dispose();
	void DisposeParent();

	long GetUnitsStartId()
	{
		return m_unitsStartId;
	}

	virtual vector<int> GetStructure()
	{
		return vector<int>(1,m_units.size());
	}

	// Values in a buffer - used for MPI-comm
	void SetValuesBuffer(vector<float> values)
	{
		m_valuesBuffer = values;
	}

	virtual vector<float> GetValuesLocal();

	// 
	vector<float> GetValuesBuffer()
	{
		if(m_valuesBuffer.size()>0 )
			return m_valuesBuffer;
		else
		{
			return GetValuesLocal();
		}
	}

	void SetTrace(bool useTrace, float traceStrength)
	{
		m_useTrace = useTrace;
		m_traceStrength = traceStrength;
	}

	void SetAttachedMeter(Meter* meter)
	{
		m_meterLayer = meter;
	}

	vector<Connectivity*> GetPreConnectivitys()
	{
		return m_preConnectivitys;
	}

	void AddUnitsProperty(UnitModifier* p)
	{
		m_unitPropertiesLayer.push_back(p);
	}

	vector<UnitModifier*>* GetUnitPropertiesLayer()
	{
		return &m_unitPropertiesLayer;
	}

	bool KeepValues()
	{
		return m_keepValues;
	}

	void SetKeepValues(bool keepValues)
	{
		m_keepValues = keepValues;
	}

	// temporary in use
	void SetUnitIdLocalId(long unitId, long localUnitId)
	{
		m_unitIdsLocalIds[unitId] = localUnitId;
	}

	// temporary in use
	long GetLocalIdUnitId(long unitId)
	{
		return m_unitIdsLocalIds[unitId];
	}
	
protected:

	int m_layerId;
	//Network* m_network; // parent network
	PopulationType m_populationType;
	long m_nrUnits;
	//NetworkType* m_networkType;
	vector<Unit*> m_units;
	vector<Population*> m_pres;
	vector<Population*> m_posts;
	vector<int> m_postNodes; // mpi nodes
	vector<int> m_preNodes; // mpi nodes
	vector<Connectivity*> m_preConnectivitys;
	vector<Connectivity*> m_postConnectivitys;

	vector<UnitModifier*> m_unitPropertiesToInitialize;
	vector<UnitModifier*> m_unitPropertiesLayer;
	
	vector<Projection*> m_projectionOutgoing;
	vector<Projection*> m_projectionIncoming;

//	map<long,Unit*> m_hashIdUnit; // currently full size on all nodes
//	map<long, vector<Unit*> > m_hashIdPreUnit;

	vector<ProjectionModifier*> m_preProjectionModifiers; // ?

	vector<long> m_localUnitIndexesInterval;
	vector<vector<long> > m_unitIndexesInterval;

	vector<PopulationModifier*> m_populationModifiers;
	
	MPIDistribution* m_mpiDistribution;
	vector<int> m_mpiProcessesUsed;

	vector<int> m_nodeLayerIndexes;

	//m_on in NetworkObject // keeps the parameters (values) of units fixed (ie during training phase)
	bool m_firstRun;
	vector<vector<float> > m_recordedValues;
	long m_unitsStartId;

	// Buffer of all values, both global and local - set i.e. in MPIDistribution::MPIMakeLayerValuesLocal
	vector<float> m_valuesBuffer;

	// trace - could be put in Populationcolumns atm.
	bool m_useTrace;
	float m_traceStrength;

	// Meter - here to be able to check which type of file-writing it is supposed to be (as that determines logic)
	Meter* m_meterLayer;

	// Parallelization strategy (normally not specified + change to more general setting)
	MPIDistribution::ParallelizationSchemeLayer m_parallelizationScheme;

	// Updating
	bool m_keepValues; // if true, units should be updated in next timestep even if they are manually set

#if USE_UNORDERED_MAP == 1
	unordered_map<long,long> m_unitIdsLocalIds;
#else
	map<long,long> m_unitIdsLocalIds;
#endif
};

class PopulationColumns : public Population
{
public:

	enum UnitType
	{
		Graded,
		GradedThresholded,
		IF,
		adEIF
	};

	PopulationColumns(Network* net, unsigned long nrHypercolumns, unsigned long nrRateUnits, UnitType unitType, MPIDistribution::ParallelizationSchemeLayer parallelizationScheme = MPIDistribution::ParallelizationDefault, bool useSilentHypercolumns = false, float silentHypercolumnsThreshold = 0);

	void Initialize();
	void InitializeProjectionsEventsAndParameters();
	void Simulate();

	vector<int> GetNrRateUnits()
	{
		return m_nrRateUnits;
	}

	long GetNrUnitsTotal() // will give total nr minicolumns in Populationcolumns
	{
		return m_nrUnits;
	}

	vector<int> GetLocalHypercolumnIndexes()
	{
		return m_localHypercolumnIndexes;
	}

	vector<Unit*> GetLocalRateUnits()
	{
		return m_localRateUnits;
	}

	vector<int> GetUnitIdLocals();

	// overridden, will only return minicolumn values
	vector<float> GetValuesLocal();

	// will only return the minicolumns
	vector<Unit*> GetUnits()
	{
	//	return m_units;
		return m_minicolumns;
	}

	// Manually setting the values
	void SetValuesLocal(vector<float> values);
	void SetValue(int localIndex,float value);
	void SetValuesAll(vector<float> values, bool alsoSimulate = false);
//	void SetValuesAll(vector<short> values);

	void Reset();

	vector<Unit*> GetUnits(string type)
	{
		if(type.compare("minicolumn")==0)
			return (vector<Unit*>)m_minicolumns;
		else if(type.compare("hypercolumn") == 0)
			return (vector<Unit*>)m_hypercolumns;
		else return m_units;
	}

	void AddUnit(Unit* unit);

	vector<RateUnit*> GetRateUnits()	{ vector<RateUnit*> mcs; for(int i=0;i<m_minicolumns.size();i++) mcs.push_back((RateUnit*)m_minicolumns[i]); return mcs;	}
	vector<Hypercolumn*> GetHypercolumns()	{ vector<Hypercolumn*> hcs; for(int i=0;i<m_hypercolumns.size();i++) hcs.push_back((Hypercolumn*)m_hypercolumns[i]); return hcs;	}


	Unit* GetHypercolumn(int index)
	{
		return (m_hypercolumns[index]);
	}

	/*vector<int> GetRateUnitsIndexes(int hcIndex)
	{
		return m_hashMcHcIndexes[hcIndex];//m_mcHcIndexes[hcIndex];
	}*/

	// assumes a linear distribution of the minicolumn ids - currently not preprocessed
	vector<long> GetRateUnitsIds(int hcIndex)
	{
		vector<long> mcIds(m_nrRateUnits[hcIndex]);
		int cId = 0;
		for(int i=0;i<hcIndex;i++)
			cId+=m_nrRateUnits[i];

		for(int i=0;i<m_nrRateUnits[hcIndex];i++)
			mcIds[i] = cId+i+m_startRateUnitId;
		
		return mcIds;
		//return m_mcHcIds[hcIndex];
	}

	vector<long> GetNodeIndexes(int hcIndex)
	{
		return m_nodeHcIndexes[hcIndex];
	}

	vector<vector<int> > GetAllNodeIndexes()
	{
		// change to default int
		vector<vector<int> > nhi(m_nodeHcIndexes.size());
		for(int i=0;i<nhi.size();i++)
		{
			vector<int> v(m_nodeHcIndexes[i].size());
			nhi[i] = v;
			for(int j=0;j<v.size();j++)
				nhi[i][j] = m_nodeHcIndexes[i][j];
		}

		return nhi;
	}

	vector<int> GetMPIDistribution();//int processId);
	vector<int> GetMPIDistributionHypercolumns(int processId);

	void ResetLocalities();
	void Dispose();
	~PopulationColumns()
	{
		Dispose();
	}
	//void SendAndReceiveVersionISend();
	//void SendAndReceiveVersionAllgather();

	void ClearEventsIncoming();

	
	vector<int> GetStructure()
	{
		return m_nrRateUnits;
	}

	void SetSilentHypercolumns(bool useSilent, float threshold);

	int GetNrHypercolumns()
	{
		return m_nrHypercolumns;
	}

	int GetNrRateUnitsInHypercolumn()
	{
		return m_nrRateUnitsInHypercolumn;
	}
	
protected:

	UnitType m_unitType;

	long m_startRateUnitId;

	vector<Unit*> m_minicolumns;
	vector<Unit*> m_hypercolumns;
	int m_nrHypercolumns;
	int m_nrRateUnitsInHypercolumn;
	vector<int> m_nrRateUnits;
	vector<int> m_localHypercolumnIndexes;
//	vector<long> m_localHypercolumnIds;

	vector<Unit*> m_localRateUnits;
	vector<vector<int> > m_mcHcIndexes;
	map<int, vector<int> > m_hashMcHcIndexes;
	vector<vector<long> > m_mcHcIds;

	vector<vector<long> > m_nodeHcIndexes; // splitted hypercolumns over several processes -> change to int
	vector<vector<int> > m_nodeHcIds; // splitted hypercolumns (ids)
	
	vector<UnitModifier*> m_allEventsIncoming;

	// Silent hypercolumns
	bool m_useSilentHypercolumns;
	float m_silentHypercolumnsThreshold;
};

#endif