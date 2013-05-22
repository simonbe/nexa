#pragma once
#ifndef NETWORKCONNS_H
#define NETWORKCONNS_H

#include <string>
#include <set>
#include "Network.h"
#include "NetworkPopulation.h"
#include "NetworkProjectionModifier.h"
#include "NetworkUnits.h"

class Network;
class Population;
class ProjectionModifier;
class Unit;

class Projection;
class ProjectionFixed;
class UnitModifier;

using namespace std;

// Note: if symmetric Projections (ie for random connectivities), they needs to be constructed specifically for each connectivity type (if not by default, ie full connectivity), not implemented in any generic way

class Connectivity : public NetworkObject
{
public:

	Connectivity()
	{
		m_symmetric = false; // can be true by default as defined by how connectivity is built
		m_name = "Connectivity";
		m_useIterationSpecificIndexes = false;
		m_useExtraConstraint = false;
	}

	virtual ~Connectivity() { }
	
	virtual void Initialize()
	{		
		InitializeWeightsAndProjectionModifiers();
	}

	virtual void InitializeWeightsAndProjectionModifiers();

	virtual void MakeProjectionsSymmetric(vector<vector<long> >* preCache,vector<long>* postCache, ProjectionFixed* conn);

	virtual bool ExtraConstraints(Unit* preUnit, Unit* postUnit)// allows constraint on unit property level
	{
		return true;
	}

	virtual bool IterationPreConstraint(Unit* postUnit,int currentIndex, int currentHcIndex)		// faster to use for random Projections etc that are only dependent on current iteration number or less
	{
		return true;
	}

	virtual int IterationAdvance()								// fast way to put constraint on Projections (but cannot use any properties)
	{
		return 0;
	}

	virtual vector<long> IterationSpecificIndexes(Unit* postUnit, vector<int> nrPreRateUnits)				// fast way to build Projections between post and only specifically generated pre indexes (but cannot use any properties)
	{
		return vector<long>(); // needs to be sorted lowest to highest!
	}
	
	virtual void GlobalSettings() { }; // can be used to initialize things such as global indexes before Projections are built
	
	//virtual bool SetProjections();

	void SetPreAndPost(Population* pre, Population* post)
	{
		m_pre = pre;
		m_post = post;
	}

	void AddProjectionsEvent(ProjectionModifier* e)
	{
		m_eventProjections.push_back(e);
	}

	void AddUnitsProperty(UnitModifier* p)
	{
		m_unitsProperties.push_back(p);
	}

	vector<UnitModifier*> GetUnitsProperties()
	{
		return m_unitsProperties;
	}

	void Dispose()
	{
	}

	virtual void SetWeightValues(long preId, long postId) {};
	//virtual void SetWeightValues(Unit* pre, Unit* post) = 0;

	virtual void ExtraPostUnit(Unit* postUnit,Projection* conn) { }

	// Multiple parameters implementation
	virtual void AddParameters() { }

	void SetUnitType(string unitType)
	{
		m_unitType = unitType;
	}

protected:

	//Network* m_network; // parent network
	Population* m_pre;
	Population* m_post;
	bool m_useIterationSpecificIndexes;
	bool m_useExtraConstraint;

	string m_unitType;

	vector<ProjectionModifier*> m_eventProjections;
	vector<UnitModifier*> m_unitsProperties;
	bool m_symmetric;
};

// can be used to initialize Projections that will be copied from another connectivity later on
class EmptyConnectivity : public virtual Connectivity
{
public:
	EmptyConnectivity()
	{
		m_unitType = "";
		m_name = "EmptyConnectivity";
		m_useIterationSpecificIndexes = true; // will return empty vector of Projections to build
	}

	void ExtraPostUnit(Unit* postUnit,Projection* conn);

private:
};

// TODO: Full connectivity within hypercolumns only

class FullConnectivity : public virtual Connectivity
{
public:
	
	FullConnectivity()
	{
		m_unitType = "";
		m_name = "FullConnectivity";
		m_allowSelfProjections = true;
		m_setRandomWeights = true;
		m_weightsMinVal = 0;
		m_weightsMaxVal = 1;
		m_setRandomDelays = false;	
		m_useExtraConstraint = false;
	}

	FullConnectivity(bool allowSelfProjections, string unitType)
	{
		m_name = "FullConnectivity";
		m_allowSelfProjections = allowSelfProjections;
		m_unitType = unitType;//m_unitType = "";
		m_setRandomWeights = false;//true;
		m_weightsMinVal = 0;
		m_weightsMaxVal = 1;
		m_setRandomDelays = false;
		m_useExtraConstraint = false;
	}

	void SetRandomWeights(float minVal, float maxVal)
	{
		m_setRandomWeights = true;
		m_weightsMinVal = minVal;
		m_weightsMaxVal = maxVal;
	}

	void SetRandomWeightsMin(float minVal)
	{
		m_setRandomWeights = true;
		m_weightsMinVal = minVal;
	}

	void SetRandomWeightsMax(float maxVal)
	{
		m_setRandomWeights = true;
		m_weightsMaxVal = maxVal;
	}

	// currently need to supply network as that may be set first when Projection made (could supply in constructor instead)
	void SetRandomWeights(vector<float> minVal, float maxVal,Network* net);
	void SetRandomWeights(float minVal, vector<float> maxVal,Network* net);

	void SetRandomDelays(float minVal, float maxVal)
	{
		m_setRandomDelays = true;
		m_delaysMinVal = minVal;
		m_delaysMaxVal = maxVal;
	}

/*	FullConnectivity(string unitType)
	{
		m_unitType = unitType;
	}
*/
	/*FullConnectivity(vector<string> unitTypes)
	{
		m_unitTypes = unitTypes;
	}*/

	//virtual void Initialize();
	virtual void SetWeightValues(long preId, long postId);
	//virtual void SetWeightValues(Unit* pre, Unit* post);
	virtual bool ExtraConstraints(Unit* preUnit, Unit* postUnit);
	//virtual void PatchyConnect(vector<long> pres,vector<long> posts); 
	//virtual void PatchyConnect(vector<long> pres,vector<long> posts,vector<float> weights); // this method should be used to set up connectivity patterns

	// Multiple parameters implementation
	void AddParameters();


protected:

	bool m_setRandomWeights;
	bool m_setRandomDelays;
	bool m_allowSelfProjections;

	float m_weightsMinVal, m_weightsMaxVal;
	float m_delaysMinVal, m_delaysMaxVal;

	vector<float> m_weightsMinValParams;
	vector<float> m_weightsMaxValParams;
};

class RandomConnectivity : public FullConnectivity
{
public:

	RandomConnectivity(float fractionConnected, bool symmetric = true)
	{
		m_fracConnected = fractionConnected;
		m_name = "RandomConnectivity";
		
		m_setRandomWeights = false;//true;
		m_symmetric = symmetric; // building symmetric will be alot slower atm - forcing to go through all units
		m_weightsMinVal = 0;
		m_weightsMaxVal = 1;
		m_setRandomDelays = false;
		
		m_useExtraConstraint = false;
		m_useIterationSpecificIndexes = true;//false;
		//m_useIterationSpecificIndexes = true;
		/*if(symmetric==true)
			m_useExtraConstraint = true;
		else
		{
			m_useExtraConstraint = false;
			//m_useIterationSpecificIndexes = true; // not forced but checked now
		}*/
	}

	// override extra constraints
	bool ExtraConstraints(Unit* preUnit, Unit* postUnit);
	//bool IterationPreConstraint(int currentIndex);
	//int IterationAdvance();

	vector<long> IterationSpecificIndexes(Unit* postUnit, vector<int> nrPreRateUnits);

private:
	float m_fracConnected;
};

class FullConnectivityNoLocalHypercolumns : public FullConnectivity
{
public:

	FullConnectivityNoLocalHypercolumns()
	{
		m_unitType = "";
		m_name = "FullConnectivityNoLocalHypercolumns";
		m_allowSelfProjections = true;
		m_setRandomWeights = true;
		m_weightsMinVal = 0;
		m_weightsMaxVal = 1;
		m_setRandomDelays = false;	
		m_useExtraConstraint = true;
	}

	// override extra constraints
	bool ExtraConstraints(Unit* preUnit, Unit* postUnit);

private:
};

class OneToOneConnectivity : public FullConnectivity
{
public:
	
	OneToOneConnectivity(bool betweenHypercolumns) // if false, between minicolumns
	{
		m_name = "OneToOneConnectivity";
		
		m_setRandomWeights = false;
		m_weightsMinVal = 1;
		m_weightsMaxVal = 1;
		m_setRandomDelays = false;
		m_betweenHypercolumns = betweenHypercolumns;
		m_useExtraConstraint = true;
		//m_useIterationSpecificIndexes = true;
	}


	// override extra constraints
	bool ExtraConstraints(Unit* preUnit, Unit* postUnit);
	bool IterationPreConstraint(Unit* postUnit,int currentIndex, int currentHcIndex);
	vector<long> IterationSpecificIndexes(Unit* postUnit, vector<int> nrPreRateUnits);	

private:

	bool m_betweenHypercolumns;
};

class FanInConnectivity : public FullConnectivity
{
public:
	
	FanInConnectivity(int nrIns)
	{
		m_name = "FanInConnectivity";
		
		m_setRandomWeights = false;
		m_weightsMinVal = 0;
		m_weightsMaxVal = 1;
		m_setRandomDelays = false;
		m_nrIns = nrIns;
		m_useExtraConstraint = true;
	}

	// override extra constraints
	bool ExtraConstraints(Unit* preUnit, Unit* postUnit);

private:

	int m_nrIns;
};

class Projection : public NetworkObject
{
public:

	Projection()
	{
		m_keepActiveBuffer = false;
	}

	~Projection();
	void SetProjection(Unit* pre, Unit* post);
	void Initialize(Population* pre, Population* post)
	{
		m_preLayer = pre;
		m_postLayer = post;
		m_storeSamplingCounter = 0;

		m_allowPreIdsChanges = false;
		m_keepActiveBuffer = false;
	}

	void SetRandomWeights(float minVal, float maxVal); // will work on post side
	void SetRandomWeightsIn(long localPostIndex, float minVal, float maxVal);
	void SetRandomWeightsOut(long localPreIndex, float minVal, float maxVal);
	
	vector<ProjectionModifier*> GetProjectionModifiers()
	{
		return m_events;
	}

	ProjectionModifier* GetEvent(string hashName)
	{
		return m_hashEvents[hashName];
	}

	void AddEvent(ProjectionModifier* e, string hashName);

	virtual void SimulateEvent(UnitModifier* e) = 0;
	virtual void ModifyProjection() = 0;
	virtual void GenerateHashTables();

	void AddProjection(long preId, long postId, long postLocalId, bool firstRun);
	void AddProjections(vector<long> preIds, long postId, long postLocalId); // faster version - will speed up further by removing vectors of post ids

	Population* PreLayer()
	{
		return m_preLayer;
	}

	Population* PostLayer()
	{
		return m_postLayer;
	}

	vector<float> GetPreValues(long unitId);
	vector<float> GetPreValuesAll(); // values for union of pre ids
	vector<long> GetPreLocalIds(long unitId);
	vector<float> GetPostValues();
	vector<long> GetPostLocalIds();
	vector<long> GetPostIds();
	vector<long>* GetPreIdsAll(); // in case of changing pre ids over time, this should be allowed to vary, otherwise keep static
	void CreatePreIdsUnion();

	void SetConnectivity(Connectivity* connectivityType)
	{
		m_connectivityType = connectivityType;
	}

	Connectivity* GetConnectivity()
	{
		return m_connectivityType;
	}

	vector<long> GetPreIds(long unitId);
	vector<pair<long,float> >* GetPreIdsActiveLocal(long unitId);
		
	void Clear();
	void Clear(long postId);

	vector<vector<float> > GetValuesToRecord();

	// Copy all weight values to another connectivity
	void CopyProjectionsOtherPost(Population* newPost);

	// Set all weights to the same value
	void SetWeightValues(float value);

	void SetRecording(bool on, int samplingRate) // sampling rate over 0 then weights evolution stored
	{
		m_recording = on;
		m_storeSamplingWeightsEvolution = samplingRate;
	}

	void AddPostUnit(long unitId)
	{
		m_postIds.push_back(unitId);
	}

	void AddUnitsProperty(UnitModifier* p)
	{
		if(p!=NULL) // can set to NULL to not use any
			m_unitPropertiesProjection.push_back(p);
	}

	vector<UnitModifier*>* GetUnitPropertiesProjection()
	{
		return &m_unitPropertiesProjection;
	}

	UnitModifier* GetUnitModifier(int id);
	UnitModifier* GetUnitModifier(string name);

	void ClearEventsIncoming();

	void AddActiveEvent(long preId, float value); // used in hashed active communication mode
	void AddActiveEvents(vector<long> preIds, vector<float> values);
	set<long>& GetPostIds(long preId);
	set<long>& GetLocalPostIds(long preId);

	void ClearActiveBuffer();

	void KeepActiveBuffer(bool keep)
	{
		m_keepActiveBuffer = keep;
	}

	bool KeepActiveBuffer()
	{
		return m_keepActiveBuffer;
	}

	void EraseSynapses(vector<pair<long,long> > vectorPostIdPreId);
	void EraseSynapse(long postId, long preId);


protected:

	Population* m_preLayer;
	Population* m_postLayer;
	Connectivity* m_connectivityType;

	bool m_allowPreIdsChanges;
	bool m_keepActiveBuffer;
	long m_totalNrLocalProjections;
	
	//ProjectionType* m_projection;
	vector<ProjectionModifier*> m_events;
	map<string,ProjectionModifier*> m_hashEvents;
	map<long,long,string> m_hashIdEvents;

	vector<UnitModifier*> m_unitPropertiesProjection; // can also be on a unit or layer level
	
	vector<long> m_postIds;
	vector<long> m_preIdsUnion; // collection of all pre ids for all post ids, should be allowed to change over time e.g. if structural plasticity is taking place

	vector<vector<pair<long,float> > > m_preIdsActive; // accessed by local unit id

#if USE_UNORDERED_MAP == 1
	unordered_map<long, vector<long> > m_preIds;
#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	unordered_map<long, set<long> > m_postIdsPre;
	unordered_map<long, set<long> > m_localPostIdsPre;
#endif

#else
	map<long, vector<long> > m_preIds; // id-based (wrt position in m_postIds) - used anymore?
#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	map<long, set<long> > m_postIdsPre;
	map<long, set<long> > m_localPostIdsPre;
#endif

#endif

	int m_storeSamplingWeightsEvolution; // 0 if not, otherwise store by this rate (1 highest resolution)
	int m_storeSamplingCounter;

	// could place in parent
	vector<vector<float> > m_recordedValues;
};

class ProjectionFixed : public Projection
{
public:

	void SimulateEvent(UnitModifier* e);
	void ModifyProjection();
	
protected:

};



#endif