#pragma once
#ifndef NETWORKCONNS_H
#define NETWORKCONNS_H

#include <string>
#include "Network.h"
#include "NetworkPopulation.h"
#include "NetworkConnectionModifier.h"
#include "NetworkUnits.h"

class Network;
class Population;
class ConnectionModifier;
class Unit;

class Connection;
class ConnectionFixed;
class UnitModifier;

using namespace std;

// Note: if symmetric connections (ie for random connectivities), they needs to be constructed specifically for each connectivity type (if not by default, ie full connectivity), not implemented in any generic way

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
		InitializeWeightsAndConnectionModifiers();
	}

	virtual void InitializeWeightsAndConnectionModifiers();

	virtual void MakeConnectionsSymmetric(vector<vector<long> >* preCache,vector<long>* postCache, ConnectionFixed* conn);

	virtual bool ExtraConstraints(Unit* preUnit, Unit* postUnit)// allows constraint on unit property level
	{
		return true;
	}

	virtual bool IterationPreConstraint(Unit* postUnit,int currentIndex, int currentHcIndex)		// faster to use for random connections etc that are only dependent on current iteration number or less
	{
		return true;
	}

	virtual int IterationAdvance()								// fast way to put constraint on connections (but cannot use any properties)
	{
		return 0;
	}

	virtual vector<long> IterationSpecificIndexes(Unit* postUnit, vector<int> nrPreRateUnits)				// fast way to build connections between post and only specifically generated pre indexes (but cannot use any properties)
	{
		return vector<long>(); // needs to be sorted lowest to highest!
	}
	
	virtual void GlobalSettings() { }; // can be used to initialize things such as global indexes before connections are built
	
	//virtual bool SetConnections();

	void SetPreAndPost(Population* pre, Population* post)
	{
		m_pre = pre;
		m_post = post;
	}

	void AddConnectionsEvent(ConnectionModifier* e)
	{
		m_eventConnections.push_back(e);
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

	virtual void ExtraPostUnit(Unit* postUnit,Connection* conn) { }

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

	vector<ConnectionModifier*> m_eventConnections;
	vector<UnitModifier*> m_unitsProperties;
	bool m_symmetric;
};

// can be used to initialize connections that will be copied from another connectivity later on
class EmptyConnectivity : public virtual Connectivity
{
public:
	EmptyConnectivity()
	{
		m_unitType = "";
		m_name = "EmptyConnectivity";
		m_useIterationSpecificIndexes = true; // will return empty vector of connections to build
	}

	void ExtraPostUnit(Unit* postUnit,Connection* conn);

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
		m_allowSelfConnections = true;
		m_setRandomWeights = true;
		m_weightsMinVal = 0;
		m_weightsMaxVal = 1;
		m_setRandomDelays = false;	
		m_useExtraConstraint = false;
	}

	FullConnectivity(bool allowSelfConnections, string unitType)
	{
		m_name = "FullConnectivity";
		m_allowSelfConnections = allowSelfConnections;
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

	// currently need to supply network as that may be set first when connection made (could supply in constructor instead)
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
	bool m_allowSelfConnections;

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
		m_allowSelfConnections = true;
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

class Connection : public NetworkObject
{
public:

	Connection()
	{
		m_keepActiveBuffer = false;
	}

	~Connection();
	void SetConnection(Unit* pre, Unit* post);
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
	
	vector<ConnectionModifier*> GetConnectionModifiers()
	{
		return m_events;
	}

	ConnectionModifier* GetEvent(string hashName)
	{
		return m_hashEvents[hashName];
	}

	void AddEvent(ConnectionModifier* e, string hashName);

	virtual void SimulateEvent(UnitModifier* e) = 0;
	virtual void ModifyConnection() = 0;
	//virtual void AddConnection(Unit* pre, Unit* post, bool firstRun) = 0;
	virtual void GenerateHashTables();

	void AddConnection(long preId, long postId, long postLocalId, bool firstRun);
	void AddConnections(vector<long> preIds, long postId, long postLocalId); // faster version - will speed up further by removing vectors of post ids

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
	vector<long>* GetPreIdsAll();//Union(); // in case of changing pre ids over time, this should be allowed to vary, otherwise keep static
	void CreatePreIdsUnion();

	void SetConnectivity(Connectivity* connectivityType)
	{
		m_connectivityType = connectivityType;
	}

	Connectivity* GetConnectivity()
	{
		return m_connectivityType;
	}

	//vector<long>* GetPostIds(long preId);

	vector<long> GetPreIds(long unitId);
	//vector<pair<long,float> >* GetPreIdsActive(long unitId); // gets only the active (from buffer) units
	vector<pair<long,float> >* GetPreIdsActiveLocal(long unitId);

	/*vector<vector<long>* > GetPreIds()
	{
		return m_listPosts;
	}*/

	/*void AddWeight(float weight, int index)
	{
		if(m_weights.size()<=index)
			m_weights.push_back(weight);
		else
			m_weights[index] = weight;
	}*/

	// will be removed
/*	void SetWeight(float weight, long preId, long postId)
	{
		((Network*)network())->SetWeight(weight, preId, postId);
		//m_hashWeights[postId][preId] = weight;
		//((this->network()->GetHashSynapses()))[preId][postId].weight = weight;//m_hashSynapses[preId][postId].weight = weight;
	}
	
	// will be removed
	void SetDelay(float delay, long preId, long postId)
	{
		network()->SetDelay(delay, preId, postId);

		//(*(this->network()->GetHashSynapses()))[preId][postId].delay = delay;//m_hashSynapses[preId][postId].delay = delay;
	}

	float GetDelay(long preId, long postId)
	{
		return network()->GetDelay(preId,postId);
		//return (*this->network()->GetHashSynapses())[preId][postId].delay;
		//return m_hashSynapses[preId][postId].delay;
		//return m_listSynapses[preId][postId].delay;
	}

	float GetWeight(long preId, long postId)
	{
		return network()->GetWeight(preId,postId);
		//return (*this->network()->GetHashSynapses())[preId][postId].weight;
		//return m_hashSynapses[preId][postId].weight;
		//return m_hashWeights[postId][preId];
		//return m_listSynapses[preId][postId].weight;
	}
	*/

	/*void SetWeight(float weight, long preId, long postId);
	void SetDelay(float delay, long preId, long postId);
	float GetWeight(long preId, long postId);
	float GetDelay(long preId, long postId);
	*/

			/*map<long,map<long,float> >::iterator itr;

		if ( (itr = m_hashWeights.find(postId)) != m_hashWeights.end())
		{
			map<long,float>::iterator itr2;

			if ( (itr2 = m_hashWeights[postId].find(preId)) != m_hashWeights[postId].end())
				return true;
			else
				return false;
		}
		else
			return false;
		*/

	// could send from which population
	/*bool ConnectionExists(long preId, long postId)
	{


		if(m_listPosts.size() == 0)
			return false;

		if(preId>m_listPosts.size()-1)
			return false;

		vector<long>* posts = m_listPosts[preId];
		
		for(unsigned int i=0;i<posts->size();i++)
		{
			if((*posts)[i] == postId)
				return true;
		}

		return false;
	}*/

	void Clear();
	void Clear(long postId);


	// will be removed anytime (!)
	/*vector<long>* PostIds() 
	{ 
		return &m_postIds;
	}*/

	/*vector<vector<long> >* PreIds() 
	{ 
		return &m_preIds;
	}*/

	std::vector<std::vector<float> > GetValuesToRecord();
	
	void CopyConnectionsOtherPost(Population* newPost);

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
			m_unitPropertiesConnection.push_back(p);
	}

	vector<UnitModifier*>* GetUnitPropertiesConnection()
	{
		return &m_unitPropertiesConnection;
	}

	UnitModifier* GetUnitModifier(int id);
	UnitModifier* GetUnitModifier(string name);

	void ClearEventsIncoming();

	void AddActiveEvent(long preId, float value); // used in hashed active communication mode
	void AddActiveEvents(vector<long> preIds, vector<float> values);
	vector<long> GetPostIds(long preId);
	vector<long> GetLocalPostIds(long preId);

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
	long m_totalNrLocalConnections;
	
	//ConnectionType* m_connection;
	vector<ConnectionModifier*> m_events;
	map<string,ConnectionModifier*> m_hashEvents;
	map<long,long,string> m_hashIdEvents;

	vector<UnitModifier*> m_unitPropertiesConnection; // can also be on a unit or layer level
	
//	map<long,map<long,float> > m_hashWeights; // <post id, <pre id, weight> >, map could be replace by faster/more memory efficient hash class or list
	
	
	//vector<vector<SynapseStandard> > m_listSynapses; // id-based list of synapses [preId, list of synapses]
	//vector<vector<long>* > m_listPosts; // id-based list of post ids [preId, list of post ids]

	vector<long> m_postIds;
	vector<long> m_preIdsUnion; // collection of all pre ids for all post ids, should be allowed to change over time e.g. if structural plasticity is taking place

	//vector<vector<long> > m_preIds; // index-based (wrt position in m_postIds) - used anymore?

	
	vector<vector<pair<long,float> > > m_preIdsActive; // accessed by local unit id

#if USE_UNORDERED_MAP == 1
	unordered_map<long, vector<long> > m_preIds;
//	unordered_map<long, vector<pair<long, float> > > m_preIdsActive;
#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	unordered_map<long, vector<long> > m_postIdsPre;
	unordered_map<long, vector<long> > m_localPostIdsPre;
#endif

#else
	map<long, vector<long> > m_preIds; // id-based (wrt position in m_postIds) - used anymore?
//	map<long, vector<pair<long,float> > > m_preIdsActive;

#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	map<long, vector<long> > m_postIdsPre;
	map<long, vector<long> > m_localPostIdsPre;
#endif

#endif

	int m_storeSamplingWeightsEvolution; // 0 if not, otherwise store by this rate (1 highest resolution)
	int m_storeSamplingCounter;
	// could place in parent
	vector<vector<float> > m_recordedValues;
};

class ConnectionFixed : public Connection
{
public:

/*	void SetWeight(float weight, long connectionId)
	{
		m_hashWeights[connectionId] = weight;
	}*/

	/*void SetWeight(float weight)
	{
		m_weight = weight;
	}

	float GetWeight()
	{
		return m_weight;
	}*/

	void SimulateEvent(UnitModifier* e);
	void ModifyConnection();

	//void AddConnection(Unit* pre, Unit* post, bool firstRun);
	
protected:

	//float m_weight;
	
};



#endif