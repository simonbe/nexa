#include <iostream>
#include <set>
#include "NetworkProjections.h"
#include "NetworkUnits.h"
#include "NetworkProjectionModifier.h"
#include "NetworkUnitModifier.h"

// Currently need to do the following if adding a new plasticity class:
// 1. Add .h-file here as they are called in the Projection initializations
// 2. Initialize them down where the others are initialized (this will be moved to inside these classes to get rid of this)

#include "NetworkBCPNN.h"
#include "NetworkCL.h"
#include "NetworkMI.h"
#include "NetworkMDS.h"
#include "NetworkVQ.h"
#include "NetworkKussul.h"
#include "NetworkFoldiak.h"
#include "NetworkTriesch.h"
#include "NetworkCorr.h"
#include "NetworkHebbSimple.h"
#include "NetworkBCM.h"
#include "NetworkSanger.h"
#include "NetworkSTDP.h"


void Projection::GenerateHashTables()
{
	// any extra hash tables can be initialized here
}

// Used to add a new Projection between units based on their ids. firstRun by default false.
void Projection::AddProjection(long preId, long postId, long postLocalId, bool firstRun)
{
	bool first = true;
	long postIndex;

	vector<long>::iterator index = std::find(m_postIds.begin(),m_postIds.end(),postId);
	if(index !=m_postIds.end())
	{
		postIndex = index - m_postIds.begin();
		first = false;
	}

	if(first == true)
	{
		postIndex = m_postIds.size();
		m_postIds.push_back(postId);
	}

	if(m_preIds.size()<postIndex+1)
	{
		for(int i=m_preIds.size();i<postIndex+1;i++)
		{
			vector<long> v;
			//m_preIds.push_back(v);
			m_preIds[postId] = v;
		}
	}

#if USE_HASHED_ACTIVE_COMMUNICATION == 1
		m_postIdsPre[preId].insert(postId);
		if(postLocalId<0)
			postLocalId = this->PostLayer()->GetLocalIdUnitId(postId);
		m_localPostIdsPre[preId].insert(postLocalId);
#endif

	m_preIds[postIndex].push_back(preId);
}


void Projection::AddProjections(vector<long> preIds, long postId, long postLocalId)
{
	bool first = true;
	int postIndex;

	postIndex = m_postIds.size();

	if(find(m_postIds.begin(),m_postIds.end(),postId) == m_postIds.end())
		m_postIds.push_back(postId);

	if(m_preIds.find(postId) == m_preIds.end())
		m_preIds[postId] = preIds;
	else
	{
		vector<long> v = m_preIds[postId];
		for(int i=0;i<preIds.size();i++)
			v.push_back(preIds[i]);
		
		// remove duplicates
		sort(v.begin(),v.end());
		v.erase(unique(v.begin(),v.end()),v.end());
		m_preIds[postId] = v;
	}

#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	for(int i=0;i<preIds.size();i++)
	{
		// always guarantee no duplicates?
		// cannot be necessary to handle both
		m_postIdsPre[preIds[i]].insert(postId);
		if(postLocalId<0)
			postLocalId = this->PostLayer()->GetLocalIdUnitId(postId);
		m_localPostIdsPre[preIds[i]].insert(postLocalId);
	}
#endif

	// for mpi_Alltoall-communication pattern
#if USE_COMMUNICATION_ALLTOALL == 1

#endif

}

void Projection::Clear()
{
	//for(int i=0;i<m_preIds.size();i++)
#if USE_UNORDERED_MAP == 1
	for(unordered_map<long, vector<long> >::iterator iter = m_preIds.begin();iter!=m_preIds.end();++iter)
#else
	for(map<long, vector<long> >::iterator iter = m_preIds.begin();iter!=m_preIds.end();++iter)
#endif
	{
		for(int j=0;j<m_preIds[iter->first].size();j++)
			m_network->EraseSynapseValues(m_preIds[iter->first][j],iter->first);

		m_preIds[iter->first].clear();
	}

#if USE_HASHED_ACTIVE_COMMUNICATION == 1
#if USE_UNORDERED_MAP == 1
	for(unordered_map<long, vector<long> >::iterator iter = m_preIds.begin();iter!=m_preIds.end();++iter)
#else
	for(map<long, vector<long> >::iterator iter = m_postIdsPre.begin();iter!=m_postIdsPre.end();++iter)
#endif
	{
		m_postIdsPre[iter->first].clear();
	}

	m_postIdsPre.clear();

#endif

	m_preIds.clear();
	m_postIds.clear();
	
// alternative way
/*	for(map<long,map<long,float> >::iterator iter =m_hashWeights.begin();iter!=m_hashWeights.end();++iter)
	{
		long postId = (*iter).first;
		Clear(postId);	
	}
	*/
}

void Projection::Clear(long postId)
{
	map<long,float> empty;

	int postIndex = -1;
	for(int i=0;i<m_postIds.size();i++)
	{
		if(postId == m_postIds[i])
		{
			postIndex = i;
			break;
		}
	}
	
	m_preIds[postIndex].clear();
}

void Projection::ClearActiveBuffer()
{
#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	m_preIdsActive.clear();
	//m_preIdsActive = vector<vector<pair<long, float> > >(this->PostLayer()->GetNrUnitsTotal());//this->PostLayer()->GetUnits().size());
#endif
}

void Projection::AddEvent(ProjectionModifier* e, string hashName)
{
	e->SetProjection(this);
	m_events.push_back(e);
	m_hashEvents[hashName] = e;
}

void Projection::SetRandomWeights(float minVal, float maxVal)
{
	vector<long> postIds = this->GetPostIds();

	for(int i=0;i<postIds.size();i++)
	{
		vector<long> preIds = this->GetPreIds(postIds[i]);

		for(int j=0;j<preIds.size();j++)
		{
			float r = rand()/(float(RAND_MAX)+1);
			float weight = (maxVal-minVal)*r + minVal;

			network()->SetWeight(weight,preIds[j],postIds[i]);
		}
	}
}

void Projection::SetRandomWeightsIn(long localPostIndex, float minVal, float maxVal)
{
	vector<long> postIds = this->GetPostIds();
	vector<long> localPostIds = this->GetPostLocalIds();

	for(int i=0;i<localPostIds.size();i++)
	{
		if(localPostIds[i] == localPostIndex)
		{
			vector<long> preIds = this->GetPreIds(postIds[i]);

			for(int j=0;j<preIds.size();j++)
			{
				float r = rand()/(float(RAND_MAX)+1);
				float weight = (maxVal-minVal)*r + minVal;

				network()->SetWeight(weight,preIds[j],postIds[i]);
			}
		}
	}
}

void Projection::SetRandomWeightsOut(long localPreIndex, float minVal, float maxVal)
{
	vector<long> postIds = this->GetPostIds();

	for(int i=0;i<postIds.size();i++)
	{
		vector<long> preIds = this->GetPreIds(postIds[i]);
		vector<long> preLocalIds = this->GetPreLocalIds(postIds[i]);

		for(int j=0;j<preLocalIds.size();j++)
		{
			if(preLocalIds[j] == localPreIndex)
			{
				float r = rand()/(float(RAND_MAX)+1);
				float weight = (maxVal-minVal)*r + minVal;

				network()->SetWeight(weight,preIds[j],postIds[i]);
			}
		}
	}
}

void ProjectionFixed::ModifyProjection()
{
	for(int i=0;i<m_events.size();i++)
	{
		if(m_events[i]->IsOn())
			m_events[i]->Modify();

	}

	// record new weights if storing on evolution of weights

	if(IsRecording() && m_storeSamplingWeightsEvolution>0)
	{
		if(m_storeSamplingCounter % m_storeSamplingWeightsEvolution == 0)
		{
			m_storeSamplingCounter=0;

			// append
			for(int i=0;i<this->GetPostIds().size();i++)
			{
				long postId = this->GetPostIds()[i];
				vector<long> preIds = this->GetPreIds(postId);
				vector<float> ws(preIds.size());

				for(int j=0;j<preIds.size();j++)
				{
					ws[j] = m_network->GetWeight(preIds[j],postId);
				}

				m_recordedValues.push_back(ws);
			}
		}

		m_storeSamplingCounter++;
	}
}

void ProjectionFixed::SimulateEvent(UnitModifier* e)
{
	for(int i=0;i<m_events.size();i++)
	{
		if(m_events[i]->GetEventId() == 1)
		{
			ProjectionModifierBcpnnOnline* b = (ProjectionModifierBcpnnOnline*)m_events[i];
			b->Simulate(e);
		}
	}
}

vector<float> Projection::GetPreValues(long unitId)
{
	vector<long> preIds = m_preIds[unitId]; // could change to pointer
	vector<float> out(preIds.size());
	
	for(int i=0;i<preIds.size();i++)
	{
		out[i] = m_network->GetPreValue(preIds[i]);
	}

	return out;
}

// re-implement
vector<long> Projection::GetPreLocalIds(long unitId)
{
	vector<long> preIds = m_preIds[unitId];
	vector<long> localPreIds(preIds.size());

	// assumes a linear distribution and startId is hypercolumn id, so one less
	long startId = m_preLayer->GetUnitsStartId() + 1;

	for(int i=0;i<preIds.size();i++)
	{
		long localIndex = preIds[i]-startId;
		//long localId = ((RateUnit*)m_preLayer->network()->GetUnitFromId(preIds[i]))->GetUnitIdLocal();
		localPreIds[i] = localIndex;
	}

	return localPreIds;
}

vector<long> Projection::GetPreIds(long unitId)
{
	return m_preIds[unitId]; // if an index-based implementation were to be used this could be changed to be retrieved by unitId-m_postIds[0];
}

void Projection::AddActiveEvent(long preId, float value)
{
	set<long> postIds;

	// get all the units this unit is connected to on this Projection
	postIds = GetPostIds(preId);
	
	// put in their active queues/buffers
	for(set<long>::iterator it = postIds.begin(); it!=postIds.end();++it)
	{
		pair<long, float> p;
		p.first = preId;
		p.second = value;
		// put in vector with position determined by hash (this is set up once)
		m_preIdsActive[*it].push_back(p);
	}
}

void Projection::AddActiveEvents(vector<long> preIds, vector<float> values)
{
	//vector<set<long> > localPostIds(preIds.size()); 
	if(m_preIdsActive.size()>0)
		m_preIdsActive.clear();
	else
		m_preIdsActive = vector<vector<pair<long,float> > >(this->PostLayer()->GetNrUnitsTotal());//this->PostLayer()->GetUnits().size());

	for(int i=0;i<preIds.size();i++)
	{
		// get all the units this unit is connected to on this Projection
		//localPostIds[i] = m_localPostIdsPre[preIds[i]];//GetLocalPostIds(preIds[i]);//GetPostIds(preIds[i]);
		for(set<long>::iterator it=m_localPostIdsPre[preIds[i]].begin();it!=m_localPostIdsPre[preIds[i]].end();++it)
		{
			m_preIdsActive[*it].reserve(preIds.size()); // optimization
		}
	}

	pair<long, float> p;
	for(int i=0;i<preIds.size();i++)
	{
		// put in their active queues/buffers

		//pair<long, float> p; // new here will take time (track in profiler to get exact measures)
		p.first = preIds[i];
		p.second = values[i];
		long localId;

		for(set<long>::iterator it=m_localPostIdsPre[preIds[i]].begin();it!=m_localPostIdsPre[preIds[i]].end();++it)//postIds.size();j++)
		{
			m_preIdsActive[*it].emplace_back(p);
		}
	}
}


vector<pair<long,float> >* Projection::GetPreIdsActiveLocal(long localUnitId)
{
	if(m_preIdsActive.size() == 0) // move (!)
		this->ClearActiveBuffer();
	// fix this -> should be index-based (?)
	//for(int i=0;i<m_postIds.size
	if(m_preIdsActive.size()==300)
		bool b = false;
	
	return &m_preIdsActive[localUnitId];//unitId-m_postIds[0]];
}

// Necessary to set USE_HASHED_ACTIVE_COMMUNICATION to 1 to be able to use this (consumes unnecessary memory otherwise (scales with nr synapses))
set<long>& Projection::GetPostIds(long preId)
{
#if	USE_HASHED_ACTIVE_COMMUNICATION == 1
	return m_postIdsPre[preId];
#else
	return vector<long>(0);
#endif
}

set<long>& Projection::GetLocalPostIds(long preId)
{
#if	USE_HASHED_ACTIVE_COMMUNICATION == 1
	return m_localPostIdsPre[preId];
#else
	return vector<long>(0);
#endif
}

void Projection::EraseSynapses(vector<pair<long,long> > vectorPostIdPreId)
{
	for(int i=0;i<vectorPostIdPreId.size();i++)
	{
		pair<long,long> p = vectorPostIdPreId[i];
		EraseSynapse(p.first,p.second);
	}

	// recreate hash (ids of pre-units) union
	CreatePreIdsUnion();
}

void Projection::EraseSynapse(long postId, long preId)
{
	vector<long> preIds = m_preIds[postId];

	// slow search, and not using same for all preIds of one postId (EraseSynapses specifies vector)
	int index = -1;
	
	for(int i=0;i<preIds.size();i++)
	{
		if(preIds[i] == preId)
		{
			index = i;
			break;
		}
	}

#if USE_HASHED_ACTIVE_COMMUNICATION

	int index2 = -1;
	set<long> postIds = m_postIdsPre[postId];
	set<long>::iterator it = postIds.find(postId);
	if(it!=postIds.end()) {
		postIds.erase(it);
		m_postIdsPre[preId] = postIds;
	}
		

#endif

	if(index != -1)
	{
		preIds.erase(preIds.begin()+index);
		m_preIds[postId] = preIds;

		// erase value from network hash
		this->network()->EraseSynapseValues(postId,preId);
	}


}

vector<long> Projection::GetPostLocalIds()
{
	//return m_postIds;//m_network->GetPostIds();
	
	vector<long> localIds;

	for(int i=0;i<m_postIds.size();i++) // all local?
	{
		RateUnit* m = (RateUnit*)m_postLayer->network()->GetUnitFromId(m_postIds[i]);
		if(m->IsLocal()==true)
			localIds.push_back(m->GetUnitIdLocal());
			//localIds.push_back(m_postIds[i]);
	}

	return localIds;
}

vector<long> Projection::GetPostIds()
{
	/*vector<long> ids(m_hashWeights.size());
	//int i=m_hashWeights.size()-1;
	int i=0;

	for(map<long,map<long,float> >::iterator iter = m_hashWeights.begin();iter!=m_hashWeights.end();++iter)
	{
		long postId = (*iter).first;

		ids[i] = postId;
		//i--;
		i++;
	}

	return ids;*/

	return m_postIds;
}

vector<long>* Projection::GetPreIdsAll()//GetPreIdsUnion()
{
	if(m_preIdsUnion.size() == 0)
	{
		CreatePreIdsUnion();

		return &m_preIdsUnion;
	}
	else return &m_preIdsUnion;

}

// does same as CreateAllPreIdsUnion but only these Projections, could merge for build performance (!)
// default: called in initialization phase
void Projection::CreatePreIdsUnion()
{
	m_preIdsUnion.clear();
	m_totalNrLocalProjections = 0; //

	map<long,int> tempMap;
	
	for(int i=0;i<m_postIds.size();i++)
	{
		if(m_network->PreSynapsesExist(m_postIds[i]))
		{

#if USE_UNORDERED_MAP == 1
			unordered_map<long,Network::SynapseStandard>* preSynapses = m_network->GetPreSynapses(m_postIds[i]);
			unordered_map<long, Network::SynapseStandard>::iterator it2;
#else
			map<long,Network::SynapseStandard>* preSynapses = m_network->GetPreSynapses(m_postIds[i]);
			map<long, Network::SynapseStandard>::iterator it2;
#endif

			for(it2 = preSynapses->begin();it2!=preSynapses->end();it2++)
			{
				// no duplicates (not using nr times)
				tempMap[it2->first] = tempMap[it2->first]+1;
				m_totalNrLocalProjections++;
			}
		}
	}

	map<long,int>::iterator it3;

	for(it3 = tempMap.begin();it3!=tempMap.end();it3++)
		m_preIdsUnion.push_back(it3->first);

	if(m_network->MPIGetNodeId() == 0)
	{
		cout<<"Total nr incoming local Projections (layer "<<m_name<<") == "<<m_totalNrLocalProjections<<"\n";
		cout.flush();
	}
}

/*vector<long>* Projection::GetPostIds(long preId)
{
	if(m_listPosts.size() == 0) // will be 0 if node has no local units in layer
		return NULL;
	else if(preId >= m_listPosts.size())
		return NULL;
	else
		return m_listPosts[preId]; // use pointer instead
}*/

vector<float> Projection::GetPostValues()
{
	/*vector<float> out(m_hashWeights.size());//network()->HashWeights()->size());

	int i=0;

	for(map<long,map<long,float> >::iterator iter = m_hashWeights.begin();iter != m_hashWeights.end();++iter)
	{
		long postId = (*iter).first;
		
		RateUnit* m = (RateUnit*)m_postLayer->network()->GetUnitFromId(postId);
		out[i] = m->GetValue();
		
		i++;
	}

	return out;*/

	vector<float> out(m_postIds.size());//network()->HashWeights()->size());

	int i=0;

	for(int i=0;i<m_postIds.size();i++)
	{
		RateUnit* m = (RateUnit*)m_postLayer->network()->GetUnitFromId(m_postIds[i]);
		out[i] = m->GetValue();
	}

	return out;
}

vector<float> Projection::GetPreValuesAll()
{
	vector<float> out(m_preIdsUnion.size());//network()->HashWeights()->size());

	int i=0;

	for(int i=0;i<m_preIdsUnion.size();i++)
	{
		out[i] = m_network->GetPreValue(m_preIdsUnion[i]);//GetIncomingBufferData(m_preIdsUnion[i]);
	}

	return out;	
}

/*float Projection::GetWeight(long preId, long postId)
{
	return (*network()->HashWeights())[postId][preId];
}

void Projection::SetWeight(float weight, long preId, long postId)
{
	(*network()->HashWeights())[postId][preId] = weight;
}*/

vector<vector<float> > Projection::GetValuesToRecord()
{
	vector<vector<float> > outData;

	if(IsRecording())
	{
		if(m_storeSamplingWeightsEvolution>0)
		{
			outData = m_recordedValues;
			m_recordedValues.clear();
		}
		else
		{
			float w;

			vector<vector<float> > data;
#if USE_UNORDERED_MAP == 1
			unordered_map<long, vector<long> >::iterator it;
#else
			map<long, vector<long> >::iterator it;
#endif
			// (?) SLOW! Could change pre/post -> post/pre in hash for direct access
			for(it = m_preIds.begin();it!=m_preIds.end();it++)
			//for(int i=0;i<m_preIds.size();i++) // m_listPosts may not exist for much longer
			{
				vector<float> d;
				for(int j=0;j<it->second.size();j++)//m_preIds[i].size();j++)
				{
					w = this->network()->GetWeight(it->second[j], it->first);//i,j);
					d.push_back(w);
				}

				if(d.size()>0)
				{
					data.push_back(d);
				}
			}

			return data;
		}

		/*for(map<long, map<long, SynapseStandard> >::iterator iter = m_hashSynapses.begin();iter!=m_hashSynapses.end();++iter)
		{
			long preId = (*iter).first;
			map<long, SynapseStandard> postUnitSynapses = m_hashSynapses[preId];

			vector<float> weightValues(postUnitSynapses.size());

			int i=0;
			for(map<long, SynapseStandard>::iterator iter = postUnitSynapses.begin();iter!=postUnitSynapses.end();++iter)
			{
				float weight = (*iter).second.weight;
				weightValues[i] = weight;
				i++;
			}

			outData.push_back(weightValues);
		}*/
	}

	return outData;
}

// has been moved to Network class
void Projection::CopyProjectionsOtherPost(Population* newPost)
{
	// Get preids and postids

	// Change list of post ids to new population

	// use connectivitytype::AddProjections

	// go through all Projections, copy weight values etc (optional)
}

void Projection::SetWeightValues(float value)
{
	vector<long> postIds = this->GetPostIds();

	for(int i=0;i<postIds.size();i++)
	{
		long postId = postIds[i];
		vector<long> preIds = this->GetPreIds(postId);//this->GetPreSynapses(postId);

		for(int j=0;j<preIds.size();j++)
		{
			this->network()->SetWeight(value,preIds[j],postId);
		}
	}
}

void FullConnectivity::SetWeightValues(long preId, long postId)//Unit* pre, Unit* post)
{
//	long preId = pre->GetUnitId();
//	long postId = post->GetUnitId();

	if(m_setRandomWeights == true)
	{
		float r = rand()/(float(RAND_MAX)+1);
		float weight = (m_weightsMaxVal-m_weightsMinVal)*r + m_weightsMinVal;

		if(weight>0.5)
			bool b = false;
		network()->SetWeight(weight,preId,postId);//preUnits[i]->GetUnitId(),postUnits[j]->GetUnitId());
	}
	else // should this be as default or selected ?
	{
		network()->SetWeight(0.0,preId,postId);
	}

	if(m_setRandomDelays == true)
	{
		float r = rand()/(float(RAND_MAX)+1);
		float delay = (m_delaysMaxVal-m_delaysMinVal)*r + m_delaysMinVal;

		network()->SetDelay(delay,preId,postId);//preUnits[i]->GetUnitId(),postUnits[j]->GetUnitId());
	}
}

void Connectivity::InitializeWeightsAndProjectionModifiers()
{
	TimingStart(m_name);

//	if(m_symmetric == true) 
//		this->m_network->SetSeed(true); // set same seed, may not stay like this

	// may change compute-wise
	if(m_initialized == false)
	{
		if(network()->MPIGetNodeId() == 0)
			cout<<m_name<<": Initialize\n";

		vector<Unit*> preUnits;
		vector<Unit*> postUnits;

		vector<vector<long> > preCache;
		vector<long> postCache;
		vector<long> postLocalCache;

		//TransferBcpnnOnline* transferBcpnn = new TransferBcpnnOnline(); // transfer function, may not be used (so should not be initialized here)

		// change from all real units to local and non-local units
		bool specificType = false;
		bool isOverHypercolumns = false;

		if(m_unitType.compare("") != 0)
		{
			preUnits = m_pre->GetUnits(m_unitType);
			postUnits = m_post->GetUnits(m_unitType);
			if(m_unitType.compare("hypercolumn") == 0)
				isOverHypercolumns = true;
			specificType = true;
		}
		else
		{
			preUnits = m_pre->GetUnits();
			postUnits = m_post->GetUnits();
		}

		ProjectionFixed* conn = new ProjectionFixed();
		conn->Initialize(m_pre,m_post);
		conn->network(this->network());
		conn->SetConnectivity(this);

		m_pre->AddOutgoingProjection(conn);
		m_post->AddIncomingProjection(conn);

		ProjectionModifierBcpnnOnline* eBcpnnOnline; // allocated upon request
		ProjectionModifierMIRateUnit* eMIRateUnits;
		ProjectionModifierMIHypercolumn* eMIHypercolumns;
		ProjectionModifierMDS* eMDS;
		ProjectionModifierVQ* eVQ;
		ProjectionModifierCL* eCL;
		ProjectionModifierKussul* eKussul;
		ProjectionModifierFoldiak* eFoldiak;
		ProjectionModifierCSL* eCSL;
		ProjectionModifierTriesch* eTriesch;
		ProjectionModifierPearson* ePearson;
		ProjectionModifierHebbSimple* eHebbSimple;
		ProjectionModifierBCM* eBCM;
		ProjectionModifierSanger* eSanger;
		ProjectionModifierSTDP* eSTDP;

		bool bcpnnInit = false, miMcInit = false, miHcInit = false, mdsInit = false, vqInit = false, clInit = false, ksInit = false, foInit = false, cslInit = false, trieschInit = false, pearsonInit = false, hebbSimpleInit = false, bcmInit = false, sangerInit = false, stdpInit = false;

		// pre och post alla local för hypercolumn! borde inte vara så, kolla resetlocalities etc.

		long ProjectionId = 0;
		bool firstRun = true;
		int cIndex = -1;

		vector<int> nrPreRateUnits = ((PopulationColumns*)m_pre)->GetStructure(); // assumes we use Populationcolumns, otherwise local id will not be correct in created pre unit (if property is to be used in ExtraConstraint)

		if(isOverHypercolumns == true)
		{
			for(int i=0;i<nrPreRateUnits.size();i++)
				nrPreRateUnits[i] = 1;
		}

		for(long j=0;j<postUnits.size();j++)
		{
			if(j%4==0 && network()->MPIGetNodeId() == 0)
			{
				cout<<".";
				cout.flush();
				if(j%100 == 0)
				{
					int totPre = 0;
					for(int kk=0;kk<preCache.size();kk++)
						totPre+=preCache[kk].size();

					cout<<100*(float)j/(float)postUnits.size()<<"% (post,pre,totPre;"<<postCache.size()<<","<<preCache.size()<<","<<totPre<<")";
				}
			}

			ExtraPostUnit(postUnits[j],conn);

			if(j>0)
				firstRun = false;

			bool addedPostCache = false;
			bool addedPreCache = false;

			//for(long i=0;i<preUnits.size();i++)
			long currentPreIndex = 0; // now fixed by startId //nrPreRateUnits.size(); // assumes first ids are reserved for the hypercolumns

			if(isOverHypercolumns == true) // build Projections over the hypercolumns
			{
				currentPreIndex = 0;
			}



			int steps = IterationAdvance();

			TimingStart("SpecificIndexes");
			vector<long> specificIndexes = IterationSpecificIndexes(postUnits[j],nrPreRateUnits); // needs to be ordered lowest to highest
			if(specificIndexes.size()>0)
				m_useIterationSpecificIndexes = true; // could also be forced
			TimingStop("SpecificIndexes");

			sort(specificIndexes.begin(),specificIndexes.end());

			long currentSpecificIndex = 0;
			long currentSpecificIteration = 0;

			//TimingStart("IterPre");

			int nrHcs = nrPreRateUnits.size();
			int nrMcs = 1;
			int nrItersForSpec = 1;
			if(specificIndexes.size()>0)
			{
				nrItersForSpec = specificIndexes.size();
				nrHcs = 1;
				nrMcs = 1;
			}

			for(int iterSpec = 0;iterSpec<nrItersForSpec;iterSpec++)
			{
				for(long m2=0;m2<nrHcs;m2++)
				{
					if(specificIndexes.size()==0)
						nrMcs = nrPreRateUnits[m2];

					for(int i2=0;i2<nrMcs;i2++)
					{
						bool continueBuild = true;

						if(specificIndexes.size()>0)
						{
							// always cont
						}
						else
						{

							if(m_useIterationSpecificIndexes == true)//specificIndexes.size()>0)
							{
								continueBuild = false;
								if(currentSpecificIndex!=specificIndexes.size())
								{
									if(specificIndexes[currentSpecificIndex] == currentSpecificIteration)
									{
										currentSpecificIndex++;
										continueBuild = true;
									}
								}

								currentSpecificIteration++;
							}
						}

						if(continueBuild == true)
						{
							long m;
							int i;

							if(specificIndexes.size()>0)
							{
								currentPreIndex = specificIndexes[iterSpec];
								//currentPreIndex += specificIndexes[iterSpec] -1 ;
								// assuming same nr mcs in each hc
								i = (specificIndexes[iterSpec]) % nrPreRateUnits[0];
								m = (specificIndexes[iterSpec]) / nrPreRateUnits[0];

							}
							else
							{
								m = m2;
								i = i2;
							}

							if(steps>0)
								steps--;
							else
							{
								//TimingStart("IterPreAdd");

								steps = IterationAdvance(); // fastest way to put constraint on Projections (but cannot use any properties)

								// create new preUnit so that also non-local units get checked (viable while using ExtraConstraint to build Projections)
								// could switch to iterated build to get rid of non-local checks, but then need to build iterators instead of using ExtraConstraint
								if(IterationPreConstraint(postUnits[j],currentPreIndex,m))
								{
									bool extraContinue = true;
									long startId = m_pre->GetUnitsStartId();
									long preUnitId = currentPreIndex+startId;

									if(m_useExtraConstraint == true)
									{
										RateUnit preUnit(false);

										preUnit.SetUnitId(preUnitId);
										preUnit.SetUnitIdLocal(i+nrPreRateUnits[0]*m); // assumes same number of mcs in each hc
										Hypercolumn h;
										h.SetUnitIdLocal(m);
										for(int mi=0;mi<nrPreRateUnits[m];mi++) // this scales bad - do not allow ExtraConstraint to have access to this info (!)
										{
											RateUnit m2(false);
											h.AddRateUnit(&m2);
										}
										preUnit.AddHypercolumn(&h);
										//								preUnit.SetHypercolumnId();

										if(specificType == true)
										{
											// implement !
										}

										extraContinue = ExtraConstraints(&preUnit,postUnits[j]);
									}


									/*if(preUnits[i]->IsLocal() == true) // add the mpi node of post unit to the list of send nodes (if not already there
									{
									preUnits[i]->AddPostNode(postUnits[j]->GetNodeId()); // notify unit
									m_pre->AddPostNode(postUnits[j]->GetNodeId()); // notify layer
									}*/

									//if(postUnits[j]->IsLocal() == true && ExtraConstraints(preUnits[i],postUnits[j]) == true)
									if(postUnits[j]->IsLocal() == true && extraContinue == true)
									{
										//m_post->AddPreNode(preUnits[i]->GetNodeId()); // notify layer

										Population::PopulationType typePost = m_post->GetPopulationType();
										Population::PopulationType typePre = m_pre->GetPopulationType();

										if(typePost == Population::FixedWeightsAndRateUnits && typePre == Population::FixedWeightsAndRateUnits) // default
										{
											//ProjectionFixed* conn = new ProjectionFixed();
											//conn->SetProjection(preUnits[i],postUnits[j]);

											// cached for performance
											if(addedPostCache == false)
											{
												postCache.push_back(postUnits[j]->GetUnitId());
												postLocalCache.push_back(postUnits[j]->GetUnitIdLocal());
												cIndex++;
												addedPostCache = true;
											}

											if(addedPreCache == false)
											{
												vector<long> vitem;
												//vitem.push_back(preUnits[i]->GetUnitId());
												vitem.push_back(preUnitId);
												preCache.push_back(vitem);
												addedPreCache = true;
											}
											else
												preCache[cIndex].push_back(preUnitId);
											//preCache[cIndex].push_back(preUnits[i]->GetUnitId());
										}
									}
								}

								//TimingStop("IterPreAdd");
							}
						}

						if(specificIndexes.size()==0)
							currentPreIndex++;
					}
				}

			}

			//TimingStop("IterPre");
			//}

			//TimingStart("AddProjections");
			//for(int m=0;m<postCache.size();m++)
			//{
			int m=preCache.size()-1;
			//m_network->AddProjections(preCache[m],postCache[m]);	// alternative: put the actual Projections in Network

			if(preCache.size()>0)
			{
				conn->AddProjections(preCache[m],postCache[m],postLocalCache[m]);			// then put the description on what is in Network in Projections


				for(int n=0;n<preCache[m].size();n++)
				{
					SetWeightValues(preCache[m][n],postCache[m]);

					
					// adds the available Projection events
					// the explicit check is temporary, will be put outside in the objects themselves

					for(int k=0;k<m_eventProjections.size();k++)
					{
						if(m_eventProjections[k]->GetEventId() == 1) // BCPNN-online
						{
							if(bcpnnInit == false)
							{	
								bcpnnInit = true;
								eBcpnnOnline = (ProjectionModifierBcpnnOnline*)m_eventProjections[k];
								conn->AddEvent(eBcpnnOnline,"bcpnn");
								conn->AddUnitsProperty(eBcpnnOnline->GetTransferFunction());
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 2) // MI mincolumns
						{
							if(miMcInit == false)
							{
								miMcInit = true;
								eMIRateUnits = (ProjectionModifierMIRateUnit*)m_eventProjections[k];
								conn->AddEvent(eMIRateUnits,"miminicolumn");
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 3) // MI hypercolumns
						{
							if(miHcInit == false)
							{
								miHcInit = true;
								eMIHypercolumns = (ProjectionModifierMIHypercolumn*)m_eventProjections[k];
								conn->AddEvent(eMIHypercolumns,"mihypercolumn");
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 4) // MDS (hypercolumns)
						{
							if(mdsInit == false)
							{
								mdsInit = true;
								eMDS = (ProjectionModifierMDS*)m_eventProjections[k];
								conn->AddEvent(eMDS,"mds");
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 5) // VQ shell class
						{
							if(vqInit == false)
							{
								vqInit = true;
								eVQ = (ProjectionModifierVQ*)m_eventProjections[k];
								conn->AddEvent(eVQ,"vq");
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 6) // Competitive Learning
						{
							if(clInit == false)
							{	
								clInit = true;
								eCL = (ProjectionModifierCL*)m_eventProjections[k];
								conn->AddEvent(eCL,"cl");
								conn->AddUnitsProperty(eCL->GetTransferFunction());
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 7) // Kussul (delta rule)
						{
							if(ksInit == false)
							{	
								ksInit = true;
								eKussul = (ProjectionModifierKussul*)m_eventProjections[k];
								conn->AddEvent(eKussul,"ks");
								conn->AddUnitsProperty(eKussul->GetTransferFunction());
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 8) // Foldiak
						{
							if(foInit == false)
							{	
								foInit = true;
								eFoldiak = (ProjectionModifierFoldiak*)m_eventProjections[k];
								conn->AddEvent(eFoldiak,"fo");
								conn->AddUnitsProperty(eFoldiak->GetTransferFunction());
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 9) // CSL used for clustering and feature extraction
						{
							if(cslInit == false)
							{	
								cslInit = true;
								eCSL = (ProjectionModifierCSL*)m_eventProjections[k];
								conn->AddEvent(eCSL,"csl");
								conn->AddUnitsProperty(eCSL->GetTransferFunction());
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 10) // Triesch
						{
							if(trieschInit == false)
							{	
								trieschInit = true;
								eTriesch = (ProjectionModifierTriesch*)m_eventProjections[k];
								conn->AddEvent(eTriesch,"tr");
								conn->AddUnitsProperty(eTriesch->GetTransferFunction());
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 11) // Pearson's correlation
						{
							if(pearsonInit == false)
							{	
								pearsonInit = true;
								ePearson = (ProjectionModifierPearson*)m_eventProjections[k];
								conn->AddEvent(ePearson,"pearson");
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 12) // standard hebb with weight normalization
						{
							if(hebbSimpleInit == false)
							{	
								hebbSimpleInit = true;
								eHebbSimple = (ProjectionModifierHebbSimple*)m_eventProjections[k];
								conn->AddEvent(eHebbSimple,"hs");
								conn->AddUnitsProperty(eHebbSimple->GetTransferFunction());
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 13) // BCM
						{
							if(bcmInit == false)
							{	
								bcmInit = true;
								eBCM = (ProjectionModifierBCM*)m_eventProjections[k];
								conn->AddEvent(eBCM,"bcm");
								conn->AddUnitsProperty(eBCM->GetTransferFunction());
							}
						}
						else if(m_eventProjections[k]->GetEventId() == 14) // Sanger / Oja
						{
							if(sangerInit == false)
							{	
								sangerInit = true;
								eSanger = (ProjectionModifierSanger*)m_eventProjections[k];
								conn->AddEvent(eSanger,"sanger");
								conn->AddUnitsProperty(eSanger->GetTransferFunction());
							}

						}
						else if(m_eventProjections[k]->GetEventId() == 15) // STDP
						{
							if(stdpInit == false)
							{	
								stdpInit = true;
								eSTDP = (ProjectionModifierSTDP*)m_eventProjections[k];
								conn->AddEvent(eSTDP,"stdp");
							}
						}
						else
						{
							cout<<"Warning: EventId not recognized." << m_eventProjections[k]->GetEventId();
						}
					}
				}
			}
		}
		
		if(network()->MPIGetNodeId() == 0)
		{
			cout<<"ok\n";
			cout.flush();
		}

		conn->network(network());
		
		// first time initializations need to be run after hash tables generated
		if(bcpnnInit == true)
			eBcpnnOnline->Initialize(conn);
		if(mdsInit == true)
			eMDS->Initialize(conn);
		if(miHcInit == true)
		{
			eMIHypercolumns->Initialize(conn);
			eMIHypercolumns->network(network());
		}
		if(miMcInit == true)
			eMIRateUnits->Initialize(conn); // mi hypercolumn must currently be initialized first (put both initializations in common function/class)
		if(clInit == true)
			eCL->Initialize(conn);
		if(ksInit == true)
			eKussul->Initialize(conn);
		if(foInit == true)
			eFoldiak->Initialize(conn);
		if(cslInit == true)
			eCSL->Initialize(conn);
		if(trieschInit == true)
			eTriesch->Initialize(conn);
		if(pearsonInit == true)
			ePearson->Initialize(conn);
		if(hebbSimpleInit == true)
			eHebbSimple->Initialize(conn);
		if(bcmInit == true)
			eBCM->Initialize(conn);
		if(sangerInit == true)
			eSanger->Initialize(conn);
		if(stdpInit == true)
			eSTDP->Initialize(conn);

		for(int i=0;i<this->m_unitsProperties.size();i++)
			conn->AddUnitsProperty(m_unitsProperties[i]);

		if(m_symmetric == true)
		{
			// also copies all synapses across nodes if needed
			if(preCache.size()>0 && postCache.size()>0)
				MakeProjectionsSymmetric(&preCache,&postCache,conn);
		}
	}

	m_initialized = true;
	

	TimingStop(m_name);
}

/// - Assumes an MPI_WORLD distribution
/// - Assumes ordered distribution of units in a population
/// Also, only legitimate to use in cases where we have recurrent connectivity
/// Quite specifically made for the random Projections atm.
void Connectivity::MakeProjectionsSymmetric(vector<vector<long> >* preCache, vector<long>* postCache, ProjectionFixed* conn)
{
	// send version
	// only gather all data send to node
	vector<int> nrPreRateUnits = ((PopulationColumns*)m_pre)->GetStructure();
	vector<long> postsToAdd;
	vector<vector<long> > presToAdd;
	vector<vector<float> > weightsToAdd;
	vector<vector<float> > delaysToAdd;

	vector<int> procsUsed = conn->PreLayer()->MPIGetProcessesUsed();

	MPI_Comm *comm = m_post->MPI()->MPIGetCommLayer();

	bool cont = false;
	if(procsUsed.size() == 0)
		cont = true;
	else if(binary_search(procsUsed.begin(),procsUsed.end(),this->network()->MPIGetNodeId()))
		cont = true;

	if(cont == true) // only continue if this process takes part in the population
	{
		// assumes an ordered distribution
		long startId = ((PopulationColumns*)m_post)->GetUnits()[0]->GetUnitId();
		long endId = ((PopulationColumns*)m_post)->GetUnits()[((PopulationColumns*)m_post)->GetUnits().size()-1]->GetUnitId();

		// assumes that it is

		int totNr = procsUsed.size();
		if(procsUsed.size()==0)
			totNr = this->network()->MPIGetNrProcs();

		for(int i=0;i<totNr;i++)//this->network()->MPIGetNrProcs();i++)
		{
			vector<long> tempPost;
			vector<long> tempPostLocal;
			vector<vector<long> > tempPre;
			vector<vector<float> > weights;
			vector<vector<float> > delays;

			int mpiRank;
			if(procsUsed.size() == 0)
				mpiRank = i;
			else
				mpiRank = procsUsed[i];

			if(procsUsed.size() == 1 || totNr == 1)		// whole population on one process
			{
				tempPost = *postCache;
				tempPre = *preCache;

				int nrPosts = (*postCache).size();
				for(int j=0;j<nrPosts;j++)
				{
					int nrPre = (*preCache)[j].size();
					vector<float> w(nrPre);
					vector<float> d(nrPre);

					for(int m=0;m<nrPre;m++)
					{
						w[m] = this->network()->GetWeight((*preCache)[j][m],(*postCache)[j]);
						d[m] = this->network()->GetDelay((*preCache)[j][m],(*postCache)[j]);
					}

					weights.push_back(w);
					delays.push_back(d);
				}
			}
			else							// broadcast the Projections and weights to other processes involved in this population
			{
				if(mpiRank == this->network()->MPIGetNodeId())
				{
					tempPost = *postCache;
					tempPre = *preCache;

					// send locally corresponding part
					int nrPosts = (*postCache).size();
					MPI_Bcast(&nrPosts,1,MPI_INT,mpiRank,*comm);//NETWORK_COMM_WORLD);
					if(nrPosts>0)
					{
						MPI_Bcast(&(*postCache)[0],nrPosts,MPI_LONG,mpiRank,*comm);//NETWORK_COMM_WORLD);

						for(int j=0;j<nrPosts;j++)
						{
							int nrPre = (*preCache)[j].size();
							MPI_Bcast(&nrPre,1,MPI_INT,mpiRank,*comm);//NETWORK_COMM_WORLD);
							vector<float> w(nrPre);
							vector<float> d(nrPre);

							for(int m=0;m<nrPre;m++)
							{
								w[m] = this->network()->GetWeight((*preCache)[j][m],(*postCache)[j]);
								d[m] = this->network()->GetDelay((*preCache)[j][m],(*postCache)[j]);
							}

							weights.push_back(w);
							delays.push_back(d);

							if(nrPre>0)
							{
								MPI_Bcast(&(*preCache)[j][0],nrPre,MPI_LONG,mpiRank,*comm);//,NETWORK_COMM_WORLD);
								MPI_Bcast(&w[0],nrPre,MPI_FLOAT,mpiRank,*comm);//,NETWORK_COMM_WORLD);
								MPI_Bcast(&d[0],nrPre,MPI_FLOAT,mpiRank,*comm);//,NETWORK_COMM_WORLD);
							}
						}
					}
				}
				else
				{
					// retrieve
					int nrPosts;
					MPI_Bcast(&nrPosts,1,MPI_INT,mpiRank,*comm);//,NETWORK_COMM_WORLD);
					if(nrPosts>0)
					{
						tempPost = vector<long>(nrPosts);
						tempPre = vector<vector<long> >(nrPosts);
						weights = vector<vector<float> >(nrPosts);
						delays = vector<vector<float> >(nrPosts);

						MPI_Bcast(&tempPost[0],nrPosts,MPI_LONG,mpiRank,*comm);//,NETWORK_COMM_WORLD);

						for(int j=0;j<nrPosts;j++)
						{
							int nrPre;
							MPI_Bcast(&nrPre,1,MPI_INT,mpiRank,*comm);//,NETWORK_COMM_WORLD);
							if(nrPre>0)
							{
								tempPre[j] = vector<long>(nrPre);
								weights[j] = vector<float>(nrPre);
								delays[j] = vector<float>(nrPre);

								MPI_Bcast(&tempPre[j][0],nrPre,MPI_LONG,mpiRank,*comm);//,NETWORK_COMM_WORLD);
								MPI_Bcast(&weights[j][0],nrPre,MPI_FLOAT,mpiRank,*comm);//,NETWORK_COMM_WORLD);
								MPI_Bcast(&delays[j][0],nrPre,MPI_FLOAT,mpiRank,*comm);//,NETWORK_COMM_WORLD);
							}
						}
					}
				}
			}

			vector<long>::iterator it;
			// put in list to build
			for(int j=0;j<tempPost.size();j++)
			{
				for(int m=0;m<tempPre[j].size();m++)
				{
					if(tempPre[j][m] >= startId && tempPre[j][m] <= endId)
					{
						it = find(postsToAdd.begin(),postsToAdd.end(),tempPre[j][m]);
						if(it == postsToAdd.end())
						{
							postsToAdd.push_back(tempPre[j][m]);
							vector<long> preAdd(1);
							vector<float> weightAdd(1);
							vector<float> delayAdd(1);
							preAdd[0] = tempPost[j];
							weightAdd[0] = weights[j][m];
							delayAdd[0] = delays[j][m];
							weightsToAdd.push_back(weightAdd);
							delaysToAdd.push_back(delayAdd);
							presToAdd.push_back(preAdd);
						}
						else
						{
							presToAdd[it-postsToAdd.begin()].push_back(tempPost[j]);
							weightsToAdd[it-postsToAdd.begin()].push_back(weights[j][m]);
							delaysToAdd[it-postsToAdd.begin()].push_back(delays[j][m]);
						}
					}
				}
			}
		}

		// build new Projections at the same time

		for(int i=0;i<postsToAdd.size();i++)
		{
			conn->AddProjections(presToAdd[i],postsToAdd[i],-1); // -1 will force a lookup of the local id value

			for(int j=0;j<presToAdd[i].size();j++)
			{

				m_network->SetWeight(weightsToAdd[i][j],presToAdd[i][j],postsToAdd[i]);
				//m_network->SetDelay(delaysToAdd[i][j],presToAdd[i][j],postsToAdd[i]);
			}
		}

		// check symmetry
		bool checkSymmetry = true;
		if(checkSymmetry)
		{
			vector<long> postIds = conn->GetPostIds();
			bool symmetric = true;
			for(int i=0;i<postIds.size();i++)
			{
				for(int j=0;j<postIds.size();j++)
				{
					if(i!=j)
					{
						vector<long> preIds1 = conn->GetPreIds(postIds[i]);
						vector<long> preIds2 = conn->GetPreIds(postIds[j]);
						if(find(preIds2.begin(),preIds2.end(),postIds[i]) != preIds2.end())
						{
							if(find(preIds1.begin(),preIds1.end(),postIds[j]) == preIds1.end())
								symmetric = false;
							//							if(find(preIds2.begin(),preIds2.end(),postIds[i]) == preIds2.end())
							//								symmetric = false;
						}
					}
				}
			}

			if(m_network->MPIGetNodeId() == 0)
			{
				cout<<"\nSymmetric == "<<symmetric<<"\n";cout.flush();
			}
		}
	}
}

bool FullConnectivityNoLocalHypercolumns::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	if(((RateUnit*)preUnit)->GetHypercolumnId() == ((RateUnit*)postUnit)->GetHypercolumnId())
		return false;
	else
		return true;
}

// fast version, but not supporting symmetric Projections
vector<long> RandomConnectivity::IterationSpecificIndexes(Unit* postUnit, vector<int> nrPreRateUnits)
{
	m_useIterationSpecificIndexes = true;
	long maxIndexes = 0;
	for(int i=0;i<nrPreRateUnits.size();i++)
		maxIndexes+=nrPreRateUnits[i];

	int nrItems = (int)(maxIndexes*this->m_fracConnected);
	vector<long> out;

	// better scaling for sparse conns
	int nrToAdd = nrItems;
	while(out.size()!=nrItems)
	{
		for(int i=0;i<nrToAdd;i++)
		{
			float r = rand()/(float(RAND_MAX));
			long loc = (long)(((float)maxIndexes)*r); // float limit of population size (*100 .. / 100)

			out.push_back(loc);
		}

		sort(out.begin(),out.end());

		// remove duplicates
		out.erase(unique(out.begin(),out.end()),out.end());

		nrToAdd = nrItems-out.size();
		//if(out.size()>nrItems)
		//	out.erase(out.begin()+nrItems,out.end());
	}

	sort(out.begin(),out.end());

	return out;
}

// slow
bool RandomConnectivity::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	// currently forced or non-forced self-Projections
	if(preUnit->GetUnitId() == postUnit->GetUnitId())
	{
		if(this->m_allowSelfProjections == true)
			return true;
		else return false;
	}

	return false;
}


bool OneToOneConnectivity::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	if(this->m_useIterationSpecificIndexes == true)
		return true;

	if(this->m_betweenHypercolumns == true)
	{
		return true; // taken care of in iterationpreconstraint instead (faster)
	}
	else
	{
		long preIdLocal = preUnit->GetUnitIdLocal();
		long postIdLocal = postUnit->GetUnitIdLocal();
		if(preIdLocal == postIdLocal)
			return true;
		else
			return false;
	}
}

vector<long> OneToOneConnectivity::IterationSpecificIndexes(Unit* postUnit, vector<int> nrPreRateUnits)
{
	vector<long> list(1);
	list[0] = postUnit->GetUnitIdLocal();// + 1;

	return list; // will be sorted lowest to highest
}

bool OneToOneConnectivity::IterationPreConstraint(Unit* postUnit,int currentIndex, int currentHcIndex)
{
	return true;

	if(this->m_useIterationSpecificIndexes == true)
		return true;

	if(this->m_betweenHypercolumns == true)
	{
		if(((RateUnit*)postUnit)->GetHypercolumns()[0]->GetUnitIdLocal() == currentHcIndex)
			return true;
		else 
			return false;
	}
	else
	{
		long unitIdLocal = postUnit->GetUnitIdLocal();
		if(unitIdLocal == currentIndex - 1)
			return true;
		else
			return false;
	}
}

bool FanInConnectivity::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	if((int)(preUnit->GetUnitIdLocal()/m_nrIns) == postUnit->GetUnitIdLocal())
		return true;
	else
		return false;
}

Projection::~Projection()
{
	for(int i=0;i<m_events.size();i++)
	{
		delete m_events[i];
	}
}

bool FullConnectivity::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	if(m_allowSelfProjections == false)
	{
		if(preUnit->GetUnitId() == postUnit->GetUnitId())
			return false;
		else
			return true;
	}
	else
		return true;
}

void EmptyConnectivity::ExtraPostUnit(Unit* postUnit,Projection* conn) 
{ 
	conn->AddPostUnit(postUnit->GetUnitId());
}

UnitModifier* Projection::GetUnitModifier(int id)
{
	for(int i=0;i<m_unitPropertiesProjection.size();i++)
	{
		if(m_unitPropertiesProjection[i]->GetId() == id)
			return m_unitPropertiesProjection[i];
	}

	return NULL;
}

UnitModifier* Projection::GetUnitModifier(string name)
{
	for(int i=0;i<m_unitPropertiesProjection.size();i++)
	{
		if(m_unitPropertiesProjection[i]->GetName().compare(name) == 0)
			return m_unitPropertiesProjection[i];
	}

	return NULL;
}


// Needed to support the use of multiple parameters
void FullConnectivity::SetRandomWeights(vector<float> minVal, float maxVal,Network* net)
{
	m_setRandomWeights = true;
	m_weightsMinVal = minVal[0];
	m_weightsMaxVal = maxVal;
	m_weightsMinValParams = minVal;
	this->network(net);
	AddParameters();
}

// Needed to support the use of multiple parameters
void FullConnectivity::SetRandomWeights(float minVal, vector<float> maxVal,Network* net)
{
	m_setRandomWeights = true;
	m_weightsMinVal = minVal;
	m_weightsMaxVal = maxVal[0];
	m_weightsMaxValParams = maxVal;
	this->network(net);
	AddParameters();
}

// Needed to support the use of multiple parameters
void FullConnectivity::AddParameters()
{
	if(m_setRandomWeights == true)
	{
		if(m_weightsMinValParams.size()>0)
			this->network()->Parameters()->AddParameters(this,&FullConnectivity::SetRandomWeightsMin,m_weightsMinValParams);

		if(m_weightsMaxValParams.size()>0)
			this->network()->Parameters()->AddParameters(this,&FullConnectivity::SetRandomWeightsMax,m_weightsMaxValParams);
	}
}