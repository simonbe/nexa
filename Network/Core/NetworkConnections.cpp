#include <iostream>
//#include <random>

#include "NetworkConnections.h"
#include "NetworkUnits.h"
#include "NetworkConnectionModifier.h"
#include "NetworkUnitModifier.h"

// need to access plasticity rules here, will be changed
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
//

//#include "NetworkTransferFunctions.h"

/*m_postIds.clear();
		m_preIds.clear();

		m_postIds = GetPostIds();
		for(int i=0;i<m_postIds.size();i++)
		{
			vector<long> pres = GetPreIds(m_postIds[i]);
			m_preIds.push_back(pres);
		}*/

/*void ConnectionFixed::AddConnection(Unit* pre, Unit* post, bool firstRun)
{
	map<long, map<long, float> >::iterator p;

	if(!firstRun)
	{
		map<long,float> w = m_hashWeights[post->GetUnitId()];
		w[pre->GetUnitId()] = 0.0;
		m_hashWeights[post->GetUnitId()] = w;
	}
	else // empty
	{
		map<long,float> w;
		w[pre->GetUnitId()] = 0.0;
		m_hashWeights[post->GetUnitId()] = w;
	}

	post->AddPre(pre); // needed now?
}*/


void Connection::GenerateHashTables()
{
/*	if(m_preLayer->network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"Generating hash tables (connections)...";
	}

	// always keep postIds sorted
	sort(m_postIds.begin(),m_postIds.end());

	// weights

	srand(this->network()->MPIGetNodeId());

	if(m_postIds.size()==0)
		cout<<"Warning, m_postIds.size() == 0;";

	for(int i=0;i<m_postIds.size();i++)
	{
		map<long,float> w;

		int count = m_hashWeights.count(m_postIds[i]);//m_hashWeights.count(m_postIds[i]);
		if(count>0)
			w = m_hashWeights[m_postIds[i]];
		
		for(int j=0;j<m_preIds[i].size();j++) // preIds is index-based order
		{
			float r = 0;//(float)rand()/(float)RAND_MAX;//0.0;//(rand()/(float(RAND_MAX)+1));
			w[m_preIds[i][j]] = r;
		}
		
		m_hashWeights[m_postIds[i]] = w;
	}
	
	map<long, map<long, float> >::iterator p;

	if(!firstRun)
	{
		map<long,float> w = m_hashWeights[postId];
		w[preId] = 0.0;
		m_hashWeights[postId] = w;
	}
	else // empty
	{
		map<long,float> w;
		w[preId] = 0.0;
		m_hashWeights[postId] = w;
	}

	//post->AddPre(pre);


	// remove unused for memory
	// not implemented/necessary yet
	
	if(m_preLayer->network()->MPIGetNodeId() == 0 && DEBUG_LEVEL > 2)
	{
		cout<<"done.\n";
	}
	*/
}

void Connection::AddConnection(long preId, long postId, long postLocalId, bool firstRun)
{
	bool first = true;
	long postIndex;

/*	for(int i=0;i<m_postIds.size();i++)
	{
		if(m_postIds[i] == postId)
		{
			postIndex = i;
			first = false;
			break;
		}
	}*/

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
		m_postIdsPre[preId].push_back(postId);
		if(postLocalId<0)
			postLocalId = this->PostLayer()->GetLocalIdUnitId(postId);

		m_localPostIdsPre[preId].push_back(postLocalId);
#endif

	m_preIds[postIndex].push_back(preId);
	
	/*if(m_listPosts.size()<preId+1) // change
	{
		for(int i=m_listPosts.size();i<preId+1;i++)
			m_listPosts.push_back(new vector<long>());
	}
	
	//used? remove
	m_listPosts[preId]->push_back(postId);*/
}


void Connection::AddConnections(vector<long> preIds, long postId, long postLocalId)
{
	//for(int i=0;i<preIds.size();i++)
	//	AddConnection(preIds[i],postId,true);

	bool first = true;
	int postIndex;

	postIndex = m_postIds.size();

	if(find(m_postIds.begin(),m_postIds.end(),postId) == m_postIds.end())
		m_postIds.push_back(postId);

	if(m_preIds.find(postId) == m_preIds.end())//find(m_preIds.begin(),m_preIds.end(),postId) == m_preIds.end())
		m_preIds[postId] = preIds; // used in GetPreValues - could be replaced by indexes to synapse hash table (which belongs to which connection)
	else
	{
		vector<long> v = m_preIds[postId];
		for(int i=0;i<preIds.size();i++)
			v.push_back(preIds[i]);
		
		// remove duplicates
		sort(v.begin(),v.end());
		//vector<long>::iterator it;
		//it = unique (v.begin(), v.end());
		//v.resize( it - v.begin() );
		v.erase(unique(v.begin(),v.end()),v.end());
		m_preIds[postId] = v;
	}

#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	for(int i=0;i<preIds.size();i++)
	{
		// always guarantee no duplicates?
		// cannot be necessary to handle both
		m_postIdsPre[preIds[i]].push_back(postId);
		if(postLocalId<0)
			postLocalId = this->PostLayer()->GetLocalIdUnitId(postId);

		m_localPostIdsPre[preIds[i]].push_back(postLocalId);
	}
#endif

	// for mpi_Alltoall-communication pattern
#if USE_COMMUNICATION_ALLTOALL == 1

#endif


	/*if(m_preIds.size()<postIndex+1)
	{
		for(int i=m_preIds.size();i<postIndex+1;i++)
		{
			vector<long> v;
			m_preIds.push_back(v);
		}
	}

	//used? remove
	for(int i=0;i<preIds.size();i++)
	{
		m_preIds[postIndex].push_back(preIds[i]);
	}*/

	//m_network->AddPostIds(preIds,postId);
	
		
		// used currently in sendreceiveallgather - change?
		//if(m_listPosts.size()<preIds[i]+1) // change
		//{
		//	for(int j=m_listPosts.size();j<preIds[i]+1;j++)
		//		m_listPosts.push_back(new vector<long>());
		//}

		//used? remove
		//m_listPosts[preIds[i]]->push_back(postId);
}

void Connection::Clear()
{
	/*
	for(int i=0;i<m_listPosts.size();i++)
	{
		m_listPosts[i]->clear();
		delete m_listPosts[i];
	}

	m_listPosts.clear();
	*/

	//for(int i=0;i<m_preIds.size();i++)
#if USE_UNORDERED_MAP == 1
	for(unordered_map<long, vector<long> >::iterator iter = m_preIds.begin();iter!=m_preIds.end();++iter)
#else
	for(map<long, vector<long> >::iterator iter = m_preIds.begin();iter!=m_preIds.end();++iter)
#endif
	{
		for(int j=0;j<m_preIds[iter->first].size();j++)
			m_network->EraseSynapseValues(m_preIds[iter->first][j],iter->first);//iter->first,m_preIds[iter->first][j]);

		//m_preIds[i].clear();
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
	
// old way
/*	for(map<long,map<long,float> >::iterator iter =m_hashWeights.begin();iter!=m_hashWeights.end();++iter)
	{
		long postId = (*iter).first;
		Clear(postId);	
	}
	*/
}

void Connection::Clear(long postId)
{
	map<long,float> empty;
//	m_hashWeights[postId] = empty;
	// also remove from preIds?

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

void Connection::ClearActiveBuffer()
{
#if USE_HASHED_ACTIVE_COMMUNICATION == 1
	m_preIdsActive.clear();
	m_preIdsActive = vector<vector<pair<long, float> > >(this->PostLayer()->GetNrUnitsTotal());//this->PostLayer()->GetUnits().size());
	//		m_postIdsPre.clear();
#endif
}

void Connection::AddEvent(ConnectionModifier* e, string hashName)
{
	e->SetConnection(this);
	m_events.push_back(e);
	m_hashEvents[hashName] = e;
}

void Connection::SetRandomWeights(float minVal, float maxVal)
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

void Connection::SetRandomWeightsIn(long localPostIndex, float minVal, float maxVal)
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

void Connection::SetRandomWeightsOut(long localPreIndex, float minVal, float maxVal)
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

void ConnectionFixed::ModifyConnection()
{
//	if(IsOn())
//	{
		for(int i=0;i<m_events.size();i++)
		{
			//ConnectionModifier* e = m_events[i];
			if(m_events[i]->IsOn())
				m_events[i]->Modify();
			/*if(m_events[i]->GetEventId() == 1)
			{
			ConnectionModifierBcpnnOnline* e = (ConnectionModifierBcpnnOnline*)m_events[i];
			e->Modify();
			}
			else if(m_events[i]->GetEventId() == 1)
			{
			ConnectionModifierBcpnnOnline* e = (ConnectionModifierBcpnnOnline*)m_events[i];
			e->Modify();
			}*/
		}
//	}

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
				vector<long> preIds = this->GetPreIds(postId); // does this function exist anymore ?
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

void ConnectionFixed::SimulateEvent(UnitModifier* e)
{
	for(int i=0;i<m_events.size();i++)
	{
		//ConnectionModifier* e = m_events[i];
		if(m_events[i]->GetEventId() == 1)
		{
			ConnectionModifierBcpnnOnline* b = (ConnectionModifierBcpnnOnline*)m_events[i];
			b->Simulate(e);

			//ConnectionModifierBcpnnOnline e2;
			//e2.Simulate();
		}
	}
}

vector<float> Connection::GetPreValues(long unitId)
{
/*	map<long, Network::SynapseStandard>* preSynapses = m_network->GetPreSynapses(unitId);	
	vector<float> out(preSynapses->size());
	map<long, Network::SynapseStandard>::iterator it;
	int i = 0;
	for(it = preSynapses->begin(); it!=preSynapses->end(); it++)
	{
		out[i] = m_network->GetPreValue(it->first);
		i++;
	}

	return out;*/
	
	vector<long> preIds = m_preIds[unitId]; // change to pointer
	vector<float> out(preIds.size());
	
	for(int i=0;i<preIds.size();i++)
	{
		out[i] = m_network->GetPreValue(preIds[i]);
	}

	return out;
}

// re-implement
vector<long> Connection::GetPreLocalIds(long unitId)
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

/*	map<long,float> preUnitWeights = m_hashWeights[unitId];
	vector<long> indexes(preUnitWeights.size());

	int i=0;
	for(map<long,float>::iterator iter = preUnitWeights.begin();iter!=preUnitWeights.end();++iter)
	{
		long preId = (*iter).first;
		long localId = ((RateUnit*)m_preLayer->network()->GetUnitFromId(preId))->GetUnitIdLocal();

		indexes[i] = localId;//preId;
		i++;
	}

	return indexes;*/

	return localPreIds;//vector<long>(0);
}

vector<long> Connection::GetPreIds(long unitId)
{
	// fix this -> should be index-based (?)
	//for(int i=0;i<m_postIds.size
	return m_preIds[unitId];//unitId-m_postIds[0]];
}


void Connection::AddActiveEvent(long preId, float value)
{
	vector<long> postIds;
	// get all the units this unit is connected to on this connection
	postIds = GetPostIds(preId);
	
	// put in their active queues/buffers
	for(int i=0;i<postIds.size();i++)
	{
		pair<long, float> p;
		p.first = preId;
		p.second = value;
		// put in vector with position determined by hash (set up once)
		m_preIdsActive[postIds[i]].push_back(p);//preId);
	}
}

void Connection::AddActiveEvents(vector<long> preIds, vector<float> values)
{
	vector<vector<long> > localPostIds(preIds.size());
	m_preIdsActive = vector<vector<pair<long,float> > >(this->PostLayer()->GetNrUnitsTotal());//this->PostLayer()->GetUnits().size());

	for(int i=0;i<preIds.size();i++)
	{
		// get all the units this unit is connected to on this connection
		localPostIds[i] = GetLocalPostIds(preIds[i]);//GetPostIds(preIds[i]);
		for(int j=0;j<localPostIds[i].size();j++)
		{
			m_preIdsActive[localPostIds[i][j]].reserve(preIds.size()); // optimization
		}
	}

	pair<long, float> p;
	for(int i=0;i<preIds.size();i++)
	{
		// put in their active queues/buffers
		//pair<long, float> p; // new here will take time (track in profiler)
		p.first = preIds[i];
		p.second = values[i];
		long localId;

		for(int j=0;j<localPostIds[i].size();j++)//postIds.size();j++)
		{
			//localId = this->PostLayer()->GetLocalIdUnitId(postIds[j]);
			// put in vector with position determined by hash (set up once) - optimize (!) (remove new + slow hash insert)
			m_preIdsActive[localPostIds[i][j]].push_back(p);
			//m_preIdsActive[localId].push_back(p);//preId);
			//m_preIdsActive[postIds[j]].push_back(p);//preId);

			// organized by local id
		}
	}
}

/*vector<pair<long,float> >* Connection::GetPreIdsActive(long unitId)
{
	// fix this -> should be index-based (?)
	//for(int i=0;i<m_postIds.size
	return &m_preIdsActive[unitId];//unitId-m_postIds[0]];
}*/

vector<pair<long,float> >* Connection::GetPreIdsActiveLocal(long localUnitId)
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
vector<long> Connection::GetPostIds(long preId)
{
#if	USE_HASHED_ACTIVE_COMMUNICATION == 1
	return m_postIdsPre[preId];
#else
	return vector<long>(0);
#endif
}

vector<long> Connection::GetLocalPostIds(long preId)
{
#if	USE_HASHED_ACTIVE_COMMUNICATION == 1
	return m_localPostIdsPre[preId];
#else
	return vector<long>(0);
#endif
}

void Connection::EraseSynapses(vector<pair<long,long> > vectorPostIdPreId)
{
	for(int i=0;i<vectorPostIdPreId.size();i++)
	{
		pair<long,long> p = vectorPostIdPreId[i];
		EraseSynapse(p.first,p.second);
	}

	// recreate hash (ids of pre-units) union
	CreatePreIdsUnion();
}

void Connection::EraseSynapse(long postId, long preId)
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
	vector<long> postIds = m_postIdsPre[postId];

	for(int i=0;i<postIds.size();i++)
	{
		if(postIds[i] == postId)
		{
			index2 = i;
			break;
		}
	}
	if(index2 !=-1)
	{
		postIds.erase(postIds.begin()+index2);
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

vector<long> Connection::GetPostLocalIds()
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

vector<long> Connection::GetPostIds()
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

vector<long>* Connection::GetPreIdsAll()//GetPreIdsUnion()
{
	if(m_preIdsUnion.size() == 0)
	{
		CreatePreIdsUnion();

		return &m_preIdsUnion;
	}
	else return &m_preIdsUnion;

}

// does same as CreateAllPreIdsUnion but only these connections, could merge for build performance (!)
// default: called in initialization phase
void Connection::CreatePreIdsUnion()
{
	m_preIdsUnion.clear();
	m_totalNrLocalConnections = 0; //

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
				m_totalNrLocalConnections++;
			}
		}
	}

	map<long,int>::iterator it3;

	for(it3 = tempMap.begin();it3!=tempMap.end();it3++)
		m_preIdsUnion.push_back(it3->first);

	if(m_network->MPIGetNodeId() == 0)
	{
		cout<<"Total nr incoming local connections (layer "<<m_name<<") == "<<m_totalNrLocalConnections<<"\n";
		cout.flush();
	}
}

/*vector<long>* Connection::GetPostIds(long preId)
{
	if(m_listPosts.size() == 0) // will be 0 if node has no local units in layer
		return NULL;
	else if(preId >= m_listPosts.size())
		return NULL;
	else
		return m_listPosts[preId]; // use pointer instead
}*/

vector<float> Connection::GetPostValues()
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

vector<float> Connection::GetPreValuesAll()
{
	vector<float> out(m_preIdsUnion.size());//network()->HashWeights()->size());

	int i=0;

	for(int i=0;i<m_preIdsUnion.size();i++)
	{
		out[i] = m_network->GetPreValue(m_preIdsUnion[i]);//GetIncomingBufferData(m_preIdsUnion[i]);
	}

	return out;	
}

/*float Connection::GetWeight(long preId, long postId)
{
	return (*network()->HashWeights())[postId][preId];
}

void Connection::SetWeight(float weight, long preId, long postId)
{
	(*network()->HashWeights())[postId][preId] = weight;
}*/

std::vector<std::vector<float> > Connection::GetValuesToRecord()
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

void Connection::CopyConnectionsOtherPost(Population* newPost)
{
	// Get preids and postids

	// Change list of post ids to new population

	// use connectivitytype::AddConnections

	// go through all connections, copy weight values etc (optional)
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

void Connectivity::InitializeWeightsAndConnectionModifiers()
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

		ConnectionFixed* conn = new ConnectionFixed();
		conn->Initialize(m_pre,m_post);
		conn->network(this->network());
		conn->SetConnectivity(this);

		m_pre->AddOutgoingConnection(conn);
		m_post->AddIncomingConnection(conn);

		ConnectionModifierBcpnnOnline* eBcpnnOnline; // allocated upon request
		ConnectionModifierMIRateUnit* eMIRateUnits;
		ConnectionModifierMIHypercolumn* eMIHypercolumns;
		ConnectionModifierMDS* eMDS;
		ConnectionModifierVQ* eVQ;
		ConnectionModifierCL* eCL;
		ConnectionModifierKussul* eKussul;
		ConnectionModifierFoldiak* eFoldiak;
		ConnectionModifierCSL* eCSL;
		ConnectionModifierTriesch* eTriesch;
		ConnectionModifierPearson* ePearson;
		ConnectionModifierHebbSimple* eHebbSimple;
		ConnectionModifierBCM* eBCM;
		ConnectionModifierSanger* eSanger;
		ConnectionModifierSTDP* eSTDP;

		bool bcpnnInit = false, miMcInit = false, miHcInit = false, mdsInit = false, vqInit = false, clInit = false, ksInit = false, foInit = false, cslInit = false, trieschInit = false, pearsonInit = false, hebbSimpleInit = false, bcmInit = false, sangerInit = false, stdpInit = false;

		// pre och post alla local för hypercolumn! borde inte vara så, kolla resetlocalities etc.

		long connectionId = 0;
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

			if(isOverHypercolumns == true) // build connections over the hypercolumns
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

								steps = IterationAdvance(); // fastest way to put constraint on connections (but cannot use any properties)

								// create new preUnit so that also non-local units get checked (viable while using ExtraConstraint to build connections)
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
											//ConnectionFixed* conn = new ConnectionFixed();
											//conn->SetConnection(preUnits[i],postUnits[j]);

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

			//TimingStart("AddConnections");
			//for(int m=0;m<postCache.size();m++)
			//{
			int m=preCache.size()-1;
			//m_network->AddConnections(preCache[m],postCache[m]);	// alternative: put the actual connections in Network

			if(preCache.size()>0)
			{
				conn->AddConnections(preCache[m],postCache[m],postLocalCache[m]);			// then put the description on what is in Network in Connections


				for(int n=0;n<preCache[m].size();n++)
				{

					//conn->AddConnection(preUnits[i]->GetUnitId(),postUnits[j]->GetUnitId(),firstRun);//, conn->GetCurrentConnectionId());
					//conn->AddConnection(preUnits[i]->GetUnitId(),postUnits[j]->GetUnitId(),firstRun);//, conn->GetCurrentConnectionId());

					SetWeightValues(preCache[m][n],postCache[m]);


					//float r = (rand()/(float(RAND_MAX)+1))/2;
					//network()->SetWeight(0.0,conn->GetCurrentConnectionId());
					//conn->NextConnectionId();

					//network()->SetWeight(0.0);//r);//1.0);

					/*					if(preUnits[i]->GetNodeId() == m_network->MPIGetNodeId())
					preUnits[i]->SetLocal(true);
					else
					preUnits[i]->SetLocal(false);

					if(postUnits[j]->GetNodeId() == m_network->MPIGetNodeId())
					postUnits[j]->SetLocal(true);
					else
					postUnits[j]->SetLocal(false);
					*/
					// connection events - will be made object neutral

					for(int k=0;k<m_eventConnections.size();k++)
					{
						if(m_eventConnections[k]->GetEventId() == 1) // BCPNN-online
						{
							if(bcpnnInit == false)
							{	
								bcpnnInit = true;
								eBcpnnOnline = (ConnectionModifierBcpnnOnline*)m_eventConnections[k];//new ConnectionModifierBcpnnOnline(*((ConnectionModifierBcpnnOnline*)m_eventConnections[k])); // deep copy
								conn->AddEvent(eBcpnnOnline,"bcpnn");
								conn->AddUnitsProperty(eBcpnnOnline->GetTransferFunction()); //conn->PostLayer()->AddUnitsProperty(eBcpnnOnline->GetTransferFunction()); // bcpnn transfer function
							}

							/*ConnectionModifierBcpnnOnline* e = new ConnectionModifierBcpnnOnline(*((ConnectionModifierBcpnnOnline*)m_eventConnections[k])); // deep copy
							e->Initialize();
							conn->AddEvent(e,"bcpnn");
							conn->Post()->AddUnitModifier(t); // bcpnn transfer function
							*/
						}
						else if(m_eventConnections[k]->GetEventId() == 2) // MI mincolumns
						{
							if(miMcInit == false)
							{
								miMcInit = true;
								eMIRateUnits = (ConnectionModifierMIRateUnit*)m_eventConnections[k];//new ConnectionModifierMIRateUnit(*((ConnectionModifierMIRateUnit*)m_eventConnections[k]));
								//eMIRateUnits->Initialize();
								conn->AddEvent(eMIRateUnits,"miminicolumn");
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 3) // MI hypercolumns
						{
							if(miHcInit == false)
							{
								miHcInit = true;
								eMIHypercolumns = (ConnectionModifierMIHypercolumn*)m_eventConnections[k];//new ConnectionModifierMIHypercolumn(*((ConnectionModifierMIHypercolumn*)m_eventConnections[k]));
								//eMIRateUnits->Initialize();
								conn->AddEvent(eMIHypercolumns,"mihypercolumn");
							}
							/*ConnectionModifierMIHypercolumn* e = new ConnectionModifierMIHypercolumn(*((ConnectionModifierMIHypercolumn*)m_eventConnections[k]));
							e->Initialize();
							conn->AddEvent(e,"mihypercolumn");*/
						}
						else if(m_eventConnections[k]->GetEventId() == 4) // MDS (hypercolumns)
						{
							if(mdsInit == false)
							{
								mdsInit = true;
								eMDS = (ConnectionModifierMDS*)m_eventConnections[k];//new ConnectionModifierMDS(*((ConnectionModifierMDS*)m_eventConnections[k]));
								//e->Initialize();
								conn->AddEvent(eMDS,"mds");
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 5) // VQ shell class
						{
							if(vqInit == false)
							{
								vqInit = true;
								eVQ = (ConnectionModifierVQ*)m_eventConnections[k];
								conn->AddEvent(eVQ,"vq");
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 6) // Competitive Learning
						{
							if(clInit == false)
							{	
								clInit = true;
								eCL = (ConnectionModifierCL*)m_eventConnections[k];
								conn->AddEvent(eCL,"cl");
								conn->AddUnitsProperty(eCL->GetTransferFunction());//conn->PostLayer()->AddUnitsProperty(eCL->GetTransferFunction());
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 7) // Kussul (delta rule)
						{
							if(ksInit == false)
							{	
								ksInit = true;
								eKussul = (ConnectionModifierKussul*)m_eventConnections[k];
								conn->AddEvent(eKussul,"ks");
								conn->AddUnitsProperty(eKussul->GetTransferFunction());//conn->PostLayer()->AddUnitsProperty(eKussul->GetTransferFunction());
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 8) // Foldiak
						{
							if(foInit == false)
							{	
								foInit = true;
								eFoldiak = (ConnectionModifierFoldiak*)m_eventConnections[k];
								conn->AddEvent(eFoldiak,"fo");
								conn->AddUnitsProperty(eFoldiak->GetTransferFunction());//conn->PostLayer()->AddUnitsProperty(eFoldiak->GetTransferFunction());
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 9) // CSL used for feature extraction
						{
							if(cslInit == false)
							{	
								cslInit = true;
								eCSL = (ConnectionModifierCSL*)m_eventConnections[k];
								conn->AddEvent(eCSL,"csl");
								conn->AddUnitsProperty(eCSL->GetTransferFunction());//conn->PostLayer()->AddUnitsProperty(eCSL->GetTransferFunction());
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 10) // Triesch
						{
							if(trieschInit == false)
							{	
								trieschInit = true;
								eTriesch = (ConnectionModifierTriesch*)m_eventConnections[k];
								conn->AddEvent(eTriesch,"tr");
								conn->AddUnitsProperty(eTriesch->GetTransferFunction());//conn->PostLayer()->AddUnitsProperty(eTriesch->GetTransferFunction());
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 11) // Pearson's correlation
						{
							if(pearsonInit == false)
							{	
								pearsonInit = true;
								ePearson = (ConnectionModifierPearson*)m_eventConnections[k];
								conn->AddEvent(ePearson,"pearson");
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 12) // standard hebb with weight normalization
						{
							if(hebbSimpleInit == false)
							{	
								hebbSimpleInit = true;
								eHebbSimple = (ConnectionModifierHebbSimple*)m_eventConnections[k];
								conn->AddEvent(eHebbSimple,"hs");
								conn->AddUnitsProperty(eHebbSimple->GetTransferFunction());//conn->PostLayer()->AddUnitsProperty(eHebbSimple->GetTransferFunction());
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 13) // BCM
						{
							if(bcmInit == false)
							{	
								bcmInit = true;
								eBCM = (ConnectionModifierBCM*)m_eventConnections[k];
								conn->AddEvent(eBCM,"bcm");
								conn->AddUnitsProperty(eBCM->GetTransferFunction());//conn->PostLayer()->AddUnitsProperty(eBCM->GetTransferFunction());
							}
						}
						else if(m_eventConnections[k]->GetEventId() == 14) // Sanger / Oja
						{
							if(sangerInit == false)
							{	
								sangerInit = true;
								eSanger = (ConnectionModifierSanger*)m_eventConnections[k];
								conn->AddEvent(eSanger,"sanger");
								conn->AddUnitsProperty(eSanger->GetTransferFunction());//conn->PostLayer()->AddUnitsProperty(eSanger->GetTransferFunction());
							}

						}
						else if(m_eventConnections[k]->GetEventId() == 15) // STDP
						{
							if(stdpInit == false)
							{	
								stdpInit = true;
								eSTDP = (ConnectionModifierSTDP*)m_eventConnections[k];
								conn->AddEvent(eSTDP,"stdp");
								//conn->PostLayer()->AddUnitsProperty(eSTDP->GetTransferFunction());
							}
						}
						else
						{
							cout<<"Warning: EventId not recognized." << m_eventConnections[k]->GetEventId();
						}
					}
				}


				//preCache[m].clear(); // lower temporary memory consumption
			}
		}

		//TimingStop("AddConnections");


		if(network()->MPIGetNodeId() == 0)
		{
			cout<<"ok\n";
			cout.flush();
		}

		conn->network(network());
		//conn->GenerateHashTables();

		// first time need to be run after hash tables generated

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
			// copy all synapses across nodes
			if(preCache.size()>0 && postCache.size()>0)
				MakeConnectionsSymmetric(&preCache,&postCache,conn);
		}
	}

	m_initialized = true;
	
//	if(m_symmetric == true)
//	{
	//	this->network()->SetSeed(false);
//	}

	TimingStop(m_name);
}

/// - Assumes an MPI_WORLD distribution
/// - Assumes ordered distribution of units in a population
/// Also, only legitimate to use in cases where we have recurrent connectivity
/// Quite specifically made for the random connections atm.
void Connectivity::MakeConnectionsSymmetric(vector<vector<long> >* preCache, vector<long>* postCache, ConnectionFixed* conn)
{
	if(true)//this->network()->MPIGetNrProcs()>1)
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
				if(procsUsed.size())
					mpiRank = i;
				else
					mpiRank = procsUsed[i];

				if(procsUsed.size() == 1)		// whole population on one process
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
				else							// broadcast the connections and weights to other processes involved in this population
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
									MPI_Bcast(&(*preCache)[j][0],nrPre,MPI_LONG,mpiRank,NETWORK_COMM_WORLD);
									MPI_Bcast(&w[0],nrPre,MPI_FLOAT,mpiRank,NETWORK_COMM_WORLD);
									MPI_Bcast(&d[0],nrPre,MPI_FLOAT,mpiRank,NETWORK_COMM_WORLD);
								}
							}
						}
					}
					else
					{
						// retrieve
						int nrPosts;
						MPI_Bcast(&nrPosts,1,MPI_INT,mpiRank,NETWORK_COMM_WORLD);
						if(nrPosts>0)
						{
							tempPost = vector<long>(nrPosts);
							tempPre = vector<vector<long> >(nrPosts);
							weights = vector<vector<float> >(nrPosts);
							delays = vector<vector<float> >(nrPosts);

							MPI_Bcast(&tempPost[0],nrPosts,MPI_LONG,mpiRank,NETWORK_COMM_WORLD);

							for(int j=0;j<nrPosts;j++)
							{
								int nrPre;
								MPI_Bcast(&nrPre,1,MPI_INT,mpiRank,NETWORK_COMM_WORLD);
								if(nrPre>0)
								{
									tempPre[j] = vector<long>(nrPre);
									weights[j] = vector<float>(nrPre);
									delays[j] = vector<float>(nrPre);

									MPI_Bcast(&tempPre[j][0],nrPre,MPI_LONG,mpiRank,NETWORK_COMM_WORLD);
									MPI_Bcast(&weights[j][0],nrPre,MPI_FLOAT,mpiRank,NETWORK_COMM_WORLD);
									MPI_Bcast(&delays[j][0],nrPre,MPI_FLOAT,mpiRank,NETWORK_COMM_WORLD);
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

			// build new connections at the same time

			for(int i=0;i<postsToAdd.size();i++)
			{
				conn->AddConnections(presToAdd[i],postsToAdd[i],-1); // -1 will force a lookup of the local id value

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
}

bool FullConnectivityNoLocalHypercolumns::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	if(((RateUnit*)preUnit)->GetHypercolumnId() == ((RateUnit*)postUnit)->GetHypercolumnId())
		return false;
	else
		return true;
}

// fast version, but not supporting symmetric connections
vector<long> RandomConnectivity::IterationSpecificIndexes(Unit* postUnit, vector<int> nrPreRateUnits)
{
	m_useIterationSpecificIndexes = true;
	//if(m_symmetric == false)
	//{
		long maxIndexes = 0;
		for(int i=0;i<nrPreRateUnits.size();i++)
			maxIndexes+=nrPreRateUnits[i];

		int nrItems = (int)(maxIndexes*this->m_fracConnected);
		vector<long> out;

		// better scaling for sparse conns
		while(out.size()!=nrItems)
		{
			for(int i=0;i<nrItems;i++)
			{
				float r = rand()/(float(RAND_MAX)+1);
				long loc = (long)(((float)maxIndexes)*r); // float limit of population size (*100 .. / 100)

				out.push_back(loc);
			}

			sort(out.begin(),out.end());
			
			// remove duplicates
			out.erase(unique(out.begin(),out.end()),out.end());

			if(out.size()>nrItems)
				out.erase(out.begin()+nrItems,out.end());
		}

		sort(out.begin(),out.end());

		/*for(int i=0;i<nrItems;i++)
		{
		float r = rand()/(float(RAND_MAX)+1);
		long loc = (long)((maxIndexes-0)*r + 0);

		while(find(out.begin(),out.end(),loc) != out.end())
		{
		r = rand()/(float(RAND_MAX)+1);
		loc = (long)((maxIndexes-0)*r + 0);
		}

		out.push_back(loc);
		}*/

		return out;
	//}
	//else return vector<long>();
}

// slow
bool RandomConnectivity::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	// currently forced or non-forced self-connections
	if(preUnit->GetUnitId() == postUnit->GetUnitId())
	{
		if(this->m_allowSelfConnections == true)
			return true;
		else return false;
	}


	/*if(m_symmetric == true)
	{
		//unsigned int seed = preUnit->GetUnitId()+postUnit->GetUnitId();
		//unsigned int seed = (float)60000*(float(preUnit->GetUnitId())/(float(m_pre->GetNrUnits())))+(float)60000*(float(postUnit->GetUnitId())/(float(m_post->GetNrUnits())));

		//srand(seed);//preUnit->GetUnitId()*postUnit->GetUnitId()+preUnit->GetUnitId()+postUnit->GetUnitId()); // guarantee symmetry
		//srand(rand());
		float r = rand()/(float(RAND_MAX)+1);
		if(r<this->m_fracConnected)
			return true;
		else
			return false;
	}
	else*/ return false;
}

/*int RandomConnectivity::IterationAdvance()
{
	// go forward X steps based on a poisson process of mean 1/fraction connected

	/*	std::tr1::poisson_distribution<int> pd(1/this->m_fracConnected); // gives compiler error vsc++
	int val = 0; pd(val);
	return val;*/

	// approximation to poisson process
/*	double x,t;
	t=(double)rand()/(double)RAND_MAX;
	int step = -log(t)*1/this->m_fracConnected;
	return step;
}*/

/*bool RandomConnectivity::IterationPreConstraint(int currentIndex) // faster
{
	if((float)rand()/(float)RAND_MAX < m_fracConnected)
		return true;
	else
		return false;
}*/

/*bool RandomConnectivity::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	if((float)rand()/(float)RAND_MAX < m_fracConnected)
		return true;
	else
		return false;
}*/

bool OneToOneConnectivity::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	if(this->m_useIterationSpecificIndexes == true)
		return true;
	//return true;

	if(this->m_betweenHypercolumns == true)
	{
	/*	if(((RateUnit*)preUnit)->GetHypercolumns()[0]->GetUnitIdLocal() == ((RateUnit*)postUnit)->GetHypercolumns()[0]->GetUnitIdLocal())
			return true;
		else
			return false;*/
		
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

	return list; // needs to be sorted lowest to highest!
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

Connection::~Connection()
{
	for(int i=0;i<m_events.size();i++) // they can not be removed here since they are used in resetting
	{
		//m_events[i]->Dispose(); 
		delete m_events[i];
	}

	/*for(int i=0;i<m_listPosts.size();i++)
	{
		m_listPosts[i]->clear();
		delete m_listPosts[i];
	}*/

}

bool FullConnectivity::ExtraConstraints(Unit* preUnit, Unit* postUnit)
{
	if(m_allowSelfConnections == false)
	{
		if(preUnit->GetUnitId() == postUnit->GetUnitId())
			return false;
		else
			return true;
	}
	else
		return true;
}

void EmptyConnectivity::ExtraPostUnit(Unit* postUnit,Connection* conn) 
{ 
	conn->AddPostUnit(postUnit->GetUnitId());
}

// will be removed
/*void Connection::SetWeight(float weight, long preId, long postId)
{
	network()->SetWeight(weight, preId, postId);
}

void Connection::SetDelay(float delay, long preId, long postId)
{
	network()->SetDelay(delay, preId, postId);
}

float Connection::GetDelay(long preId, long postId)
{
	return network()->GetDelay(preId,postId);
}

float Connection::GetWeight(long preId, long postId)
{
	return network()->GetWeight(preId,postId);
}*/


UnitModifier* Connection::GetUnitModifier(int id)
{
	for(int i=0;i<m_unitPropertiesConnection.size();i++)
	{
		if(m_unitPropertiesConnection[i]->GetId() == id)
			return m_unitPropertiesConnection[i];
	}

	return NULL;
}

UnitModifier* Connection::GetUnitModifier(string name)
{
	for(int i=0;i<m_unitPropertiesConnection.size();i++)
	{
		if(m_unitPropertiesConnection[i]->GetName().compare(name) == 0)
			return m_unitPropertiesConnection[i];
	}

	return NULL;
}

void FullConnectivity::SetRandomWeights(vector<float> minVal, float maxVal,Network* net)
{
	m_setRandomWeights = true;
	m_weightsMinVal = minVal[0];
	m_weightsMaxVal = maxVal;
	m_weightsMinValParams = minVal;
	this->network(net);
	AddParameters();
	//this->network()->Parameters()->AddParameters(this,&FullConnectivity::SetRandomWeightsMin,minVal);
}

void FullConnectivity::SetRandomWeights(float minVal, vector<float> maxVal,Network* net)
{
	m_setRandomWeights = true;
	m_weightsMinVal = minVal;
	m_weightsMaxVal = maxVal[0];
	m_weightsMaxValParams = maxVal;
	this->network(net);
	AddParameters();
	//this->network()->Parameters()->AddParameters(this,&FullConnectivity::SetRandomWeightsMax,maxVal);
}

// Multiple parameters implementation
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