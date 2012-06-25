#include <math.h>
#include "NetworkUnits.h"
#include <deque>

// will be used across all incoming connections - put in a specific connection instead to make it specific
void Unit::AddUnitModifier(UnitModifier* p)
{
	m_unitProperties.insert(m_unitProperties.begin(),p);
	//m_unitProperties.push_back(p);
	p->SetUnit(this);
}

// not used atm.
void Unit::SendReceive()
{
	// send created events (unit)
	/*UnitModifier* e = this->CreateEvent();
	this->SetEventOutgoing(e);

	m_maxNode = -1;
	for(int i=0;i<m_posts.size();i++)
	{
		if(m_posts[i]->IsLocal())
			m_posts[i]->AddEventIncoming(e);
		else
		{
			int nodeId = m_posts[i]->GetNodeId();
			if(AlreadySentEventToNode(nodeId) == false)
			{
				e->Send(m_posts[i]->GetNodeId(),this->GetUnitId());
				m_nrEventsSent++;
			}
		}
			//this->SendEvent(&e,,m_posts[i]->GetUnitIdLocal()); // MPI specific
	}

	m_sentToNode.clear();

	cout<<GetUnitId()<<": events sent = "<< m_nrEventsSent<<", ";

	// retrieve events (unit)
	ReceiveEvents();*/
}

// not used atm.
void Unit::SendReceiveNext()
{
	// send created events (unit)
	/*UnitModifier* e = this->CreateEvent();
	this->SetEventOutgoing(e);

	m_maxNode = -1;
	for(int i=0;i<m_posts.size();i++)
	{
		if(m_posts[i]->IsLocal())
			m_posts[i]->AddEventIncoming(e);
		else
		{
			int nodeId = m_posts[i]->GetNodeId();
			if(AlreadySentEventToNode(nodeId) == false)
			{
				e->Send(m_posts[i]->GetNodeId(),this->GetUnitId());
				m_nrEventsSent++;
			}
		}
			//this->SendEvent(&e,,m_posts[i]->GetUnitIdLocal()); // MPI specific
	}

	m_sentToNode.clear();
	cout<<GetUnitId()<<": events sent = "<< m_nrEventsSent<<", ";

	// retrieve events (unit)
	ReceiveEvents();*/
}

// Loop around all units in population ends up here, overridable
void Unit::Simulate()
{
	// work on the event queues for unit and connections
	SimulateEventQueue();

	/*int nr = m_eventsIncoming.size();
	for(int i=0;i<nr;i++)
	{
		if(m_eventsIncoming[i]!=NULL)
		{
			delete m_eventsIncoming[i];
			m_eventsIncoming[i] = NULL;
		}
	}

	m_eventsIncoming.clear();*/

	SimulateMisc();

	this->ClearEventsIncoming();

	//m_eventsIncoming.clear();
}

/*void RateUnit::AddEvent(UnitModifier* e)
{
	// do nothing
}*/

// not used atm.
bool Unit::AlreadySentEventToNode(int nodeId)
{
	bool alreadySent = false;

	if(nodeId == m_maxNode)
	{
		alreadySent = true;
	}
	else if(nodeId < m_maxNode)
	{
		for(int i=0;i<m_sentToNode.size();i++)
			if(nodeId == m_sentToNode[i])
			{
				alreadySent = true;
				break;
			}
	}
	else
	{
		m_maxNode = nodeId;
	}

	if(alreadySent == false)
	{
		m_sentToNode.push_back(nodeId);
	}

	return alreadySent;
}

// not used atm.
void Unit::ReceiveEvents()
{
	for(int i=0;i<m_pres.size();i++)
	{
		if(m_pres[i]->IsLocal() == true)
		{
			//this->AddEventIncoming(m_pres[i]->GetEventOutgoing());
		}
		else
		{
			UnitModifier* e = this->CreateEvent();
			int nodeId = m_pres[i]->GetNodeId();
			e->Receive(m_pres[i]->GetNodeId());
			m_pres[i]->SetEventOutgoing(e);
			
			// ok to remove? only place local set
			//m_pres[i]->SetLocal(true); // fetch locally from now on

			m_nrEventsReceived++;
		}
	}

	cout<<" received = "<< m_nrEventsReceived<<"\n";

	/*for(int i=0;i<m_pres.size();i++)
	{
		if(m_pres[i]->IsLocal() == false)
		{
			MPI_Status status;
			long unitIdLocal;
			float value;

			int tagAct = 0;
			MPI_Recv(&unitIdLocal,1,MPI_LONG,m_pres[i]->GetNodeId(),tagAct,NETWORK_COMM_WORLD,&status);
			tagAct = 1;
			MPI_Recv(&value,1,MPI_FLOAT,m_pres[i]->GetNodeId(),tagAct,NETWORK_COMM_WORLD,&status);

			UnitModifier* eg = new UnitModifierGraded;
			eg->SetValue(value);
			eg->SetPreLocalIndex(unitIdLocal);
		}
		else // check if event has been sent
		{

		}

		// if 
	}*/
}

// Simulates and handles all incoming data (from other units most of the time) to unit
// currently:
// - have some separate cases for certain transfer functions (such as bcpnn which has log around hypercolumns)
// - not data-driven as loops would need to change their order in that case
// - main optimize place (check timing results)
void RateUnit::SimulateEventQueue()
{
//	this->Population()->TimingStart(m_name);
	// let connections modify incoming
	int index=0;
	/*for(int i=0;i<m_preConnections.size();i++)
	{
		if(m_preConnections[i]->Pre()->IsLocal() == true)
		{
			Connection* c = m_preConnections[i];
			c->SimulateEvent(m_eventsIncoming[index]);
			index++;
		}
	}*/

	/*for(int i=0;i<m_eventsIncoming.size();i++)
	{
		Connection* c = m_hashIdConnection[m_eventsIncoming[i]->GetFromUnitId()];
		c->SimulateEvent(m_eventsIncoming[index]);
	}*/

	// currently creating and adding new events here, this should be changed for performance
	// change to iterate over synapses instead
	
	//int nrProperties = m_unitProperties.size();
	
	//if(nrProperties>1)
	//	nrProperties = 1; // debug setting

	bool doNoUpdating = false;

	if(m_noUpdatingCurrentTimeStep == true)
	{
		if(m_value!=0.0)
			m_isNewEvent = true;
		else
			m_isNewEvent = false;

		doNoUpdating = true;
		m_noUpdatingCurrentTimeStep = false;
	}
	else if(this->GetPopulation()->KeepValues() == true)
	{
		bool b = false;
		if(m_value > 0.0)
			b = true;
		// values have been manually set but unit (m_value) should be updated as well
	}
	else
	{
		m_value = 0.0;

		/*if(m_useThreshold == false)
			m_value = 0.0;
		else
		{
			if(m_value == 1.0)	// right after send
				//m_value = 0.0;
				doNoUpdating = true; // no updating until reset called
		}*/
	}

	if(doNoUpdating == false)
	{
		vector<Connection*> conns = this->GetPopulation()->GetIncomingConnections(); // change to pointer
		
		/*if(nrProperties > conns.size()) // currently assuming max one unit property per connection
			if(m_network->MPIGetNodeId()==0)
			cout<<"Error: More transfer functions used than currently possible.";*/
		
//		vector<UnitModifier*> allEus;
		vector<float> allValues;
		vector<float> allWeights;
		vector<long> allHypercolumnIds;
		vector<UnitModifier*>* layerUnitProperties = m_population->GetUnitPropertiesLayer();

		vector<long> preIds;
		vector<pair<long, float> >* preIdsValues;
		vector<UnitModifier*>* unitProperties;

		// may be more allocation and weight fetching if connection level check here, could move it down (but then check if preId exist in connection is needed)
		for(int m=0;m<conns.size();m++)
		{
			if(conns[m]->IsOn())
			{
				unitProperties = conns[m]->GetUnitPropertiesConnection();

				// move so this is not in inner loop
				bool isTransferCSL = false;
				for(int n=0;n<unitProperties->size();n++)
				{
					if((*unitProperties)[n]->GetId() == 4)
						isTransferCSL = true;
				}

				//map<long, Network::SynapseStandard>* preSynapses = m_network->GetPreSynapses(this->GetUnitId());
				//this->Population()->TimingStart("GetPreIds");

				if(isTransferCSL == true)
					preIds = conns[m]->GetPreIds(this->GetUnitId());
				else
#if USE_HASHED_ACTIVE_COMMUNICATION == 1
					//preIdsValues = conns[m]->GetPreIdsActive(this->GetUnitId()); // optimization
					preIdsValues = conns[m]->GetPreIdsActiveLocal(this->GetUnitIdLocal()); // optimization
#else
					preIds = conns[m]->GetPreIds(this->GetUnitId());
#endif

				//this->Population()->TimingStop("GetPreIds");
				//map<long, Network::SynapseStandard>::iterator it;
				// iterate

	//			vector<UnitModifier*> eus;
				vector<float> weights;
				vector<float> values;
				vector<long> hypercolumnIds; // only filled if hypercolumn ids tracked
				float incomingBufferData;
				long preIdIter;
				float w;
				long hId = -1;

				bool cont = false;

//#if USE_UNORDERED_MAP == 0
				//				map<long,float>::iterator iter;
				//#endif
#if USE_HASHED_ACTIVE_COMMUNICATION == 1

				if(isTransferCSL)
				{
					if(preIds.size()>0)
					{
						weights.reserve(preIds.size()); // allocates but does not construct (fastest?)
						values.reserve(preIds.size());
						allWeights.reserve(preIds.size()); // allocates but does not construct (fastest?)
						allValues.reserve(preIds.size());

						if(this->network()->IsTrackingHypercolumnIds())
						{
							hypercolumnIds.reserve(preIds.size());
							allHypercolumnIds.reserve(preIds.size());
						}

						cont = true;
					}
				}
				else if(preIdsValues->size()>0)
				{
					weights.reserve(preIdsValues->size()); // allocates but does not construct (fastest?)
					values.reserve(preIdsValues->size());
					allWeights.reserve(preIdsValues->size()); // allocates but does not construct (fastest?)
					allValues.reserve(preIdsValues->size());

					if(this->network()->IsTrackingHypercolumnIds())
					{
						hypercolumnIds.reserve(preIdsValues->size());
						allHypercolumnIds.reserve(preIdsValues->size());
					}

					cont = true;
				}
#else
				if(preIds.size()>0)
				{
					weights.reserve(preIds.size()); // allocates but does not construct (fastest?)
					values.reserve(preIds.size());
					allWeights.reserve(preIds.size()); // allocates but does not construct (fastest?)
					allValues.reserve(preIds.size());

					if(this->network()->IsTrackingHypercolumnIds())
					{
						hypercolumnIds.reserve(preIds.size());
						allHypercolumnIds.reserve(preIds.size());
					}

					cont = true;
				}
#endif
				//this->Population()->TimingStart("GetBufferHashData");

				if(cont == true)
				{
					int nrIter;

					if(preIds.size()>0)
						nrIter = preIds.size();
					else
						nrIter = preIdsValues->size();

					for(int j=0;j<nrIter;j++)
						/*#if USE_HASHED_ACTIVE_COMMUNICATION == 1
						for(int j=0;j<preIdsValues->size();j++)
						#else
						for(int j=0;j<preIds.size();j++)
						#endif*/
					{
						//vector<float>* incomingBufferData = m_network->GetIncomingBufferData(it->first);

						//this->Population()->TimingStart("GetPreValue");
#if USE_HASHED_ACTIVE_COMMUNICATION == 1

						if(isTransferCSL)
						{
							preIdIter = preIds[j];
							incomingBufferData = m_network->GetPreValue(preIdIter);
						}
						else
						{
							incomingBufferData = preIdsValues->at(j).second;
							preIdIter = preIdsValues->at(j).first;
						}
#else
						preIdIter = preIds[j];
						incomingBufferData = m_network->GetPreValue(preIdIter);//GetIncomingBufferData(preIds[j]);//it->first);
#endif
						//this->Population()->TimingStop("GetPreValue");
						//if((*incomingBufferData).size()>0)

						//this->Population()->TimingStart("GetWeights");
						if(incomingBufferData!=0 || isTransferCSL)
						{
							if(this->network()->IsTrackingHypercolumnIds())
							{
								hId = m_network->GetIncomingBufferHypercolumnIds(preIdIter);//it->first);
								hypercolumnIds.push_back(hId);
								if(m_unitProperties.size()>0 || layerUnitProperties->size()>0) // if global unit properties exist
								{
									allHypercolumnIds.push_back(hId);
								}
							}

//							UnitModifier* eu = new UnitModifierGraded(preIdIter,hId,incomingBufferData);//it->first,hId,incomingBufferData);//(*incomingBufferData)[0]); // assumes no delays
//							eus.push_back(eu);
							
							w = m_network->GetWeight(preIdIter,this->GetUnitId());

							// speed up (!)
							weights.push_back(w);//it->first,this->GetUnitId()));
							values.push_back(incomingBufferData);
							if(m_unitProperties.size()>0 || layerUnitProperties->size()>0) // if global unit properties exist
							{ 
								//allEus.push_back(eu);
								allWeights.push_back(w);
								allValues.push_back(incomingBufferData);
							}
						}
						//this->Population()->TimingStop("GetWeights");
					}
					//this->Population()->TimingStop("GetBufferHashData");

					/*for(int j=0;j<1;j++)//conns.size();j++)
					{
					if(m_eventsIncoming.size()>0)
					for(int i=0;i<m_eventsIncoming[0].size();i++)//m_eventsIncoming.size();i++)
					{
					UnitModifier* eu = network()->GetUnitModifierIncoming(m_eventsIncoming[0][i]);

					//m_value += m_eventsIncoming[i]->GetValue() * conns[j]->GetWeight(m_eventsIncoming[i]->GetFromUnitId(),this->GetUnitId());
					bool exists = true;

					if(conns.size()>1) // if more than one incoming connection we need to know which connection it is coming from
					{
					exists = true;//network()->ConnectionExists(eu->GetFromUnitId(),this->GetUnitId());//conns[j]->ConnectionExists(eu->GetFromUnitId(),this->GetUnitId());
					}

					if(exists)
					{
					float weight = network()->GetWeight(eu->GetFromUnitId(),this->GetUnitId());//conns[j]->GetWeight(eu->GetFromUnitId(),this->GetUnitId());
					weights.push_back(weight);
					eus.push_back(eu);
					//eu->SetValue(eu->GetValue() * weight);
					}
					else
					{
					float test = network()->GetWeight(eu->GetFromUnitId(),this->GetUnitId());
					bool b=false;
					}
					}
					}
					*/
					//for(int i=0;i<nrProperties;i++) // will do it twice atm for transfer bcpnn if rec+ff
					//{

					//int i = m;

					// unit properties on a connection level
					//SimulateUnitProperties(unitProperties,eus,weights);
					SimulateUnitPropertiesV2(unitProperties,&values,&weights,&hypercolumnIds);
					//}

					//m_beta = 0; // reset, not used

					// change to not use allocated vector !
					/*if(m_unitProperties.size()==0 && layerUnitProperties.size() == 0) // if no global unit properties exist / deleted afterwards otherwise	
						for(int i=0;i<eus.size();i++)
							delete eus[i];*/
				}
			}
		}

		if(layerUnitProperties->size()>0)
		{
			//SimulateUnitProperties(layerUnitProperties,allEus,allWeights);
			SimulateUnitPropertiesV2(layerUnitProperties,&allValues,&allWeights,&allHypercolumnIds);

			/*if(m_unitProperties.size() == 0)
			{
				for(int i=0;i<allEus.size();i++)
					delete allEus[i];
			}*/
		}

		if(m_unitProperties.size()>0) // if no global unit properties exist / deleted afterwards otherwise	
		{	
			SimulateUnitPropertiesV2(layerUnitProperties,&allValues,&allWeights,&allHypercolumnIds);
			/*SimulateUnitProperties(m_unitProperties,allEus,allWeights);
			for(int i=0;i<allEus.size();i++)
				delete allEus[i];*/
		}

		if(m_value!=0.0)
			m_isNewEvent = true;
		else
			m_isNewEvent = false;
	}

	if(m_useThreshold == true)
	{
		float tau = 40;
		//m_value+=m_beta;
		m_subThreshValue = m_subThreshValue + (-m_subThreshValue + m_value)/tau;//m_tau; 

		float threshold = 1;
		if(m_subThreshValue > threshold)//m_threshold)
		{
			m_value = m_subThreshValue-threshold;//threshold;//mx-m_threshold;
			//m_subThreshValue = threshold;
			m_isNewEvent = true;
		}
		else
		{
			m_value = 0;
			m_isNewEvent = false;
		}

		/*
		// above threshold dynamics determined by transfer fcn and (wta) eventlayer
		if(m_value == 1.0)
			m_isNewEvent = true;
		else
			m_isNewEvent = false;
		*/
	}
	else m_subThreshValue = m_value;

	if(m_isNewEvent == true)
	{
		//m_eventsOutgoing.push_back(this->CreateEvent(m_value));
	}

//	this->Population()->TimingStop(m_name);
}

// unit properties and transfer functions
// special cases for certain transfer functions
void RateUnit::SimulateUnitProperties(vector<UnitModifier*> unitProperties, vector<UnitModifier*> eus, vector<float> weights)
{
	//this->Population()->TimingStart("SimulateUnitProperties");
	bool bcpnnRun = false;

	for (int j=0;j<unitProperties.size();j++)
	{
		if(unitProperties[j]->IsOn())
		{
			if(unitProperties[j]->GetId() == 2)//m_unitProperties[i]->GetId() == 2) // transfer bcpnn
			{
				if(bcpnnRun == false)
				{
					((TransferBcpnnOnline*)unitProperties[j])->SetValue(0.0);//((TransferBcpnnOnline*)m_unitProperties[i])->SetValue(0.0);
					if(m_beta != 0.0) // fix ?
						((TransferBcpnnOnline*)unitProperties[j])->SetBeta(m_beta);//((TransferBcpnnOnline*)m_unitProperties[i])->SetBeta(m_beta);
					unitProperties[j]->Simulate(eus,weights, this);//m_unitProperties[i]->Simulate(eus,weights, this);
					m_value += ((TransferBcpnnOnline*)unitProperties[j])->GetValue();//m_value += ((TransferBcpnnOnline*)m_unitProperties[i])->GetValue();
					bcpnnRun = true;
				}
			}
			else// if(m_unitProperties[i]->GetId() == 1) // transfer linear (+ nonlinearity) (CSL, foldiak, triesch)
			{
				unitProperties[j]->SetValue(m_value);//m_unitProperties[i]->SetValue(m_value);
				unitProperties[j]->Simulate(eus,weights, this);//m_unitProperties[i]->Simulate(eus,weights, this);
				m_value += unitProperties[j]->GetValue();//m_value += m_unitProperties[i]->GetValue();//((TransferBcpnnOnline*)m_unitProperties[i])->GetValue();

				// fix (!) - should unit properties instead be able to either add or set value (?)
				//if(this->Population()->KeepValues() == true)
				//	m_value = unitProperties[j]->GetValue();
			}

			if(m_value != 1.0)
				bool b = false;
		}
	}

	//this->Population()->TimingStop("SimulateUnitProperties");
}

void RateUnit::SimulateUnitPropertiesV2(vector<UnitModifier*>* unitProperties, vector<float>* values, vector<float>* weights, vector<long>* hypercolumnIds)
{
	//this->Population()->TimingStart("SimulateUnitProperties");
	bool bcpnnRun = false;

	for (int j=0;j<unitProperties->size();j++)
	{
		if((*unitProperties)[j]->IsOn())
		{
			if((*unitProperties)[j]->GetId() == 2)//m_unitProperties[i]->GetId() == 2) // transfer bcpnn
			{
				if(bcpnnRun == false)
				{
					((TransferBcpnnOnline*)(*unitProperties)[j])->SetValue(0.0);//((TransferBcpnnOnline*)m_unitProperties[i])->SetValue(0.0);
					if(m_beta != 0.0) // fix ?
						((TransferBcpnnOnline*)(*unitProperties)[j])->SetBeta(m_beta);//((TransferBcpnnOnline*)m_unitProperties[i])->SetBeta(m_beta);
					(*unitProperties)[j]->SimulateV2HypercolumnIds(*values,*weights,*hypercolumnIds,this);
					//unitProperties[j]->Simulate(eus,weights, this);//m_unitProperties[i]->Simulate(eus,weights, this);
					m_value += ((TransferBcpnnOnline*)(*unitProperties)[j])->GetValue();//m_value += ((TransferBcpnnOnline*)m_unitProperties[i])->GetValue();
					bcpnnRun = true;
				}
			}
			else// if(m_unitProperties[i]->GetId() == 1) // transfer linear (+ nonlinearity) (CSL, foldiak, triesch)
			{
				(*unitProperties)[j]->SetValue(m_value);//m_unitProperties[i]->SetValue(m_value);
				(*unitProperties)[j]->SimulateV2(values,weights,this);
				//unitProperties[j]->Simulate(eus,weights, this);//m_unitProperties[i]->Simulate(eus,weights, this);
				m_value += (*unitProperties)[j]->GetValue();//m_value += m_unitProperties[i]->GetValue();//((TransferBcpnnOnline*)m_unitProperties[i])->GetValue();

				// fix (!) - should unit properties instead be able to either add or set value (?)
				//if(this->Population()->KeepValues() == true)
				//	m_value = unitProperties[j]->GetValue();
			}

			if(m_value != 1.0)
				bool b = false;
		}
	}

	//this->Population()->TimingStop("SimulateUnitProperties");
}

// used for some models
void RateUnit::AddGains()
{
	m_value += m_inhibBeta;
	m_inhibBeta = 0;
}

void Hypercolumn::SimulateEventQueue()
{
	// determine if it should be silent or not. If so, set all minicolumn values to 0.
	if(m_useSilent == true)
	{
		float maxVal = -1e9;
		float data;
		for(int i=0;i<this->GetRateUnits().size();i++)
		{
			data = this->GetRateUnits()[i]->GetValue();
			if(data>maxVal)
				maxVal = data;
		}

		if(maxVal>this->m_silentThreshold)
			m_isSilent = false;
		else
		{
			m_isSilent = true;
			for(int i=0;i<this->GetRateUnits().size();i++)
			{
				this->GetRateUnits()[i]->SetValue(0.0);
			}
		}
	}

	// work on the modifying event queues for the connections
	/*for(int i=0;i<m_preConnections.size();i++)
	{
		Connection* c = m_preConnections[i];
		c->ModifyConnection();
	}*/
}

// 
float RateUnit::GetValueFromGroup(vector<int> nodeIndexes)
{
	int thisNodeId = this->GetPopulation()->network()->MPIGetNodeId();

	if(this->IsLocal() == true)
	{
		for(int i=0;i<nodeIndexes.size();i++)
		{
			if(nodeIndexes[i] != thisNodeId)
			{
				int tagAct = thisNodeId*10 + 0;
				int t = m_localIdHypercolumn+1;
				MPI_Send(&t,1,MPI_INT,nodeIndexes[i],tagAct,NETWORK_COMM_WORLD); // id of sender (unit, not node)
				tagAct = thisNodeId*10 + 1;
				float value = this->GetValue()*10;
				MPI_Send(&value,1,MPI_FLOAT,nodeIndexes[i],tagAct,NETWORK_COMM_WORLD);
			}
		}
	}
	else
	{
		int unitIdLocalHypercolumn = 0;
		int nodeId = this->GetNodeId();
		float value;
		int tagAct = nodeId*10 + 0;
		MPI_Status status;
		MPI_Recv(&unitIdLocalHypercolumn,1,MPI_INT,nodeId,tagAct,NETWORK_COMM_WORLD,&status);

		unitIdLocalHypercolumn=unitIdLocalHypercolumn-1;
		tagAct = nodeId*10 + 1;
		MPI_Recv(&value,1,MPI_FLOAT,nodeId,tagAct,NETWORK_COMM_WORLD,&status);

		m_value = value;
		m_isLocal = true;
	}

	return m_value;
}

// not used atm?
UnitModifier* RateUnit::CreateEvent(float value)
{
	m_value = value;

	UnitModifierGraded* e = new UnitModifierGraded(m_unitId,m_hypercolumnId,m_value);
	//m_isNewEvent = true;
	//e->SetValue(m_value);
	//e->SetFromUnitId(m_unitId);
	//e->SetFromHypercolumnId(m_hypercolumnId);

	return e;
}


UnitModifier* RateUnit::CreateEvent()
{
	UnitModifierGraded* e = new UnitModifierGraded(m_unitId,m_hypercolumnId,m_value);
	/*e->SetValue(m_value);
	e->SetFromUnitId(m_unitId);
	e->SetFromHypercolumnId(m_hypercolumnId);*/

	return e;
}

UnitModifier* Hypercolumn::CreateEvent(float value)
{
	UnitModifierGraded* e = new UnitModifierGraded(this->GetUnitId(),this->GetUnitId(),value);
	return e;
}

// not used atm.
UnitModifier* Hypercolumn::CreateEvent()
{
	UnitModifierGraded* e = NULL;// = new UnitModifierGraded(this->GetUnitId(),this->GetUnitId(),value);

	return e;
}

/*void RateUnit::SetTraceTimeSteps(int traceTimeSteps)
{
	m_traceTimeSteps = traceTimeSteps;
	if(m_traceTimeSteps>0)
		m_useTrace = true;
	else
		m_useTrace = false;
}*/

// move
int poisson(double c) { // c is the intensity
	int x = 0;
	double t = 0.0;
	while(true)
	{
		float r = (float)rand()/(float)RAND_MAX;
	    t -= log(r)/c;//log(r01())/c;
	    if (t > 1.0)
			return x;
	    x++;
	}
}

// overridden to implement trace capability for minicolumns
void RateUnit::SimulateMisc()
{
	// update adaptation
	// impl as an eventlayer atm
/*	if(m_useAdaptation == true)
	{
		m_adaptationValue += (m_adaptationAmpl*m_value - m_adaptationValue) * m_adaptationTau;
		m_value -= m_adaptationValue;
	}
	*/

	// update trace
	// (timestep is now 1/timestep)
	if(m_useTrace == true)
	{
		m_value = m_value+m_traceValueLast*m_traceStrength;
		/*for(map<int,float>::const_iterator it = m_trace.begin(); it != m_trace.end(); ++it)
		{
			m_trace[it->first] += (m_value - it->second) * 1.0/(float)it->first;
			//trgtrc[j] += (trgact[j] - trgtrc[j]) * zjtaudt;
			//trgctrc[j] += (trgcact[j] - trgctrc[j]) * zjtaudt;
		}*/
	}
}

/*if (m_noiseampl<=0)
			m_dsup[j] += (gain*sup[j] - dsup[j])* taumdt;
		else {
			dsup[j] += (gain*(sup[j] +
				m_noiseampl*(poisson(noiseintens) - noiseintens)) -
				dsup[j])* taumdt;
		}*/

// Recording return function
vector<vector<float> > RateUnit::GetValuesToRecord()
{
	vector<vector<float> > f(1);
	vector<float> f2(1);
	f2[0] = m_value;
	f[0] = f2;

	return f;
}

// Recording return function
vector<vector<float> > Hypercolumn::GetValuesToRecord()
{
	vector<vector<float> > f(1);
	f[0] = m_values;
	return f;
}

/*vector<vector<float> > UnitIF::GetValuesToRecord()
{
	vector<vector<float> > f2(2);
	vector<float> f(4);
	f[0] = S_.y0_;
	f[1] = S_.y1_;
	f[2] = S_.y2_;
	f[3] = S_.y3_; // membrane potential RELATIVE TO RESTING POTENTIAL.

	f2[0] = f;

	return f2;
}*/

void Unit::AddEventOutgoing()
{
	m_eventsOutgoing.push_back(this->CreateEvent(m_value));
}

// in changes
void Unit::AddEventIncoming(long eventIndex)
{
	// eventIndex assumed to be preId here
	int delay = 0;

	// could specify global delay flag if it is on, instead of checking in hash even if not used - large speedup

	if(m_network->IsUsingDelays())
	{
		delay = network()->GetDelay(eventIndex,this->GetUnitId());//conns[j]->GetDelay(eventIndex,this->GetUnitId());

		float timeSteps = 1.0/this->network()->GetTimeResolution();
		delay = delay*(int)timeSteps; // [ms]
	}

/*	if(m_eventsIncoming.size()<delay+1) // will change in size during runtime
		for(int i=m_eventsIncoming.size();i<delay+1;i++)
		{
			//vector<long> v(0);
			vector<long> v(0);
			m_eventsIncoming.push_back(v);
		}

	// speed up by using a double-ended que, deque, instead according to http://www.codeproject.com/KB/stl/vector_vs_deque.aspx
	// deques should grow more efficiently
	// (this function currently called in communication loop)

	m_eventsIncoming[delay].push_back(eventIndex);
	*/
}

/*UnitModifier* Unit::GetUnitModifier(int id)
{
	for(int i=0;i<m_unitProperties.size();i++)
	{
		if(m_unitProperties[i]->GetId() == id)
			return m_unitProperties[i];
	}

	return NULL;
}

UnitModifier* Unit::GetUnitModifier(string name)
{
	for(int i=0;i<m_unitProperties.size();i++)
	{
		if(m_unitProperties[i]->GetName().compare(name) == 0)
			return m_unitProperties[i];
	}

	return NULL;
}*/