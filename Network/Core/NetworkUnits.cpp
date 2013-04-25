#include <math.h>
#include "NetworkUnits.h"
#include <deque>

/// <summary>	Adds a unit modifier (e.g. transfer function) on a unit
/// 			This will be used across all incoming projections - can put on a specific projection instead to make it specific </summary>
///
/// <param name="p">	[in,out] If non-null, the UnitModifier* to process. </param>

void Unit::AddUnitModifier(UnitModifier* p)
{
	m_unitProperties.insert(m_unitProperties.begin(),p);
	p->SetUnit(this);
}

/// <summary>	Loop around all units in population gets here. </summary>

void Unit::Simulate()
{
	// work on the event queues for added properties on units and projections (e.g. transfer functions etc)
	SimulateEventQueue();

	// Simulate any additional tasks for unit (typically not used but can be overridden)
	SimulateMisc();

	this->ClearEventsIncoming();
}


/// <summary>	Simulates and handles all incoming data (from other units most of the time) to unit
/// currently:
///  - have some separate cases for certain transfer functions (such as bcpnn which has log around hypercolumns), this may be moved
///  - not data-driven as loops would need to change their order in that case (could get better performance if changed)
///  - has been optimized, hence the reserves and all different cases </summary>

void RateUnit::SimulateEventQueue()
{
	int index=0;

	bool doNoUpdating = false;

	if(m_noUpdatingCurrentTimeStep == true)
	{
		if(!(fabs(m_value)<EPS))
			m_isNewEvent = true;
		else
			m_isNewEvent = false;

		doNoUpdating = true;
		m_noUpdatingCurrentTimeStep = false;
	}
	else if(this->GetPopulation()->KeepValues() == true)
	{
		// do nothing
	}
	else
	{
		m_value = 0.0;
	}

	if(doNoUpdating == false)
	{
		vector<Projection*> conns = this->GetPopulation()->GetIncomingProjections(); // optionally change to pointer

		vector<float> allValues;
		vector<float> allWeights;
		vector<long> allHypercolumnIds;
		vector<UnitModifier*>* layerUnitProperties = m_population->GetUnitPropertiesLayer();

		vector<long> preIds;
		vector<pair<long, float> >* preIdsValues;
		vector<UnitModifier*>* unitProperties;

		// may be more allocation and weight fetching if Projection level check here, could move it down (but then check if preId exist in Projection is needed)
		for(int m=0;m<conns.size();m++)
		{
			preIds.clear();

			if(conns[m]->IsOn())
			{
				unitProperties = conns[m]->GetUnitPropertiesProjection();

				// could move so this is not in inner loop, special case for CSL-type of classes
				bool isTransferCSL = false;
				for(int n=0;n<unitProperties->size();n++)
				{
					if((*unitProperties)[n]->GetId() == 4)
						isTransferCSL = true;
				}

				if(isTransferCSL == true)
					preIds = conns[m]->GetPreIds(this->GetUnitId());
				else
#if USE_HASHED_ACTIVE_COMMUNICATION == 1
					//preIdsValues = conns[m]->GetPreIdsActive(this->GetUnitId()); // optimization, turned off
					preIdsValues = conns[m]->GetPreIdsActiveLocal(this->GetUnitIdLocal()); // optimization
#else
					preIds = conns[m]->GetPreIds(this->GetUnitId());
#endif

				vector<float> weights;
				vector<float> values;
				vector<long> hypercolumnIds; // only filled if hypercolumn ids tracked
				float incomingBufferData;
				long preIdIter;
				float w;
				long hId = -1;

				bool cont = false;

#if USE_HASHED_ACTIVE_COMMUNICATION == 1
				if(isTransferCSL)
				{
					if(preIds.size()>0)
					{
						// reserve allocates but does not construct
						weights.reserve(preIds.size()); 
						values.reserve(preIds.size());
						allWeights.reserve(preIds.size());
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
					weights.reserve(preIdsValues->size());
					values.reserve(preIdsValues->size());
					allWeights.reserve(preIdsValues->size());
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
					weights.reserve(preIds.size());
					values.reserve(preIds.size());
					allWeights.reserve(preIds.size());
					allValues.reserve(preIds.size());

					if(this->network()->IsTrackingHypercolumnIds())
					{
						hypercolumnIds.reserve(preIds.size());
						allHypercolumnIds.reserve(preIds.size());
					}

					cont = true;
				}
#endif

				if(cont == true)
				{
					int nrIterations;

					if(preIds.size()>0)
						nrIterations = preIds.size();
					else
						nrIterations = preIdsValues->size();

					for(int j=0;j<nrIterations;j++)
					{
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
						if(!(fabs(incomingBufferData)<EPS) || isTransferCSL)
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

							w = m_network->GetWeight(preIdIter,this->GetUnitId());

							weights.push_back(w);//it->first,this->GetUnitId()));
							values.push_back(incomingBufferData);
							if(m_unitProperties.size()>0 || layerUnitProperties->size()>0) // if global unit properties exist
							{ 
								allWeights.push_back(w);
								allValues.push_back(incomingBufferData);
							}
						}
					}

					// calls e.g. transfer functions
					SimulateUnitPropertiesV2(unitProperties,&values,&weights,&hypercolumnIds);
				}
			}
		}

		// if properties have been added to populations and not units
		if(layerUnitProperties->size()>0)
		{
			SimulateUnitPropertiesV2(layerUnitProperties,&allValues,&allWeights,&allHypercolumnIds);
		}

		// if no global unit properties exist / deleted afterwards otherwise	
		if(m_unitProperties.size()>0)
		{	
			SimulateUnitPropertiesV2(layerUnitProperties,&allValues,&allWeights,&allHypercolumnIds);
		}

		if(!(fabs(m_value)<EPS))
			m_isNewEvent = true;
		else
			m_isNewEvent = false;
	}
}

// unit properties and transfer functions
// special cases for certain transfer functions which do need global information at the moment
// these are separated for performance.
void RateUnit::SimulateUnitProperties(vector<UnitModifier*> unitProperties, vector<UnitModifier*> eus, vector<float> weights)
{
	bool bcpnnRun = false;

	for (int j=0;j<unitProperties.size();j++)
	{
		if(unitProperties[j]->IsOn())
		{
			if(unitProperties[j]->GetId() == 2) // transfer bcpnn
			{
				if(bcpnnRun == false)
				{
					((TransferBcpnnOnline*)unitProperties[j])->SetValue(0.0);
					if(!(fabs(m_beta)<EPS)) // unnecessary (?)
						((TransferBcpnnOnline*)unitProperties[j])->SetBeta(m_beta);
					unitProperties[j]->Simulate(eus,weights, this);
					m_value += ((TransferBcpnnOnline*)unitProperties[j])->GetValue();
					bcpnnRun = true;
				}
			}
			else // e.g. linear transfer functions (+ nonlinearity) (CSL, foldiak, triesch)
			{
				unitProperties[j]->SetValue(m_value);
				unitProperties[j]->Simulate(eus,weights, this);
				m_value += unitProperties[j]->GetValue();
			}

			if(m_value != 1.0)
				bool b = false;
		}
	}
}

// unit properties and transfer functions
// special cases for certain transfer functions which do need global information at the moment
// these are separated for performance.
void RateUnit::SimulateUnitPropertiesV2(vector<UnitModifier*>* unitProperties, vector<float>* values, vector<float>* weights, vector<long>* hypercolumnIds)
{
	bool bcpnnRun = false;
	m_subThreshDrive = 0;

	for (int j=0;j<unitProperties->size();j++)
	{
		if((*unitProperties)[j]->IsOn())
		{
			if((*unitProperties)[j]->GetId() == 2) // special case for bcpnn
			{
				if(bcpnnRun == false)
				{
					((TransferBcpnnOnline*)(*unitProperties)[j])->SetValue(0.0);
					if(!(fabs(m_beta)<EPS))
						((TransferBcpnnOnline*)(*unitProperties)[j])->SetBeta(m_beta);
					(*unitProperties)[j]->SimulateV2HypercolumnIds(*values,*weights,*hypercolumnIds,this);
					m_value += ((TransferBcpnnOnline*)(*unitProperties)[j])->GetValue();
					bcpnnRun = true;
				}
			}
			else // e.g. linear transfer functions (+ nonlinearity) (CSL, foldiak, triesch)
			{
				(*unitProperties)[j]->SetSubThresholdValue(m_subThreshValue);
				(*unitProperties)[j]->SetValue(m_value);
				(*unitProperties)[j]->SimulateV2(values,weights,this);
				m_value += (*unitProperties)[j]->GetValue();
				
				for(int i=0;i<values->size();i++)
					if(weights->at(i)>0)
						m_subThreshDrive += values->at(i)*weights->at(i);

				m_subThreshValue = (*unitProperties)[j]->GetSubThresholdValue();
			}

		}
	}
}

// used for some plasticity models (bcpnn)
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
}

// TODO: Check if used any longer
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
		int processId = this->GetNodeId();
		float value;
		int tagAct = processId*10 + 0;
		MPI_Status status;
		MPI_Recv(&unitIdLocalHypercolumn,1,MPI_INT,processId,tagAct,NETWORK_COMM_WORLD,&status);

		unitIdLocalHypercolumn=unitIdLocalHypercolumn-1;
		tagAct = processId*10 + 1;
		MPI_Recv(&value,1,MPI_FLOAT,processId,tagAct,NETWORK_COMM_WORLD,&status);

		m_value = value;
		m_isLocal = true;
	}

	return m_value;
}

// TODO: since createEvents not used any longer, check dependencies and remove
UnitModifier* RateUnit::CreateEvent(float value)
{
	m_value = value;

	UnitModifierGraded* e = new UnitModifierGraded(m_unitId,m_hypercolumnId,m_value);
	return e;
}


UnitModifier* RateUnit::CreateEvent()
{
	UnitModifierGraded* e = new UnitModifierGraded(m_unitId,m_hypercolumnId,m_value);
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


// TODO: move
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

/// <summary>	Implements trace capability for RateUnits. </summary>

void RateUnit::SimulateMisc()
{
	// update trace
	// (timestep is now 1/timestep)
	if(m_useTrace == true)
	{
		m_value = m_value+m_traceValueLast*m_traceStrength;
	}
}

/// <summary>	Returns what by default should be recorded (if recording on) of a RateUnit. </summary>
///
/// <returns>	Activity value. </returns>

vector<vector<float> > RateUnit::GetValuesToRecord()
{
	vector<vector<float> > f(1);
	vector<float> f2(1);
	f2[0] = m_value;
	f[0] = f2;

	return f;
}

/// <summary>	Returns what by default should be recorded (if recording on) of a RateUnit. </summary>
///
/// <returns>	Values of all units in hypercolumn. </returns>

vector<vector<float> > Hypercolumn::GetValuesToRecord()
{
	vector<vector<float> > f(1);
	f[0] = m_values;
	return f;
}

/// <summary>	Creates what should be distributed from Unit by global communication. </summary>

void Unit::AddEventOutgoing()
{
	m_eventsOutgoing.push_back(this->CreateEvent(m_value));
}

/// <summary>	Adds incoming data at the correct delay.
/// 			TODO: Create test to check if this is working as it should after optimizations. </summary>
///
/// <param name="eventIndex">	Zero-based index of the event. </param>

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

	if(m_eventsIncoming.size()<delay+1) // will change in size during runtime
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
}