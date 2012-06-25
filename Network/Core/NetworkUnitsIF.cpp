#include "NetworkUnitsIF.h"
#include "NetworkUnits.h"
#include <iostream>
#include <cmath>
#include <deque>


UnitIF::UnitIF()
{
	//m_unitType = "minicolumn";
	Km = 0.0041;
	Vth = 60000;
	Vahp = 0;
	Vk = 0;
	Ik = 0;
}

void UnitIF::SimulateEventQueue()
{
//	if(m_eventsIncoming.size()>0)
//	{
		//vector<long> currentIncomingSpikes = m_eventsIncoming[0];
		vector<long> currentIncomingSpikes;// = m_eventsIncoming[0];

		//if(currentIncomingSpikes.size()>0)
		//{
			float totWeights = 0;

			//			m_eventsIncoming[0].clear();

			vector<Connection*> conns = this->GetPopulation()->GetIncomingConnections();
			vector<UnitModifier*> eus;
			vector<float> weights;

			// put connection vars in network instead of connection to avoid check

			for(int m=0;m<conns.size();m++)
			{
				//for(int j=0;j<1;j++)//conns.size();j++)
				//{
				if(conns[m]->IsOn())
				{

					// currently fix after eventsIncoming removed, used as in graded units, even though incoming non-spikes should not need to be checked
					// change into storing incoming spike-ids in bufferSpikes or something similar

					vector<long> preIds = conns[m]->GetPreIds(this->GetUnitId());

					if(preIds.size()>0)
					{
						for(int j=0;j<preIds.size();j++)
						{

							//for(int i=0;i<currentIncomingSpikes.size();i++)//m_eventsIncoming.size();i++)
							//{
							float incomingBufferData = m_network->GetPreValue(preIds[j]);//GetIncomingBufferData(preIds[j]);//it->first);

							if(incomingBufferData!=0)
							{
								//UnitModifier* eu = network()->GetUnitModifierIncoming(currentIncomingSpikes[i]);

								//m_value += m_eventsIncoming[i]->GetValue() * conns[j]->GetWeight(m_eventsIncoming[i]->GetFromUnitId(),this->GetUnitId());
								bool exists = true;

								if(conns.size()>1) // if more than one incoming connection we need to know which connection it is coming from
								{
									exists = true;//network()->ConnectionExists(eu->GetFromUnitId(),this->GetUnitId());//conns[j]->ConnectionExists(eu->GetFromUnitId(),this->GetUnitId());
								}

								if(exists)
								{
									//totWeights += conns[j]->GetWeight(eu->GetFromUnitId(),this->GetUnitId());
									float weight = network()->GetWeight(preIds[j],this->GetUnitId());//eu->GetFromUnitId(),this->GetUnitId());//conns[j]->GetWeight(eu->GetFromUnitId(),this->GetUnitId());
									//if(weight<0) // inhibitory
									//{
									Ik+=weight*incomingBufferData;//weight*eu->GetEventData()[0];
									//if(eu->GetEventTypeId()==2)
//										bool b=false;
									//	S_.y_[DG_IN] += -weight * V_.g0_in_;
									//}
									//else // excitatory
									//	;//S_.y_[DG_EX] += weight * V_.g0_ex_;

									//weights.push_back(weight);
									//eus.push_back(eu);
									//eu->SetValue(eu->GetValue() * weight);
								}
							}
						}
					}
				}
			}

			//if(S_.y_[DG_EX] > 300 * V_.g0_ex_)
			//	S_.y_[DG_EX] = 300 * V_.g0_ex_;
			//}

			// elsewhere now
			//m_eventsIncoming.erase(m_eventsIncoming.begin(),m_eventsIncoming.begin()+1);
			//	}

	// matlab code
	/*	% output: spikes
	% Vk1: updated potential
	% Vk: previous potential
	% Ik: input current

	Km = 0.0041; 

	% km=0.14;
	Vth = 60000;  

	Vahp = 0;             

	Vk1=Vk+Ik-Vk*Km;

	if Vk1 >= Vth                                      
	Vk1 = Vahp;
	output=1;
	else 
	output=0;
	end
	*/
	m_value = 0.0; // for Population recording

	float Ke = 0.27;
	Ik = Ik -Ik*Ke;

	float Vk1 = Vk+Ik-Vk*Km;

	if(IsRecording())
	{
		if(Vk1>=Vth)
			m_storedValues.push_back(Vth);
		else
			m_storedValues.push_back(Vk1);
	}

	if(Vk1 >= Vth)
	{
		// spike

		Vk1 = Vahp;

		// send spike
		m_isNewEvent = true;
		float t = 0.0; // not used atm
		m_eventsOutgoing.push_back(this->CreateEvent(t));

		if(this->GetPopulation()->IsRecording() == true)
		{
			m_value = 1;
		}
	}

	Vk = Vk1;

	if(IsRecording())
	{
		m_storedValues.insert(m_storedValues.begin(),this->network()->GetCurrentTimeStep());
		m_recordedValues.push_back(m_storedValues);
		m_storedValues.clear();
	}
}


UnitModifier* UnitIF::CreateEvent(float time)
{
	// m_hypercolumnId ?
	UnitModifierSpike* e = new UnitModifierSpike(m_unitId,m_hypercolumnId,time);
	return e;
}


#if GSL_AVAILABLE==1

#endif