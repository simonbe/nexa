#include "NetworkUnitsIF.h"
#include "NetworkUnits.h"
#include <iostream>
#include <cmath>
#include <deque>


UnitIF::UnitIF(float km, float vth, float vahp, float vk, float ik)
{
	Km = km;
	Vth = vth;
	Vahp = vahp;
	Vk = vk;
	Ik = ik;
}

void UnitIF::SimulateEventQueue()
{
		vector<long> currentIncomingSpikes;

			float totWeights = 0;

			vector<Projection*> conns = this->GetPopulation()->GetIncomingProjections();
			vector<UnitModifier*> eus;
			vector<float> weights;

			for(int m=0;m<conns.size();m++)
			{

				if(conns[m]->IsOn())
				{
					// currently fix after eventsIncoming removed, used as in rate units, even though incoming non-spikes should not need to be checked
					// TODO: change into storing incoming spike-ids in bufferSpikes or something similar

					vector<long> preIds = conns[m]->GetPreIds(this->GetUnitId());

					if(preIds.size()>0)
					{
						for(int j=0;j<preIds.size();j++)
						{

							float incomingBufferData = m_network->GetPreValue(preIds[j]);

							if(!(fabs(incomingBufferData)<EPS))
							{
								bool exists = true;

								if(conns.size()>1) // if more than one incoming Projection we need to know which Projection it is coming from
								{
									exists = true;
								}

								if(exists)
								{
									float weight = network()->GetWeight(preIds[j],this->GetUnitId());

									Ik+=weight*incomingBufferData;
								}
							}
						}
					}
				}
			}

	// original matlab code
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
	UnitModifierSpike* e = new UnitModifierSpike(m_unitId,m_hypercolumnId,time);
	return e;
}


#if GSL_AVAILABLE==1
// implementation excluded from this release
#endif