#include "NetworkUnitModifier.h"




void GeometryUnit::SetPosition(float x, float y, float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
	if(this->m_unit->GetUnitType().compare("hypercolumn") == 0)
	{
		// position minicolumns automatically
		Hypercolumn* hc = (Hypercolumn*)this->m_unit;
		vector<RateUnit*> mcs = hc->GetRateUnits();

		for(int i=0;i<mcs.size();i++)
		{
			float xd = 0.5 - (float)rand()/(float)RAND_MAX;
			float yd = 0.5 - (float)rand()/(float)RAND_MAX;
			float zd = 0.5 - (float)rand()/(float)RAND_MAX;

			// replace to get from Projection instead
//			((GeometryUnit*)mcs[i]->GetUnitModifier(this->m_id))->SetPosition(x+xd,y+yd,z+zd); // change to have spread
		}
	}
}

TransferReTIDe::TransferReTIDe(Population* layer, float threshold, bool useDivNormalization, float maxValue, float tau)
{
	m_id = 8;
	m_threshold = threshold;//1;
	m_maxValue = maxValue; //-1
	m_tau = tau;//10;
	m_firstRun = true;
	m_useHashImpl = false;
	m_useDivNormalization = useDivNormalization;
	m_nrUnits = 0;

	if(m_useHashImpl == false)
	{
		vector<Unit*> units = layer->GetUnits(); // assumes ordered units
		if(units.size()>0)
		{
			m_lowestUnitLocalId = units[0]->GetUnitIdLocal();
			m_nrUnits = units.size(); // can change to make only local list of values
			m_xtList = vector<float>(m_nrUnits);
		}	
	}
}

void TransferReTIDe::Clear()
{
	m_xt.clear();
	m_xtList.clear();
	m_xtList = vector<float>(m_nrUnits);
}
// replace map by vector (local id gives rel index) to increase speed
void TransferReTIDe::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	float xt = 0;
	
//	TimingStart("TransferReTIDe_summing");

	for(int i=0;i<events.size();i++)
	{
		xt += events[i]->GetEventData()[0] * weights[i];
	}

//	TimingStop("TransferReTIDe_summing");

//	TimingStart("TransferReTIDe_calc");
	// hash impl
	long id;
	float mx;

	float totActivity = 0;

	if(m_useDivNormalization)
	{
		totActivity = 0;
		// no need to calc every time
		#if USE_UNORDERED_MAP == 1
			unordered_map<long,float>::iterator itr;
		#else
			map<long,float>::iterator itr;
		#endif

		for(itr = m_xt.begin(); itr != m_xt.end(); ++itr){
			totActivity += (*itr).second;
		}

		//if(totActivity == 0)
		//	totActivity = 1;
	}

	if(m_useHashImpl == true)
	{
		id = unit->GetUnitId();
		mx =  m_xt[id];
		if(m_useDivNormalization)
			m_xt[id] = mx + ((-mx+xt)/(1+totActivity))/m_tau;
		else
			m_xt[id] = mx + (-mx+xt)/m_tau; 
		
		mx = m_xt[id];
	}
	else
	{
		id = unit->GetUnitIdLocal() - m_lowestUnitLocalId; // assumes ordered local unit ids
		if(m_useDivNormalization)
			m_xtList[id] = m_xtList[id] + ((-m_xtList[id]+xt)/(1+totActivity))/m_tau; 
		else
			m_xtList[id] = m_xtList[id] + (-m_xtList[id]+xt)/m_tau; 

		mx = m_xtList[id];
	}
	
	if(mx>m_threshold)
	{
		m_value = mx-m_threshold;
		/*if(m_value>10)
		{
			bool b=false;
			cout<<"|";cout.flush();
		}*/
		if(m_maxValue > 0)
		{
			if(m_value > m_maxValue)
				m_value = m_maxValue;
		}
	}
	else m_value = 0;

//	TimingStop("TransferReTIDe_calc");

	// index-based impl
	/*
	if(m_firstRun == true)
	{
		vector<Unit*> localUnits = unit->Population()->GetLocalUnits();
		m_localFirstId = localUnits[1]->GetUnitIdLocal(); // assumes we use Populationcolumns and first unit is a hypercolumn (not used)
		m_firstRun = false;
	}

	int index = unit->GetUnitIdLocal() - m_localFirstId;
	while(index>=m_xt.size())
	{
		m_xt.push_back(0.0);
	}

	float mx =  m_xt[index];
	m_xt[index] = mx + (-mx+xt)/m_tau; 
	mx =  m_xt[index];

	if(mx>m_threshold)
		m_value = mx-m_threshold;
	else m_value = 0;
	*/
}

void TransferReTIDe::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	float xt = 0;

	for(int i=0;i<values->size();i++)
	{
		xt += values->at(i) * weights->at(i);
	}

	long id;
	float mx;

	float totActivity = 0;

	if(m_useHashImpl == true)
		id = unit->GetUnitId();
	else
		id = unit->GetUnitIdLocal() - m_lowestUnitLocalId; // assumes ordered local unit ids

	/*if(m_useDivNormalization)
	{
		totActivity = 0;
		// no need to calc every time
		#if USE_UNORDERED_MAP == 1
			unordered_map<long,float>::iterator itr;
		#else
			map<long,float>::iterator itr;
		#endif

		if(m_useHashImpl == true)
		{
			if(m_xt[id]>m_threshold)
			{
				for(itr = m_xt.begin(); itr != m_xt.end(); ++itr)
					totActivity += (*itr).second-m_threshold;
			}
		}
		else
		{
			if(m_xtList[id]>m_threshold)
			{
				for(int i=0;i<m_xtList.size();i++)
					totActivity+=m_xtList[i]-m_threshold;
			}
		}

		//if(totActivity == 0)
		//	totActivity = 1;
	}*/


	if(m_useHashImpl == true)
	{
		mx = m_subThresholdValue;//m_xt[id];
		
		/*if(m_useDivNormalization)
		{
			m_xt[id] = mx + ((-mx+xt)/(1+totActivity))/m_tau;
		}
		else*/
			m_xt[id] = mx + (-mx+xt)/m_tau;

		//m_xt[id] = mx + (-mx+xt)/m_tau; 
		mx = m_xt[id];
	}
	else
	{
		/*if(m_useDivNormalization)
		{
			m_xtList[id] = m_xtList[id] + (-m_xtList[id]+xt)/m_tau; 
			
			if(m_xtList[id]>m_threshold)
			{
			//	m_xtList[id] = m_xtList[id]+(100-totActivity)/m_tau;
				//m_xtList[id] = (m_xtList[id] + ((m_xtList[id])/(1+totActivity)))/2;
			}
		}
		else*/
		
		mx = m_subThresholdValue;
		mx = mx + (-mx+xt)/m_tau;
		//m_xtList[id] = m_xtList[id] + (-m_xtList[id]+xt)/m_tau; 
	
		//mx = m_xtList[id];
	}

	m_subThresholdValue = mx;
	
	if(mx>m_threshold)
	{
		m_value = mx-m_threshold;

		if(m_maxValue > 0)
		{
			if(m_value > m_maxValue)
				m_value = m_maxValue;
		}
	}
	else m_value = 0;
}


// replace map by vector (local id gives rel index) to increase speed
void TransferThreshold::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	float xt = unit->GetValue();//0;
	/*for(int i=0;i<events.size();i++)
	{
		xt += events[i]->GetEventData()[0] * weights[i];
	}*/

	// hash impl
	long id = unit->GetUnitId();
	float mx =  m_xt[id];
	m_xt[unit->GetUnitId()] = mx + (-mx+xt)/m_tau; 
	mx = m_xt[unit->GetUnitId()];

	if(mx>m_threshold)
		m_value = m_threshold;//mx-m_threshold;
	else m_value = 0;
}

void TransferThreshold::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	float xt = unit->GetValue();//0;
	/*for(int i=0;i<events.size();i++)
	{
		xt += events[i]->GetEventData()[0] * weights[i];
	}*/

	// hash impl
	long id = unit->GetUnitId();
	float mx =  m_xt[id];
	m_xt[unit->GetUnitId()] = mx + (-mx+xt)/m_tau; 
	mx = m_xt[unit->GetUnitId()];

	if(mx>m_threshold)
		m_value = m_threshold;//mx-m_threshold;
	else m_value = 0;
}

void TransferBcpnnOnline::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	
	m_value = 0.0;
	map<long,float> hashHcIdSum;

	// initial hash values
	for(int i=0;i<events.size();i++)
		hashHcIdSum[((UnitModifierGraded*)events[i])->GetFromHypercolumnId()] = 0.0;

	for(int i=0;i<events.size();i++)
	{
		//float weight = m_unit->Population()->GetWeights GetWeight( events[i]->GetFromUnitId(), m_unit->GetUnitId() );
		hashHcIdSum[((UnitModifierGraded*)events[i])->GetFromHypercolumnId()] += events[i]->GetEventData()[0] * weights[i];
	}

	for(map<long,float>::const_iterator it = hashHcIdSum.begin(); it != hashHcIdSum.end(); ++it)
	{
		if(it->second>0.00001)
			m_value += log(it->second);
		else if(it->second<-0.00001) // put this case here for the inhibitory recurrent Projections (defined by impact weights -1.0 etc)
			m_value -= log(-it->second);
	}

	m_value += m_beta;

	// (currently used only for top-down context dependence in olfactory cortex->bulb)
	if(m_maxValue>-0.09)
		if(m_value>m_maxValue)
			m_value = m_maxValue;

	if(m_useThreshold)
		if(m_value < m_threshold)
			m_value = 0;
		else
			m_value = 1.0;
}

// need hypercolumn ids, so implemented below
void TransferBcpnnOnline::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	bool b = false;
}

void TransferBcpnnOnline::SimulateV2HypercolumnIds(vector<float> values, vector<float> weights, vector<long> hypercolumnIds, Unit* unit)
{
	m_value = 0.0;
	map<long,float> hashHcIdSum;

	// initial hash values
	for(int i=0;i<weights.size();i++)
		hashHcIdSum[hypercolumnIds[i]] = 0.0;

	for(int i=0;i<weights.size();i++)
	{
		//float weight = m_unit->Population()->GetWeights GetWeight( events[i]->GetFromUnitId(), m_unit->GetUnitId() );
		hashHcIdSum[hypercolumnIds[i]] += values[i] * weights[i];
	}

	for(map<long,float>::const_iterator it = hashHcIdSum.begin(); it != hashHcIdSum.end(); ++it)
	{
		if(it->second>0.00001)
			m_value += log(it->second);
		else if(it->second<-0.00001) // put this case here for the inhibitory recurrent Projections (defined by impact weights -1.0 etc)
			m_value -= log(-it->second);
	}

	m_value += m_beta;

	// (currently used only for top-down context dependence in olfactory cortex->bulb)
	if(m_maxValue>-0.09)
		if(m_value>m_maxValue)
			m_value = m_maxValue;

	if(m_useThreshold)
		if(m_value < m_threshold)
			m_value = 0;
		else
			m_value = 1.0;
}

void TransferLinear::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	//m_value = 0.0;

	if(unit==NULL)
		cerr << "TransferLinear::Simulate. Null Unit." << endl;

	if(weights.size()!=events.size())
		cerr << "TransferLinear::Simulate. Unequal vector sizes: events and weights." << endl;

	float y = 0;
	for(int i=0;i<events.size();i++)
	{
		if(events[i]==NULL)
			cerr << "TransferLinear::Simulate. Null Event: "<< i << endl;
		y += events[i]->GetEventData()[0] * weights[i];
	}

	if(y>0)
		bool b=false;
	if(m_useThreshold == true) // debug test, not changeable parameters atm.
	{
		float r = 1.0/4.0 * (float)rand()/(float)RAND_MAX;// - 0.5;//0.0;//(rand()/(float(RAND_MAX)+1));
		m_value+=m_eta*(1+r)*y;//m_eta*1.0/ ( 1.0 + exp(-y));;

		//m_value+=s[unit->GetUnitId()];
		// s[unit->GetUnitId()] += m_eta3*(m_alpha-m_value); // f: -gamma (+s instead of -f)
		float threshold = 0.99;

		if(m_value<0.0) m_value = 0.0;

		if(m_value + s[unit->GetUnitId()] > threshold)
		{
			m_value = 1.0;
			long id = unit->GetUnitId();
			s[id] -= m_eta3;//*m_alpha;//(m_alpha-m_value);
			
			for(map<long,float>::const_iterator it = s.begin(); it != s.end(); ++it)
			{
				if(it->first!=id)
					s[it->first] += m_eta3*m_alpha;
			}

			if(s[id]<0.0)
				s[id] = 0.0;
		}
	}
	else m_value = y;

	if(m_useSign == true)
	{
		/*
		if(m_value>0.5)
			m_value = 1;
		else
			m_value = 0.0;*/
		if(m_value>1)
			m_value = 1;
	}
	if(m_value!=m_value)
		cerr << "TransferFunctions::TransferLinear::Simulate(). Yields NaN value." << endl;
}


void TransferLinear::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	//m_value = 0.0;

	if(unit==NULL)
		cerr << "TransferLinear::Simulate. Null Unit." << endl;

	if(weights->size()!=values->size())
		cerr << "TransferLinear::Simulate. Unequal vector sizes: events and weights." << endl;

	float y = 0;
	for(int i=0;i<values->size();i++)
	{
		y += values->at(i) * weights->at(i);
	}

	if(y>0)
		bool b=false;
	if(m_useThreshold == true) // debug test, not changeable parameters atm.
	{
		float r = 1.0/4.0 * (float)rand()/(float)RAND_MAX;// - 0.5;//0.0;//(rand()/(float(RAND_MAX)+1));
		m_value+=m_eta*(1+r)*y;//m_eta*1.0/ ( 1.0 + exp(-y));;

		//m_value+=s[unit->GetUnitId()];
		// s[unit->GetUnitId()] += m_eta3*(m_alpha-m_value); // f: -gamma (+s instead of -f)
		float threshold = 0.99;

		if(m_value<0.0) m_value = 0.0;

		if(m_value + s[unit->GetUnitId()] > threshold)
		{
			m_value = 1.0;
			long id = unit->GetUnitId();
			s[id] -= m_eta3;//*m_alpha;//(m_alpha-m_value);
			
			for(map<long,float>::const_iterator it = s.begin(); it != s.end(); ++it)
			{
				if(it->first!=id)
					s[it->first] += m_eta3*m_alpha;
			}

			if(s[id]<0.0)
				s[id] = 0.0;
		}
	}
	else m_value = y;

	if(m_useSign == true)
	{
		/*
		if(m_value>0.5)
			m_value = 1;
		else
			m_value = 0.0;*/
		if(m_value>1)
			m_value = 1;
	}
	if(m_value!=m_value)
		cerr << "TransferFunctions::TransferLinear::Simulate(). Yields NaN value." << endl;
}


void TransferFoldiak::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	m_value = 0.0;

	for(int i=0;i<events.size();i++)
	{
		m_value += events[i]->GetEventData()[0] * weights[i]; // both Wx and Uy / ff and rec?
	}

	// equilbrium...
	m_value += s[unit->GetUnitId()];

	m_value = 1.0/ ( 1.0 + exp(-m_beta*m_value));
	s[unit->GetUnitId()] += m_eta3*(m_alpha-m_value); // f: -gamma (+s instead of -f)
}

void TransferFoldiak::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	m_value = 0.0;

	for(int i=0;i<values->size();i++)
	{
		m_value += values->at(i) * weights->at(i); // both Wx and Uy / ff and rec?
	}

	// equilbrium...
	m_value += s[unit->GetUnitId()];

	m_value = 1.0/ ( 1.0 + exp(-m_beta*m_value));
	s[unit->GetUnitId()] += m_eta3*(m_alpha-m_value); // f: -gamma (+s instead of -f)
}

float TransferCSL::RRValue(const vector<float>& x1, const vector<float>& x2)
{
  float val = 0;
  
  for(int i=0;i<x1.size();i++) {
    val+=(x1[i]-x2[i])*(x1[i]-x2[i]);//pow(x1[i]-x2[i],2);
  }

  return sqrt(val);
}

void TransferCSL::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	m_value = 0.0;

	vector<float> values(events.size());

	for(int i=0;i<events.size();i++)
	{
		values[i] = events[i]->GetEventData()[0];
//		m_value += events[i]->GetValue() * weights[i];
	}

	m_value = 1-RRValue(values,weights);
}

void TransferCSL::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	m_value = 1-RRValue(*values,*weights);
}

void TransferTriesch::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	// from Butko and Triesch (2007)
	
	float h = 0.0;

	for(int i=0;i<events.size();i++)
	{
		h += events[i]->GetEventData()[0] * weights[i];
	}

	float y = 1.0/(1.0+exp(-a[unit->GetUnitId()]*h-b[unit->GetUnitId()]));
	float y2 = y*y;//pow(y,2.0f);

	if(a[unit->GetUnitId()] == 0)
		a[unit->GetUnitId()] = 0.01;
	if(b[unit->GetUnitId()] == 0)
		b[unit->GetUnitId()] = 0.01;
	
	if(m_isThresholded == true)
	{
		if(m_value == 0) // last was spike
		{
			m_value += m_tau*y;

			if(m_thresholdLastSpike == true)
			{
				m_thresholdLastSpike = false;
				y=1;//m_thresholdLastY;//1;
				h = m_thresholdLastH;
				y2 = y*y;//pow(y,2.0f);
			}
			else
			{
				y = 0.0;//m_thresholdLastY;//0.5;
				h = m_thresholdLastH;
				y2 = y*y;//pow(y,2.0f);
			}

			a[unit->GetUnitId()] = a[unit->GetUnitId()] + eta_ip*(1.0/a[unit->GetUnitId()] + h - (2.0+1.0/mu)*h*y + 1.0/mu * h * y2);
			b[unit->GetUnitId()] = b[unit->GetUnitId()] + eta_ip*(1.0 - (2.0 + 1.0/mu)*y + 1.0/mu * y2);
		}
		else
		{
			m_value += m_tau*y;
		}

		float threshold = 0.5;
		
		if(m_value>threshold)
		{
			m_thresholdLastSpike = true;
			m_value = 1.0;
		}

		m_thresholdLastH = h;
		m_thresholdLastY = y;
	}
	else
	{
		m_value = y;
		a[unit->GetUnitId()] = a[unit->GetUnitId()] + eta_ip*(1.0/a[unit->GetUnitId()] + h - (2.0+1.0/mu)*h*y + 1.0/mu * h * y2);
		b[unit->GetUnitId()] = b[unit->GetUnitId()] + eta_ip*(1.0 - (2.0 + 1.0/mu)*y + 1.0/mu * y2);
	}

//	if(y<0.1)
//		y = 0.0;


}


void TransferTriesch::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	// from Butko and Triesch (2007)
	
	float h = 0.0;

	for(int i=0;i<values->size();i++)
	{
		h += values->at(i) * weights->at(i);
	}

	float y = 1.0/(1.0+exp(-a[unit->GetUnitId()]*h-b[unit->GetUnitId()]));
	float y2 = y*y;//pow(y,2.0f);

	if(a[unit->GetUnitId()] == 0)
		a[unit->GetUnitId()] = 0.01;
	if(b[unit->GetUnitId()] == 0)
		b[unit->GetUnitId()] = 0.01;
	
	if(m_isThresholded == true)
	{
		if(m_value == 0) // last was spike
		{
			m_value += m_tau*y;

			if(m_thresholdLastSpike == true)
			{
				m_thresholdLastSpike = false;
				y=1;//m_thresholdLastY;//1;
				h = m_thresholdLastH;
				y2 = y*y;//pow(y,2.0f);
			}
			else
			{
				y = 0.0;//m_thresholdLastY;//0.5;
				h = m_thresholdLastH;
				y2 = y*y;//pow(y,2.0f);
			}

			a[unit->GetUnitId()] = a[unit->GetUnitId()] + eta_ip*(1.0/a[unit->GetUnitId()] + h - (2.0+1.0/mu)*h*y + 1.0/mu * h * y2);
			b[unit->GetUnitId()] = b[unit->GetUnitId()] + eta_ip*(1.0 - (2.0 + 1.0/mu)*y + 1.0/mu * y2);
		}
		else
		{
			m_value += m_tau*y;
		}

		float threshold = 0.5;
		
		if(m_value>threshold)
		{
			m_thresholdLastSpike = true;
			m_value = 1.0;
		}

		m_thresholdLastH = h;
		m_thresholdLastY = y;
	}
	else
	{
		m_value = y;
		a[unit->GetUnitId()] = a[unit->GetUnitId()] + eta_ip*(1.0/a[unit->GetUnitId()] + h - (2.0+1.0/mu)*h*y + 1.0/mu * h * y2);
		b[unit->GetUnitId()] = b[unit->GetUnitId()] + eta_ip*(1.0 - (2.0 + 1.0/mu)*y + 1.0/mu * y2);
	}

//	if(y<0.1)
//		y = 0.0;


}

void TransferTriesch::Clear()
{
	a.clear();
	b.clear();
}


void TransferSigmoid::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	m_value = 0.0;

	for(int i=0;i<events.size();i++)
	{
		m_value += events[i]->GetEventData()[0] * weights[i];
	}
	
	// tanh sigmoid
	// [-1,1]

	if(m_value > 0.1)
		bool b = false;
	
	m_value = 2.0/(1.0+exp(-2.0*m_value)) - 1.0;


	if(m_useThreshold == true)
	{
		if(m_value >0.6)
			m_value = 1.0;
		else
			m_value = 0.0;
	}
}


void TransferSigmoid::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	m_value = 0.0;

	for(int i=0;i<values->size();i++)
	{
		m_value += values->at(i) * weights->at(i);
	}
	
	// tanh sigmoid
	// [-1,1]

	if(m_value > 0.1)
		bool b = false;
	
	m_value = 2.0/(1.0+exp(-2.0*m_value)) - 1.0;


	if(m_useThreshold == true)
	{
		if(m_value >0.6)
			m_value = 1.0;
		else
			m_value = 0.0;
	}
}

void TransferSpikesToRate::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	// can change this entire fcn to use the time stamps in the incoming events instead (calculate from last 5 spikes etc), makes more sense for instantaneous firing rate

	m_value = 0.0;

	for(int i=0;i<events.size();i++)
	{
		//m_value += events[i]->GetEventData()[0] * weights[i];
		// or just
		m_value += events[i]->GetEventData()[0];
		// if weights not set (could make optional)
	}

	(m_spikesHistory[unit->GetUnitId()]).push_back(m_value);
	if(m_spikesHistory[unit->GetUnitId()].size()>m_timeStepsWindow)
	{
		m_spikesHistory[unit->GetUnitId()].erase(m_spikesHistory[unit->GetUnitId()].begin(),m_spikesHistory[unit->GetUnitId()].begin()+1);
	}

	m_value = 0;
	
	// calculate rate
	for(int i=0;i<m_spikesHistory[unit->GetUnitId()].size();i++)
	{
		m_value+=m_spikesHistory[unit->GetUnitId()][i];
	}

	m_value = m_value/(float)(m_spikesHistory[unit->GetUnitId()].size());

	// relate to time step
	//m_value = m_value/m_network->GetTimeStep();
}


void TransferSpikesToRate::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	// can change this entire fcn to use the time stamps in the incoming events instead (calculate from last 5 spikes etc), makes more sense for instantaneous firing rate

	m_value = 0.0;

	for(int i=0;i<values->size();i++)
	{
		//m_value += events[i]->GetEventData()[0] * weights[i];
		// or just
		m_value += values->at(i);
		// if weights not set (could make optional)
	}

	(m_spikesHistory[unit->GetUnitId()]).push_back(m_value);
	if(m_spikesHistory[unit->GetUnitId()].size()>m_timeStepsWindow)
	{
		m_spikesHistory[unit->GetUnitId()].erase(m_spikesHistory[unit->GetUnitId()].begin(),m_spikesHistory[unit->GetUnitId()].begin()+1);
	}

	m_value = 0;
	
	// calculate rate
	for(int i=0;i<m_spikesHistory[unit->GetUnitId()].size();i++)
	{
		m_value+=m_spikesHistory[unit->GetUnitId()][i];
	}

	m_value = m_value/(float)(m_spikesHistory[unit->GetUnitId()].size());

	// relate to time step
	//m_value = m_value/m_network->GetTimeStep();
}


void TransferPositive::Simulate(vector<UnitModifier*> events, vector<float> weights, Unit* unit)
{
	if(unit->GetValue() < 0)
	{
		m_value = 0;
		//unit->SetValue(0.0);
	}
}

void TransferPositive::SimulateV2(vector<float>* values, vector<float>* weights, Unit* unit)
{
	if(unit->GetValue() < 0)
	{
		m_value = 0;
		//unit->SetValue(0.0);
	}
}

void UnitModifierGraded::Send(int processId, long unitId)
{
	int tagAct = 0;
	MPI_Send(&unitId,1,MPI_LONG,processId,tagAct,NETWORK_COMM_WORLD); // id of sender (unit, not node)
	tagAct = 1;
	float value = this->GetEventData()[0];
	MPI_Send(&value,1,MPI_FLOAT,processId,tagAct,NETWORK_COMM_WORLD);
	tagAct = 2;
	MPI_Send(&m_fromHypercolumnId,1,MPI_INT,processId,tagAct,NETWORK_COMM_WORLD); // id of sender (hypercolumn)
}

void UnitModifierGraded::Receive(int processId)
{
	long unitId;
	float value;
	int tagAct = 0;
	MPI_Status status;

	MPI_Recv(&unitId,1,MPI_LONG,processId,tagAct,NETWORK_COMM_WORLD,&status);
	tagAct = 1;
	MPI_Recv(&value,1,MPI_FLOAT,processId,tagAct,NETWORK_COMM_WORLD,&status);

	tagAct = 2;
	MPI_Recv(&m_fromHypercolumnId,1,MPI_INT,processId,tagAct,NETWORK_COMM_WORLD,&status);

	//m_value = value;
}