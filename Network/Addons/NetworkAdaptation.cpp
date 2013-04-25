#include "NetworkAdaptation.h"

void PopulationModifierAdaptation::Simulate()
{
	if(!IsOn())
		return;

	if(fabs(m_adampl) < EPS || fabs(m_adtaudt) < EPS)
		return;

	PopulationColumns* layer = (PopulationColumns*)m_population;
	vector<Unit*> minicolumns = layer->GetLocalRateUnits();

	if(m_adaptValues.size()!=minicolumns.size())
		m_adaptValues = vector<float>(minicolumns.size());

	vector<float> data = vector<float>(minicolumns.size());

	for(int i=0;i<minicolumns.size();i++)
	{
		data[i] = ((RateUnit*)minicolumns[i])->GetValue();
		
//		m_adaptValues[i] += (m_adampl*data[i] - m_adaptValues[i]) * m_adtaudt;

		if(m_type == this->StandardMinusBeta)
		{
			float beta = ((RateUnit*)minicolumns[i])->GetBeta();
			((RateUnit*)minicolumns[i])->SetValue(data[i]-m_adaptValues[i]-beta);
		}
		else 
		{
			((RateUnit*)minicolumns[i])->SetValue(data[i]-m_adaptValues[i]);
			((RateUnit*)minicolumns[i])->SetSubThresholdValue(data[i]-m_adaptValues[i]); // for debugging
		}
	}
}

void PopulationModifierAdaptation::Modify()
{
	if(!IsOn())
		return;

	PopulationColumns* layer = (PopulationColumns*)m_population;
	vector<Unit*> minicolumns = layer->GetLocalRateUnits();
	vector<float> data = vector<float>(minicolumns.size());

	for(int i=0;i<minicolumns.size();i++)
	{
		data[i] = ((RateUnit*)minicolumns[i])->GetValue();
		if(m_type == this->FeatureExtraction)
			m_adaptValues[i] += (m_adampl*data[i]) - (1-data[i])*m_adampl*m_adtaudt;//m_adaptValues[i] = (1-m_adtaudt)*m_adaptValues[i];
		else
		{
			if(!(fabs(data[i])<EPS))
				m_adaptValues[i] += ((m_adampl*1) - m_adaptValues[i]) * m_adtaudt;//(m_adampl*data[i]*fabs((((RateUnit*)minicolumns[i])->GetSubThresholdValue())) - m_adaptValues[i]) * m_adtaudt;
			else
				m_adaptValues[i] += ( -m_adaptValues[i]) * m_adtaudt*5;
		}
	}
}



void PopulationModifierTrace::Simulate()
{
	PopulationColumns* layer = (PopulationColumns*)m_population;
	vector<Unit*> minicolumns = layer->GetLocalRateUnits();

	if(m_lastValues.size() == 0)
		m_lastValues = vector<float>(minicolumns.size(),0.0);

	for(int i=0;i<minicolumns.size();i++)
	{
		minicolumns[i]->SetValue(minicolumns[i]->GetValue()+m_lastValues[i]*m_traceStrength);

		if(minicolumns[i]->GetValue()>1.0)
			minicolumns[i]->SetValue(1.0);

		m_lastValues[i] = minicolumns[i]->GetValue();
	}
}

void PopulationModifierTrace::Modify()
{
}