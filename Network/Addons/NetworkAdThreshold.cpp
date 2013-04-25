#include "NetworkAdThreshold.h"

void AdaptiveThreshold::Modify()
{
	PopulationColumns* layer = (PopulationColumns*)m_population;
	vector<Unit*> units = layer->GetUnits();

	float h = 0.0;

	if(m_thresholdType == AdaptiveThreshold::Average) // exponential moving average
	{
		if(m_average.size()==0)
		{
			m_average = vector<float>(units.size(),0.01);
			m_firstRun = false;
		}

		for(int j=0;j<units.size();j++)
		{
			RateUnit* mc = (RateUnit*)units[j];
			h = mc->GetValue();

			m_average[j] = m_average[j] + m_alpha*(((1-m_lambda0)*h + m_lambda0)-m_average[j]);
		}
	}
}

void AdaptiveThreshold::Simulate()
{
	float h = 0.0;

	PopulationColumns* layer = (PopulationColumns*)m_population;

	vector<int> localHypercolumns = layer->GetLocalHypercolumnIndexes();
	vector<Unit*> units = layer->GetUnits();

	if(m_thresholdType == AdaptiveThreshold::Triesch)
	{
		for(int i=0;i<units.size();i++)
		{
			RateUnit* unit = (RateUnit*)units[i];
			h = unit->GetValue();

			float y = 1.0/(1.0+exp(-a[unit->GetUnitId()]*h-b[unit->GetUnitId()]));
			float y2 = y*y;//pow(y,2.0f);

			if(fabs(a[unit->GetUnitId()]) < EPS)
				a[unit->GetUnitId()] = 0.01;
			if(fabs(b[unit->GetUnitId()]) < EPS)
				b[unit->GetUnitId()] = 0.01;

			if(m_isThresholded == true)
			{
				/*			if(m_value == 0) // last was spike
				{
				m_value += m_tau*y;

				if(m_thresholdLastSpike == true)
				{
				m_thresholdLastSpike = false;
				y=1;//m_thresholdLastY;//1;
				h = m_thresholdLastH;
				y2 = pow(y,2.0f);
				}
				else
				{
				y = 0.0;//m_thresholdLastY;//0.5;
				h = m_thresholdLastH;
				y2 = pow(y,2.0f);
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
				m_thresholdLastY = y;*/
			}
			else
			{
				unit->SetValue(y);
				a[unit->GetUnitId()] = a[unit->GetUnitId()] + eta_ip*(1.0/a[unit->GetUnitId()] + h - (2.0+1.0/mu)*h*y + 1.0/mu * h * y2);
				b[unit->GetUnitId()] = b[unit->GetUnitId()] + eta_ip*(1.0 - (2.0 + 1.0/mu)*y + 1.0/mu * y2);
			}
		}
	}
	else if(m_thresholdType == AdaptiveThreshold::BCPNNBeta)
	{
		// not correct impl - can get neg log if done after transf func.

		if(m_firstRun == true)
		{
			m_Aj = vector<float>(units.size(),m_alpha);
			m_firstRun = false;
		}
		
		for(int j=0;j<units.size();j++)
		{
			RateUnit* mc = (RateUnit*)units[j];
			h = mc->GetValue();

			m_Aj[j] = m_Aj[j] + m_alpha*(( (1-m_lambda0)*h+ m_lambda0)-m_Aj[j]);
			float beta = m_impactBeta * log(m_Aj[j]);
			h+=beta;

			mc->SetValue(h);
		}
	}
	else if(m_thresholdType == AdaptiveThreshold::Average) // exponential moving average
	{
		if(m_averageGain.size() == 0)
		{
			m_averageGain = vector<float>(units.size());
		}
		else
		{
			for(int j=0;j<units.size();j++)
			{
				RateUnit* mc = (RateUnit*)units[j];
				h = mc->GetValue();
				//m_averageGain[j]+=m_impactAverage*m_average[j];
				//-0.01*m_averageGain[j]
				m_averageGain[j] +=  + m_impactAverage*(m_averageAim - m_average[j]);//(1-m_impactAverage)*m_averageGain[j]+10*m_impactAverage*(m_averageAim - m_average[j]);
				h += m_averageGain[j];//m_impactAverage*m_average[j];

				mc->SetValue(h);
			}
		}
	}
}