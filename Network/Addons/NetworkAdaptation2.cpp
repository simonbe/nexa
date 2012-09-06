#include "NetworkAdaptation2.h"
// #include "NetworkBCPNN.h"

void PopulationModifierAdaptation2::Initm_Aj(float aj) // m_aj was set to 0.01 in BCPNN initialization
	// unlike the class NetworkAdaption, here m_Aj should be initialized after construction, 
	// although a condition in Simulates checks for initialization and initalizes later if necessary
{
	int n=((PopulationColumns*)m_population)->GetLocalRateUnits().size();
	m_Aj = vector<float>(n); 
	for(int i=0;i<n;i++)
		m_Aj[i]=aj;
}

void PopulationModifierAdaptation2::Simulate()
{

	PopulationColumns* layer = (PopulationColumns*)m_population;

	vector<Unit*> minicolumns = layer->GetLocalRateUnits();

	if(m_adaptValues.size()!=minicolumns.size())
		m_adaptValues = vector<float>(minicolumns.size());  // initialization to 0

	if(m_Aj.size()!=minicolumns.size())  // initialize m_Aj late if necessary
	{
		m_Aj = vector<float>(minicolumns.size());
		for(int i=0;i<minicolumns.size();i++)
			m_Aj[i]=1.0;
	}

	if(m_adampl.size()!=minicolumns.size())
		m_adampl = vector<float>(minicolumns.size()); // initialize

	vector<float> data = vector<float>(minicolumns.size());
	// pulling the values from BCPNN? Or having it independent of BCPNN? 
	// ProjectionModifierBcpnnOnline* bcpnn = (ProjectionModifierBcpnnOnline*)this->Layer()->GetIncomingProjections()[0]->GetEvent("bcpnn"); // assumes bcpnn is in incoming Projection 0 and that bcpnn is used
	// m_adampl = bcpnn->GetBeta();

	for(int i=0;i<minicolumns.size();i++){
		data[i] = ((RateUnit*)minicolumns[i])->GetValue();  
	}


	if(IsRecording()){			
		for(int i=0;i<minicolumns.size();i++) 	
			switchdata.push_back(data[i]); // copy the data for the output meter
	}

	//////////////////////////////////////////////////////////////////
	// normalization of activations: start
	//////////////////////////////////////////////////////////////////
	// data should be positive

	float maxD=*(std::max_element( data.begin(), data.end() ));
	float datasum=0;
	if(maxD>0)
		for(int i=0;i<minicolumns.size();i++){
			data[i] = exp(data[i]-maxD);
			datasum+=data[i];
		}
	else{
		for(int i=0;i<minicolumns.size();i++){
			data[i] = exp(data[i]);
			datasum+=data[i];
		}
	}

	for(int i=0;i<minicolumns.size();i++){
		data[i] = data[i]/datasum;		
	}
	
	//////////////////////////////////////////////////////////////////
	// normalization of activations: end
	//////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////
	// WTA as per NetworkPopulationModifier.cpp
	//////////////////////////////////////////////////////////////////
	/*
	float maxValue = 0.2; // arbitrary threshold
	vector<int> maxIndex(0);

	for(int j=0;j<minicolumns.size();j++)
	{
		if(data[j]>maxValue)
		{
			//maxValue = data[j];
			maxIndex.push_back(j);
		}
	}

	if(maxValue > -1e10)//> m_threshold)
	{
		for(int j=0;j<minicolumns.size();j++)
		{
			if(j == maxIndex)
			{
				data[j] = 1.0;
			}
			else
				data[j] = 0;
		}
	}
	
	for(int j=0;j<minicolumns.size();j++)
			data[j] = 0;
	for(int i=0;i<maxIndex.size();i++)
		data[maxIndex[i]]=1;
	*/

	//////////////////////////////////////////////////////////////////
	// end WTA
	//////////////////////////////////////////////////////////////////

	//if(m_adtaudt == 0)
	//{	
		// m_Aj should be calculated during Network training... adaptation only when switched on...

		//////////////////////////////////////////////////////////////////
		// calculation of activation statistics (Pj)
		//////////////////////////////////////////////////////////////////
		for(int i=0;i<minicolumns.size();i++) // calculate adaptation values m_adampl[]
		//for all cell activations (postValues), code taken from BCPNN to include Anders' switching formulas
		{
			//data[i] = -((RateUnit*)minicolumns[i])->GetValue();  // in Anders' formula the activations are positive???
			m_Aj[i] = m_Aj[i] + m_alpha*(( (1-m_lambda0)*data[i] + m_lambda0)-m_Aj[i]);
			m_adampl[i] = - m_impactBeta * log(m_Aj[i]);  // Anders' first formula for adaptation. 
			// Note that the variable in the log term, m_Aj, should be positive (and possibly between 0 and 1)
			// I changed it to -
			if(m_adampl[i]!=m_adampl[i]) // if NaN
				m_adampl[i]=0;
		}
		//////////////////////////////////////////////////////////////////
		// end calculation of activation statistics
		//////////////////////////////////////////////////////////////////

	//}

	if(m_adtaudt > 0) // when switching turned on
	{


		// only m_Aj (beta) update, while switching is not on
		//////////////////////////////////////////////////////////////////
		// now adapt the activations according to adaptation
		//////////////////////////////////////////////////////////////////
		for(int i=0;i<minicolumns.size();i++) 
		{
			m_adaptValues[i] += (m_adampl[i]*data[i] - m_adaptValues[i]) * m_adtaudt; // Anders' second formula for adaptation
			// m_adaptValues[i] += (10*data[i] - m_adaptValues[i]) * m_adtaudt; // old formula, m_adaptValues oscillate roughly between -60 and -120
			((RateUnit*)minicolumns[i])->SetValue(data[i]-m_adaptValues[i]);
		}
	}

}

vector<vector<float> > PopulationModifierAdaptation2::GetValuesToRecord()
{	
	PopulationColumns* layer = (PopulationColumns*)m_population;
	vector<Unit*> minicolumns = layer->GetLocalRateUnits();
	int iters=switchdata.size()/minicolumns.size();
	vector<vector<float> > f(iters);	
	int index=0;
	for(int j=0;j<iters;j++)		
		for(int i=0;i<minicolumns.size();i++){
			f[j].push_back(switchdata[index]);
			index++;
		}

	return f;
}
