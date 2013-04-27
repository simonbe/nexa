#include <math.h>
#include "NetworkPopulationModifier.h"
#include "MPIDistribution.h"
#include <algorithm>
#include <functional>

vector<double> SoftMax::ksoftwinners(vector<double> data,int k){	
	vector<double> sorted(data); 
	vector<double> output(data.size(),0); 
	sort(sorted.begin(), sorted.end(),greater<double>()); // descending order	

	int ks=0;
 	for(int i=data.size()-1;i>=0;i--){
		if( (data[i]>=sorted[k-1]) ){	
			output[i]=exp(m_G*data[i]);
			if(output[i]!=output[i])
				cerr<<"NAN("<<data[i]<<")";

			ks++;
			if(ks>=k) 
				break;
		}
		if(data[i]!=data[i])
			cerr << "KSOFT::Nan encountered" << endl;
	}

	return output;
}

vector<double> SoftMax::tempWTAFunction(vector<double> data)
{
	vector<double> output(data.size());
	double sum = 0;

	double maxValue = -1e10;
	int maxIndex = -1;
	
	for(int i=0;i<data.size();i++)
	{
		if(data[i]>maxValue)
		{
			maxValue = data[i];
			maxIndex = i;
		}
	}

	output[maxIndex] = 1.0;

	return output;
}

/// <summary>	Performs a soft winner-take-all operation </summary>

void SoftMax::Simulate()
{
	PopulationColumns* layer = (PopulationColumns*)m_population;

	layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // MPI distributor

	vector<int> localHypercolumns = layer->GetLocalHypercolumnIndexes();
	vector<Unit*> units = layer->GetUnits();

	// no need to let all nodes do this since values are distributed anyway now

	if(m_type == SoftMax::WTAThresholded)
	{
		for(int i=0;i<localHypercolumns.size();i++)
		{
			Hypercolumn* h = (Hypercolumn*)layer->GetHypercolumn(i);
			vector<RateUnit*> mcs = h->GetRateUnits();

			vector<float> data = h->GetValues();

			float maxValue = -1e10;
			int maxIndex = -1;

			for(int j=0;j<data.size();j++)
			{
				if(data[j]>maxValue)
				{
					maxValue = data[j];
					maxIndex = j;
				}
			}

			if(maxValue == 1.0)
			{

				for(int j=0;j<data.size();j++)
				{
					if(j == maxIndex)
					{
						data[j] = 1.0;
					}
					else
						data[j] = 0;
				}
			}

			for(int j=0;j<mcs.size();j++)
				mcs[j]->SetValue(data[mcs[j]->GetUnitIdLocalInHypercolumn()]);

		}
	}
	else if(m_type==KSOFT){
		for(int i=0;i<localHypercolumns.size();i++)
		{
			Hypercolumn* h = (Hypercolumn*)layer->GetHypercolumn(i);
			vector<RateUnit*> mcs = h->GetRateUnits();
			vector<double> data(mcs.size());
			for(int j=mcs.size()-1;j>-1;j--)
			{
				data[j] = mcs[j]->GetValue();
			}
			
			int k=int(floor(data.size()*0.5));//Hsum/data.size())); // regulate sparseness
			data = ksoftwinners(data,k);  
			
			for(int j=0;j<mcs.size();j++){
				mcs[j]->SetValue(data[j]);
			}
		}
	}
	else  // standard; wta, probwta
	{
		for(int i=0;i<localHypercolumns.size();i++)
		{
			Hypercolumn* h = (Hypercolumn*)layer->GetHypercolumn(i);
			vector<RateUnit*> mcs = h->GetRateUnits();

			vector<float> dt = h->GetValues();
			vector<double> data(dt.size());
			for(int j=0;j<dt.size();j++) data[j] = dt[j];

			data = Function(data,m_G);

			for(int j=0;j<mcs.size();j++)
			{
				mcs[j]->SetValue(data[mcs[j]->GetUnitIdLocalInHypercolumn()]);
				mcs[j]->AddGains();
			}
		}

		// distribute all added gains
		layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // MPI distributor
	}
}

vector<double> SoftMax::WTAProb(vector<double> data)
{
	double sum = 0;
	for(int i=0;i<data.size();i++)
		sum+=data[i];

	float r = (float)rand()/(float)RAND_MAX;
	r=r*sum;

	int win = data.size()-1;
	float currentVal = data[0];

	for(int i=1;i<data.size();i++)
	{
		if(r<currentVal)
		{
			win = i-1;
			i = data.size()-1;
			break;
		}
		else
		{
			currentVal += data[i];
		}
	}

	vector<double> output(data.size());
	output[win] = 1.0;
	return output;
}

/// <summary>	Softmax function. </summary>
///
/// <param name="data">	Column activities. </param>
/// <param name="G"> Sharpness of softmax.  </param>
///
/// <returns> Column activities. </returns>

vector<double> SoftMax::Function(vector<double> data, float G)
{
	vector<double> output;
	double sum = 0;

	bool lower = false;
	for(int i=0;i<data.size();i++)
	{
		if(data[i]>300)
		{
			data[i] = 200;
			lower = true;
		}
	}

	if(lower == true)
	{
		if(m_population->network()->MPIGetNodeId() == 0)
			cout<<"<";//cout<<"Warning: Large in softmax. Forcing lowering.";
	}

	for(int i=0;i<data.size();i++)
		sum += exp(G*data[i]);

	if(fabs(sum)<EPS)
		cout<<"Error: Sum==0 in softmax.";

	for(int i=0;i<data.size();i++)
		output.push_back(exp(G*data[i]) / sum);

	return output;
}

/// <summary>	Winner-take-all operation. </summary>

void WTA::Simulate()
{
	if(!IsOn()) return; // TODO: move outside

	TimingStart(m_name);

	PopulationColumns* population = (PopulationColumns*)m_population;

//	layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // MPI distributor

	vector<int> localHypercolumns = population->GetLocalHypercolumnIndexes();
	vector<Unit*> units = population->GetUnits();

	// no need to let all nodes do this since values are distributed anyway now

	vector<int> mcsIndexes;

	////////////////////////////////////////////////////////////////
	// Compared to earlier version: replaced the added gain method by maps in transfer function
	////////////////////////////////////////

	// distribute all added gains
	// this will often not do anything as all minicolumns are local (unless many processes compared to network size are used)
	population->MPI()->MakeActivityValuesLocal();//MPIMakeHypercolumnsValuesLocal(); // MPI distributor

	// redo wta including added gains
	for(int i=0;i<localHypercolumns.size();i++)
	{
		Hypercolumn* h = (Hypercolumn*)population->GetHypercolumn(i);
		
		if(h->IsSilent() == false) // special case if true (not used in any standard model)
		{
			//Hypercolumn* h = (Hypercolumn*)m_network->GetUnitFromId(localHypercolumns[i]);
			//vector<RateUnit*> mcs = h->GetRateUnits();
			//vector<int> mcsIndexes = layer->GetRateUnitsIndexes(localHypercolumns[i]);

			vector<float> data = h->GetValues();
			data = wta(data);

			h->SetValues(data);

			/*for(int j=0;j<mcs.size();j++)
			{
				mcs[j]->SetValue(data[mcs[j]->GetUnitIdLocalInHypercolumn()]);//data[j]);
			}*/
		}
	}

	TimingStop(m_name);
}

/// <summary>	Winner-take-all. </summary>
///
/// <param name="data">	Column/population activities. </param>
///
/// <returns>	Column/population activities. </returns>

vector<float> WTA::wta(vector<float> data)
{
	vector<float> output(data.size());
	double sum = 0;

	double maxValue = -1e10;
	int maxIndex = -1;
	float r;

	for(int i=0;i<data.size();i++)
	{
		if(data[i]>maxValue)
		{
			maxValue = data[i];
			maxIndex = i;
		}
		else if(data[i]==maxValue) // randomly switch
		{
			r = (float)rand()/(float)RAND_MAX;
			if(r>0.5)
				maxIndex = i;
		}
	}

	if(maxIndex!=-1)
		output[maxIndex] = 1.0;

	return output;
}


/// <summary>	Winner-take-all-operation with a threshold. </summary>

void WTAThreshold::Simulate()
{
	if(!IsOn()) return;

	TimingStart(m_name);

	PopulationColumns* layer = (PopulationColumns*)m_population;

	vector<int> localHypercolumns = layer->GetLocalHypercolumnIndexes();
	vector<Unit*> units = layer->GetUnits();
	vector<int> mcsIndexes;

	layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // MPI distributor

	for(int i=0;i<localHypercolumns.size();i++)
	{
		Hypercolumn* h = (Hypercolumn*)layer->GetHypercolumn(i);
		vector<RateUnit*> mcs = h->GetRateUnits();
		vector<float> data = h->GetValues();

		bool oneAbove = false;
		for(int j=0;j<data.size();j++)
		{
			if(data[j]>0.0000001)
			{
				oneAbove=true;
				break;
			}
		}

		if(oneAbove == true)
		{
			data = m_wta.wta(data);

			for(int j=0;j<mcs.size();j++)
			{
				if(data[mcs[j]->GetUnitIdLocalInHypercolumn()] == 1.0)
					bool b=false;//mcs[j]->SetValue(data[mcs[j]->GetUnitIdLocalInHypercolumn()]);//data[j]);
				else
				{
					mcs[j]->SetValue(data[mcs[j]->GetUnitIdLocalInHypercolumn()]);//data[j]);
					mcs[j]->SetSubThresholdValue(data[mcs[j]->GetUnitIdLocalInHypercolumn()]);//data[j]);
				}
			}
		}
	}

	TimingStop(m_name);
}


/// <summary>	Threshold-operation. </summary>

void Threshold::Simulate()
{
	if(!IsOn()) return;

	PopulationColumns* layer = (PopulationColumns*)m_population;
	vector<RateUnit*> minicolumns = layer->GetRateUnits();

	m_numActiveUnits = 0;
	for(int i=0;i<minicolumns.size();i++)
	{
	   if(minicolumns[i]->GetValue()>=m_threshold){
		   minicolumns[i]->SetValue(1.0);
		   m_numActiveUnits++;
	   } else
		   minicolumns[i]->SetValue(0.0);
	}
}

/// <summary>	Divise normalization-operation. </summary>

void DivisiveNormalization::Simulate()
{
	if(!IsOn()) return;

	TimingStart(m_name);

	PopulationColumns* layer = (PopulationColumns*)m_population;

	vector<int> localHypercolumns = layer->GetLocalHypercolumnIndexes();
	vector<Unit*> units = layer->GetUnits();

	vector<int> mcsIndexes;
	layer->MPI()->MakeActivityValuesLocal(true);//MPIMakeHypercolumnsValuesLocal(); // MPI distributor, not needed anymore?

	vector<float> valuesLocal = layer->GetValuesLocal();
	vector<float> allValues = layer->GetValuesBuffer();
	float totValues = 0;
	for(int i=0;i<allValues.size();i++)
		totValues += allValues[i]*allValues[i];//pow(allValues[i],2);

	if(totValues>0)
	{
		vector<RateUnit*> mcs = layer->GetRateUnits();

		float value, subThreshValue;
		for(int j=0;j<mcs.size();j++)
		{
			value = mcs[j]->GetValue();
			if(value >0)
			{
				value = sqrt(value*value/(1.0+totValues));
				mcs[j]->SetValue(value);
			}
			else // subthreshold normalization
			{
				subThreshValue = mcs[j]->GetSubThresholdValue();
				subThreshValue = sqrt(subThreshValue*subThreshValue/(1.0+totValues));
				mcs[j]->SetSubThresholdValue(subThreshValue);
			}

		}
	}

	TimingStop(m_name);
}
