#include <math.h>
#include "NetworkPopulationModifier.h"
#include "MPIDistribution.h"
#include <algorithm>


vector<double> SoftMax::ksoftwinners(vector<double> data,int k){	
	vector<double> sorted(data); 
	vector<double> output(data.size(),0); 
	sort(sorted.begin(), sorted.end(),greater<double>()); // descending order	
	//double sum = 0;	
	int ks=0;
 	for(int i=data.size()-1;i>=0;i--){
		if( (data[i]>=sorted[k-1]) ){	
			output[i]=exp(m_G*data[i]);
			if(output[i]!=output[i])
				cerr<<"NAN("<<data[i]<<")";

	//		if(output[i]>sum)
	//			sum=output[i];
			ks++;
			if(ks>=k) 
				break;
		}
		if(data[i]!=data[i])
			cerr << "KSOFT::Nan encountered" << endl;
	}
	//if(sum>0)  // avoid nan
	//	for(int i=0;i<output.size();i++){
	//		output[i]=output[i]/sum; 
	//	}	
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

void SoftMax::Simulate()
{
	PopulationColumns* layer = (PopulationColumns*)m_population;

	/*if(layer->IsFirstRun() == true)
	{
		layer->MPI()->MPICreateCommsHypercolumns(); // could be initialized sooner
		layer->SetFirstRun(false);
	}*/

	layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // MPI distributor

	vector<int> localHypercolumns = layer->GetLocalHypercolumnIndexes();
	vector<Unit*> units = layer->GetUnits();

	// no need to let all nodes do this since values are distributed anyway now

	if(m_type == SoftMax::WTAThresholded)
	{
		for(int i=0;i<localHypercolumns.size();i++)
		{
			//vector<int> mcsIndexes = layer->GetRateUnitsIndexes(localHypercolumns[i]);
			//vector<double> data(mcsIndexes.size());
			Hypercolumn* h = (Hypercolumn*)layer->GetHypercolumn(i);
			vector<RateUnit*> mcs = h->GetRateUnits();

			/*for(int j=mcsIndexes.size()-1;j>-1;j--)
			{
				data[j] = ((RateUnit*)units[mcsIndexes[j]])->GetValue();
			}*/
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

			if(maxValue == 1.0)//> m_threshold)
			{

				for(int j=0;j<data.size();j++)
				{
					if(j == maxIndex)
					{
						data[j] = 1.0;
					}
					else
						data[j] = 0;

					//					((RateUnit*)units[mcsIndexes[j]])->SetValue(data[j]);
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
			for(int j=mcs.size()-1;j>-1;j--){
				data[j] = mcs[j]->GetValue();

				//vector<int> mcsIndexes = layer->GetRateUnitsIndexes(localHypercolumns[i]);
				//vector<double> data(mcsIndexes.size());
				//int Hsum=0;

				//for(int j=mcsIndexes.size()-1;j>-1;j--){
				//	data[j] = ((RateUnit*)units[mcsIndexes[j]])->GetValue();
				//Hsum+=data[j];
			}
			//if(Hsum>0){
			int k=int(floor(data.size()*0.5));//Hsum/data.size())); // regulate sparseness here				
			data = ksoftwinners(data,k);  
			//}
			//for(int j=0;j<mcsIndexes.size();j++){
			//	((RateUnit*)units[mcsIndexes[j]])->SetValue(data[j]);
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
			//vector<int> mcsIndexes = layer->GetRateUnitsIndexes(localHypercolumns[i]);
			//vector<double> data(mcsIndexes.size());

			/*for(int j=mcsIndexes.size()-1;j>-1;j--)
			{
				data[j] = ((RateUnit*)units[mcsIndexes[j]])->GetValue();
			}*/
			vector<float> dt = h->GetValues();
			vector<double> data(dt.size());
			for(int j=0;j<dt.size();j++) data[j] = dt[j];

			data = Function(data,m_G);

			/*for(int j=0;j<mcsIndexes.size();j++)
			{
				((RateUnit*)units[mcsIndexes[j]])->SetValue(data[j]);
				//if(m_probabilisticWTA == true)
				((RateUnit*)units[mcsIndexes[j]])->AddGains(); // only local units will !=0.0
			}*/

			for(int j=0;j<mcs.size();j++)
			{
				mcs[j]->SetValue(data[mcs[j]->GetUnitIdLocalInHypercolumn()]);
				mcs[j]->AddGains();
			}
		}

		// distribute all added gains
		layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // MPI distributor

		// this work anymore?
		if(m_type == WTA || m_type == ProbWTA)
		{
			for(int i=0;i<localHypercolumns.size();i++)
			{
				//vector<int> mcsIndexes = layer->GetRateUnitsIndexes(localHypercolumns[i]);
				/*Hypercolumn* h = layer->GetHypercolumn(i);
				vector<double> data(h->GetTotalNrRateUnits());//mcsIndexes.size());
				vector<RateUnit*> mcs = h->GetRateUnits();

				for(int j=0;j<//mcsIndexes.size()-1;j>-1;j--)
				{
					data[j] = ((RateUnit*)units[mcsIndexes[j]])->GetValue();
				}

				switch(m_type){
				case ProbWTA: 
					data = WTAProb(data);//tempWTAFunction(data);//WTAProb(data);//tempWTAFunction(data);//WTAProb(data);// rand run on all nodes - same seed
					break;
				case WTA:
					data = tempWTAFunction(data);
					break;
				}	
				//		data = tempWTAFunction(data);//Function(data,m_G);

				for(int j=0;j<mcsIndexes.size();j++)
				{
					((RateUnit*)units[mcsIndexes[j]])->SetValue(data[j]);
				}*/
			}
		}
	}
}

std::vector<double> SoftMax::WTAProb(std::vector<double> data)
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

// In hypercolumn setting: currently performed on all nodes involved in hypercolumn

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

		/*for(int i=0;i<data.size();i++)
		{
			if(G>0)
				data[i] = data[i]/(G*4);
			else
				data[i] = data[i]/(4);
		}*/
	}

	for(int i=0;i<data.size();i++)
		sum += exp(G*data[i]);

	if(sum==0.0)
		cout<<"Error: Sum==0 in softmax.";

	for(int i=0;i<data.size();i++)
		output.push_back(exp(G*data[i]) / sum);

	return output;
}


void WTA::Simulate()
{
	if(!IsOn()) return;

	TimingStart(m_name);

	PopulationColumns* layer = (PopulationColumns*)m_population;

	/*if(layer->IsFirstRun() == true)
	{
		layer->MPI()->MPICreateCommsHypercolumns(); // could be initialized sooner
		layer->SetFirstRun(false);
	}*/

//	layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // MPI distributor

	vector<int> localHypercolumns = layer->GetLocalHypercolumnIndexes();
	vector<Unit*> units = layer->GetUnits();

	// no need to let all nodes do this since values are distributed anyway now

	vector<int> mcsIndexes;

	////////////////////////////////////////////////////////////////
	// OBS: replaced the added gain method by maps in transfer function
	////////////////////////////////////////


/*	for(int i=0;i<localHypercolumns.size();i++)
	{
		mcsIndexes = layer->GetRateUnitsIndexes(localHypercolumns[i]);
		vector<double> data(mcsIndexes.size());

		for(int j=mcsIndexes.size()-1;j>-1;j--)
		{
			data[j] = ((RateUnit*)units[mcsIndexes[j]])->GetValue();
		}

		//data = Function(data);

		for(int j=0;j<mcsIndexes.size();j++)
		{
			((RateUnit*)units[mcsIndexes[j]])->SetValue(data[j]);
			((RateUnit*)units[mcsIndexes[j]])->AddGains(); // only local units will !=0.0
		}
	}
*/
	// distribute all added gains
	// this will usually not do anything as all minicolumns are usually local (unless very many processes cmp to network size are used)
	layer->MPI()->MPIMakeHypercolumnsValuesLocal(); // MPI distributor

	// redo wta including added gains
	for(int i=0;i<localHypercolumns.size();i++)
	{
		Hypercolumn* h = (Hypercolumn*)layer->GetHypercolumn(i);
		
		if(h->IsSilent() == false)
		{
			//Hypercolumn* h = (Hypercolumn*)m_network->GetUnitFromId(localHypercolumns[i]);
			vector<RateUnit*> mcs = h->GetRateUnits();
			//vector<int> mcsIndexes = layer->GetRateUnitsIndexes(localHypercolumns[i]);

			vector<float> data = h->GetValues();


			/*		vector<double> data(mcs.size());//(mcs.size());

			for(int j=0;j<mcs.size();j++)
			{
			data[j] = mcs[j]->GetValue();
			}
			*/
			data = Function(data);

			for(int j=0;j<mcs.size();j++)
			{
				//	cout<<mcs[j]->GetUnitIdLocalInHypercolumn()<<" ";cout.flush();
				mcs[j]->SetValue(data[mcs[j]->GetUnitIdLocalInHypercolumn()]);//data[j]);
			}
		}

		/*
		vector<int> mcsIndexes = layer->GetRateUnitsIndexes(localHypercolumns[i]);
		vector<double> data(mcsIndexes.size());

		for(int j=mcsIndexes.size()-1;j>-1;j--)
		{
			data[j] = ((RateUnit*)m_network->GetUnitFromId(mcsIndexes[j]))->GetValue();//((RateUnit*)units[mcsIndexes[j]])->GetValue();
		}

		data = Function(data);

		for(int j=0;j<mcsIndexes.size();j++)
		{
			((RateUnit*)m_network->GetUnitFromId(mcsIndexes[j]))->SetValue(data[j]);//((RateUnit*)units[mcsIndexes[j]])->SetValue(data[j]);
		}
		*/
	}

	TimingStop(m_name);
}

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
			data = m_wta.Function(data);

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

vector<float> WTA::Function(vector<float> data)
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
