#include "NetworkDemoVis0.h"

// Defines network structure
void NetworkDemoVis0::NetworkSetupStructure()
{
	// total nrHypercolumns*nrRateUnits neural units
	int nrHypercolumns = 1;
	int nrRateUnits = 128*128;

	PopulationColumns* layer1 = new PopulationColumns(this,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	this->AddLayer(layer1); // this population/layer will have index 0

	// Population fully connected to itself by random synaptic weights
	FullConnectivity* full1 = new FullConnectivity();
	full1->SetRandomWeights(0,10);
	layer1->AddPre(layer1,full1);

	// Winner-take-all operation
	//WTA* wta = new WTA();
	//layer1->AddLayerEvent(wta);

	this->SetSeed(10);
}

// Setup outputs
void NetworkDemoVis0::NetworkSetupMeters()
{
	// example of file output
	Meter* layerMeter = new Meter("activity.csv",Storage::CSV,Storage::Standard);
	layerMeter->AttachLayer(this->GetLayer(0));
	this->AddMeter(layerMeter);

	// also check simulation times
	this->SetTiming(true,this);
}

void NetworkDemoVis0::NetworkSetupParameters()
{
	// not used
}

void NetworkDemoVis0::NetworkRun()
{
	int nrTimeSteps = 100;

	// local minicolumns/neural units on this process for population 0
	vector<RateUnit*> minicolumns = ((PopulationColumns*)this->GetLayer(0))->GetRateUnits();
	((PopulationColumns*)this->GetLayer(0))->GetNrRateUnits();
	
//	for(int i=0;i<minicolumns.size();i++)
//		minicolumns[i]->GetUnitIdLocal();

	// simulate
	for(int i=0;i<nrTimeSteps;i++)
	{
		// inject some random data every nth iteration
		if(i%10 == 0)
		{
			vector<float> randData(minicolumns.size());
			for(int m=0;m<randData.size();m++) randData[m] = rand()/(float(RAND_MAX)+1);
			this->GetLayer(0)->SetValuesAll(randData, true);
		}

		// check buffer
		// insert into network

		this->Simulate(); // 1 timestep

		// tic

		// example of how to retrieve current activity (local for this process)
		vector<float> values(minicolumns.size());
		for(int j=0;j<minicolumns.size();j++)
			values[j] = minicolumns[j]->GetValue();
	}

	// store data afterwards for this case
	this->RecordAll();
	this->StoreTimings();
}