#include "NetworkDemo0.h"

// Defines network structure
void NetworkDemo0::NetworkSetupStructure()
{
	// total nrHypercolumns*nrRateUnits neural units
	int nrHypercolumns = 10;
	int nrRateUnits = 100;

	// define population, sett units as rate units (PopulationColumns::Graded)
	PopulationColumns* pop1 = new PopulationColumns(this,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	this->AddPopulation(pop1); // this population/layer will have index 0

	// Population fully connected to itself by random synaptic weights
	FullConnectivity* full1 = new FullConnectivity();
	full1->SetRandomWeights(0,10);
	pop1->AddPre(pop1,full1);

	// use a linear transfer function
	TransferLinear* transferLin = new TransferLinear();
	pop1->AddUnitsModifier(transferLin);
	
	// Winner-take-all operation on the values in each column
	WTA* wta = new WTA();
	pop1->AddPopulationModifier(wta);
}

// Setup outputs
void NetworkDemo0::NetworkSetupMeters()
{
	// example of file output
	Meter* layerMeter = new Meter("activity.csv",Storage::CSV,Storage::Standard);
	layerMeter->AttachPopulation(this->GetLayer(0));
	this->AddMeter(layerMeter);

	// also check simulation times
	this->SetTiming(true,this);
}

void NetworkDemo0::NetworkSetupParameters()
{
	// not used
}

void NetworkDemo0::NetworkRun()
{
	int nrTimeSteps = 100; // total nr of simulation time steps
	int n = 10; // inject data every nth iteration
	
	// simulate loop
	for(int i=0;i<nrTimeSteps;i++)
	{
		if(i % n == 0)
		{
			// inject some random data
			vector<float> randData(minicolumns.size());
			for(int m=0;m<randData.size();m++) randData[m] = rand()/(float(RAND_MAX)+1);
			this->GetPopulation(0)->SetValuesAll(randData, true);
		}

		this->Simulate(); // simulates 1 timestep
	}
}