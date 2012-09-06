#include "NetworkUnitTests.h"
#include "NetworkBCPNN.h"

bool statusRun;

void SetRunStatus(bool b)
{
	statusRun = b;
}

void UnitTests::Run()
{
	int mpiSize,mpiRank;
	MPI_Comm_size(NETWORK_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(NETWORK_COMM_WORLD, &mpiRank);

	bool ProjectionsOk = true;
	
	// Projections
	NetworkProjectionTests testProjections;
	testProjections.Run();
	
	if(!statusRun)
		ProjectionsOk = false;

	// ...
	if(mpiRank == 0)
	{
		cout<<"\n\nStatus:\n-------------\n";
		if(ProjectionsOk) cout<<"\nProjections test ok.\n"; else cout<<"Projections test error.\n";
	}
}

void NetworkProjectionTests::NetworkSetupStructure()
{
	//PopulationColumns* layer1 = new PopulationColumns(this,10,10,PopulationColumns::Graded);
	//this->AddPopulation(layer1);

	// total nrHypercolumns*nrRateUnits neural units
	m_nrHypercolumns = 128*1;//128;//*128;
	m_nrRateUnits = 1;

	PopulationColumns* layer0 = new PopulationColumns(this, m_nrHypercolumns,
		m_nrRateUnits, PopulationColumns::Graded);
	this->AddPopulation(layer0); // this population/layer will have index 0

	PopulationColumns* layer1 = new PopulationColumns(this, m_nrHypercolumns,
		m_nrRateUnits, PopulationColumns::Graded);
	this->AddPopulation(layer1); // this population/layer will have index 1
	layer1->AddUnitsProperty(new TransferLinear(false));

	PopulationColumns* layer2 = new PopulationColumns(this, m_nrHypercolumns,
		m_nrRateUnits, PopulationColumns::Graded);
	this->AddPopulation(layer2); // this population/layer will have index 2
	
	OneToOneConnectivity* conn = new OneToOneConnectivity(false);
	RandomConnectivity* rndConn = new RandomConnectivity(0.5,false); // correct ??? 
	FullConnectivity* fullConn = new FullConnectivity();
	fullConn->SetRandomWeights(0,1);
	conn->SetRandomWeights(1,1);
	rndConn->SetRandomWeights(1,1);
	layer1->AddPre(layer0,conn);
	layer2->AddPre(layer0,fullConn);

	float lambda0 = 0.0001;
	float alpha = 0.01;
	ProjectionModifierBcpnnOnline* bcpnn = new ProjectionModifierBcpnnOnline(alpha,lambda0);

	fullConn->AddProjectionsEvent(bcpnn);

	//conn->SetUnitType("hypercolumn");

	//layer1->AddPre(layer0,conn);

	// Population fully connected to itself by random synaptic weights
	//FullConnectivity* full1 = new FullConnectivity();
	//full1->SetRandomWeights(0,10);
	//layer1->AddPre(layer1,full1);

	// Winner-take-all operation
	//WTA* wta = new WTA();
	//layer1->AddPopulationModifier(wta);

	this->SetSeed(10);
}

void NetworkProjectionTests::NetworkSetupMeters()
{
	Meter* layerMeterOut = new Meter("TestProjectionsLayer1.csv", Storage::CSV);
	layerMeterOut->AttachPopulation(this->GetLayer(0));
	Meter* layerMeterOut2 = new Meter("TestProjectionsLayer2.csv", Storage::CSV);
	layerMeterOut2->AttachPopulation(this->GetLayer(1));

	this->AddMeter(layerMeterOut);
	this->AddMeter(layerMeterOut2);
}

void NetworkProjectionTests::NetworkSetupParameters()
{

}

void NetworkProjectionTests::NetworkRun()
{
	this->GetLayer(0)->SwitchOnOff(false);
	
	vector<float> r(128);
	r[2] = 1;
	r[3] = 1;
	r[126] = 1;

	bool oneToOneOk = true;

	vector<float> prevValues;

	for(int i=0;i<100;i++)
	{
		this->GetLayer(0)->SetValuesAll(r);
		//this->GetLayer(2)->SetValuesAll(r);

		this->Simulate();
		this->GetLayer(0)->MPI()->MPIMakeLayerValuesLocal();
		this->GetLayer(1)->MPI()->MPIMakeLayerValuesLocal();

		// Tests one-to-one Projection
		if(i>0)
		{
			vector<float> values = this->GetLayer(1)->GetValuesBuffer();
			for(int j=0;j<values.size();j++)
				if(values[j] != prevValues[j])
					oneToOneOk = false;
		}

		prevValues = this->GetLayer(0)->GetValuesBuffer();
	}

	// updates for all processes
	SetRunStatus(true);
	if(oneToOneOk == false) SetRunStatus(false);

	this->RecordAll();
}




void NetworkSynapsesTests::NetworkSetupStructure()
{
	int nrHypercolumns = 3000;//*128;
	int nrRateUnits = 1;

	PopulationColumns* layer0 = new PopulationColumns(this, nrHypercolumns,
		nrRateUnits, PopulationColumns::Graded);
	this->AddPopulation(layer0); // this population/layer will have index 0
	
	//layer0->AddUnitsProperty(new TransferLinear(false));

	FullConnectivity* conn = new FullConnectivity();

	m_bcpnn = new ProjectionModifierBcpnnOnline();
	conn->AddProjectionsEvent(m_bcpnn);

	layer0->AddPre(layer0,conn);
	this->AddTiming(this);
	this->AddTiming(layer0);
	vector<RateUnit*> mcs = layer0->GetRateUnits();
	for(int i=0;i<mcs.size();i++)
		this->AddTiming(mcs[i]);
}

void NetworkSynapsesTests::NetworkRun()
{
	bool synapsesOk = true;
	bool eraseAll = false;

	m_bcpnn->SwitchOnOff(false);

	///////////////////////////////////////////////////
	///// Check deletion of some synapses 
	///////////////////////////////////////////////////

	vector<long> postIds = this->GetLayer(0)->GetIncomingProjections()[0]->GetPostIds();

	vector<pair<long,long> > ids;
	for(int i=0;i<postIds.size();i++)
	{
		vector<long> preIds = this->GetLayer(0)->GetIncomingProjections()[0]->GetPreIds(postIds[i]);
		
		for(int j=0;j<preIds.size();j++)
		{
			pair<long,long> p;
			p.first = postIds[i];
			p.second = preIds[j];
			if(eraseAll)
				ids.push_back(p);
			else
			{
				if(rand()%2 == 0)
					ids.push_back(p);
			}
		}
	}

	// check
	vector<long>* preIdsAll = this->GetLayer(0)->GetIncomingProjections()[0]->GetPreIdsAll();
	int totalNr = preIdsAll->size();
	if(preIdsAll->size() == 0)
		synapsesOk = false;

	TimingStart("RunAllSynapses");
	this->Simulate(30);
	TimingStop("RunAllSynapses");

	// erase some/all synapses
	this->GetLayer(0)->GetIncomingProjections()[0]->EraseSynapses(ids);

	// check (this also checks CreatePreIdsUnion work)
	preIdsAll = this->GetLayer(0)->GetIncomingProjections()[0]->GetPreIdsAll();
	
	if(eraseAll == true && preIdsAll->size() != 0)
		synapsesOk = false;

	MPI_Barrier(MPI_COMM_WORLD);

	TimingStart("RunSomeSynapses");
	this->Simulate(30);
	TimingStop("RunSomeSynapses");

	// updates for all processes
	SetRunStatus(true);
	if(synapsesOk == false) SetRunStatus(false);

	this->RecordAll();
	this->StoreAnalysis();
}