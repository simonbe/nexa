#include "Network.h"
#include "Meter.h"
#include "Storage.h"
#include "DataSources.h"

#include "StructureMIMDSVQ.h"
#include "NetworkBCPNN.h"
#include "NetworkMDS.h"
#include "NetworkMI.h"
#include "NetworkVQ.h"
//#include "EarlyOlfSys.h"
#include "NetworkAdaptation2.h"

#include "NetworkScalingDemos.h"

void NetworkScalingStrong::NetworkSetupStructure()
{
	m_run = 2;
	int nrPop2;
	bool doAnalysis = false;

	if(m_run == 0) // model 1, small run
	{
		m_sizePopulation1 = 24e2;
		m_sizePopulation2 = 24e2;
		m_nrColumns = 192*5;//24e2;
		m_sizePopulation3 = 24e2/m_nrColumns;
		m_probRecurr = 0.1;
		m_probRecurr2 = 0.01;
		m_probForward = 0.2;
		nrPop2 = 5;
	}
	else if(m_run == 1) // weak scaling run
	{
		m_sizePopulation1 =  3e2;//this->MPIGetNrProcs()*1e2;//this->MPIGetNrProcs()*10;//24e3;//125*this->MPIGetNrProcs();//24e3;
		m_sizePopulation2 = 3e2;//3e3;//24e3;
		m_nrColumns = this->MPIGetNrProcs();//192*10;//24e2;
		m_sizePopulation3 = 100;//24e3/m_nrColumns;
		m_probRecurr = 0.2;
		m_probRecurr2 = 0.1;
		m_probForward = 0.2;
		nrPop2 = this->MPIGetNrProcs();//4;//4*this->MPIGetNrProcs()/192;
	}
	else if(m_run == 2) // strong scaling run
	{
		m_sizePopulation1 =  3e2;//this->MPIGetNrProcs()*1e2;//this->MPIGetNrProcs()*10;//24e3;//125*this->MPIGetNrProcs();//24e3;
		m_sizePopulation2 = 3e3;//24e3;
		m_nrColumns = 0;//24*10;//192*10;//24e2; // not used
		m_sizePopulation3 = 100;//24e3/m_nrColumns;
		m_probRecurr = 0.2;
		m_probRecurr2 = 0.1;
		m_probForward = 0.2;
		nrPop2 = 384;//this->MPIGetNrProcs();//4;//4*this->MPIGetNrProcs()/192;
	}

	m_layer1 = new PopulationColumns(this,1,m_sizePopulation1,PopulationColumns::Graded,MPIDistribution::ParallelizationHypercolumns);
	//m_layer3 = new PopulationColumns(this,m_nrColumns,m_sizePopulation3,PopulationColumns::Graded,MPIDistribution::ParallelizationHypercolumns);

	OneToOneConnectivity* oneToOne = new OneToOneConnectivity(false);
	RandomConnectivity* randConn3 = new RandomConnectivity(m_probRecurr2, false);

	//m_layer3->AddPre(m_layer3,randConn3);//fullConn);
	//TransferLinear* transferLin = new TransferLinear(false);
	//m_layer3->AddUnitsProperty(transferLin);
	
	float lambda0 = 10e-6;
	float alpha = 0.05;
	//m_bcpnn= new ProjectionModifierBcpnnOnline(alpha,lambda0);
	
	//randConn3->AddProjectionsEvent(m_bcpnn);

	// population operations
	WTA* wta = new WTA();
	//m_layer3->AddPopulationModifier(wta);

	this->AddPopulation(m_layer1);
	//this->AddPopulation(m_layer3);

	vector<RandomConnectivity*> randConn1;
	vector<RandomConnectivity*> randConn2;

	// middle (recurrent) populations with linear-threshold units (implemented in transfer function)
	for(int i=0;i<nrPop2;i++)
	{
		PopulationColumns* layer = new PopulationColumns(this,1,m_sizePopulation2,PopulationColumns::Graded);//,MPIDistribution::ParallelizationHypercolumns);
		m_layer2.push_back(layer);

		// sets up pre-defined recurrent structure
		StructureReTIDe* reTIDe = new StructureReTIDe(m_probRecurr,-1,1,true);
		reTIDe->Initialize(this,layer);
		m_reTIDe.push_back(reTIDe);
		
		RandomConnectivity* randConnForward = new RandomConnectivity(m_probForward,false);
		FullConnectivity* full1 = new FullConnectivity();
		OneToOneConnectivity* oneToOne = new OneToOneConnectivity(false);

		randConn1.push_back(randConnForward);
		if(i==-1)
		{
			randConnForward->SetRandomWeights(0.1,0.1);
			full1->SetRandomWeights(0.1,0.1);
			oneToOne->SetRandomWeights(1,1);
		}
		else
		{
			randConnForward->SetRandomWeights(0.08,0.08);//10.0/float(this->MPIGetNrProcs()));//0.01);
			full1->SetRandomWeights(0.0,0.1);
		}

		RandomConnectivity* randConnForward2 = new RandomConnectivity(m_probForward,false);
		randConn2.push_back(randConnForward2);
		randConnForward2->SetRandomWeights(0,1);

		// connect
		layer->AddPre(m_layer1,randConnForward);//full1);//oneToOne);//randConnForward);//full1);//randConnForward);
		//m_layer3->AddPre(layer,randConnForward2);
		// add populations to network
		this->AddPopulation(layer);

		// also check simulation times
		AddTiming(layer);
	}

	// Create some long-range Projections (not symmetric)
	this->SetSeed(true); // guarantee same description on all

	srand(8);
	for(int i=0;i<nrPop2;i++)
	{
		for(int j=0;j<nrPop2;j++)
		{
			if(i!=j)
			{
				if(i==j+1 || i==j-1 || (i==0 && j==nrPop2-1) || (i==nrPop2-1 && j == 0))//rand()%(this->MPIGetNrProcs()) == 0) // 10% chance populations are connected
				{
					RandomConnectivity* randRecurr2 = new RandomConnectivity(0.0005,false);//0.0001,false);
					randRecurr2->SetRandomWeights(-0.1,-0.1); // weak long-range inhibitory Projections

					m_layer2[i]->AddPre(m_layer2[j],randRecurr2);
				}
			}
		}
	}

	this->SetSeed(); // go back to default

	// Set up analysis (optional)
	if(doAnalysis == true)
	{
		// calculates distances in all populations on all processes
		AnalysisDistance* analysisDistance = new AnalysisDistance(m_layer2[this->MPIGetNodeId()],AnalysisDistance::Euclidean);
		this->AddAnalysis(analysisDistance);
	}

	AddTiming(this);
	AddTiming(m_layer1);
	//AddTiming(m_layer3);

	this->SetSeed(true); // to guarantee same connectivities on each process
}

void NetworkScalingStrong::NetworkSetupMeters()
{
	// change (addunitsproperty) so that they can be placed in networkstructure instead
	for(int i=0;i<m_reTIDe.size();i++)
		m_reTIDe[i]->SetupStructure(this,m_layer2[i]);
}

void NetworkScalingStrong::NetworkSetupParameters()
{
	// not used here
}

// resets activities in parts of network
void NetworkScalingStrong::ClearActivities()
{
	// x time steps
	MPI_Barrier(MPI_COMM_WORLD); // to get accurate timing

	TimingStart("ClearActivities1");
	this->ClearEventsIncoming();
	for(int j=0;j<m_layer2.size();j++)
	{
		if(m_layer2[j]->GetUnits().size()>0)
		{
			(*m_layer2[j]->GetUnitPropertiesLayer())[0]->Clear(); // clear transfer fcn history

			vector<RateUnit*> mcs = m_layer2[j]->GetRateUnits();
			for(int m=0;m<mcs.size();m++)
			{
				mcs[m]->SetValue(0.0);
				mcs[m]->SetSubThresholdValue(0.0);
			}
		}
	}

	vector<float> emptyIn(m_sizePopulation1);
	m_layer1->SetValuesAll(emptyIn);
	TimingStop("ClearActivities1");
	MPI_Barrier(MPI_COMM_WORLD); // to get accurate timing

	for(int i=0;i<1;i++)
	{
		this->Simulate();
	}
}

void NetworkScalingStrong::NetworkRun()
{
	int nrPatterns = 10; // 20;
	int nrTimeStepsPerPattern = 300;
	bool lowCommunication = true;
	DataSources dataSource;
	vector<vector<float> > inputData;
	
//	if(m_run == 1)
//	{
//		inputData = dataSource.GetRandomBinary(m_sizePopulation1,0.05*(192.0/(float)this->MPIGetNrProcs()),nrPatterns);
//	}
//	else
//	{
		inputData = dataSource.GetRandomBinary(m_sizePopulation1,0.05,nrPatterns);
		for(int i=0;i<inputData.size();i++)
		{
			for(int j=0;j<inputData[i].size();j++)
			{
				if((j+i)%5 == 0)
					inputData[i][j] = 1.0;
				else
					inputData[i][j] = 0.0;
			}
		}
//	}

	m_layer1->SwitchOnOff(false); // turn off since we set this manually
//	m_bcpnn->SwitchOnOff(false); // plasticity turned on only on last time steps for each pattern

	for(int i=0;i<nrPatterns;i++)
	{
		if(this->MPIGetNodeId()==0)
		{
			cout<<"\n"<<i<<":\n";
			cout.flush();
		}

		m_layer1->SetValuesAll(inputData[i]);

//		m_bcpnn->SwitchOnOff(false);

		for(int j=0;j<nrTimeStepsPerPattern;j++)
		{
//			if(j==nrTimeStepsPerPattern-20)
//				m_bcpnn->SwitchOnOff(true);

			this->Simulate(); // simulate one time step

			// same input data for following time steps, so can freeze this part to lower the communication
			if(lowCommunication)
			{
				for(int k=0;k<m_layer2.size();k++)
				{
					m_layer2[k]->GetIncomingProjections()[1]->KeepActiveBuffer(true);
					vector<float> emptyIn(m_sizePopulation1);
					m_layer1->SetValuesAll(emptyIn);
				}
			}
		}
		
		if(lowCommunication)
		{
			// stop freezing input
			for(int k=0;k<m_layer2.size();k++)
			{
				m_layer2[k]->GetIncomingProjections()[0]->KeepActiveBuffer(false);
			}
		}
		
//		m_bcpnn->SwitchOnOff(false);

		MPI_Barrier(MPI_COMM_WORLD);
		TimingStart("ClearActivities");
		ClearActivities();
		ClearActivities();
		TimingStop("ClearActivities");
	}

	// store timings
	this->RecordAll();
	this->StoreAnalysis();
}

/////////////////////////////////////////////
/// Old scaling runs
/////////////////////////////////////////////

NetworkScalingDemos::NetworkScalingDemos()
{
	m_architecture = JUGENE;//CRAY//JUGENE;//PC;
	m_useBinaryFileWrite = true;
	m_allowFileWrite = true;
}

void NetworkScalingDemos::RunAll()
{
	//m_extraFilenameString = "rec1";
	//NetworkScalingRunRecurrent(false); // no storage of data
	m_extraFilenameString = "rec22";
	NetworkScalingRunRecurrentBCPNN(true); // store data

	m_extraFilenameString = "spik1";
	NetworkScalingRunSpiking(); // storing
	m_extraFilenameString = "mimdsvq1";
	NetworkScalingRunMIMDSVQ();
}

void NetworkScalingDemos::NetworkScalingRunRecurrentBCPNN(bool storeData)
{
	int nrHypercolumns;
	int nrRateUnits;

	if(m_allowFileWrite==false)
		storeData = false;

	float probConnectivity = 0.0003;//0.00005;//0.000025;//0.00003;//0.000025;//0.001;//0.00044;//0.00011;//0.00024;//0.00006;//0.0004;//0.00001;
	//if(m_architecture == PC)
	//{

	if(m_architecture == CRAY)
	{
		nrHypercolumns = 6000*4;//1024*72*4;//1024*72*4;//200;//1024*200;//1024*200;//65536;//1024*10;//65536*4;
		nrRateUnits = 100;//200;//100;
		probConnectivity = 8000/(float)(nrHypercolumns*nrRateUnits);
	}
	else if(m_architecture == JUGENE)
	{
		// medium
		/*nrHypercolumns = 1024*72*4;//1024*72*4;//200;//1024*200;//1024*200;//65536;//1024*10;//65536*4;
		nrRateUnits = 500;//200;//100;
		probConnectivity = 3000/(float)(nrHypercolumns*nrRateUnits);
		*/
		// large
		//3e3 per mc
		nrHypercolumns = 1024*64*4;//1024*64*4;//1024*72*4;//200;//1024*200;//1024*200;//65536;//1024*10;//65536*4;
		nrRateUnits = 250;//200;//100;
		probConnectivity = 6000/(float)(nrHypercolumns*nrRateUnits);
		
		/*nrHypercolumns = 1024*72*4;//1024*72*4;//200;//1024*200;//1024*200;//65536;//1024*10;//65536*4;
		nrRateUnits = 200;//200;//100;
		probConnectivity = 8000/(float)(nrHypercolumns*nrRateUnits);*/
	}
	else if(m_architecture == BGL)
	{
		nrHypercolumns = 1024*1.8*4;//1024*72*4;//128*100;//65536;//65536;//128*16*10;//128*16*10;//128*16*10;//128*16;
		nrRateUnits = 250;//500;//50;//100;//100;//40;//10
		probConnectivity = 6000/(float)(nrHypercolumns*nrRateUnits);//20/(float)(nrHypercolumns*nrRateUnits); // -> 5e9 synapses
	}
	else
	{
		nrHypercolumns = 1024*1;//72*2;//128*20;
		nrRateUnits = 5;
		probConnectivity = 200/(float)(nrHypercolumns*nrRateUnits);//0;//0.01;//100/(float)(nrHypercolumns*nrRateUnits);
	}

	//float activity = 0.05;
	int nrItems = 3;
	int nrSameStimuli = 5;
//	bool storeData = false;

	//}

	bool doTesting = true; // false for some scaling tests

	// network construction
	Network* network = new Network();
	network->SetExtraFilenameString((char*)(m_extraFilenameString.c_str()));
	network->AddTiming(network);

	PopulationColumns* layer1 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	network->AddPopulation(layer1);
	
	FullConnectivity* full = new FullConnectivity();//false,"minicolumns");
	full->SetRandomWeights(0,0);
	RandomConnectivity* randConn = new RandomConnectivity(probConnectivity);//0.1);
	randConn->SetRandomWeights(0,0);	
	//network->AddTiming(randConn);
	
	layer1->AddPre(layer1,randConn); // recurrent

	// Add Projection changes
	float lambda0 = 10e-6;
	float alpha = 0.05;
	ProjectionModifierBcpnnOnline* bStandard = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	
	//full->AddProjectionsEvent(bStandard);
	randConn->AddProjectionsEvent(bStandard);

	WTA* wta = new WTA();
	layer1->AddPopulationModifier(wta);

	network->AddTiming(randConn);
	network->AddTiming(layer1);

	// Construct initial network
	network->Initialize();
	vector<int> partsOfDataToUseAsInput = layer1->GetMPIDistribution();

	// Specify input data
	// - change to only create local part
	DataSources source;
	// not correct right now as SetValuesAll working locally and this is global
	//vector<vector<float> > data = source.GetRandomHCsOrthogonal(nrHypercolumns/100,nrRateUnits,nrItems);
	// correct
	vector<vector<float> > data = source.GetRandomHCsOrthogonalLocal(nrHypercolumns,nrRateUnits,partsOfDataToUseAsInput,nrItems);

	//cout<<"{"<<data[0].size();cout.flush();
	// Meters
	Meter* l1meter;
	Meter* c1meter;
	AnalysisDistance* analysisDistance1;
	string filename1;
	string filename2;
	string filename3;

	if(m_useBinaryFileWrite == true)
	{
		filename1 = "layer1_" + m_extraFilenameString + ".dat";
		l1meter = new Meter((char*)(filename1.c_str()), Storage::MPI_Binary);
		filename2 = "Projections1_" + m_extraFilenameString + ".dat";
		c1meter = new Meter((char*)(filename2.c_str()),Storage::MPI_Binary);
	}
	else
	{
		filename1 = "layer1_" + m_extraFilenameString + ".csv";
		l1meter = new Meter((char*)(filename1.c_str()), Storage::CSV);
		filename2 = "Projections1_" + m_extraFilenameString + ".csv";
		c1meter = new Meter((char*)(filename2.c_str()),Storage::CSV);
	}

	if(storeData)
	{
		l1meter->AttachPopulation(layer1);
		//network->AddMeter(l1meter);
		c1meter->AttachProjection(layer1->GetIncomingProjections()[0],0);
		//network->AddMeter(c1meter);
		analysisDistance1 = new AnalysisDistance(layer1,AnalysisDistance::Euclidean);
		filename3 = "Distances_"+ m_extraFilenameString + ".csv";
		analysisDistance1->SetOwnFileWriting(true,(char*)(filename3.c_str()));
		//analysisDistance1->SetCollectRepresentations(true);
		network->AddAnalysis(analysisDistance1);
	}
	
	// Timings
	network->AddTiming(bStandard);
	network->AddTiming(full);
	

	// need to access after it has been built
	network->AddTiming(layer1->GetIncomingProjections()[0]);

	//cout<<layer1->GetRateUnits().size()<<" ";cout.flush();

	// Training
	// set fixed pattern
	layer1->SwitchOnOff(false);
	if(storeData)
	{
		 analysisDistance1->SwitchOnOff(false);
	}

	int trainIterations = 1;
	int testIterations = 1;

	for(int i=0;i<trainIterations;i++)
	{
		//cout<<i<<"\n";

		for(int j=0;j<data.size();j++)
		{
			for(int k=0;k<nrSameStimuli;k++)
			{
				//layer1->SetValuesAll(data[j]);
				layer1->SetValuesAll(data[j]);
				if(k==0 && storeData==true) analysisDistance1->AddSpecificRepresentation(layer1->GetValuesLocal());

				// next time step
				network->Simulate(1);
			}
		}
	}

	// Testing
	if(storeData)
	{
		analysisDistance1->SwitchOnOff(true);
	}

	if(doTesting == true)
	{
		layer1->SwitchOnOff(true);
		bStandard->SwitchOnOff(false);

		for(int i=0;i<testIterations;i++)
		{
			for(int j=0;j<data.size();j++)
			{
					// clear all events before switching so no disturbance, can remove if moving average activity etc.
					network->ClearEventsIncoming();

					layer1->SetValuesAll(data[j]);
					for(int k=0;k<nrSameStimuli;k++)
					{
						if(k==1)
							layer1->SetValuesAll(data[j]);

						network->Simulate(1);

						//for(int i=0;i<testIterations;i++)
						// next time step
						//			network->Simulate();
					}
			}
		}
	}
	
	network->RecordAll();
	network->StoreAnalysis();
	delete network;
}

void NetworkScalingDemos::NetworkScalingRunRecurrentSimple(bool storeData)
{
	int nrHypercolumns;
	int nrRateUnits;

	if(m_allowFileWrite==false)
		storeData = false;

		float probConnectivity = 0.0003;//0.00005;//0.000025;//0.00003;//0.000025;//0.001;//0.00044;//0.00011;//0.00024;//0.00006;//0.0004;//0.00001;
	//if(m_architecture == PC)
	//{
	

	if(m_architecture == JUGENE)
	{
		nrHypercolumns = 1024*100;//65536;//1024*10;//65536*4;
		nrRateUnits = 100;
		probConnectivity = 10000/(float)(nrHypercolumns*nrRateUnits);
	}
	else if(m_architecture == BGL)
	{
		nrHypercolumns = 1024*10;//128*100;//65536;//65536;//128*16*10;//128*16*10;//128*16*10;//128*16;
		nrRateUnits = 100;//50;//100;//100;//40;//10
		probConnectivity = 5000/(float)(nrHypercolumns*nrRateUnits); // -> 5e9 synapses
	}
	else
	{
		nrHypercolumns = 128*4;
		nrRateUnits = 5;
		probConnectivity = 0.01;//100/(float)(nrHypercolumns*nrRateUnits);
	}

	//float activity = 0.05;
	int nrItems = 5;
	int nrSameStimuli = 20;
//	bool storeData = false;

	//}

	bool doTesting = true; // false for some scaling tests

	// network construction
	Network* network = new Network();
	network->SetExtraFilenameString((char*)(m_extraFilenameString.c_str()));
	network->AddTiming(network);

	PopulationColumns* layer1 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	network->AddPopulation(layer1);
	
	FullConnectivity* full = new FullConnectivity();//false,"minicolumns");
	full->SetRandomWeights(0,0);
	RandomConnectivity* randConn = new RandomConnectivity(probConnectivity);//0.1);
	randConn->SetRandomWeights(0,0);	
	//network->AddTiming(randConn);
	
	layer1->AddPre(layer1,randConn); // recurrent

	// Add Projection changes
	float lambda0 = 10e-6;
	float alpha = 0.05;
	ProjectionModifierBcpnnOnline* bStandard = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	
	//full->AddProjectionsEvent(bStandard);
	randConn->AddProjectionsEvent(bStandard);

	WTA* wta = new WTA();
	layer1->AddPopulationModifier(wta);

	// Construct initial network
	network->Initialize();
	vector<int> partsOfDataToUseAsInput = layer1->GetMPIDistribution();

	// Specify input data
	// - change to only create local part
	DataSources source;
	// not correct right now as SetValuesAll working locally and this is global
	//vector<vector<float> > data = source.GetRandomHCsOrthogonal(nrHypercolumns/100,nrRateUnits,nrItems);
	// correct
	vector<vector<float> > data = source.GetRandomHCsOrthogonalLocal(nrHypercolumns,nrRateUnits,partsOfDataToUseAsInput,nrItems);

	//cout<<"{"<<data[0].size();cout.flush();
	// Meters
	Meter* l1meter;
	Meter* c1meter;
	AnalysisDistance* analysisDistance1;
	string filename1;
	string filename2;
	string filename3;

	if(m_useBinaryFileWrite == true)
	{
		filename1 = "layer1_" + m_extraFilenameString + ".dat";
		l1meter = new Meter((char*)(filename1.c_str()), Storage::MPI_Binary);
		filename2 = "Projections1_" + m_extraFilenameString + ".dat";
		c1meter = new Meter((char*)(filename2.c_str()),Storage::MPI_Binary);
	}
	else
	{
		filename1 = "layer1_" + m_extraFilenameString + ".csv";
		l1meter = new Meter((char*)(filename1.c_str()), Storage::CSV);
		filename2 = "Projections1_" + m_extraFilenameString + ".csv";
		c1meter = new Meter((char*)(filename2.c_str()),Storage::CSV);
	}

	if(storeData)
	{
		l1meter->AttachPopulation(layer1);
		//network->AddMeter(l1meter);
		c1meter->AttachProjection(layer1->GetIncomingProjections()[0],0);
		//network->AddMeter(c1meter);
		analysisDistance1 = new AnalysisDistance(layer1,AnalysisDistance::Euclidean);
		filename3 = "Distances_"+ m_extraFilenameString + ".csv";
		analysisDistance1->SetOwnFileWriting(true,(char*)(filename3.c_str()));
		//analysisDistance1->SetCollectRepresentations(true);
		network->AddAnalysis(analysisDistance1);
	}
	
	// Timings
	network->AddTiming(bStandard);
	network->AddTiming(layer1);
	network->AddTiming(full);

	// need to access after it has been built
	network->AddTiming(layer1->GetIncomingProjections()[0]);

	//cout<<layer1->GetRateUnits().size()<<" ";cout.flush();

	// Training
	// set fixed pattern
	layer1->SwitchOnOff(false);
	if(storeData)
	{
		 analysisDistance1->SwitchOnOff(false);
	}

	int trainIterations = 1;
	int testIterations = 5;

	for(int i=0;i<trainIterations;i++)
	{
		//cout<<i<<"\n";

		for(int j=0;j<data.size();j++)
		{
			for(int k=0;k<nrSameStimuli;k++)
			{
				//layer1->SetValuesAll(data[j]);
				layer1->SetValuesAll(data[j]);
				if(k==0 && storeData==true) analysisDistance1->AddSpecificRepresentation(layer1->GetValuesLocal());

				// next time step
				network->Simulate(1);
			}
		}
	}

	// Testing
	if(storeData)
	{
		analysisDistance1->SwitchOnOff(true);
	}

	if(doTesting == true)
	{
		layer1->SwitchOnOff(true);
		bStandard->SwitchOnOff(false);

		for(int i=0;i<testIterations;i++)
		{
			for(int j=0;j<data.size();j++)
			{
					// clear all events before switching so no disturbance, can remove if moving average activity etc.
					network->ClearEventsIncoming();

					layer1->SetValuesAll(data[j]);
					for(int k=0;k<nrSameStimuli;k++)
					{
						if(k==1)
							layer1->SetValuesAll(data[j]);

						network->Simulate(1);

						//for(int i=0;i<testIterations;i++)
						// next time step
						//			network->Simulate();
					}
			}
		}
	}
	
	network->RecordAll();
	network->StoreAnalysis();
	delete network;
}

void NetworkScalingDemos::NetworkScalingRunSpiking()
{
	
}

void NetworkScalingDemos::NetworkScalingRunMIMDSVQ()
{
}