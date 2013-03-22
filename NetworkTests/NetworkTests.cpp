#include <string>
#include <sstream>

#include "NetworkBCPNN.h"
#include "NetworkKussul.h"
#include "NetworkCL.h"
#include "NetworkMDS.h"
#include "NetworkMI.h"
#include "NetworkVQ.h"
#include "NetworkCorr.h"
#include "NetworkTriesch.h"
#include "NetworkFoldiak.h"
//#include "EarlyOlfSys.h"
#include "NetworkAdaptation2.h"
#include "NetworkHebbSimple.h"
#include "NetworkBCM.h"
#include "NetworkSanger.h"
#include "Network.h"
#include "Meter.h"
#include "Storage.h"
#include "DataSources.h"
#include "NetworkTests.h"
#include <math.h>
//#include "StimulusPatterns.h"
//#include "KernelConnectivity.h"
//#include "ParamReader.h"
#include "StructureMIMDSVQ.h"

using namespace std;

NetworkTests::NetworkTests()
{
	m_architecture = PC;
}

void NetworkTests::NetworkTestBCPNNRecurrent()
{
	// network construction
	Network* network = new Network();

	PopulationColumns* layer1 = new PopulationColumns(network,256,1,PopulationColumns::Graded);
	PopulationColumns* layer2 = new PopulationColumns(network,1,5, PopulationColumns::Graded);

	FullConnectivity* full = new FullConnectivity();
	FullConnectivity* full2 = new FullConnectivity();
	FullConnectivity* full3 = new FullConnectivity();
	FullConnectivity* full4 = new FullConnectivity();

	network->AddPopulation(layer1);
	network->AddPopulation(layer2);

	layer1->AddPost(layer2,full); // Feedforward, with BCPNN
	layer2->AddPre(layer1,full2); // Feedforward

	layer2->AddPost(layer2,full3); // Recurrent, with inhibitory BCPNN
	layer2->AddPre(layer2,full4); // Recurrent

	SoftMax* softmax = new SoftMax(1.0, SoftMax::Standard);
	layer2->AddPopulationModifier(softmax);

	// Add Projection changes
	float lambda0 = 0.0001;
	float alpha = 0.01;
	ProjectionModifierBcpnnOnline* bStandard = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	ProjectionModifierBcpnnOnline* bInhib = new ProjectionModifierBcpnnOnline(alpha,lambda0);

	full2->AddProjectionsEvent(bStandard);
	full4->AddProjectionsEvent(bInhib);

	// Construct initial network
	network->Initialize();
	vector<int> partsOfDataToUseAsInput = layer1->GetUnitIdLocals();//GetMPIDistribution(network->MPIGetNodeId());

	// Specify input data

	vector<vector<float> > iData2;// = storageH5.LoadDataFloatHDF5("C:\\Projects\\Databases\\Bars\\bars1616.h5","dataset1",partsOfDataToUseAsInput,0,5);

	// Simulate
	int iterations = 20;

	for(int i=0;i<iterations;i++)
	{
		cout<<i<<"\n";

		for(int j=0;j<iData2.size();j++)
		{
			layer1->SetValuesLocal(iData2[j]);

			// next time step
			network->Simulate();
		}

		// save weights
	}
}

void NetworkTests::NetworkTestMDSVQ(int mpiRank, int mpiSize)
{
	// network construction
	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	int nrInputHypercolumns = 24;
	int nrInputRateUnits = 25;
	int nrOutputHypercolumns = 5;
	int nrOutputRateUnits = 5;

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded);
	PopulationColumns* layer2 = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::Graded);

	FullConnectivity* full = new FullConnectivity();
	FullConnectivity* full2 = new FullConnectivity();
	FullConnectivity* full3 = new FullConnectivity();
	FullConnectivity* full4 = new FullConnectivity();
	FullConnectivity* full5 = new FullConnectivity();
	FullConnectivity* full6 = new FullConnectivity();
	FullConnectivity* full7 = new FullConnectivity(true,"hypercolumn");
	FullConnectivity* full8 = new FullConnectivity(true,"hypercolumn");

	network->AddPopulation(layer1);
	network->AddPopulation(layer2);

	layer1->AddPost(layer2,full); // Feedforward, modified by VQ calculations
	layer2->AddPre(layer1,full2); // Feedforward

	layer1->AddPost(layer1,full7); // Recurrent, with MI hypercolumn calculations + MDS
	layer1->AddPre(layer1,full8); // Recurrent

	layer1->AddPost(layer1,full3); // Recurrent, with MI minicolumn calculations
	layer1->AddPre(layer1,full4); // Recurrent

	layer2->AddPost(layer2,full5); // Recurrent, with inhib bcpnn
	layer2->AddPre(layer2,full6); // Recurrent

	int mdsDimension = 3;
	int miDimension = nrInputHypercolumns;
	
	// MI
	ProjectionModifierMIHypercolumn* miHypercolumns = new ProjectionModifierMIHypercolumn();
	ProjectionModifierMIRateUnit* miRateUnits = new ProjectionModifierMIRateUnit(miHypercolumns);
	full8->AddProjectionsEvent(miHypercolumns);
	full4->AddProjectionsEvent(miRateUnits);
	//miRateUnits->AddParentProjectionModifier(miHypercolumns); // allows mi hypercolumns to have access to the belonging mi minicolumns (set as default?)

	// MDS
	LayerMDS* MDS = new LayerMDS(miDimension,mdsDimension, network);
	ProjectionModifierMDS* mdsHypercolumns = new ProjectionModifierMDS();
	layer1->AddPopulationModifier(MDS);
	mdsHypercolumns->AddParentPopulationModifier(MDS); // allows MDS to have access to the hypercolumn event Projections (will be set as default)
	full8->AddProjectionsEvent(mdsHypercolumns);

	// VQ
	int nrGroups = 10;
	LayerVQ* VQ = new LayerVQ(nrGroups, LayerVQ::VQStandard);
	layer1->AddPopulationModifier(VQ);
	VQ->AddChildPopulationModifier(MDS); // Allow VQ to have access to MDS output (m_Xi)
	//full2->AddProjectionsEvent(VQconn);

	// Inhibitory bcpnn
	float lambda0 = 0.0001;
	float alpha = 0.01;
	ProjectionModifierBcpnnOnline* bInhib = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	bInhib->SetImpactBeta(-1);
	bInhib->SetImpactWeights(-1);
	full6->AddProjectionsEvent(bInhib);

	// Construct initial network
	network->Initialize();
	vector<int> partsOfDataToUseAsInput = layer1->GetUnitIdLocals();//GetMPIDistribution(mpiRank);
	vector<int> allData;
	for(int i=0;i<nrInputHypercolumns*nrInputRateUnits;i++)
		allData.push_back(i);

	// Specify input data
	Storage storageH5;
	storageH5.SetMPIParameters(mpiRank,mpiSize);

	vector<vector<float> > iData2;// = storageH5.LoadDataFloatHDF5("C:\\Projects\\Network\\Databases\\HWBCPNN\\mtact.h5","dataset1",allData,0,25);//partsOfDataToUseAsInput,0,25);

	// Simulate
	int iterations = 10;

	for(int i=0;i<iterations;i++)
	{
		if(mpiRank == 1) cout<<i<<"\n";

		for(int j=3;j<iData2.size();j++)
		{
			if(mpiRank == 1) cout<<".";

			for(int k=0;k<3;k++)
			{
				layer1->SetValuesLocal(iData2[j]);

				// next time step
				network->Simulate();
			}
		}

		// save weights

		// non-mpi version
/*		for(int i=0;i<mpiSize;i++)
		{
			if(i==mpiRank)
			{
				vector<vector<float> > weights = layer2->GetLocalWeights();
				long localId = layer2->GetLocalUnits()[0]->GetUnitIdLocal();

				storageH5.SaveDataFloatHDF5("C:\\Projects\\Databases\\Bars\\outputWeightsLayer2.h5",localId,weights);
			}

			MPI_Barrier(NETWORK_COMM_WORLD);
		}*/
	}
}

vector<float> NetworkTests::toBinary(int nr, int total)
{
	vector<float> out(total,0.0);
	out[nr] = 1.0;

	return out;
}

vector<float> NetworkTests::toBinary(vector<float> data, int nrHc, int nrMc)
{
	vector<float> out(nrHc*nrMc);

	int currentI = 0;
	for(int i=0;i<data.size();i++)
	{
		out[currentI+data[i]] = 1.0;
		currentI+=nrMc;		
	}

	return out;
}

void NetworkTests::NetworkTestMNISTRecurrent(int mpiRank, int mpiSize)
{
	// Set up network
	int nrColors = 2;
	char* filename;

	int nrInputHypercolumns = 28*5;//28;
	int nrInputRateUnits = nrColors;
	int nrMiddleHypercolumns = 5;
	int nrMiddleRateUnits = 5;
	int nrOutputHypercolumns = 1;
	int nrOutputRateUnits = 10;

	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded); // input
	PopulationColumns* layer2 = new PopulationColumns(network,nrMiddleHypercolumns,nrMiddleRateUnits,PopulationColumns::Graded); // middle
	PopulationColumns* layer3 = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::Graded); // output

	FullConnectivity* full = new FullConnectivity();
	FullConnectivity* full2 = new FullConnectivity();
	FullConnectivity* full3 = new FullConnectivity();
	FullConnectivity* full4 = new FullConnectivity();
	FullConnectivity* full7 = new FullConnectivity(true,"hypercolumn");
	FullConnectivity* full8 = new FullConnectivity(true,"hypercolumn");
	FullConnectivity* full9 = new FullConnectivity();
	FullConnectivity* full10 = new FullConnectivity();

	network->AddPopulation(layer1);

//	layer1->AddPost(layer1,full3); // Recurrent, with MI minicolumn calculations
	layer1->AddPre(layer1,full4); // Recurrent

	int mdsDimension = 5;
	int miDimension = nrInputHypercolumns;
	
	// Recurrent bcpnn for memory test
	float lambda0 = 0.0001;
	float alpha = 0.01;
	ProjectionModifierBcpnnOnline* bRecTest = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	full4->AddProjectionsEvent(bRecTest);

	SoftMax* softmax = new SoftMax(1.0,SoftMax::Standard);	
	layer1->AddPopulationModifier(softmax);

	// Construct initial network
	network->Initialize();

	// Specify input data
	Storage storageH5;
	storageH5.SetMPIParameters(mpiRank,mpiSize);

	if(m_architecture == BGL)
	{
		if(nrColors == 2)
			filename = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_2colors.h5";
		else if(nrColors <= 0)
			filename = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST.h5";
	}
	else if(m_architecture == PC)
	{
		if(nrColors == 2)
			filename = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_2colors.h5";
		else if(nrColors <= 0)
			filename = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST.h5";
	}

	vector<int> partsOfDataToUseAsInput = layer1->GetMPIDistributionHypercolumns(mpiRank);
	vector<int> partsOfDataToUseAsOutput = vector<int>();//layer1->GetMPIDistributionHypercolumns(mpiRank);

	// Training phase

	int nrTrainImages = 60000;
	int nrTestImages = 1000;

	vector<vector<float> > trainingData;

	int iterations = 10;
	int stepsStimuliOn = 1;

	// Training phase
	layer1->SwitchOnOff(false); // fixed during training phase

	for(int j=0;j<iterations;j++)
	for(int i=0;i<nrTrainImages;i++)
	{
		if(mpiRank == 0)
		{
			cout<<i;
			cout.flush();
		}
		vector<float> binData = toBinary(trainingData[i],trainingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
		layer1->SetValuesAll(binData);

		for(int t=0;t<stepsStimuliOn;t++)
		{
			// next time step
			network->Simulate();
		}
	}
}

void NetworkTests::NetworkTestMNISTClassification(int mpiRank, int mpiSize)
{
	// Set up network
	int nrColors = 2;
	char* filenameTraining, *filenameTrainingLabels, *filenameTesting, *filenameTestingLabels;

	int nrInputHypercolumns;

	if( m_architecture == PC)
		nrInputHypercolumns = 28*5;//28;
	else if( m_architecture == BGL)
		nrInputHypercolumns = 28*28;

	int nrInputRateUnits = nrColors;
	int nrMiddleHypercolumns = 16;
	int nrMiddleRateUnits; // in each hypercolumn
	int nrOutputHypercolumns = 1;
	int nrOutputRateUnits;

	if( m_architecture == PC)
	{
		nrMiddleRateUnits = 10;
		nrOutputRateUnits = 10;
	}
	else if( m_architecture == BGL)
	{
		nrMiddleRateUnits = 256;
		nrOutputRateUnits = 256;
	}

	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded); // input
	PopulationColumns* layer2 = new PopulationColumns(network,nrMiddleHypercolumns,nrMiddleRateUnits,PopulationColumns::Graded); // middle
	PopulationColumns* layer3 = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::Graded); // output

	FullConnectivity* full2 = new FullConnectivity();
	FullConnectivity* full4 = new FullConnectivity();
	FullConnectivity* full5 = new FullConnectivity();
	FullConnectivity* full6 = new FullConnectivity();
	FullConnectivity* full8 = new FullConnectivity(true,"hypercolumn");
	FullConnectivity* full10 = new FullConnectivity();

	network->AddPopulation(layer1);
	network->AddPopulation(layer2);
	network->AddPopulation(layer3);

	//layer1->AddPost(layer1,full3); // Recurrent, with MI minicolumn calculations
	layer1->AddPre(layer1,full4); // Recurrent

	//layer1->AddPost(layer1,full7); // Recurrent, with MI hypercolumn calculations + MDS
	layer1->AddPre(layer1,full8); // Recurrent

	//layer1->AddPost(layer2,full); // Feedforward, modified by VQ calculations and BCPNN
	layer2->AddPre(layer1,full2); // Feedforward

	//layer2->AddPost(layer2,full5); // Recurrent, with inhib BCPNN
	layer2->AddPre(layer2,full6); // Recurrent

	//layer2->AddPost(layer3,full9); // Feedforward, modified by Kussul or BCPNN calculations (supervised learning for classification) 
	layer3->AddPre(layer2,full10); // Feedforward

	int mdsDimension = 5;
	int miDimension = nrInputHypercolumns;

	// MI
	ProjectionModifierMIHypercolumn* miHypercolumns = new ProjectionModifierMIHypercolumn();
	ProjectionModifierMIRateUnit* miRateUnits = new ProjectionModifierMIRateUnit(miHypercolumns);
	full8->AddProjectionsEvent(miHypercolumns);
	full4->AddProjectionsEvent(miRateUnits);
	miRateUnits->AddParentProjectionModifier(miHypercolumns); // allows mi hypercolumns to have access to the belonging mi minicolumns (set as default?)

	// MDS
	LayerMDS* MDS = new LayerMDS(miDimension,mdsDimension, network);
	ProjectionModifierMDS* mdsHypercolumns = new ProjectionModifierMDS();
	layer1->AddPopulationModifier(MDS);
	mdsHypercolumns->AddParentPopulationModifier(MDS); // allows MDS to have access to the hypercolumn event Projections (will be set as default)
	full8->AddProjectionsEvent(mdsHypercolumns);

	// VQ
	int nrGroups = nrMiddleHypercolumns;
	LayerVQ* VQ = new LayerVQ(nrGroups, LayerVQ::VQCSL);
	full2->AddProjectionsEvent(VQ->GetProjectionModifier()); // Feedforward modified by VQ calculations
	layer1->AddPopulationModifier(VQ);
	VQ->AddChildPopulationModifier(MDS); // Allow VQ to have access to MDS output (m_Xi)

	// Inhibitory bcpnn + feedforward bcpnn + softmax
	float lambda0 = 10e-8;
	float alpha = 0.05;//0.01;
	float impactBeta = -0.1;//-0.03;//-0.01/nrMiddleRateUnits;
	ProjectionModifierBcpnnOnline* bInhib = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	bInhib->SetImpactWeights(0.0);
	bInhib->SetImpactBeta(impactBeta);

	//// BGL - 40 patterns (32 nodes)
	//float lambda0 = 10e-8;
	//float alpha = 0.05;
	//ProjectionModifierBcpnnOnline* bInhib = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	//bInhib->SetImpactWeights(0.0);
	//bInhib->SetImpactBeta(-0.00001);

	// PC - 10, 10 patterns
	//	float lambda0 = 10e-8;
	//float alpha = 0.05;
	//ProjectionModifierBcpnnOnline* bInhib = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	//bInhib->SetImpactWeights(0.0);
	//bInhib->SetImpactBeta(-0.0001);

	ProjectionModifierBcpnnOnline* bFeedforward = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	//bFeedforward->SetImpactBeta(0.0);

//	full6->AddProjectionsEvent(bInhib);
//	full2->AddProjectionsEvent(bFeedforward);

	float clC = nrMiddleRateUnits*20;//*2;
	ProjectionModifierCL* compLearn = new ProjectionModifierCL(nrMiddleRateUnits,0.0001,0.0005,clC);//nrMiddleRateUnits*100);

	full2->AddProjectionsEvent(compLearn);

	SoftMax* softmax = new SoftMax(1.0, SoftMax::WTA);
	WTA* wta = new WTA();
	layer2->AddPopulationModifier(wta);//wta);//softmax);

	// Classification
	SoftMax* softmaxOutput = new SoftMax(1.0, SoftMax::WTA);
	ProjectionModifierBcpnnOnline* bClassification = new ProjectionModifierBcpnnOnline(0.05, 10e-3);
	ProjectionModifierKussul* kClassification = new ProjectionModifierKussul();
	//full10->AddProjectionsEvent(bClassification);
	full10->AddProjectionsEvent(kClassification);

	//layer3->AddPopulationModifier(softmaxOutput);

	// Construct initial network
	network->Initialize();

	// Specify input data
	Storage storageH5;
	Storage storageLabels;
	storageH5.SetMPIParameters(mpiRank,mpiSize);
	storageLabels.SetMPIParameters(mpiRank,mpiSize);

	if(m_architecture == BGL)
	{
		if(nrColors == 2)
		{
			filenameTraining = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
			filenameTesting = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingData_2colors.csv";
		}
		else if(nrColors <= 0)
			filenameTraining = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST.h5";

		filenameTrainingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingLabels.csv";
		filenameTestingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingLabels.csv";
	}
	else if(m_architecture == PC)
	{
		if(nrColors == 2)
		{
			filenameTraining = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
			filenameTesting = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingData_2colors.csv";
		}
		else if(nrColors <= 0)
			filenameTraining = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingData.csv";//MNIST.h5";

		filenameTrainingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingLabels.csv";
		filenameTestingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingLabels.csv";
	}

	vector<int> partsOfDataToUseAsInput = layer1->GetMPIDistributionHypercolumns(mpiRank);
	vector<int> partsOfDataToUseAsOutput = vector<int>();//layer1->GetMPIDistributionHypercolumns(mpiRank);

	// Training phase
	int nrTrainImages;
	int nrTestImages;

	if( m_architecture == PC)
	{
		nrTrainImages = 10;
		nrTestImages = 10;
	}
	else if( m_architecture == BGL)
	{
		nrTrainImages = 512;//10000;
		nrTestImages = 10000;
	}

	vector<vector<float> > trainingData = storageH5.LoadDataFloatCSV(filenameTraining,nrTrainImages,true);//storageH5.LoadDataFloatHDF5(filename,"trainingData",0,nrTrainImages);//storageH5.LoadDataFloatHDF5(filename,"trainingData",partsOfDataToUseAsInput,0,nrTrainImages);
	vector<vector<float> > trainingLabels = storageLabels.LoadDataFloatCSV(filenameTrainingLabels,nrTrainImages,true);

	// Recordings

	char name1[20];
	sprintf(name1,"Projection%d_2.csv",mpiRank);
	Meter connMeter(name1, Storage::CSV);
	connMeter.AttachProjection(layer2->GetIncomingProjections()[0],0);//AttachUnit(layer1->GetRateUnits()[0]);
	network->AddMeter(&connMeter);

	char* name2 = "Layer2Activity_2.csv";
	char* name3 = "Layer3Activity.csv";
	char* nameTest2 = "Layer2Activity_test.csv";
	char* nameTest3 = "Layer3Activity_test.csv";

	Meter layerMeter(name2, Storage::CSV);
	Meter layer3Meter(name3, Storage::CSV);

	layerMeter.AttachPopulation(layer2);
	network->AddMeter(&layerMeter);

	layer3Meter.AttachPopulation(layer3);
	network->AddMeter(&layer3Meter);

	Meter vqMeter("vqGroups_2.csv", Storage::CSV);
	vqMeter.AttachPopulationModifier((PopulationModifier*)VQ);//AttachUnit(layer1->GetRateUnits()[0]);

	Meter miMeter("mds.csv",Storage::CSV);
	miMeter.AttachObject((NetworkObject*)miHypercolumns, Meter::MeterPopulationModifier);

	network->AddMeter(&vqMeter);
	network->AddMeter(&miMeter);

	int iterations = 30;
	int stepsStimuliOn = 1;

	// Training phase
	// fixed will not allow RateUnits->Simulate in the layer
	// will allow PopulationModifiers(!)
	layer1->SwitchOnOff(false);	// fixed during training phase
	layer3->SwitchOnOff(false); // fixed during training phase (supervised)

	int run = 3;
	// switch off for some timesteps
	mdsHypercolumns->SwitchOnOff(false);
	MDS->SwitchOnOff(false);
	VQ->SwitchOnOff(false);
	bool areOff = true;

	if(run == 2)
	{
		miHypercolumns->SwitchOnOff(false);
		miRateUnits->SwitchOnOff(false);
	}
	else if(run == 1)
	{
		bInhib->SwitchOnOff(false);
		bFeedforward->SwitchOnOff(false);
	}
	else if(run == 3)
	{
		kClassification->SwitchOnOff(false);
		compLearn->SwitchOnOff(false);
		VQ->SwitchOnOff(false);
	}

	// Semi-sequential version

	// 1. Training phase
	// 1A. Feature extraction

	vector<float> binData;
	vector<float> binDataOut;

	for(int  j=0;j<iterations;j++)
	{
		for(int i=0;i<nrTrainImages;i++)
		{
			if(run == 3 && i == 0)
			{
				if(j == (int)(iterations*0.3))
				{
					mdsHypercolumns->SwitchOnOff(true);
					MDS->SwitchOnOff(true);

					if(mpiRank == 0)
						miMeter.RecordAll(0);
				}
				if(j == (int)(iterations*0.5))
				{
					miHypercolumns->SwitchOnOff(false);
					miRateUnits->SwitchOnOff(false);
					mdsHypercolumns->SwitchOnOff(false);
					MDS->SwitchOnOff(false);
				}
				if(j == (int)(iterations*0.5+2))
				{
					VQ->SwitchOnOff(true);

					kClassification->SwitchOnOff(true);
					compLearn->SwitchOnOff(true);
					layer3->SwitchOnOff(true);
				}
			}

			if(j == 35 && run == 1)
			{
				mdsHypercolumns->SwitchOnOff(true);
				MDS->SwitchOnOff(true);

				areOff = false;
			}

			if(j == 50 && run == 1)//50)
			{
				VQ->SwitchOnOff(true);
			}

			if(mpiRank == 0)
			{
				cout<<i<<"("<<j<<") ";
				cout.flush();
			}

			binData = toBinary(trainingData[i],trainingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
			layer1->SetValuesAll(binData);

			binDataOut = toBinary(trainingLabels[i],trainingLabels[i].size(), nrOutputRateUnits);
			layer3->SetValuesAll(binDataOut);

			if(run == 2 || run == 3)
			{
				int modSize;
				int plastStopIter = -1;
				if(m_architecture == PC)
				{
					modSize = 100;//10;
					plastStopIter = int(iterations*0.8f);
				}
				else
				{
					modSize = 100;
					plastStopIter = iterations*1;//*0.7;
				}

				if(j%modSize == 0 && j!=0 && compLearn->IsOn())
				{
					alpha = alpha/2;
					impactBeta = impactBeta/2;
					bFeedforward->SetAlpha(alpha);
					bInhib->SetImpactBeta(impactBeta);

					clC = clC*0.7f;
					compLearn->SetC(clC);

					//bInhib->SetImpactBeta(0.0);
					//softmax->SetType(false);
				}

				if(j==plastStopIter)
				{
					bFeedforward->SwitchOnOff(false);//SetAlpha(0.0);
					bInhib->SwitchOnOff(false);//SetImpactBeta(0.0);
					softmax->SetType(SoftMax::WTA);
				}
			}

			// next time step
			network->Simulate();

			cout.flush();
		}

		connMeter.RecordAll(0);

		if(mpiRank == 0)
		{
			cout<<"\n";
			cout.flush();
			layerMeter.RecordAll(0);
			layer3Meter.RecordAll(0);
		}
	}

	// save data
	if(mpiRank == 0)
		vqMeter.RecordAll(0);

	// 2. Testing phase
	//network->Simulate();
	network->Reset(); // Clears values
	//network->Simulate();
	//network->Simulate();

	layer1->SwitchOnOff(false);
	layer3->SwitchOnOff(true);
	
	bFeedforward->SwitchOnOff(false);//->SetAlpha(0.0);
	bClassification->SwitchOnOff(false);
	kClassification->SwitchOnOff(false);
	bInhib->SwitchOnOff(false);//SetImpactBeta(0.0);
	softmax->SetType(SoftMax::WTA);

	vector<vector<float> > testingData = storageH5.LoadDataFloatCSV(filenameTesting,nrTestImages,true);//storageH5.LoadDataFloatHDF5(filename,"trainingData",0,nrTrainImages);//storageH5.LoadDataFloatHDF5(filename,"trainingData",partsOfDataToUseAsInput,0,nrTrainImages);
	vector<vector<float> > testingLabels = storageLabels.LoadDataFloatCSV(filenameTestingLabels,nrTestImages,true);

	// A. Training data
	for(int i=0;i<nrTrainImages;i++)
	{
		if(mpiRank == 0)
		{
			cout<<i;
			cout.flush();
		}

		binData = toBinary(trainingData[i],trainingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
		binDataOut = toBinary(trainingLabels[i],trainingLabels[i].size(), nrOutputRateUnits);
		
		layer1->SetValuesAll(binData);
		
		network->Simulate();
	}

	network->Reset(); // Clears values

	if(mpiRank == 0)
	{
		layerMeter.RecordAll(0);
		layer3Meter.RecordAll(0);
	}

	layerMeter.SetFilename(nameTest2);
	layer3Meter.SetFilename(nameTest3);

	// B. Testing data
	for(int i=0;i<nrTestImages;i++)
	{
		if(mpiRank == 0)
		{
			cout<<i;
			cout.flush();
		}

		binData = toBinary(testingData[i],testingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
		binDataOut = toBinary(testingLabels[i],testingLabels[i].size(), nrOutputRateUnits);
		
		layer1->SetValuesAll(binData);
		
		network->Simulate();
	}

	if(mpiRank == 0)
	{
		layerMeter.SetFilename(nameTest2);
		layer3Meter.SetFilename(nameTest3);

		layerMeter.RecordAll(0);
		layer3Meter.RecordAll(0);
	}
}

void NetworkTests::NetworkTestHierarchyLayers(int mpiRank, int mpiSize)
{
	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	int nrHypercolumns = 100;
	int nrRateUnits = 100;

	PopulationColumns* layer1 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	PopulationColumns* layer2 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	PopulationColumns* layer3 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	PopulationColumns* layer4 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	PopulationColumns* layer5 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
}

void NetworkTests::NetworkTestIF(int mpiRank,int mpiSize)
{
	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	int nrHypercolumns = 2;
	int nrRateUnits = 2;

	PopulationColumns* layer1 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::IF);
	PopulationColumns* layer2 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::IF);

	FullConnectivity* full = new FullConnectivity();
	FullConnectivity* full2 = new FullConnectivity();

	FullConnectivity* full3 = new FullConnectivity();
	FullConnectivity* full4 = new FullConnectivity();
	
	network->AddPopulation(layer1);
	network->AddPopulation(layer2);

	//layer1->AddPost(layer2,full); // Feedforward layer 1 -> layer 2
	layer2->AddPre(layer1,full2);

	//layer1->AddPost(layer1,full3); // Recurrent layer 1 -> layer 1
	layer1->AddPre(layer1,full4);

	float lambda0 = 0.0001f;
	float alpha = 0.01f;
	ProjectionModifierBcpnnOnline* bRecL1 = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	bRecL1->SetImpactBeta(-1);
	bRecL1->SetImpactWeights(-1);
	full4->AddProjectionsEvent(bRecL1);

	// Construct initial network
	network->Initialize();

	// Recordings
	Meter unitMeter("unit1.csv", Storage::CSV);
	unitMeter.AttachUnit(layer1->GetRateUnits()[0]);
	network->AddMeter(&unitMeter);

	int timesteps = int(1e3);
	int nrPatterns = 4;
	vector<vector<float> > patterns(nrPatterns);
	vector<float> p1(8), p2(8), p3(8), p4(8);
	p1[0] = 1.0; p1[1] = 1.0; p1[2] = 0.0; p1[3] = 0.0; p1[4] = 1.0; p1[5] = 1.0; p1[6] = 0.0; p1[7] = 0.0; // 1100 1100
	p2[0] = 0.0; p2[1] = 1.0; p2[2] = 1.0; p2[3] = 0.0; p2[4] = 0.0; p2[5] = 1.0; p2[6] = 1.0; p2[7] = 0.0; // 0110 0110
	p3[0] = 0.0; p3[1] = 0.0; p3[2] = 1.0; p3[3] = 1.0; p3[4] = 0.0; p3[5] = 0.0; p3[6] = 1.0; p3[7] = 1.0; // 0011 0011
	p4[0] = 1.0; p4[1] = 0.0; p4[2] = 0.0; p4[3] = 1.0; p4[4] = 1.0; p4[5] = 0.0; p4[6] = 0.0; p4[7] = 1.0; // 1001 1001

	patterns[0] = p1;
	patterns[1] = p2;
	patterns[3] = p3;
	patterns[4] = p4;

	// initial training phase (separated atm)
	layer1->SwitchOnOff(false);

	int iterations = 100;
	int stepsEach = 5;
	for(int i=0;i<iterations;i++)
	{
		for(int j=0;j<stepsEach;j++)
		{
			layer1->SetValuesAll(patterns[i%nrPatterns]);
		}
	}

	// free running
	layer1->SwitchOnOff(true);
	
	for(int i=0;i<timesteps;i++)
	{
		if(i%10 == 0)
		{
			layer1->SetValuesAll(patterns[0]);
		}
		else if(i%5 == 0)
		{
			layer1->SetValuesAll(patterns[1]);
		}

		network->Simulate();
	}
}

void NetworkTests::NetworkTestPearson(int mpiRank, int mpiSize)
{
	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	int nrHypercolumns = 1;
	int nrRateUnits = 2;

	PopulationColumns* layer1 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	FullConnectivity* full = new FullConnectivity();

	network->AddPopulation(layer1);
	
	layer1->AddPre(layer1,full);

	ProjectionModifierPearson* ePearson = new ProjectionModifierPearson();
	full->AddProjectionsEvent(ePearson);

	network->Initialize();

	vector<vector<float> > trainData;

	vector<float> a1; a1.push_back(1);a1.push_back(-1);
	vector<float> a2; a2.push_back(0);a2.push_back(-1);
	vector<float> a3; a3.push_back(1);a3.push_back(-1);
	vector<float> a4; a4.push_back(0);a4.push_back(0);
	vector<float> a5; a5.push_back(2);a5.push_back(2);
	vector<float> a6; a6.push_back(3);a6.push_back(-2);

	trainData.push_back(a1);
	trainData.push_back(a2);
	trainData.push_back(a3);
	trainData.push_back(a4);
	trainData.push_back(a5);
	trainData.push_back(a6);

	unsigned int iterations = 1000;

	layer1->SwitchOnOff(false);	// fixed during training phase

	for(unsigned int j=0;j<iterations;j++)
	{
		for(unsigned int i=0;i<trainData.size();i++)
		{
			layer1->SetValuesAll(trainData[i]);
			network->Simulate();
		}
	}
}

void NetworkTests::NetworkTestTrieschAndFoldiak(int mpiRank, int mpiSize)
{
	DataSources dataSources;

	int sizeX = 5;//10;
	int sizeY = 5;//10;
	int nrItems = 2500;

	bool isTriesch = true;

	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	int nrInputHypercolumns = 1;
	int nrInputRateUnits = sizeX*sizeY;
	int nrOutputHypercolumns = 2;
	int nrOutputRateUnits = 5;//sizeX+sizeY;

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::GradedThresholded);
	PopulationColumns* layer2 = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::GradedThresholded);

	network->AddPopulation(layer1);
	network->AddPopulation(layer2);

	FullConnectivity* full = new FullConnectivity();
	FullConnectivity* full2;
	FullConnectivityNoLocalHypercolumns* full3NoLocal;

	layer2->AddPre(layer1,full);

	bool thresholded = true;
	ProjectionModifierTriesch* eTriesch = new ProjectionModifierTriesch(0.002f,0.2f,0.05f,1.0f/float(nrOutputRateUnits), thresholded);//0.05,0.2,0.005,1.0/(float)nrOutputRateUnits, thresholded);

	if(isTriesch)
		full->AddProjectionsEvent(eTriesch);

	//float eta1 = 3, eta2= 2.4, eta3 = 1.5, alpha = 0.005, beta = 200;
	float eta1 = 0.5, eta2= 0.02, eta3 = 0.02, alpha = 0.0005, beta = 10;//alpha = 1.0/8.0, beta = 10;
	bool lateral = false;

	ProjectionModifierFoldiak* eFoldiak = new ProjectionModifierFoldiak(eta1, eta2, eta3, alpha, beta, lateral);
	lateral = true;
	alpha = 0.75;
	ProjectionModifierFoldiak* eFoldiakLateral = new ProjectionModifierFoldiak(eta1, eta2, eta3, alpha, beta, lateral);
	//ProjectionModifierBCM* eBCM = new ProjectionModifierBCM(0.1,0.05,20);

	if(!isTriesch)
	{
		full2 = new FullConnectivity();
		layer2->AddPre(layer2,full2);
		full->AddProjectionsEvent(eFoldiak);
		full2->AddProjectionsEvent(eFoldiakLateral);
	}
	else
	{
		full3NoLocal = new FullConnectivityNoLocalHypercolumns();
		//full3NoLocal->AddProjectionsEvent(eBCM);
		full3NoLocal->AddProjectionsEvent(eFoldiakLateral);
		layer2->AddPre(layer2,full3NoLocal);
	}

	// implements N here
	SoftMax* softmax = new SoftMax(SoftMax::WTAThresholded,0.5);//(10.0, SoftMax::ProbWTA);
	WTA* wta = new WTA();
	//layer2->AddPopulationModifier(wta);
	layer2->AddPopulationModifier(softmax);

	network->Initialize();

	//////////////////////////////
	// Meters
	char* name1 = new char[50];
	char* name2 = new char[50];
	sprintf(name1,"Projection_triesch_n%d.csv",mpiRank);
	Meter* connMeter = new Meter(name1, Storage::CSV);
	connMeter->AttachProjection(layer2->GetIncomingProjections()[0],0);
	network->AddMeter(connMeter);

	sprintf(name2,"Layer2Activity_triesch.csv");

	Meter* layerMeter = new Meter(name2, Storage::CSV);
	layerMeter->AttachPopulation(layer2);
	network->AddMeter(layerMeter);
	// end Meters
	//////////////////////////////

	vector<vector<float> > trainData = dataSources.GetBars(sizeX,sizeY, nrItems);

	int iterations = 1;
	int iterSameStimuli = 100;

	if(!isTriesch)
		iterSameStimuli = 10;

	layer1->SwitchOnOff(false);	// fixed during training phase

	for(int j=0;j<iterations;j++)
	{
		for(int i=0;i<trainData.size();i++)
		{
			/*if(!isTriesch)
			{
				// in order to settle recurrent activity
				eFoldiak->SwitchOnOff(false);
				eFoldiakLateral->SwitchOnOff(false);
			}*/

			for(int k=0;k<iterSameStimuli;k++)
			{
			/*	if(!isTriesch && k==iterSameStimuli-1)
				{
					eFoldiak->SwitchOnOff(true);
					eFoldiakLateral->SwitchOnOff(true);
				}
*/
				for(int m=0;m<1;m++)
				{
					layer1->SetValuesAll(trainData[i]);
					//for(int n=0;n<3;n++)
					network->Simulate();
				}
			}

			// allow units to reset
			network->Reset();

			/*if(i%50 == 0)
			{
				network->RecordAll();
				if(mpiRank == 0)
					cout<<"Storing.";
			}*/
		}	
	}

	network->RecordAll();
}

// Switching
void NetworkTests::NetworkTestSwitching(int mpiRank, int mpiSize)
{
	int nrHypercolumns = 5;
	int nrRateUnits = 10;
	int nrItems = 2;

	DataSources sources;
	srand(2);
	vector<vector<float> > data = sources.GetRandomHCs(nrHypercolumns,nrRateUnits,nrItems);//sources.GetRandom(size,0.1,nrItems);
	
	// setup recurrent network

	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	PopulationColumns* layer1 = new PopulationColumns(network,nrHypercolumns,nrRateUnits,PopulationColumns::Graded);
	FullConnectivity* full = new FullConnectivity();//FullConnectivity(false,"");

	layer1->AddPre(layer1,full);
	network->AddPopulation(layer1);

	ProjectionModifierBcpnnOnline* eBcpnn = new ProjectionModifierBcpnnOnline();
	ProjectionModifierTriesch* eTriesch = new ProjectionModifierTriesch();
	ProjectionModifierHebbSimple* eHebb = new ProjectionModifierHebbSimple();
	ProjectionModifierBCM* eBCM = new ProjectionModifierBCM();

	full->AddProjectionsEvent(eBcpnn);		// incl adding transfer fcn
	//full->AddProjectionsEvent(eTriesch);		// incl adding transfer fcn
	//full->AddProjectionsEvent(eHebb);
	//full->AddProjectionsEvent(eBCM);

	PopulationModifierAdaptation2* eAdaptation = new PopulationModifierAdaptation2();
	//eAdaptation->SetParameters(0,0); 		// adaptation off initially	
	eAdaptation->SetParameters(0); 		// adaptation off initially
	layer1->AddPopulationModifier(eAdaptation);
	
	WTA* wta = new WTA();
	layer1->AddPopulationModifier(wta);//wta);//softmax);

	network->Initialize();
	eAdaptation->Initm_Aj(1); // initialize m_Aj vector

	// set up meters
	char* name1 = new char[30];
	char* name2 = new char[30];
	char* name3 = new char[30];
	sprintf(name1,"Projections_n%d.csv",mpiRank);
	sprintf(name2,"Layer1ActivityWTA.csv");
	sprintf(name3,"Layer1Activity.csv");

	Meter* connMeter = new Meter(name1, Storage::CSV);
	connMeter->AttachProjection(layer1->GetIncomingProjections()[0],0);
	network->AddMeter(connMeter);

	Meter* layerMeter = new Meter(name3, Storage::CSV);
	layerMeter->AttachPopulation(layer1);
	network->AddMeter(layerMeter);

	Meter* eventLayerMeter=new Meter(name2, Storage::CSV);
	eventLayerMeter->AttachPopulationModifier(eAdaptation);
	network->AddMeter(eventLayerMeter);

	int nrIters = 10;
	int stimuliOn = 10;

	layer1->SwitchOnOff(false); // fixed input

	// store patterns
	for(unsigned int i=0;i<nrIters;i++)
	{
		for(unsigned int j=0;j<data.size();j++)
		{
			for(unsigned int k=0;k<stimuliOn; k++)
			{
				layer1->SetValuesAll(data[j]);
				network->Simulate();
			}
		}
	}


	// random stimulation
	vector<float> randVec(data[0].size());
	for(unsigned int i=0;i<randVec.size();i++)
		randVec[i] = 0.5f*float(rand()/RAND_MAX);

	// mixture
	vector<float> mixVec(data[0].size());
	for(unsigned int i=0;i<mixVec.size();i++)
		mixVec[i] = 1*(data[0][i] + data[1][i]);

	layer1->SetValuesAll(mixVec);//randVec);

	// Test without adaptation turned on

	layer1->SwitchOnOff(true);
	//eHebb->SetEtaHebb(0.0);
	eBCM->SwitchOnOff(false);
	eBcpnn->SwitchOnOff(false);

	for(int i=0;i<nrIters;i++)
	{
		for(unsigned int j=0;j<data.size();j++)
		{
			layer1->SetValuesAll(mixVec);//data[j]);
			for(int k=0;k<stimuliOn; k++)
			{	
				network->Simulate();
			}
		}
	}

	// Turn on adaptation
	//eAdaptation->SetParameters(10,0.2);
	eAdaptation->SetParameters(0.2f);

	for(int i=0;i<nrIters;i++)
	{
		for(unsigned int j=0;j<data.size();j++)
		{
			layer1->SetValuesAll(mixVec);
			for(int k=0;k<stimuliOn; k++)
			{	
				network->Simulate();
			}
		}
	}

	network->RecordAll();

	// check switching
}

void NetworkTests::NetworkTestSanger(int mpiRank, int mpiSize)
{
	int nrInputHypercolumns = 1;
	int nrInputRateUnits = 2; // 1 input (number)
	int nrMiddleHypercolumns = 1;
	int nrMiddleRateUnits = 2; // 2 code vectors
	int nrOutputHypercolumns = 1;
	int nrOutputRateUnits = 2; // 2 classes

	// Set up network
	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	PopulationColumns* layer1InclTime = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded);
	FullConnectivity* fullMidSanger = new FullConnectivity();
	
	PopulationColumns* layerSanger = new PopulationColumns(network,nrMiddleHypercolumns,nrMiddleRateUnits,PopulationColumns::Graded);
	layerSanger->AddPre(layer1InclTime,fullMidSanger);

	// PCA/Sanger extraction
	ProjectionModifierSanger* sangerLearn = new ProjectionModifierSanger();
	fullMidSanger->AddProjectionsEvent(sangerLearn);

	network->AddPopulation(layer1InclTime);
	network->AddPopulation(layerSanger);

	network->Initialize();

	///////////////////////////////////
	//////// Data generation

	vector<vector<float> > dataIn;
	int nrTrainItems = 100000;

	for(int i=0;i<nrTrainItems;i++)
	{
		float t = (float)rand()/(float)(RAND_MAX+1);
		float y = (float)rand()/(float)(RAND_MAX+1);
		
		vector<float> f(2); f[0] = t; f[1] = y;
		dataIn.push_back(f);
	}

	////////////////////////////////////
	////// Simulate

	layer1InclTime->SwitchOnOff(false); // fixed
	for(int i=0;i<nrTrainItems;i++)
	{
		layer1InclTime->SetValuesAll(dataIn[i]);
		network->Simulate();
	}
}

std::string IntToStr( int n )
  {
  std::ostringstream result;
  result << n;
  return result.str();
  }

void NetworkTests::NetworkTestOR2ORN2MT(int mpiRank, int mpiSize)
{
	/*ParamReader* params = new ParamReader("orn_network.cfg");
	// network construction
		
	Network* network = new Network();
	int nortypes=params->getval<int>("nortypes");//20;
	//int nligands=params->getval<int>("nligands");//10;
	int or_redundancy=params->getval<int>("or_redundancy");//1; // 200 or/20 types=10 per type
	int orn_redundancy=params->getval<int>("orn_redundancy");//20; 
	int mt_redundancy=params->getval<int>("mt_redundancy");//orn_redundancy;
	int nmt=nortypes*mt_redundancy;
	int nor=nortypes*or_redundancy;
	int norn=nortypes*orn_redundancy;
	int iterations = params->getval<int>("iterations"); //20;
	float min_sens=params->getval<float>("min_sensitivity"),max_sens=params->getval<float>("max_sensitivity");
	ComplexPattern* ligands=new ComplexPattern(network,new RampPattern(network,min_sens,max_sens,int(iterations/2)));	
	ligands->addPattern(new SilentPattern(network,iterations),params->getval<int>("nligands")-1);

//	network->AddPopulation(ligands);	
	//ligands->SetNormalization(params->getval<int>( "sensor_normalization"));
	int nligands=ligands->getnligands();
	/* please use EarlyOlfSys constructor 
	ORlayerColumn* orlayer = new ORlayerColumn(network,ligands,nortypes,or_redundancy);  
	ORNlayerColumn* ornlayer = new ORNlayerColumn(network,nortypes,orn_redundancy);
	MTlayerColumn* mtlayer = new MTlayerColumn(network,nortypes,mt_redundancy,PopulationColumns::Graded);

	// connectivity
	ProjectionOlfSys* or2orn = new ProjectionOlfSys(1);	
	//ornlayer->AddPre(ligands,or2orn); 
	ornlayer->AddPre(orlayer,or2orn); 
	network->AddPopulation(orlayer);
	network->AddPopulation(ornlayer);
	
	KernelConnectivity* orn2mt = new KernelConnectivity();
	//ProjectionOlfSys* orn2mt = new ProjectionOlfSys(2);
	//FullConnectivity* orn2mt = new FullConnectivity();
	mtlayer->AddPre(ornlayer,orn2mt); 
	network->AddPopulation(mtlayer);	

	network->Initialize();

	// set up meters
	Meter* ligandMeter = new Meter("ligands.txt", Storage::CSV);
	ligandMeter->AttachPopulation(ligands);network->AddMeter(ligandMeter);
	Meter* orMeter = new Meter("receptors.txt", Storage::CSV);
	orMeter->AttachPopulation(orlayer);network->AddMeter(orMeter);
	Meter* ornMeter = new Meter("orn_activations.txt", Storage::CSV);
	ornMeter->AttachPopulation(ornlayer);network->AddMeter(ornMeter);
	Meter* mtMeter = new Meter("mt_activations.txt", Storage::CSV);
	mtMeter->AttachPopulation(mtlayer);network->AddMeter(mtMeter);	

	for(int i=0;i<iterations;i++){
		cout<<i<<"\n";
		if(i==int(iterations/2)+1)
			ligands->reset();
		network->Simulate();		
	}
	network->RecordAll();
	*/
}

void NetworkTests::NetworkTestOlfCortex(int mpiRank, int mpiSize)
{
/*	Network* network = new Network();

	ParamReader* params = new ParamReader("orn_network.cfg");
	int mt_redundancy=params->getval<int>("mt_redundancy");//orn_redundancy;

	Storage storageH5;
	vector<vector<float> > data = storageH5.LoadDataFloatCSV("n-buthanol0.1%.csv",100,true);
	vector<vector<float> > training;vector<vector<float> > test;
	storageH5.SeparateTrainingTest(data,training,test);
	Patterns* ligands = new Patterns(network,training);
	network->AddPopulation(ligands); 

	int nrligands=ligands->getnligands();
	int nmt=nrligands*mt_redundancy;
	int iterations = params->getval<int>("iterations"); 

	DataLayer* dataLayer = new DataLayer(network,nrligands,1,ligands,PopulationColumns::Graded);
	network->AddPopulation(dataLayer); 

	MTlayerColumn* mtlayer = new MTlayerColumn(network,nrligands,mt_redundancy,PopulationColumns::Graded);
	// here comes the data Projection: 
	KernelConnectivity* data2mt = new KernelConnectivity();
	//FullConnectivity* data2mt = new FullConnectivity();
	mtlayer->AddPre(dataLayer,data2mt);
	network->AddPopulation(mtlayer);

	SoftMax* softmax = new SoftMax(1.0, SoftMax::KSOFT);//KSOFT);
	mtlayer->AddPopulationModifier(softmax);

	int nrInputHypercolumns=params->getval<int>("nrInputHypercolumns");
	int nrInputRateUnits=params->getval<int>("nrInputRateUnits");
	int nrMiddleHypercolumns=params->getval<int>("nrMiddleHypercolumns");
	int nrMiddleRateUnits=params->getval<int>("nrMiddleRateUnits");

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded); // input
	StructureMIMDSVQ* structurePop1 = new StructureMIMDSVQ();
	structurePop1->SetupStructure(network,layer1,nrMiddleHypercolumns,nrMiddleRateUnits, true);

    KernelConnectivity* mt2layer1 = new KernelConnectivity();
	//FullConnectivity* data2mt = new FullConnectivity();
	mtlayer->AddPre(mtlayer,mt2layer1);
	network->AddPopulation(layer1);

	network->Initialize();

	Meter* mtMeter = new Meter("mt_activations.txt", Storage::CSV);
	mtMeter->AttachPopulation(mtlayer);network->AddMeter(mtMeter);

	Meter* dataMeter = new Meter("ligands.txt", Storage::CSV);
	dataMeter->AttachPopulation(dataLayer);network->AddMeter(dataMeter);

	structurePop1->SetupMeters(mpiRank,mpiSize);

	structurePop1->CSLLearn()->SetMaxPatterns(training.size());
	structurePop1->CSLLearn()->SwitchOnOff(false);

	for(int i=0;i<iterations;i++){
		cout<<i<<"\n";
				if(i == (int)(iterations*0.7)){
					structurePop1->MIHypercolumns()->SwitchOnOff(false);
					structurePop1->MIRateUnits()->SwitchOnOff(false);
					structurePop1->MDSHypercolumns()->SwitchOnOff(true);
					structurePop1->MDS()->SwitchOnOff(true);
				}
				if(i == (int)(iterations*0.9)){
					structurePop1->MDSHypercolumns()->SwitchOnOff(false);
					structurePop1->MDS()->SwitchOnOff(false);
					structurePop1->VQ()->SwitchOnOff(true); // moved
				}

				network->Simulate();
	}

	network->RecordAll();
	*/
}

void NetworkTests::NetworkTestInclTimingBCPNNRecurrent()
{
	int nrHypercolumns = 65536;//65536;//128*16*10;//128*16*10;//128*16*10;//128*16;
	int nrRateUnits = 100;//50;//100;//100;//40;//10
	//float activity = 0.05;
	int nrItems = 5;
	bool storeData = false;
	float probConnectivity = 0.0003;//0.00005;//0.000025;//0.00003;//0.000025;//0.001;//0.00044;//0.00011;//0.00024;//0.00006;//0.0004;//0.00001;
	bool doTesting = true; // false for some scaling tests

	// network construction
	Network* network = new Network();
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
	//vector<int> partsOfDataToUseAsInput = layer1->GetMPIDistribution(network->MPIGetNodeId());

	// Specify input data
	// - change to only create local part
	DataSources source;
	// not correct right now as SetValuesAll working locally and this is global
	vector<vector<float> > data = source.GetRandomHCsOrthogonal(nrHypercolumns/100,nrRateUnits,nrItems);

	// Meters
	Meter* l1meter = new Meter("layer1.csv", Storage::CSV);
	if(storeData)
	{
		l1meter->AttachPopulation(layer1);
		network->AddMeter(l1meter);
	}

	Meter* c1meter = new Meter("Projections1.csv",Storage::CSV);
	if(storeData)
	{
		c1meter->AttachProjection(layer1->GetIncomingProjections()[0],0);
		network->AddMeter(c1meter);
	}

	// Timings
	network->AddTiming(bStandard);
	network->AddTiming(layer1);
	network->AddTiming(full);

	// need to access after it has been built
	network->AddTiming(layer1->GetIncomingProjections()[0]);

	// Training
	// set fixed pattern
	layer1->SwitchOnOff(false);

	int trainIterations = 1;
	int testIterations = 5;

	for(int i=0;i<trainIterations;i++)
	{
		//cout<<i<<"\n";

		for(int j=0;j<data.size();j++)
		{
			layer1->SetValuesAll(data[j]);

			// next time step
			network->Simulate(10);
		}
	}

	// Testing
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
				network->Simulate(10);

				//for(int i=0;i<testIterations;i++)
				// next time step
				//			network->Simulate();
			}
		}
	}
	
	network->RecordAll();
	network->StoreAnalysis();
}

void NetworkTests::NetworkTestInclTimingMIMDSVQVisual()
{
	// Testing timing on MI+MDS+VQ on specified loaded data (same as in OlfactionScenariosClassification)
}