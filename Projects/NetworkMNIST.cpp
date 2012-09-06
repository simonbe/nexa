#include <string>
#include <sstream>
#include <iostream>

#include "NetworkBCPNN.h"
#include "NetworkKussul.h"
#include "NetworkCL.h"
#include "NetworkMDS.h"
#include "NetworkMI.h"
#include "NetworkVQ.h"
#include "Network.h"

#include "NetworkMNIST.h"
#include "NetworkTriesch.h"
#include "NetworkFoldiak.h"

//using namespace std;

void NetworkMNIST2::NetworkSetupStructure()
{
	m_architecture = CRAY;

	bool scalingRun = false;//true;
	//int nrColors = 2;

	// layer 1
	// 784 pixels in each data item
	m_nrInputHypercolumns = 784;//784;//768;//720;//784;//784;//50;//784;
	m_nrInputRateUnits = 1;//784;//nrColors;

	// layer 2
	if(m_architecture == PC)
	{
		m_nrInputHypercolumns = 50;
		m_nrHypercolumns2 = 5;//96;//48;			// number clusters/patches
		m_nrRateUnits2 = 5;//100;//200;	// number code vectors/units in each cluster
	}
	else
	{
		m_nrHypercolumns2 = 8;//24;//96;
		m_nrRateUnits2 = 50;//100;
	}

	if(scalingRun)
	{
		m_nrHypercolumns2 = 192;//360;//392;
		m_nrRateUnits2 = 20;//10;
	}

	// layer 3
	m_nrHypercolumns3 = 48;//9;		// 2nd layer cluster
	m_nrRateUnits3 = 10;//20;	// 2nd layer code vectors

	if(scalingRun)
	{
		m_nrHypercolumns3 = 192;//360;//392;
		m_nrRateUnits3 = 20;//10;//20;
	}

//	int nrPatches = 10;			// nr patches
	//m_sizePopulation2 = 100;	// nr code vectors per patch

	// Structures setting up everything to calculate mutual information (MI)/correlations (Pearson) on input data + multi-dimensional scaling (MDS) + vector quantization (VQ)
	m_structureInput = new StructureMIMDSVQ(); 
	m_structureLayer2 = new StructureMIMDSVQ();

	m_layerInput = new PopulationColumns(this,m_nrInputHypercolumns,m_nrInputRateUnits,PopulationColumns::Graded,MPIDistribution::ParallelizationHypercolumns); // input
	m_structureInput->SetMDSDimension(40); // Number of dimensions in multi-dimensional scaling matrix (Size input x dimensions)
	m_structureInput->SetMDSMeasure(true);
	m_structureLayer2->SetMDSDimension(40);
	m_structureLayer2->SetMDSMeasure(true);
	
	m_structureInput->SetupStructure(this,m_layerInput,m_nrHypercolumns2,m_nrRateUnits2, true);
	m_layer2 = m_structureInput->GetLayer(1);

	m_structureLayer2->SetupStructure(this,m_layer2,m_nrHypercolumns3,m_nrRateUnits3, false);
	m_layer3 = m_structureLayer2->GetLayer(1);

	//m_structureInput->MDSHypercolumns()->SetUseThreshold(true,0.9);
	m_structureInput->VQ()->GetCSL()->SetEta(0.001,true);
	m_structureInput->MDSHypercolumns()->SetUpdateSize(2e-3);
	m_structureLayer2->VQ()->GetCSL()->SetEta(0.001,true);
	m_structureLayer2->MDSHypercolumns()->SetUpdateSize(2e-3);

	m_structureInput->MDSHypercolumns()->SwitchOnOff(false);
	m_structureInput->MDS()->SwitchOnOff(false);
	m_structureInput->VQ()->SwitchOnOff(false);
	m_structureInput->CSLLearn()->SwitchOnOff(false);
	m_structureLayer2->MDSHypercolumns()->SwitchOnOff(false);
	m_structureLayer2->MDS()->SwitchOnOff(false);
	m_structureLayer2->VQ()->SwitchOnOff(false);
	m_structureLayer2->CSLLearn()->SwitchOnOff(false);

	AddTiming(this);
	AddTiming(m_structureInput->VQ());
	AddTiming(m_structureInput->MDSHypercolumns());
	AddTiming(m_structureInput->MDS());
	AddTiming(m_layerInput);
	AddTiming(m_layer2);
	AddTiming(m_layer3);
	AddTiming(m_structureInput->CSLLearn());
}

void NetworkMNIST2::NetworkSetupMeters()
{
	m_structureInput->SetupMeters(this->MPIGetNodeId(),this->MPIGetNrProcs());
	m_structureLayer2->SetIndex(1);
	m_structureLayer2->SetupMeters(this->MPIGetNodeId(),this->MPIGetNrProcs());

	m_inputMeter = new Meter("inputLayer.csv", Storage::CSV);
	m_inputMeter->AttachPopulation(m_layerInput);
	AddMeter(m_inputMeter);

	m_layer1Meter = new Meter("layer1.csv", Storage::CSV);
	m_layer1Meter->AttachPopulation(m_layer2);
	AddMeter(m_layer1Meter);
	
	m_layer2Meter = new Meter("layer2.csv", Storage::CSV);
	m_layer2Meter->AttachPopulation(m_layer3);
	AddMeter(m_layer2Meter);
}

void NetworkMNIST2::NetworkSetupParameters()
{
	// not used
}

/*vector<float> toBinary(int nr, int total)
{
	vector<float> out(total,0.0);
	out[nr] = 1.0;

	return out;
}*/

vector<float> NetworkMNIST2::toBinary(vector<float> data, int nrHc, int nrMc)
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

void NetworkMNIST2::TrainLayer(vector<vector<float> > trainingData, PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iterationsPatches, int iterationsFeatures)
{
	// Training phase

	int nrTrainImages = trainingData.size();

	// fixed will not allow RateUnits->Simulate in the layer
	// will allow PopulationModifiers(!)
	//structure->Layers()[0]->SwitchOnOff(false);	// fixed during training phase
	//structureInput->Layers()[2]->SwitchOnOff(false); // fixed during training phase (supervised)

	int run = 3;
	// switch off for some timesteps
	structure->MDSHypercolumns()->SwitchOnOff(false);
	structure->MDS()->SwitchOnOff(false);
	structure->VQ()->SwitchOnOff(false);
	bool areOff = true;

	//kClassification->SwitchOnOff(false);
//	structure->CompLearn()->SwitchOnOff(false);
	structure->CSLLearn()->SwitchOnOff(false);
	
	structure->VQ()->SwitchOnOff(false);

	// Semi-sequential version

	// 1. Training phase
	// 1A. Patches creation
	structure->CSLLearn()->SetMaxPatterns(nrTrainImages);
	structure->CSLLearn()->SetEta(0.001);
	
	// turn of response in 2nd layer during initial training phase for speed
	structure->GetLayer(1)->SwitchOnOff(false);
	
	vector<float> binData;
	vector<float> binDataOut;

	structure->SetRecording(false);

	for(int  j=0;j<iterationsPatches;j++)
	{
		for(int i=0;i<nrTrainImages;i++)
		{
			if(i == nrTrainImages/2)
			{
				if(j == (int)(iterationsPatches*0.6))
				{
					if(structure->UsePearson() == false)
					{
						structure->MIHypercolumns()->SwitchOnOff(false);
						structure->MIRateUnits()->SwitchOnOff(false);
					}
					else
						structure->Pearson()->SwitchOnOff(false);

					structure->MDSHypercolumns()->SwitchOnOff(true);
					structure->MDS()->SwitchOnOff(true);
					//if(mpiRank == 0)
					//	miMeter.RecordAll();
				}
				//if(j ==(int)(iterationsPatches*0.7))
					
				if(j == (int)(iterationsPatches*0.7))
				{
					structure->GetLayer(1)->SwitchOnOff(true);
					//structureInput->MIHypercolumns()->SwitchOnOff(false);
					//structureInput->MIRateUnits()->SwitchOnOff(false);
					structure->MDSHypercolumns()->SwitchOnOff(false);
					structure->MDS()->SwitchOnOff(false);

					structure->VQ()->SwitchOnOff(true); // moved
				}
				/*if(j == (int)(iterations*0.5+2))
				{
				structureInput->VQ()->SwitchOnOff(true);

				kClassification->SwitchOnOff(true);
				structureInput->CompLearn()->SwitchOnOff(true);
				layerOutput->SwitchOnOff(true);
				}*/
			}

			if(this->MPIGetNodeId() == 0)
			{
				cout<<i<<"("<<j<<") ";
				cout.flush();
			}

			//binData = toBinary(trainingData[i],trainingData[i].size(), m_nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
			//inputLayer->SetValuesAll(binData);
			inputLayer->SetValuesAll(trainingData[i]);

			//binDataOut = toBinary(trainingLabels[i],trainingLabels[i].size(), nrOutputRateUnits);
			//layerOutput->SetValuesAll(binDataOut);
			
			// next time step
			this->Simulate();
			cout.flush();
		}

		if(this->MPIGetNodeId() == 0)
			cout<<"\n";

//		this->RecordAll();

		/*connMeter.RecordAll();

		if(mpiRank == 0)
		{
		cout<<"\n";
		cout.flush();
		layerMeter.RecordAll();
		layer3Meter.RecordAll();
		}*/

		//if(j==iterationsPatches-2)
		//	structure->SetRecording(true);
	}


	this->RecordAll();

	// 1B. Feature extraction

	structure->VQ()->SwitchOnOff(false);
//	structure->CompLearn()->SwitchOnOff(true);
	

	structure->SetRecording(false);

	float clC = m_nrRateUnits2*2;
	int index=0;
	for(int  j=0;j<iterationsFeatures;j++)
	{
		//for(int i=0;i<nrTrainImages;i++)
		//{
			//binData = toBinary(trainingData[i],trainingData[i].size(), m_nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
			inputLayer->SetValuesAll(trainingData[index]);//binData);
			
			index++;
			if(index==nrTrainImages)
				index = 0;

			// next time step
			this->Simulate();

			cout.flush();

			structure->CSLLearn()->SwitchOnOff(true); // turn on first time
		//}

		int modSize;
		int plastStopIter = -1;
		if(m_architecture == PC)
		{
			modSize = 100;//10;
			plastStopIter = iterationsPatches*0.8;
		}
		else
		{
			modSize = 2;
			plastStopIter = iterationsPatches*1;//*0.7;
		}

		if(j%modSize == 0 && j!=0 && structure->CompLearn()->IsOn())
		{
			//alpha = alpha/2;
			//impactBeta = impactBeta/2;
			//bFeedforward->SetAlpha(alpha);
			//bInhib->SetImpactBeta(impactBeta);

			clC = clC*0.7;
			//structure->CompLearn()->SetC(clC);

			//bInhib->SetImpactBeta(0.0);
			//softmax->SetType(false);
		}

		if(j==plastStopIter)
		{
			//bFeedforward->SwitchOnOff(false);//SetAlpha(0.0);
			//bInhib->SwitchOnOff(false);//SetImpactBeta(0.0);
			//softmax->SetType(SoftMax::WTA);
		}

		//if(this->MPIGetNodeId()== 0)
		//	cout<<"\n";
		
		//this->RecordAll();

		//if(j==iterationsFeatures-2)
		//	structure->SetRecording(true);
	}

	this->Simulate(); // one extra time to get all training examples transmitted to feature layer

//	structure->CompLearn()->SwitchOnOff(false);
// switch off pop 1 plasticity
	structure->CSLLearn()->SwitchOnOff(false);

	this->RecordAll();
}

void NetworkMNIST2::NetworkRun()
{
	int nrColors = 2;
	char* filenameTrain, *filenameTest, *filenameTrainLabels, *filenameTestLabels;
	int nrTrainImages = 1000;//2000;//1000;
	int nrTestImages = 1000;

	int iterationsPatches = 10;//300;
	
	int iterationsFeatures = nrTrainImages + nrTrainImages;//+ 300;

	if(m_architecture == PC)
	{
		nrTrainImages = 20;
		iterationsPatches = 10;
	}

	int stepsStimuliOn = 1;

	if(m_architecture == CRAY)
	{
		if(nrColors == 2)
		{
			filenameTrain = "/cfs/klemming/nobackup/s/simonbe/Databases/MNIST/MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
			filenameTest = "/cfs/klemming/nobackup/s/simonbe/Databases/MNIST/MNIST_testingData_2colors.csv";
		}
		else if(nrColors <= 0)
			filenameTrain = "/cfs/klemming/nobackup/s/simonbe/Databases/MNIST/MNIST_trainingData.csv";//MNIST.h5";

		filenameTrainLabels = "/cfs/klemming/nobackup/s/simonbe/Databases/MNIST/MNIST_trainingLabels.csv";
		filenameTestLabels = "/cfs/klemming/nobackup/s/simonbe/Databases/MNIST/MNIST_testingLabels.csv";
	}
	else if(m_architecture == PC)
	{
		if(nrColors == 2)
		{
			filenameTrain = "D:\\Databases\\MNIST\\MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
			filenameTest = "D:\\Databases\\MNIST\\MNIST_testingData_2colors.csv";
		}
		else if(nrColors <= 0)
			filenameTrain = "D:\\Databases\\MNIST\\MNIST_trainingData.csv";//MNIST.h5";

		filenameTrainLabels = "D:\\Databases\\MNIST\\MNIST_trainingLabels.csv";
		filenameTestLabels = "D:\\Databases\\MNIST\\MNIST_testingLabels.csv";
	}

	// Specify input data
	Storage storage;
	Storage storageLabels;
	storage.SetMPIParameters(this->MPIGetNodeId(),this->MPIGetNrProcs());

	vector<vector<float> > trainingData = storage.LoadDataFloatCSV(filenameTrain,nrTrainImages,true);
	vector<vector<float> > trainingLabels = storageLabels.LoadDataFloatCSV(filenameTrainLabels,nrTrainImages,true);
	vector<vector<float> > testData = storage.LoadDataFloatCSV(filenameTest,nrTestImages,true);
	vector<vector<float> > testLabels = storageLabels.LoadDataFloatCSV(filenameTestLabels,nrTestImages,true);

	vector<int> partsOfDataToUseAsInput = m_layerInput->GetMPIDistributionHypercolumns(this->MPIGetNodeId());
	vector<int> partsOfDataToUseAsOutput = vector<int>();//layer1->GetMPIDistributionHypercolumns(mpiRank);

	m_layerInput->SwitchOnOff(false);
	m_structureInput->CSLLearn()->SetMaxPatterns(trainingData.size());
	m_structureLayer2->CSLLearn()->SetMaxPatterns(trainingData.size());
		
	// Turn off recording during training
	m_structureInput->SetRecording(false);
	m_structureLayer2->SetRecording(false);
	m_layerInput->SetRecording(false);
	m_layer2->SetRecording(false);
	m_layer3->SetRecording(false);
	
	// Train 1st layer
	m_layer2->GetIncomingProjections()[1]->SwitchOnOff(false); // turn off the Projections calculating correlations in layer 2 so (change to switch off in structure by default instead) (!)
	m_layer2->GetIncomingProjections()[2]->SwitchOnOff(false);
	m_layer3->SwitchOnOff(false);

	TrainLayer(trainingData,m_layerInput,m_structureInput,iterationsPatches,iterationsFeatures);
	
	// Run through all training and test data
	m_layer2->SwitchOnOff(true);

	// turn on all recording
	m_layerInput->SetRecording(true);
	m_layer2->SetRecording(true);
	m_layer3->SetRecording(true);

	TimingStart("RunThrough");
	for(int i=0;i<trainingData.size();i++)
	{
		m_layerInput->SetValuesAll(trainingData[i]);
		this->Simulate();
	}

	for(int i=0;i<testData.size();i++)
	{
		m_layerInput->SetValuesAll(testData[i]);
		this->Simulate();
	}

	TimingStop("RunThrough");
	// Train 2nd layer
	m_layer3->SwitchOnOff(true);
	//TrainLayer(trainingData,m_layerInput,m_structureLayer2,iterationsPatches,iterationsFeatures);

	// save data
	//if(mpiRank == 0)
	//	vqMeter.RecordAll();

	//network->Simulate();
	//this->Reset(); // Clears values
	//network->Simulate();
	//network->Simulate();

	// clear activities

	// run training and test sets

}

/*
NetworkMNIST::NetworkMNIST()
{
	m_architecture = PC;
}

vector<float> NetworkMNIST::toBinary(int nr, int total)
{
	vector<float> out(total,0.0);
	out[nr] = 1.0;

	return out;
}

vector<float> NetworkMNIST::toBinary(vector<float> data, int nrHc, int nrMc)
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

void NetworkMNIST::NetworkMNISTRun1(int mpiRank, int mpiSize)
{
	int runNrLayers = 3;
	// Set up network
	int nrColors = 2;
	char* filenameTraining, *filenameTrainingLabels, *filenameTesting, *filenameTestingLabels;

	// Input structure
	StructureMIMDSVQ* structureInput = new StructureMIMDSVQ();
	StructureMIMDSVQ* structurePop2 = new StructureMIMDSVQ();
	StructureMIMDSVQ* structurePop3 = new StructureMIMDSVQ();

	int nrInputHypercolumns;

	if( m_architecture == PC)
		nrInputHypercolumns = 28*10;//28;
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
		nrMiddleRateUnits = 90; // check so nr of nodes not too many (too few mcs since it's spread even atm).
		nrOutputRateUnits = 90;
	}

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

	// Training phase
	int nrTrainImagesPatches;
	int nrTrainImagesFeatures;
	int nrTestImages;

	if( m_architecture == PC)
	{
		nrTrainImagesPatches = 15;
		nrTrainImagesFeatures = nrTrainImagesPatches;//nrMiddleRateUnits;//20;
		nrTestImages = 30;
	}
	else if( m_architecture == BGL)
	{
		nrTrainImagesPatches = 1000;//10000;
		nrTrainImagesFeatures = nrTrainImagesPatches/3;//nrMiddleRateUnits;//256;
		nrTestImages = 10000;
	}

	vector<vector<float> > trainingData = storageH5.LoadDataFloatCSV(filenameTraining,nrTrainImagesPatches,true);//storageH5.LoadDataFloatHDF5(filename,"trainingData",0,nrTrainImages);//storageH5.LoadDataFloatHDF5(filename,"trainingData",partsOfDataToUseAsInput,0,nrTrainImages);
	vector<vector<float> > trainingLabels = storageLabels.LoadDataFloatCSV(filenameTrainingLabels,nrTrainImagesPatches,true);

	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded); // input
	structureInput->SetupStructure(network,layer1,nrMiddleHypercolumns,nrMiddleRateUnits, true);//,nrOutputHypercolumns,nrOutputRateUnits);
	
	network->Initialize();
	structureInput->SetupMeters(mpiRank,mpiSize);

	// Classification
	PopulationColumns* layerOutput = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::Graded); // output
	bool useClassification = false;
	FullConnectivity* full10 = new FullConnectivity();
	//layer2->AddPost(layer3,full9); // Feedforward, modified by Kussul or BCPNN calculations (supervised learning for classification) 
	if(useClassification)
		layerOutput->AddPre(structureInput->Layers()[1],full10); // Feedforward

	SoftMax* softmaxOutput = new SoftMax(1.0, SoftMax::WTA);
	ProjectionModifierBcpnnOnline* bClassification = new ProjectionModifierBcpnnOnline(0.05, 10e-3);
	ProjectionModifierKussul* kClassification = new ProjectionModifierKussul();
	
	if(useClassification)
		//full10->AddProjectionsEvent(bClassification);
		full10->AddProjectionsEvent(kClassification);
	// Recordings

	int iterationsPatches = 15;
	int iterationsFeatures = 10;
	int stepsStimuliOn = 1;

	vector<int> partsOfDataToUseAsInput = layer1->GetMPIDistributionHypercolumns(mpiRank);
	vector<int> partsOfDataToUseAsOutput = vector<int>();//layer1->GetMPIDistributionHypercolumns(mpiRank);

	// Training phase
	// fixed will not allow RateUnits->Simulate in the layer
	// will allow PopulationModifiers(!)
	structureInput->Layers()[0]->SwitchOnOff(false);	// fixed during training phase
	//structureInput->Layers()[2]->SwitchOnOff(false); // fixed during training phase (supervised)

	int run = 3;
	// switch off for some timesteps
	structureInput->MDSHypercolumns()->SwitchOnOff(false);
	structureInput->MDS()->SwitchOnOff(false);
	structureInput->VQ()->SwitchOnOff(false);
	bool areOff = true;

	kClassification->SwitchOnOff(false);
	structureInput->CompLearn()->SwitchOnOff(false);
	structureInput->CSLLearn()->SwitchOnOff(false);
	
	structureInput->VQ()->SwitchOnOff(false);

	// Semi-sequential version

	// 1. Training phase
	// 1A. Patches creation
	structureInput->CSLLearn()->SetMaxPatterns(nrTrainImagesPatches);

	vector<float> binData;
	vector<float> binDataOut;

	structureInput->SetRecording(false);

	for(int  j=0;j<iterationsPatches;j++)
	{
		for(int i=0;i<nrTrainImagesPatches;i++)
		{
			if(i == 0)
			{
				if(j == (int)(iterationsPatches*0.7))
				{
					structureInput->MIHypercolumns()->SwitchOnOff(false);
					structureInput->MIRateUnits()->SwitchOnOff(false);
					structureInput->MDSHypercolumns()->SwitchOnOff(true);
					structureInput->MDS()->SwitchOnOff(true);

					//if(mpiRank == 0)
					//	miMeter.RecordAll();
				}
				if(j == (int)(iterationsPatches*0.9))
				{
					//structureInput->MIHypercolumns()->SwitchOnOff(false);
					//structureInput->MIRateUnits()->SwitchOnOff(false);
					structureInput->MDSHypercolumns()->SwitchOnOff(false);
					structureInput->MDS()->SwitchOnOff(false);

					structureInput->VQ()->SwitchOnOff(true); // moved
				}
				//if(j == (int)(iterations*0.5+2))
				//{
				//structureInput->VQ()->SwitchOnOff(true);
				//
//				kClassification->SwitchOnOff(true);
//				structureInput->CompLearn()->SwitchOnOff(true);
//				layerOutput->SwitchOnOff(true);
//				}
			}

			if(mpiRank == 0)
			{
				cout<<i<<"("<<j<<") ";
				cout.flush();
			}

			binData = toBinary(trainingData[i],trainingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
			layer1->SetValuesAll(binData);

			//binDataOut = toBinary(trainingLabels[i],trainingLabels[i].size(), nrOutputRateUnits);
			//layerOutput->SetValuesAll(binDataOut);
			
			// next time step
			network->Simulate();
			cout.flush();
		}

		if(mpiRank == 0)
			cout<<"\n";

		network->RecordAll();

		//connMeter.RecordAll();
		//
		//if(mpiRank == 0)
		//{
		//cout<<"\n";
		//cout.flush();
		//layerMeter.RecordAll();
		//layer3Meter.RecordAll();
		//}

		if(j==iterationsPatches-2)
			structureInput->SetRecording(true);
	}

	// 1B. Feature extraction

	structureInput->VQ()->SwitchOnOff(false);
	structureInput->CompLearn()->SwitchOnOff(true);
	structureInput->CSLLearn()->SwitchOnOff(true);

	structureInput->SetRecording(false);

	float clC = nrMiddleRateUnits*2;

	for(int  j=0;j<iterationsFeatures;j++)
	{
		for(int i=0;i<nrTrainImagesFeatures;i++)
		{
			binData = toBinary(trainingData[i],trainingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
			layer1->SetValuesAll(binData);

			// next time step
			network->Simulate();

			cout.flush();
		}

		int modSize;
		int plastStopIter = -1;
		if(m_architecture == PC)
		{
			modSize = 100;//10;
			plastStopIter = iterationsPatches*0.8;
		}
		else
		{
			modSize = 2;
			plastStopIter = iterationsPatches*1;//*0.7;
		}

		if(j%modSize == 0 && j!=0 && structureInput->CompLearn()->IsOn())
		{
			//alpha = alpha/2;
			//impactBeta = impactBeta/2;
			//bFeedforward->SetAlpha(alpha);
			//bInhib->SetImpactBeta(impactBeta);

			clC = clC*0.7;
			structureInput->CompLearn()->SetC(clC);

			//bInhib->SetImpactBeta(0.0);
			//softmax->SetType(false);
		}

		if(j==plastStopIter)
		{
			//bFeedforward->SwitchOnOff(false);//SetAlpha(0.0);
			//bInhib->SwitchOnOff(false);//SetImpactBeta(0.0);
			//softmax->SetType(SoftMax::WTA);
		}

		if(mpiRank == 0)
			cout<<"\n";
		
		network->RecordAll();

		if(j==iterationsFeatures-2)
			structureInput->SetRecording(true);
	}

	structureInput->CompLearn()->SwitchOnOff(false);
	structureInput->CSLLearn()->SwitchOnOff(false);

	// switch off pop 1 plasticity

	// save data
	//if(mpiRank == 0)
	//	vqMeter.RecordAll();

	//network->Simulate();
	network->Reset(); // Clears values
	//network->Simulate();
	//network->Simulate();

	// 2A. setup layer 2

	if(runNrLayers>1)
	{
		PopulationColumns* layerOut1 = structureInput->Layers()[1];
		structurePop2->SetupStructure(network,layerOut1,nrMiddleHypercolumns/2,nrMiddleRateUnits, false);//,nrOutputHypercolumns/2,nrOutputRateUnits*2);

		network->Initialize(); // initializes new parts
		structurePop2->SetIndex(1); // changes output filenames
		structurePop2->SetupMeters(mpiRank,mpiSize);

		// run with changes to layer 2

		// switch off for some timesteps
		structurePop2->MDSHypercolumns()->SwitchOnOff(false);
		structurePop2->MDS()->SwitchOnOff(false);
		structurePop2->VQ()->SwitchOnOff(false);

		kClassification->SwitchOnOff(false);
		structurePop2->CompLearn()->SwitchOnOff(false);
		structurePop2->CSLLearn()->SwitchOnOff(false);
		structurePop2->VQ()->SwitchOnOff(false);

		structureInput->SetRecording(false);
		structurePop2->SetRecording(false);

		structurePop2->CSLLearn()->SetMaxPatterns(nrTrainImagesFeatures);

		for(int  j=0;j<iterationsPatches;j++)
		{
			for(int i=0;i<nrTrainImagesPatches;i++)
			{
				if(i == 0)
				{
					if(j == (int)(iterationsPatches*0.7))
					{
						structurePop2->MIHypercolumns()->SwitchOnOff(false);
						structurePop2->MIRateUnits()->SwitchOnOff(false);
						structurePop2->MDSHypercolumns()->SwitchOnOff(true);
						structurePop2->MDS()->SwitchOnOff(true);

						//if(mpiRank == 0)
						//	miMeter.RecordAll();
					}
					if(j == (int)(iterationsPatches*0.9))
					{
						//structurePop2->MIHypercolumns()->SwitchOnOff(false);
						//structurePop2->MIRateUnits()->SwitchOnOff(false);
						structurePop2->MDSHypercolumns()->SwitchOnOff(false);
						structurePop2->MDS()->SwitchOnOff(false);

						structurePop2->VQ()->SwitchOnOff(true); // moved
					}
					//if(j == (int)(iterations*0.5+2))
					//{
					//structureInput->VQ()->SwitchOnOff(true);
					//
					//kClassification->SwitchOnOff(true);
					//structureInput->CompLearn()->SwitchOnOff(true);
					//layerOutput->SwitchOnOff(true);
					//}
				}

				if(mpiRank == 0)
				{
					cout<<i<<"("<<j<<") ";
					cout.flush();
				}

				binData = toBinary(trainingData[i],trainingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
				layer1->SetValuesAll(binData);

				//binDataOut = toBinary(trainingLabels[i],trainingLabels[i].size(), nrOutputRateUnits);
				//layerOutput->SetValuesAll(binDataOut);

				int modSize;
				int plastStopIter = -1;
				if(m_architecture == PC)
				{
					modSize = 100;//10;
					plastStopIter = iterationsPatches*0.8;
				}
				else
				{
					modSize = 100;
					plastStopIter = iterationsPatches*1;//*0.7;
				}

				if(j%modSize == 0 && j!=0 && structureInput->CompLearn()->IsOn())
				{
					//alpha = alpha/2;
					//impactBeta = impactBeta/2;
					//bFeedforward->SetAlpha(alpha);
					//bInhib->SetImpactBeta(impactBeta);

					clC = clC*0.7;
					structurePop2->CompLearn()->SetC(clC);

					//bInhib->SetImpactBeta(0.0);
					//softmax->SetType(false);
				}

				if(j==plastStopIter)
				{
					//bFeedforward->SwitchOnOff(false);//SetAlpha(0.0);
					//bInhib->SwitchOnOff(false);//SetImpactBeta(0.0);
					//softmax->SetType(SoftMax::WTA);
				}

				// next time step
				network->Simulate();
				cout.flush();
			}

			if(mpiRank == 0)
				cout<<"\n";

			network->RecordAll();

			if(j == iterationsPatches-2)
			{
				structureInput->SetRecording(true);
				structurePop2->SetRecording(true);
			}

			//connMeter.RecordAll();
//
	//		if(mpiRank == 0)
	//		{
	//		cout<<"\n";
	//		cout.flush();
	//		layerMeter.RecordAll();
	//		layer3Meter.RecordAll();
	//		}
		}

		// 2B. Feature extraction

		structureInput->SetRecording(false);
		structurePop2->SetRecording(false);

		structurePop2->VQ()->SwitchOnOff(false);
		structurePop2->CompLearn()->SwitchOnOff(true);
		structurePop2->CSLLearn()->SwitchOnOff(true);

		for(int  j=0;j<iterationsFeatures;j++)
		{
			for(int i=0;i<nrTrainImagesFeatures;i++)
			{
				binData = toBinary(trainingData[i],trainingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
				layer1->SetValuesAll(binData);

				// next time step
				network->Simulate();

				cout.flush();
			}

			if(mpiRank == 0)
				cout<<"\n";

			network->RecordAll();

			if(j == iterationsFeatures-2)
			{
				structureInput->SetRecording(true);
				structurePop2->SetRecording(true);
			}
		}

		structureInput->CompLearn()->SwitchOnOff(false);
		structureInput->CSLLearn()->SwitchOnOff(false);

		structureInput->SetRecording(false);
		structurePop2->SetRecording(false);

		// 3A. Setup layer 3
		if(runNrLayers>2)
		{
			PopulationColumns* layerOut2 = structurePop2->Layers()[1];
			structurePop3->SetupStructure(network,layerOut2,nrMiddleHypercolumns/4,nrMiddleRateUnits, false);//,nrOutputHypercolumns/2,nrOutputRateUnits*2);

			network->Initialize(); // initializes new parts
			structurePop3->SetIndex(2); // changes output filenames
			structurePop3->SetupMeters(mpiRank,mpiSize);

			// run with changes to layer 3

			// switch off for some timesteps
			structurePop3->MDSHypercolumns()->SwitchOnOff(false);
			structurePop3->MDS()->SwitchOnOff(false);
			structurePop3->VQ()->SwitchOnOff(false);

			kClassification->SwitchOnOff(false);
			structurePop3->CompLearn()->SwitchOnOff(false);
			structurePop3->CSLLearn()->SwitchOnOff(false);
			structurePop3->VQ()->SwitchOnOff(false);

			structureInput->SetRecording(false);
			structurePop2->SetRecording(false);
			structurePop3->SetRecording(false);

			structurePop3->CSLLearn()->SetMaxPatterns(nrTrainImagesFeatures);

			for(int  j=0;j<iterationsPatches;j++)
			{
				for(int i=0;i<nrTrainImagesPatches;i++)
				{
					if(i == 0)
					{
						if(j == (int)(iterationsPatches*0.7))
						{
							structurePop3->MIHypercolumns()->SwitchOnOff(false);
							structurePop3->MIRateUnits()->SwitchOnOff(false);
							structurePop3->MDSHypercolumns()->SwitchOnOff(true);
							structurePop3->MDS()->SwitchOnOff(true);

							//if(mpiRank == 0)
							//	miMeter.RecordAll();
						}
						if(j == (int)(iterationsPatches*0.9))
						{
							//structurePop2->MIHypercolumns()->SwitchOnOff(false);
							//structurePop2->MIRateUnits()->SwitchOnOff(false);
							structurePop3->MDSHypercolumns()->SwitchOnOff(false);
							structurePop3->MDS()->SwitchOnOff(false);

							structurePop3->VQ()->SwitchOnOff(true); // moved
						}
						//if(j == (int)(iterations*0.5+2))
						//{
						//structureInput->VQ()->SwitchOnOff(true);
						//
						//kClassification->SwitchOnOff(true);
						//structureInput->CompLearn()->SwitchOnOff(true);
						//layerOutput->SwitchOnOff(true);
						//}
					}

					if(mpiRank == 0)
					{
						cout<<i<<"("<<j<<") ";
						cout.flush();
					}

					binData = toBinary(trainingData[i],trainingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
					layer1->SetValuesAll(binData);

					//binDataOut = toBinary(trainingLabels[i],trainingLabels[i].size(), nrOutputRateUnits);
					//layerOutput->SetValuesAll(binDataOut);

					int modSize;
					int plastStopIter = -1;
					if(m_architecture == PC)
					{
						modSize = 100;//10;
						plastStopIter = iterationsPatches*0.8;
					}
					else
					{
						modSize = 100;
						plastStopIter = iterationsPatches*1;//*0.7;
					}

					if(j%modSize == 0 && j!=0 && structureInput->CompLearn()->IsOn())
					{
						//alpha = alpha/2;
						//impactBeta = impactBeta/2;
						//bFeedforward->SetAlpha(alpha);
						//bInhib->SetImpactBeta(impactBeta);

						clC = clC*0.7;
						structurePop3->CompLearn()->SetC(clC);

						//bInhib->SetImpactBeta(0.0);
						//softmax->SetType(false);
					}

					if(j==plastStopIter)
					{
						//bFeedforward->SwitchOnOff(false);//SetAlpha(0.0);
						//bInhib->SwitchOnOff(false);//SetImpactBeta(0.0);
						//softmax->SetType(SoftMax::WTA);
					}

					// next time step
					network->Simulate();
					cout.flush();
				}

				if(mpiRank == 0)
					cout<<"\n";

				network->RecordAll();

				if(j == iterationsPatches-2)
				{
					structureInput->SetRecording(true);
					structurePop2->SetRecording(true);
					structurePop3->SetRecording(true);
				}

				//connMeter.RecordAll();
				//
				//if(mpiRank == 0)
				//{
				//cout<<"\n";
				//cout.flush();
				//layerMeter.RecordAll();
				//layer3Meter.RecordAll();
				//}
			}

			// 3B. Feature extraction

			structureInput->SetRecording(false);
			structurePop2->SetRecording(false);
			structurePop3->SetRecording(false);

			structurePop3->VQ()->SwitchOnOff(false);
			structurePop3->CompLearn()->SwitchOnOff(true);
			structurePop3->CSLLearn()->SwitchOnOff(true);

			for(int  j=0;j<iterationsFeatures;j++)
			{
				for(int i=0;i<nrTrainImagesFeatures;i++)
				{
					binData = toBinary(trainingData[i],trainingData[i].size(), nrInputRateUnits);//vector<float> binData = toBinary(currentTrainingData[0],currentTrainingData[0].size(),nrInputRateUnits);//trainingData[i],trainingData[i].size(), nrInputRateUnits);
					layer1->SetValuesAll(binData);

					// next time step
					network->Simulate();

					cout.flush();
				}

				if(mpiRank == 0)
					cout<<"\n";

				network->RecordAll();

				if(j == iterationsFeatures-2)
				{
					structureInput->SetRecording(true);
					structurePop2->SetRecording(true);
					structurePop3->SetRecording(true);
				}
			}

			structureInput->CompLearn()->SwitchOnOff(false);
			structureInput->CSLLearn()->SwitchOnOff(false);

			structureInput->SetRecording(false);
			structurePop2->SetRecording(false);
			structurePop3->SetRecording(false);
		}

		// 4. Testing set
		// perform pattern completion on partial input data
	}
}

void NetworkMNIST::NetworkMNISTRun2(int mpiRank, int mpiSize)
{
	int sizeX = 28;
	int sizeY = 28;
	int nrItems = 30000;

	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	int nrInputHypercolumns = 1;
	int nrInputRateUnits = sizeX*sizeY;
	int nrOutputHypercolumns = 1;
	int nrOutputRateUnits = 128;//sizeX+sizeY;

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded);
	PopulationColumns* layer2 = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::Graded);

	network->AddPopulation(layer1);
	network->AddPopulation(layer2);

	FullConnectivity* full = new FullConnectivity();
	layer2->AddPre(layer1,full);

	ProjectionModifierTriesch* eTriesch = new ProjectionModifierTriesch(0.05,0.2,0.005,1.0/(float)nrOutputRateUnits, false);//0.1);//(0.01,0.0,0.01,0.1);//1/500);//(float)nrOutputRateUnits);//0.01);
	full->AddProjectionsEvent(eTriesch);

	// implements N here
	WTA* wta = new WTA();
	layer2->AddPopulationModifier(wta);

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


	//////////////////////////////
	///// Input data

	char* filenameTraining, *filenameTesting, *filenameTrainingLabels, *filenameTestingLabels;

	if(m_architecture == BGL)
	{
		filenameTraining = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
		filenameTesting = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingData_2colors.csv";

		filenameTrainingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingLabels.csv";
		filenameTestingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingLabels.csv";
	}
	else if(m_architecture == PC)
	{
		filenameTraining = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
		filenameTesting = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingData_2colors.csv";

		filenameTrainingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingLabels.csv";
		filenameTestingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingLabels.csv";
	}

	int iterations = 1;

	if( m_architecture == PC)
	{
		iterations = 1;
	}
	else if( m_architecture == BGL)
	{
		iterations = 10;
	}

	Storage storage;
	DataSources dataSources;

	vector<vector<float> > trainingData = storage.LoadDataFloatCSV(filenameTraining,nrItems,false); //dataSources.GetBars(sizeX,sizeY, nrItems);
	vector<vector<float> > trainingLabels = storage.LoadDataFloatCSV(filenameTrainingLabels,nrItems,false);

	//float beta = 0;
	//float mu = 0.2;
	float betaMin = 0.2;//0;
	float betaMax = 0.25, betaIt = 0.05;
	float muMin = 0.1;//0.005;
	float muMax = 0.11;//0.05;
	float muIt = 0.01;

	layer1->SwitchOnOff(false);	// fixed during training phase

	for(float beta = betaMin; beta<betaMax;beta += betaIt)
	{
		eTriesch->SetBeta(beta);

		for(float mu=muMin; mu<muMax; mu+= muIt)
		{
			eTriesch->Clear();
			eTriesch->SetMu(mu);

			for(int j=0;j<iterations;j++)
			{
				for(int i=0;i<trainingData.size();i++)
				{
					layer1->SetValuesAll(trainingData[i]);
					network->Simulate();

					if(mpiRank == 0)
					if(i%100 == 0)
					{
						cout<<"\n"<<i<<" ";
						cout.flush();
					}
				}
			}	
		}

		if(mpiRank == 0)
		{
			cout<<"\n("<<beta<<") "; cout.flush();
		}
	}

	//if(j%30 == 0)	
	network->RecordAll();
}

void NetworkMNIST::NetworkMNISTRun3(int mpiRank, int mpiSize)
{
	int sizeX = 28;
	int sizeY = 28;
	int nrItems = 2000;

	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	int nrInputHypercolumns = 1;
	int nrInputRateUnits = sizeX*sizeY;
	int nrOutputHypercolumns = 2;
	int nrOutputRateUnits = 128;//sizeX+sizeY;

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded);
	PopulationColumns* layer2 = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::Graded);

	network->AddPopulation(layer1);
	network->AddPopulation(layer2);

	// feedforward feature extraction
	FullConnectivity* full = new FullConnectivity();
	layer2->AddPre(layer1,full);

	ProjectionModifierTriesch* eTriesch = new ProjectionModifierTriesch(0.05,0.2,0.005,1.0/(float)nrOutputRateUnits, false);//0.1);//(0.01,0.0,0.01,0.1);//1/500);//(float)nrOutputRateUnits);//0.01);
	full->AddProjectionsEvent(eTriesch);

	// inhibitory Projections
	FullConnectivityNoLocalHypercolumns* full2 = new FullConnectivityNoLocalHypercolumns();
	layer2->AddPre(layer2,full2);

	float lambda0 = 0.0001;
	float alpha = 0.01;
	ProjectionModifierBcpnnOnline* bInhib = new ProjectionModifierBcpnnOnline(alpha,lambda0);
	bInhib->SetImpactBeta(-1);
	bInhib->SetImpactWeights(-1);
	full2->AddProjectionsEvent(bInhib);

	// implements N here
	WTA* wta = new WTA();
	layer2->AddPopulationModifier(wta);

	network->Initialize();

	//////////////////////////////
	// Meters
	char* name1 = new char[50];
	char* name2 = new char[50];
	char* name3 = new char[80];
	sprintf(name1,"Projection_triesch_n%d.csv",mpiRank);
	Meter* connMeter = new Meter(name1, Storage::CSV);
	connMeter->AttachProjection(layer2->GetIncomingProjections()[0],0);
	network->AddMeter(connMeter);

	sprintf(name2,"Layer2Activity_triesch.csv");

	Meter* layerMeter = new Meter(name2, Storage::CSV);
	layerMeter->AttachPopulation(layer2);
	network->AddMeter(layerMeter);

	sprintf(name3,"Projection_triesch_inhib_n%d.csv",mpiRank);
	Meter* connMeter2 = new Meter(name3, Storage::CSV);
	connMeter2->AttachProjection(layer2->GetIncomingProjections()[1],0);
	network->AddMeter(connMeter2);
	// end Meters
	//////////////////////////////


	//////////////////////////////
	///// Input data

	char* filenameTraining, *filenameTesting, *filenameTrainingLabels, *filenameTestingLabels;

	if(m_architecture == BGL)
	{
		filenameTraining = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
		filenameTesting = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingData_2colors.csv";

		filenameTrainingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingLabels.csv";
		filenameTestingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingLabels.csv";
	}
	else if(m_architecture == PC)
	{
		filenameTraining = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
		filenameTesting = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingData_2colors.csv";

		filenameTrainingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingLabels.csv";
		filenameTestingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingLabels.csv";
	}

	int iterations = 1;

	if( m_architecture == PC)
	{
		iterations = 1;
	}
	else if( m_architecture == BGL)
	{
		iterations = 1;
	}

	Storage storage;
	DataSources dataSources;

	vector<vector<float> > trainingData = storage.LoadDataFloatCSV(filenameTraining,nrItems,false); //dataSources.GetBars(sizeX,sizeY, nrItems);
	vector<vector<float> > trainingLabels = storage.LoadDataFloatCSV(filenameTrainingLabels,nrItems,false);

	//float beta = 0;
	//float mu = 0.2;
	float betaMin = 0.2;//0;
	float betaMax = 0.25, betaIt = 0.05;
	float muMin = 0.1;//0.005;
	float muMax = 0.11;//0.05;
	float muIt = 0.01;

	layer1->SwitchOnOff(false);	// fixed during training phase

	for(float beta = betaMin; beta<betaMax;beta += betaIt)
	{
		eTriesch->SetBeta(beta);

		for(float mu=muMin; mu<muMax; mu+= muIt)
		{
			eTriesch->Clear();
			eTriesch->SetMu(mu);

			for(int j=0;j<iterations;j++)
			{
				for(int i=0;i<trainingData.size();i++)
				{
					wta->SwitchOnOff(false);

					for(int k=0;k<40;k++)
					{
						layer1->SetValuesAll(trainingData[i]);
						
						if(k==39)
							wta->SwitchOnOff(true);

						network->Simulate();
					}

					if(mpiRank == 0)
					if(i%100 == 0)
					{
						cout<<"\n"<<i<<" ";
						cout.flush();
					}
				}
			}	
		}

		if(mpiRank == 0)
		{
			cout<<"\n("<<beta<<") "; cout.flush();
		}
	}

	//if(j%30 == 0)	
	network->RecordAll();
}

void NetworkMNIST::NetworkMNISTRunLateralMaps(int mpiRank, int mpiSize)
{
	int sizeX = 28;
	int sizeY = 28;
	int nrItems = 10;

	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	int nrInputHypercolumns = 1;
	int nrInputRateUnits = sizeX*sizeY;
	int nrOutputHypercolumns = 1;
	int nrOutputRateUnits = 128;//sizeX+sizeY;

	if(m_architecture == PC)
	{
		nrOutputRateUnits = 10;
	}

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded);
	PopulationColumns* layer2 = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::Graded);

	network->AddPopulation(layer1);
	network->AddPopulation(layer2);

	// feedforward feature extraction
	FullConnectivity* full = new FullConnectivity();
	layer2->AddPre(layer1,full);

	float eta = 0.8, decay = 0.1, tau = 100;
	ProjectionModifierBCM* eBCM = new ProjectionModifierBCM(eta,decay,tau);
	full->AddProjectionsEvent(eBCM);

	// inhibitory Projections
	FullConnectivityNoLocalHypercolumns* full2 = new FullConnectivityNoLocalHypercolumns();
	layer2->AddPre(layer2,full2);

	//float eta2 = 0.4, decay2 = 0.05, tau2 = 10;
	//ProjectionModifierBCM* eBCM2 = new ProjectionModifierBCM();//(eta,decay,tau);
	//full2->AddProjectionsEvent(eBCM2);

	// implements N here
	WTA* wta = new WTA();
	SoftMax* soft = new SoftMax(2.0,SoftMax::ProbWTA);
	layer2->AddPopulationModifier(soft);
	//layer2->AddPopulationModifier(wta);

	network->Initialize();

	//////////////////////////////
	// Meters
	char* name1 = new char[50];
	char* name2 = new char[50];
	char* name3 = new char[80];
	sprintf(name1,"Projection_bcm_n%d.csv",mpiRank);
	Meter* connMeter = new Meter(name1, Storage::CSV);
	connMeter->AttachProjection(layer2->GetIncomingProjections()[0],0);
	network->AddMeter(connMeter);

	sprintf(name2,"Layer2Activity_bcm.csv");

	Meter* layerMeter = new Meter(name2, Storage::CSV);
	layerMeter->AttachPopulation(layer2);
	network->AddMeter(layerMeter);

	sprintf(name3,"Projection_bcm_inhib_n%d.csv",mpiRank);
	Meter* connMeter2 = new Meter(name3, Storage::CSV);
	connMeter2->AttachProjection(layer2->GetIncomingProjections()[1],0);
	network->AddMeter(connMeter2);
	// end Meters
	//////////////////////////////


	//////////////////////////////
	///// Input data

	char* filenameTraining, *filenameTesting, *filenameTrainingLabels, *filenameTestingLabels;

	if(m_architecture == BGL)
	{
		filenameTraining = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
		filenameTesting = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingData_2colors.csv";

		filenameTrainingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingLabels.csv";
		filenameTestingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingLabels.csv";
	}
	else if(m_architecture == PC)
	{
		filenameTraining = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
		filenameTesting = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingData_2colors.csv";

		filenameTrainingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingLabels.csv";
		filenameTestingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingLabels.csv";
	}

	int iterations = 1;

	if( m_architecture == PC)
	{
		iterations = 50;
	}
	else if( m_architecture == BGL)
	{
		iterations = 1;
	}

	Storage storage;
	DataSources dataSources;

	vector<vector<float> > trainingData = storage.LoadDataFloatCSV(filenameTraining,nrItems,false); //dataSources.GetBars(sizeX,sizeY, nrItems);
	vector<vector<float> > trainingLabels = storage.LoadDataFloatCSV(filenameTrainingLabels,nrItems,false);

	layer1->SwitchOnOff(false);	// fixed during training phase

	for(int j=0;j<iterations;j++)
	{
		for(int i=0;i<trainingData.size();i++)
		{
			//for(int k=0;k<40;k++)
			//{
			if(i==9)
				bool b = false;
				layer1->SetValuesAll(trainingData[i]);

				network->Simulate();
			//}

			if(mpiRank == 0)
				if(i%100 == 0)
				{
					cout<<"\n"<<i<<" ";
					cout.flush();
				}
		}

		if(mpiRank == 0)
			cout<<"["<<j<<"]";
	}

	//if(j%30 == 0)	
	network->RecordAll();	
}

void NetworkMNIST::NetworkMNISTRunLateralMaps2(int mpiRank, int mpiSize)
{
	DataSources dataSources;

	int sizeX = 28;
	int sizeY = 28;
	int nrItems = 4;

	bool isTriesch = true;

	Network* network = new Network();
	network->SetMPIParameters(mpiRank,mpiSize);

	int nrInputHypercolumns = 1;
	int nrInputRateUnits = sizeX*sizeY;
	int nrOutputHypercolumns = 2;
	int nrOutputRateUnits = 2;//sizeX+sizeY;

	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::GradedThresholded);
	PopulationColumns* layer2 = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::GradedThresholded);

	network->AddPopulation(layer1);
	network->AddPopulation(layer2);

	FullConnectivity* full = new FullConnectivity();
	FullConnectivity* full2;
	FullConnectivityNoLocalHypercolumns* full3NoLocal;

	layer2->AddPre(layer1,full);

	bool thresholded = true;
	ProjectionModifierTriesch* eTriesch = new ProjectionModifierTriesch(0.0005,0.2,0.05,1.0/(float)nrOutputRateUnits, thresholded);//0.05,0.2,0.005,1.0/(float)nrOutputRateUnits, thresholded);

	if(isTriesch)
		full->AddProjectionsEvent(eTriesch);

	//float eta1 = 3, eta2= 2.4, eta3 = 1.5, alpha = 0.005, beta = 200;
	float eta1 = 0.5, eta2= 0.02, eta3 = 0.02, alpha = 0.0005, beta = 10;//alpha = 1.0/8.0, beta = 10;
	bool lateral = false;

	ProjectionModifierFoldiak* eFoldiak = new ProjectionModifierFoldiak(eta1, eta2, eta3, alpha, beta, lateral);
	lateral = true;
	alpha = 0;//0.1;//0.75;
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
		//layer2->AddPre(layer2,full3NoLocal);
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

	//vector<vector<float> > trainData = dataSources.GetBars(sizeX,sizeY, nrItems);
	char* filenameTraining, *filenameTesting, *filenameTrainingLabels, *filenameTestingLabels;

	if(m_architecture == BGL)
	{
		filenameTraining = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
		filenameTesting = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingData_2colors.csv";

		filenameTrainingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_trainingLabels.csv";
		filenameTestingLabels = "/gpfs/scratch/s/simonbe/Databases/MNIST/MNIST_testingLabels.csv";
	}
	else if(m_architecture == PC)
	{
		filenameTraining = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingData_2colors.csv";//MNIST_2colors.h5";
		filenameTesting = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingData_2colors.csv";

		filenameTrainingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingLabels.csv";
		filenameTestingLabels = "C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingLabels.csv";
	}
	
	Storage storage;
	vector<vector<float> > trainData = storage.LoadDataFloatCSV(filenameTraining,nrItems,false); //dataSources.GetBars(sizeX,sizeY, nrItems);
	vector<vector<float> > trainLabels = storage.LoadDataFloatCSV(filenameTrainingLabels,nrItems,false);

	int iterations = 15;
	int iterSameStimuli = 100;

	if(!isTriesch)
		iterSameStimuli = 10;

	layer1->SwitchOnOff(false);	// fixed during training phase

	for(int j=0;j<iterations;j++)
	{
		for(int i=0;i<trainData.size();i++)
		{
			if(mpiRank == 0 && i%2==0)
				cout<<i;
			//if(!isTriesch)
			//{
			//	// in order to settle recurrent activity
			//	eFoldiak->SwitchOnOff(false);
			//	eFoldiakLateral->SwitchOnOff(false);
			//}

			for(int k=0;k<iterSameStimuli;k++)
			{
			//	if(!isTriesch && k==iterSameStimuli-1)
			//	{
			//		eFoldiak->SwitchOnOff(true);
			//		eFoldiakLateral->SwitchOnOff(true);
			//	}

				for(int m=0;m<1;m++)
				{
					layer1->SetValuesAll(trainData[i]);
					//for(int n=0;n<3;n++)
					network->Simulate();
				}
			}

			// allow units to reset
			network->Reset();

			//if(i%50 == 0)
			//{
			//	network->RecordAll();
			//	if(mpiRank == 0)
			//		cout<<"Storing.";
			//}
		}	
	}

	network->RecordAll();
}*/