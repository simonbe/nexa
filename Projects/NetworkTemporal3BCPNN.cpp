#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "NetworkBCPNN.h"
#include "NetworkKussul.h"
#include "NetworkCL.h"
#include "NetworkMDS.h"
#include "NetworkMI.h"
#include "NetworkVQ.h"
#include "Network.h"

#include "NetworkTemporal3BCPNN.h"
#include "NetworkTriesch.h"
#include "NetworkFoldiak.h"

//using namespace std;

void NetworkTemporal3BCPNN::NetworkSetupStructure()
{
	ofstream wout;
	wout.open("time_start.txt",ios::trunc);
	wout << "Timing started!";
	wout.close();
	m_dataAssigner=new SequentialAssignmentStrategy();
	m_architecture = CRAY;
	m_nrInputHypercolumns = 13;
	m_nrInputRateUnits = 33;
	m_verbose=true;
	// layer 2
	m_nrHypercolumns2 = 20;//24;//96;
	m_nrRateUnits2 = 33;//100;

	m_nrHypercolumns3 = 1;//24;//96;
	m_nrRateUnits3 = m_nrHypercolumns2*m_nrRateUnits2;//100;
	int m_mdsSize=15;

	// Structures setting up everything to calculate mutual information (MI)/correlations (Pearson) on input data + multi-dimensional scaling (MDS) + vector quantization (VQ)
	m_structureInput = new StructureMIMDSVQ(); 
	m_layerInput = new PopulationColumns(this,m_nrInputHypercolumns,m_nrInputRateUnits,PopulationColumns::Graded,MPIDistribution::ParallelizationHypercolumns); // input
	m_structureInput->SetMDSDimension(m_mdsSize); // Number of dimensions in multi-dimensional scaling matrix (Size input x dimensions)
	m_structureInput->SetMDSMeasure(true);
	m_structureInput->SetupStructure(this,m_layerInput,m_nrHypercolumns2,m_nrRateUnits2, true);
	m_layer2 = m_structureInput->GetLayer(1);
	m_layer3 = new PopulationColumns(this,m_nrHypercolumns3,m_nrRateUnits3,PopulationColumns::Graded,MPIDistribution::ParallelizationHypercolumns); 
	PopulationModifier* wta = new WTA();
	m_layer3->AddPopulationModifier(wta);
	FullConnectivity* full10 = new FullConnectivity();
	m_layer3->AddPre(m_layer2,full10); // Feedforward
	m_bcpnn = new ProjectionModifierBcpnnOnline(0.0001, 10e-6); //real value for 11000: 0.00001
	full10->AddProjectionsEvent(m_bcpnn);
	m_bcpnn->SwitchOnOff(false);
	AddPopulation(m_layer3);
	//m_structureInput->MDSHypercolumns()->SetUseThreshold(true,0.9);
	m_structureInput->VQ()->GetCSL()->SetEta(0.001,true);
	m_structureInput->MDSHypercolumns()->SetUpdateSize(2e-3);
	m_structureInput->MDSHypercolumns()->SwitchOnOff(false);
	m_structureInput->MDS()->SwitchOnOff(false);
	m_structureInput->VQ()->SwitchOnOff(false);
	m_structureInput->CSLLearn()->SwitchOnOff(false);
	AddTiming(this);
	AddTiming(m_structureInput->VQ());
	AddTiming(m_structureInput->MDSHypercolumns());
	AddTiming(m_structureInput->MDS());
	AddTiming(m_layerInput);
	AddTiming(m_layer2);
	AddTiming(m_layer3);
	AddTiming(m_structureInput->CSLLearn());
}

void NetworkTemporal3BCPNN::NetworkSetupMeters()
{
	m_structureInput->SetupMeters(this->MPIGetNodeId(),this->MPIGetNrProcs());
	m_activityMeter = new Meter("ClassOutput.csv",Storage::CSV);
	m_activityMeter->AttachPopulation(m_layer3);
	//m_activityMeter->SwitchOnOff(false);
	AddMeter(m_activityMeter);
	m_inputMeter = new Meter("inputLayer.csv", Storage::CSV);
	m_inputMeter->AttachPopulation(m_layerInput);
	AddMeter(m_inputMeter);

	m_layer1Meter = new Meter("layer1.csv", Storage::CSV);
	m_layer1Meter->AttachPopulation(m_layer2);
	AddMeter(m_layer1Meter);
}

void NetworkTemporal3BCPNN::NetworkSetupParameters()
{
	// not used
}


vector<float> NetworkTemporal3BCPNN::toBinary(const vector<float>& data, int nrHc, int nrMc)
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

void NetworkTemporal3BCPNN::ComputeCorrelation(const vector<vector<float> >& trainingData,  PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iter) {
	structure->Pearson()->SwitchOnOff(true);
	for(int i=0;i<trainingData.size();i++) {
		if(m_verbose && this->MPIGetNodeId() == 0) {
			cout<<i<<"("<<iter<<") ";
			cout.flush();
		}
		inputLayer->SetValuesAll(m_dataAssigner->prepareValues(i,trainingData));
		this->Simulate();
		cout.flush();
	}
	if(m_verbose &&this->MPIGetNodeId() == 0)
		cout<<"\n";
	structure->Pearson()->SwitchOnOff(false);
}

bool NetworkTemporal3BCPNN::ComputeMDS(const vector<vector<float> >& trainingData,  PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iter) {
	structure->MDSHypercolumns()->SwitchOnOff(true);
	structure->MDS()->SwitchOnOff(true);
	for(int i=0;i<trainingData.size();i++) {
		if(m_verbose && this->MPIGetNodeId() == 0) {
			cout<<i<<"("<<iter<<") ";
			cout.flush();
		}
		inputLayer->SetValuesAll(m_dataAssigner->prepareValues(i,trainingData));
		this->Simulate();
		cout.flush();
		if(!structure->MDS()->IsOn() || !structure->MDSHypercolumns()->IsOn()) {
			structure->MDS()->SwitchOnOff(false);
			structure->MDSHypercolumns()->SwitchOnOff(false);
		}
	}
	if(m_verbose && this->MPIGetNodeId() == 0)
		cout<<"\n";
	if(!structure->MDS()->IsOn() || !structure->MDSHypercolumns()->IsOn()) {
			return false;
	}
	structure->MDSHypercolumns()->SwitchOnOff(false);
	structure->MDS()->SwitchOnOff(false);
	return true;
}

void NetworkTemporal3BCPNN::ComputeVQ(const vector<vector<float> >& trainingData,  PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iter) {
	structure->GetLayer(1)->SwitchOnOff(true);
	structure->VQ()->SwitchOnOff(true);
	for(int i=0;i<trainingData.size();i++) {
		if(m_verbose && this->MPIGetNodeId() == 0) {
			cout<<i<<"("<<iter<<") ";
			cout.flush();
		}
		inputLayer->SetValuesAll(m_dataAssigner->prepareValues(i,trainingData));
		this->Simulate();
		cout.flush();
	}
	if(m_verbose && this->MPIGetNodeId() == 0)
		cout<<"\n";
	structure->VQ()->SwitchOnOff(false);
}

bool NetworkTemporal3BCPNN::ExtractFeatures(const vector<vector<float> >& trainingData,  PopulationColumns* inputLayer, StructureMIMDSVQ* structure) {
	structure->CSLLearn()->SwitchOnOff(true);
	structure->SetRecording(false);
	float clC = m_nrRateUnits2*2;
	int idx=0;
	for(int i=0;i<trainingData.size();i++)
	{
		inputLayer->SetValuesAll(m_dataAssigner->prepareValues(i,trainingData));
		// next time step
		this->Simulate();
		cout.flush();
	}
	if(m_verbose && this->MPIGetNodeId()== 0)
		cout<<"\n";
	if(!structure->CSLLearn()->IsOn())
			return false;
	//this->RecordAll();
	structure->CSLLearn()->SwitchOnOff(false);
	return true;
}

void NetworkTemporal3BCPNN::TrainBCPNN(const vector<vector<float> >& trainingData, const vector<vector<float> >& trainingLabels, PopulationColumns* inputLayer, int iterationsBCPNN) {
	int prevLabel=-1;
	for (int k=0;k<iterationsBCPNN;k++) {
		for(int j=0;j<trainingData.size();j++) {
			if(prevLabel<0&&j==0) {
				m_bcpnn->SwitchOnOff(false);
				prevLabel=j;
				m_layerInput->SetValuesAll(m_dataAssigner->prepareValues(j,trainingData));
				this->Simulate();
				m_bcpnn->SwitchOnOff(true);
			}
			vector<float> v(m_nrRateUnits3,0.0);
			v[prevLabel]=1.0;
			prevLabel=(int)trainingLabels[j][0]-1;
			m_layerInput->SetValuesAll(m_dataAssigner->prepareValues(j,trainingData));
			inputLayer->SetValuesAll(v);
			this->Simulate();
			if(m_verbose && this->MPIGetNodeId() == 0)
				cout << j << "(" << k << ") ";
		}
		if(m_verbose && this->MPIGetNodeId() == 0)
				cout <<endl;
	}
	this->Simulate();
	m_bcpnn->SwitchOnOff(false);
}
void NetworkTemporal3BCPNN::TrainLayer(const vector<vector<float> >& trainingData, PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iterationsCorrs, int iterationsMDS, int iterationsVQ, int iterationsFeatures)
{
	// Training phase

	int nrTrainImages = trainingData.size();
	structure->MDSHypercolumns()->SwitchOnOff(false);
	structure->MDS()->SwitchOnOff(false);
	structure->VQ()->SwitchOnOff(false);
	structure->CSLLearn()->SwitchOnOff(false);
	structure->GetLayer(1)->SwitchOnOff(false);
	structure->SetRecording(false);

	// Semi-sequential version

	// 1. Training phase
	// 1A. Patches creation
	structure->CSLLearn()->SetMaxPatterns(nrTrainImages);
	structure->CSLLearn()->SetEta(0.001);
	
	// turn of response in 2nd layer during initial training phase for speed
	structure->GetLayer(1)->SwitchOnOff(false);
	
	int j=0;

	while(j<iterationsCorrs) {
		//if(j==(int)iterationsPatches*0.8)
		//	break;
		//for(int i=0;i<trainingData.size();i++) {
		ComputeCorrelation(trainingData,inputLayer,structure,j);
		j++;
	}
	
	j=0;
	while(j<iterationsMDS) {
		//if(j==(int)iterationsPatches*0.9)
		//	break;
		//for(int i=0;i<trainingData.size();i++) {
		if(!ComputeMDS(trainingData,inputLayer,structure,j))
			break;
		j++;
	}
	structure->GetLayer(1)->SwitchOnOff(true);
	j=0;
	while(j<iterationsVQ) {
		//for(int i=0;i<trainingData.size();i++) {
			if(m_verbose && this->MPIGetNodeId() == 0)
				cout << "DataPoint (patches): " << j <<endl;
		ComputeVQ(trainingData,inputLayer,structure,j);
		j++;
	}

	structure->CSLLearn()->SetMaxPatterns(nrTrainImages);
	j=0;
	while(j<iterationsFeatures) {
		//for(int i=0;i<trainingData.size();i++) {
			if(m_verbose && this->MPIGetNodeId() == 0)
				cout << "DataPoint (features): " << j <<endl;
			if(!ExtractFeatures(trainingData,inputLayer,structure))
				break;
		//}
			j++;
	}
	this->RecordAll();
}

void NetworkTemporal3BCPNN::NetworkRun()
{
	int nrColors = 2;
	char* filenameTrain, *filenameTest, *filenameTrainLabels;
	int nrTrainImages = 1000;
	int nrTestImages = 1000;

	int iterationsPatches = 10;

	if(m_architecture == PC)
	{
		nrTrainImages = 150;
		nrTestImages = 150;
		//iterationsPatches = 10;
		filenameTrain="D:\\Databases\\tidigits\\trainData.csv";
		filenameTrainLabels="D:\\Databases\\tidigits\\trainDataLabels.csv";
		filenameTest="D:\\Databases\\tidigits\\trainData.csv";
	}
	else {
		nrTrainImages = 25000;
		nrTestImages = 25000;
		filenameTrain="/cfs/klemming/nobackup/p/paherman/Tin/data/trainData_np.csv";
		filenameTrainLabels="/cfs/klemming/nobackup/p/paherman/Tin/data/trainDataLabels_np.csv";
		filenameTest="/cfs/klemming/nobackup/p/paherman/Tin/data/testData_np.csv";
	}
	int stepsStimuliOn = 1;
	// Specify input data
	Storage storage;
	Storage storageLabels;
	storage.SetMPIParameters(this->MPIGetNodeId(),this->MPIGetNrProcs());

	vector<vector<float> > trainingData = storage.LoadDataFloatCSV(filenameTrain,nrTrainImages,true);
	vector<vector<float> > trainingLabels = storageLabels.LoadDataFloatCSV(filenameTrainLabels,nrTrainImages,true);
	vector<vector<float> > testData = storage.LoadDataFloatCSV(filenameTest,nrTestImages,true);

	vector<int> partsOfDataToUseAsInput = m_layerInput->GetMPIDistributionHypercolumns(this->MPIGetNodeId());
	vector<int> partsOfDataToUseAsOutput = vector<int>();//layer1->GetMPIDistributionHypercolumns(mpiRank);

	m_layerInput->SwitchOnOff(false);
	m_structureInput->CSLLearn()->SetMaxPatterns(trainingData.size());
	// Turn off recording during training
	m_structureInput->SetRecording(false);
	m_layerInput->SetRecording(false);
	m_layer2->SetRecording(false);
	m_layer3->SwitchOnOff(false);
	m_layer3->SetRecording(false);
	// Train 1st layer
	//m_layer2->GetIncomingProjections()[1]->SwitchOnOff(false); // turn off the Projections calculating correlations in layer 2 so (change to switch off in structure by default instead) (!)
	//m_layer2->GetIncomingProjections()[2]->SwitchOnOff(false);

	int iterationsCorrs = 10; // run through x times, 5
	int iterationsMDS = 10;//200
	int iterationsVQ = 2;//200
	int iterationsFeatures = 2; //+nrtraindata the timesteps to run CSL (first totalNrTrainingData steps are to collect the data to run on)
	int iterationsBCPNN=3;
	TrainLayer(trainingData,m_layerInput,m_structureInput,iterationsCorrs,iterationsMDS,iterationsVQ,iterationsFeatures);
	// Run through all training and test data
	m_layer2->SwitchOnOff(true);
	// turn on all recording
	m_layerInput->SetRecording(true);
	//m_layer2->SetRecording(true);

	TimingStart("RunThrough");
	for(int i=0;i<trainingData.size();i++)
	{
		if(i==1)
			m_layer2->SetRecording(true);
		m_layerInput->SetValuesAll(m_dataAssigner->prepareValues(i,trainingData));
		this->Simulate();
	}

	for(int i=0;i<testData.size();i++)
	{
		m_layerInput->SetValuesAll(m_dataAssigner->prepareValues(i,testData));
		this->Simulate();
	}
	m_layerInput->SetRecording(false);
	this->Simulate();
	m_layer2->SetRecording(false);
	TimingStop("RunThrough");
	PopulationModifier* wta = new WTA();
	m_layer2->AddPopulationModifier(wta);
	TrainBCPNN(trainingData, trainingLabels, m_layer3, iterationsBCPNN);
	m_layer3->SwitchOnOff(true);
	for(int i=0;i<trainingData.size();i++)
	{
		if(i==2)
			m_layer3->SetRecording(true);
		m_layerInput->SetValuesAll(m_dataAssigner->prepareValues(i,trainingData));
		this->Simulate();
	}
	for(int i=0;i<testData.size();i++)
	{
		m_layerInput->SetValuesAll(m_dataAssigner->prepareValues(i,testData));
		this->Simulate();
	}
	this->Simulate();
	this->Simulate();
	TimingStop("RunThrough");
}