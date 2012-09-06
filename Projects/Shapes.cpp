#include "Shapes.h"
#include <fstream>
#include <time.h>

#define DEFAULT_TIMESTEP 1e-3
#define RETINA_WIDTH 128
#define RETINA_HEIGHT 128
#define RETINA_SIZE RETINA_HEIGHT*RETINA_WIDTH

#define BCPNN_THRESHOLD 5
#define BCPNN_LEARNING_RATE 0.1 //0.1

#define WEIGHTS_TO_LAYER2 0.01

#define ACTIVATION_THRESHOLD 1
#define INPUT_AMPLITUDE 1

#define ADAPTATION_AMP 0.05
#define ADAPTATION_TAU 0.01
#define DEPRESSION_STRENGTH 0.1

#define TRAINING_STEPS 10
#define TEST_STEPS 150

#define PATTERN_NOISE_FACTOR 0.6 //0.6
#define BACKGROUND_NOISE_FACTOR 0.002 //0.005

void ShapesNetwork::Train(){
	int x,y,i,j,rowStart;
	int trainingSteps = TRAINING_STEPS;

	int numPatterns = 3;
	const char* patternFiles[] = {"patterns/plus-128.csv",
				"patterns/circle-128.csv",
				"patterns/triangle-128.csv"};

	vector<float> empty(RETINA_SIZE,0);

	//Load training patterns
	vector< vector<float> > patterns(numPatterns);
	vector< vector<float> > sequence;
	for(i=0;i<numPatterns;i++){
		 sequence = this->ReadPatternSequence(patternFiles[i]);
		 patterns[i] = sequence[0]; //just one frame
	}

	if(rank==0)
	  	std::cout << "[" << rank << "]" << "Training patterns loaded" << std::endl;

	recConn->SwitchOnOff(false);
	recConnBcpnn->SwitchOnOff(true);
	recConnBcpnn->GetTransferFunction()->SwitchOnOff(false);
	recConnAdaptation->SwitchOnOff(false);
	recConnDepression->SwitchOnOff(false);
	this->GetLayer(2)->SwitchOnOff(false);

	if(rank==0)
		std::cout << "[" << rank << "]" << "Training" << std::endl;

	for(i=0;i<trainingSteps*numPatterns;i++){
		this->GetLayer(0)->SetValuesAll(patterns[i%numPatterns]);
		this->Simulate();
		this->GetLayer(0)->SetValuesAll(empty); //clear input
		this->Simulate(1);
	
	}

	//Disable training, enable testing
	recConnBcpnn->SwitchOnOff(false);
	recConnBcpnn->GetTransferFunction()->SwitchOnOff(true);
	recConn->SwitchOnOff(true);
	recConnAdaptation->SwitchOnOff(true);
	recConnDepression->SwitchOnOff(true);
	this->GetLayer(2)->SwitchOnOff(true);

	this->GetLayer(0)->SetValuesAll(empty); //clear input
	this->GetLayer(1)->SetValuesAll(empty); //clear input
	this->Simulate(1);


	if(rank==0)
		std::cout << "[" << rank << "]" << "Training finished" << std::endl;

	if(rank==0)
	 	std::cout << "[" << rank << "]" << "Removing synapses" << std::endl;

	this->RemoveUnusedSynapses();

	if(rank==0)
		std::cout << "[" << rank << "]" << "Testing noisy training patterns" << std::endl;
	
	TimingStart("Test");
	
	float rnd;
	srand( (unsigned)time(0) );

	vector<int> run1(3);
	vector<int> run2(3);
	vector<int> run3(3);
	vector<int> run4(3);
	vector<vector<int> > runs;
	run1[0] = 0;
	run1[1] = 1;
	run1[2] = 2;
	run2 = run1;
	run3[0] = 0;
	run3[1] = 2;
	run3[2] = 1;
	run4[0] = 2;
	run4[1] = 1;
	run4[2] = 0;

	runs.push_back(run1);
	runs.push_back(run2);
	runs.push_back(run3);
	runs.push_back(run4);
	
	for(int iter=0;iter<runs.size();iter++)
	{
		vector<int> currentPatterns = runs[iter];
		// run 0: independently stimulated
		// run 1: sequentially stimulated

		for(i=0;i<numPatterns;i++){
			//add noise
			if(iter == 0)
			{
				for(y=0;y<RETINA_HEIGHT;y++){
					rowStart = y*RETINA_HEIGHT;
					for(x=0;x<RETINA_WIDTH;x++){
						rnd = (float)rand()/(float)RAND_MAX;
						if(patterns[i][rowStart+x] > 0){
							if(rnd < PATTERN_NOISE_FACTOR){
								patterns[i][rowStart+x] = 0;
							}
						} else{
							if(rnd < BACKGROUND_NOISE_FACTOR){
								patterns[i][rowStart+x] = 1;
							}
						}
					}
				}
			}

			for(j=0;j<TEST_STEPS;j++){
				this->GetLayer(0)->SetValuesAll(patterns[currentPatterns[i]]);//patterns[i]);
				this->Simulate();
				std::cout.flush();
			}

			//this->GetLayer(0)->SetValuesAll(empty); //clear input

			if(iter==0 || i == numPatterns-1)
			{
				for(int i=0;i<5;i++)
				{
					this->GetLayer(0)->SetValuesAll(empty); //clear input
					this->GetLayer(1)->SetValuesAll(empty); //clear recurrent layer
					this->Simulate();

					(*((PopulationColumns*)this->GetLayer(2))->GetUnitPropertiesLayer())[0]->Clear(); // clear transfer fcn history
					this->ClearEventsIncoming();
					vector<RateUnit*> mcs = ((PopulationColumns*)this->GetLayer(2))->GetRateUnits();
					for(int m=0;m<mcs.size();m++)
					{
						mcs[m]->SetValue(0.0);
						mcs[m]->SetSubThresholdValue(0.0);
					}
				}
			}
			else
			{
				for(int m=0;m<10;m++)
				{
					this->GetLayer(0)->SetValuesAll(empty); //clear input
					this->GetLayer(1)->SetValuesAll(empty); //clear input
					this->ClearEventsIncoming();
					this->Simulate(1);
				}
			}
		}

		this->RecordAll();
	}
	TimingStop("Test");
}

void ShapesNetwork::Test(){
	int i,j;

	int numPatterns = 3;
	const char* patternFiles[] = {"patterns/plus-128-test.csv",
							"patterns/circle-128-test.csv",
							"patterns/triangle-128-test.csv"};

	vector< vector<float> > pattern;
	vector<float> empty(RETINA_SIZE,0);

	if(rank==0)
		std::cout << "[" << rank << "]" << "Testing sequences" << std::endl;

	for(i=0;i<numPatterns;i++){
		pattern = ReadPatternSequence(patternFiles[i]);

		for(j=0;j<pattern.size();j++){ //for every pattern/frame in sequence
			this->GetLayer(0)->SetValuesAll(pattern[j]);
			this->Simulate();
		}

		this->GetLayer(0)->SetValuesAll(empty); //clear input
		this->Simulate();
	}

	if(rank==0)
		std::cout << "[" << rank << "]" << "Testing finished" << std::endl;

}

void ShapesNetwork::RemoveUnusedSynapses(){

	Projection* conn = this->GetLayer(1)->GetIncomingProjections()[1];
	vector<long> postIds = conn->GetPostIds();
	vector<long> preIds;
	int num = 0;

	for(int i=0; i<postIds.size();i++){
		preIds = conn->GetPreIds(postIds[i]);
		
		for(int j=0; j<preIds.size();j++){
		  if(this->GetWeight(preIds[j], postIds[i]) < 1.0001 && 
		     this->GetWeight(preIds[j], postIds[i]) > 0.999 ){

		    conn->EraseSynapse(preIds[j],postIds[i]);
		    num++;
		  }
		}
	}
	if(rank==0)
	  std::cout << "[" << rank << "]" << "Synapses removed: " << num << std::endl;
}
	


void ShapesNetwork::NormalizeWeights(){

	Projection* conn = this->GetLayer(1)->GetIncomingProjections()[1];
	vector<long> postIds = conn->GetPostIds();
	vector<long> preIds;
	double prod,sum;
	double normFactor;
	float e = exp(1.0);

	for(int i=0; i<postIds.size();i++){
		preIds = conn->GetPreIds(postIds[i]);
		prod = 1;
		for(int j=0; j<preIds.size();j++){
		  if (this->GetWeight(preIds[j], postIds[i]) > 0.01)
			prod *= this->GetWeight(preIds[j], postIds[i]);
		}
		
		normFactor = pow(e/prod,1.0/(float)postIds.size()); 

		for(int j=0; j<preIds.size();j++){
			this->SetWeight(this->GetWeight(preIds[j], postIds[i])*normFactor, preIds[j], postIds[i]);
		}

		if(rank == 0 && i < 10){

		  sum=0;
		  for(int j=0; j<preIds.size();j++){
		    if(this->GetWeight(preIds[j], postIds[i]) != 0)
		      sum += log(this->GetWeight(preIds[j], postIds[i]));
		  }
		  std::cout << prod << ", " << e << ", " << 1.0/(float)postIds.size() << std::endl;
		  std::cout << "[" << rank << "]" << normFactor << " Sum log(w): " << sum << std::endl;
	        }
	}
}	



// Defines network structure
void ShapesNetwork::NetworkSetupStructure() {

	// total nrHypercolumns*nrRateUnits neural units
	int nrHypercolumns = RETINA_SIZE;
	int nrRateUnits = 1;

	MPI_Comm_rank(NETWORK_COMM_WORLD,&(this->rank));

	PopulationColumns* layer0 = new PopulationColumns(this, nrHypercolumns,
			nrRateUnits, PopulationColumns::Graded);
	this->AddPopulation(layer0); // this population/layer will have index 0
	//layer0->MPI()->MPIMakeLayerValuesLocal();

	PopulationColumns* layer1 = new PopulationColumns(this, nrHypercolumns,
			nrRateUnits, PopulationColumns::Graded);
	this->AddPopulation(layer1); // this population/layer will have index 1

	OneToOneConnectivity* ffConn = new OneToOneConnectivity(false);
	ffConn->SetRandomWeights(1,1);
	ffConn->AddUnitsProperty(new TransferLinear(false));
	layer1->AddPre(layer0,ffConn);

	//recConn = new RandomConnectivity(0.01,true);
	recConn = new FullConnectivity();
	recConn->SetRandomWeights(0,1);

	recConnBcpnn = new ProjectionModifierBcpnnOnline(BCPNN_LEARNING_RATE,0.001);
	((TransferBcpnnOnline*)recConnBcpnn->GetTransferFunction())->SetThreshold(BCPNN_THRESHOLD);
	recConn->AddProjectionsEvent(recConnBcpnn);

	recConnDepression = new ProjectionModifierDepression(DEPRESSION_STRENGTH);
	//recConn->AddProjectionsEvent(recConnDepression);

	recConnAdaptation = new PopulationModifierAdaptation(ADAPTATION_AMP,ADAPTATION_TAU,PopulationModifierAdaptation::Standard);
	layer1Threshold = new Threshold(ACTIVATION_THRESHOLD);
	//	layer1->AddPopulationModifier(recConnAdaptation);
	layer1->AddPopulationModifier(layer1Threshold);
	layer1->AddPre(layer1,recConn);

	
	PopulationColumns* layer2 = new PopulationColumns(this,1,60*40,PopulationColumns::Graded);
	this->AddPopulation(layer2);
	float probRecurr = 0.2;
	reTIDe = new StructureReTIDe(probRecurr,-1);
	
	reTIDe->Initialize(this,layer2);


	//FullConnectivity* full1_2 = new FullConnectivity();
	//full1_2->SetRandomWeights(0,WEIGHTS_TO_LAYER2);
	//layer2->AddPre(layer1,full1_2);
	
	RandomConnectivity* rnd1_2 = new RandomConnectivity(0.5,true);
	vector<float> weightsL2(5);
	//weightsL2[0] = 0.01;
	//weightsL2[1] = 0.009;
	weightsL2[0] = 0.008;
	weightsL2[1] = 0.007;
	weightsL2[2] = 0.006;
	weightsL2[3] = 0.005;
	weightsL2[4] = 0.004;
	
	if(weightsL2.size()>0)
		rnd1_2->SetRandomWeights(0,weightsL2,this);
	else
		rnd1_2->SetRandomWeights(0,WEIGHTS_TO_LAYER2);
	layer2->AddPre(layer1,rnd1_2);


	this->SetSeed(10);
}


vector< vector<float> > ShapesNetwork::ReadPatternSequence(const char* filename){
	vector< vector<float> > sequence;
	ifstream infile( filename );

	while (infile){
		string s;
		if (!getline( infile, s )) break;
		istringstream ss( s );
		vector <float> frame;
		int i=0;
		while (ss){
			string pixel;
			if (!getline( ss, pixel, ',' )) break;
			if(i>0){ //jump first column
				frame.push_back( atof(pixel.c_str())*INPUT_AMPLITUDE );
			}
			i++;
		}
		sequence.push_back( frame );

		//if(rank==0)
		//	cout << "Read frame of size: " <<  frame.size()  << std::endl;
	}

	return sequence;
}


// Setup outputs
void ShapesNetwork::NetworkSetupMeters() {
	// example of file output
	Meter* layer0Meter = new Meter("activity_0.csv", Storage::CSV,
			Storage::Append);
	layer0Meter->AttachPopulation(this->GetLayer(0));
	this->AddMeter(layer0Meter);

	Meter* layer1Meter = new Meter("activity_1.csv", Storage::CSV,
			Storage::Append);
	layer1Meter->AttachPopulation(this->GetLayer(1));
	this->AddMeter(layer1Meter);

	Meter* layer2Meter = new Meter("activity_2.csv", Storage::CSV,
			Storage::Append);
	layer2Meter->AttachPopulation(this->GetLayer(2));
	this->AddMeter(layer2Meter);

	// also check simulation times
	this->AddTiming(this);
	this->AddTiming(this->GetLayer(0));
	this->AddTiming(this->GetLayer(1));
	this->AddTiming(this->GetLayer(2));

	if(rank == 0){
		char* readout_conn_fn = new char[30];
		sprintf(readout_conn_fn, "readout_Projections_%d.txt", rank);
		Meter* readoutProjectionMeter = new Meter(readout_conn_fn, Storage::CSV);
		Projection* conn = this->GetLayer(1)->GetIncomingProjections()[1];
		readoutProjectionMeter->AttachProjection(conn,	0);
		this->AddMeter(readoutProjectionMeter);
	}

	reTIDe->SetupStructure(this,(PopulationColumns*)this->GetLayer(2));
}

void ShapesNetwork::NetworkSetupParameters() {
	// not used
}


void ShapesNetwork::NetworkRun() {
	//Store local minicolumns
	if(rank==0)
		std::cout << "[" << rank << "]" << "Starting simulation" << std::endl;

	minicolumns = ((PopulationColumns*) this->GetLayer(0))->GetRateUnits();
	this->GetLayer(0)->SwitchOnOff(false);
	this->Train();
	//this->Test();

	if(rank==0)
		std::cout << "[" << rank << "]" << "Finished, storing data" << std::endl;

	this->RecordAll();
	this->StoreAnalysis();

	if(rank==0)
		std::cout << "[" << rank << "]" << "Finished" << std::endl;
}


/*
void ShapesNetwork::NormalizeWeights(){

	Projection* conn = this->GetLayer(1)->GetIncomingProjections()[1];
	vector<long> postIds = conn->GetPostIds();
	vector<long> preIds;
	float sum;

	for(int i=0; i<postIds.size();i++){
		preIds = conn->GetPreIds(postIds[i]);
		sum = 0;
		for(int j=0; j<preIds.size();j++){
			sum += this->GetWeight(preIds[j], postIds[i]);
		}

		for(int j=0; j<preIds.size();j++){
			this->SetWeight(this->GetWeight(preIds[j], postIds[i])/sum, preIds[j], postIds[i]);
		}

		sum=0;
		for(int j=0; j<preIds.size();j++){
			sum += this->GetWeight(preIds[j], postIds[i]);
		}

		if(rank == 0 && i < 10)
			std::cout << "[" << rank << "]" << "Sum: " << sum << std::endl;
	}
	}*/
