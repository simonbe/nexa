#include "StructureMIMDSVQ.h"

void StructureMIMDSVQ::SetupStructure(Network* network, int nrInputHypercolumns, int nrInputRateUnits, int nrMiddleHypercolumns, int nrMiddleRateUnits, bool addInputLayerToNetwork)//, int nrOutputHypercolumns, int nrOutputRateUnits)
{
	PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded); // input
	SetupStructure(network,layer1,nrMiddleHypercolumns,nrMiddleRateUnits,addInputLayerToNetwork);//,nrOutputHypercolumns,nrOutputRateUnits);
}

void StructureMIMDSVQ::SetupStructure(Network* network, PopulationColumns* layerInput,int nrMiddleHypercolumns, int nrMiddleRateUnits, bool addInputLayerToNetwork, bool useSilentHypercolumns, float silentHypercolumnsThreshold)//, int nrOutputHypercolumns, int nrOutputRateUnits)
{
	m_network = network;
	//PopulationColumns* layerOutput;

	//PopulationColumns* layer1 = new PopulationColumns(network,nrInputHypercolumns,nrInputRateUnits,PopulationColumns::Graded); // input
	PopulationColumns* layer2;
	if(m_featureExtraction == StructureMIMDSVQ::Recurr)
		layer2 = new PopulationColumns(network,nrMiddleHypercolumns,nrMiddleRateUnits,PopulationColumns::GradedThresholded,MPIDistribution::ParallelizationHypercolumns, useSilentHypercolumns,silentHypercolumnsThreshold); // middle
	else
		layer2 = new PopulationColumns(network,nrMiddleHypercolumns,nrMiddleRateUnits,PopulationColumns::Graded,MPIDistribution::ParallelizationHypercolumns, useSilentHypercolumns,silentHypercolumnsThreshold); // middle
	//layerOutput = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::Graded); // output

	FullConnectivity* full2 = new FullConnectivity();
	FullConnectivity* full4 = new FullConnectivity();
	FullConnectivity* full5 = new FullConnectivity();
	FullConnectivity* full6 = new FullConnectivity();
	FullConnectivity* full8 = new FullConnectivity(true,"hypercolumn");
//	FullConnectivity* full10 = new FullConnectivity();

	if(addInputLayerToNetwork == true)
		network->AddPopulation(layerInput);

	network->AddPopulation(layer2);
//	network->AddPopulation(layerOutput);

	layerInput->AddPre(layerInput,full4); // Recurrent, with MI minicolumn calculations
	layerInput->AddPre(layerInput,full8); // Recurrent, with MI hypercolumn calculations + MDS
	layer2->AddPre(layerInput,full2); // Feedforward, modified by VQ calculations and BCPNN
//	layer2->AddPre(layer2,full6); // Recurrent, with inhib BCPNN (if used)
//	layerOutput->AddPre(layer2,full10); // Feedforward, modified by Kussul or BCPNN calculations (supervised learning for classification) (if used)

	int mdsDimension = m_mdsDimension;
	int miDimension = layerInput->GetNrHypercolumns();//layerInput->GetNrRateUnits().size();//nrInputHypercolumns;

	if(m_mdsUsePearson == false) // default
	{
		// MI
		m_miHypercolumns = new ProjectionModifierMIHypercolumn();
		m_miRateUnits = new ProjectionModifierMIRateUnit(m_miHypercolumns);
		full8->AddProjectionsEvent(m_miHypercolumns);
		full4->AddProjectionsEvent(m_miRateUnits);
		m_miRateUnits->AddParentProjectionModifier(m_miHypercolumns); // allows mi hypercolumns to have access to the belonging mi minicolumns (set as default?)
	}
	else
	{
		m_pearson = new ProjectionModifierPearson();
		full4->AddProjectionsEvent(m_pearson);
		/*vector<int> mcs = layerInput->GetNrRateUnits();
		int sum = 0;
		for(int i=0;i<mcs.size();i++) sum += mcs[i];
			miDimension = sum;*/
		miDimension = layerInput->GetNrHypercolumns()*layerInput->GetNrRateUnitsInHypercolumn();
	}

	// MDS
	m_MDS = new LayerMDS(miDimension,mdsDimension, this->m_network);
	//m_network->SetSeed(); // go back to default (in network initialize now, important for building random networks etc, so unnecessary here)
	m_mdsHypercolumns = new ProjectionModifierMDS(m_mdsUsePearson);
	layerInput->AddPopulationModifier(m_MDS);
	
	//if(m_mdsUsePearson == false)
	m_mdsHypercolumns->AddParentPopulationModifier(m_MDS); // allows MDS to have access to the hypercolumn event Projections (will be set as default)
	
	if(m_mdsUsePearson == true)
		full4->AddProjectionsEvent(m_mdsHypercolumns);
	else
		full8->AddProjectionsEvent(m_mdsHypercolumns);

	// VQ
	int nrGroups = nrMiddleHypercolumns;
	m_VQ = new LayerVQ(nrGroups, LayerVQ::VQCSL);
	full2->AddProjectionsEvent(m_VQ->GetProjectionModifier()); // Feedforward modified by VQ calculations
	layerInput->AddPopulationModifier(m_VQ);
	m_VQ->AddChildPopulationModifier(m_MDS); // Allow VQ to have access to MDS output (m_Xi)

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
	m_cslLearn = NULL;

	if(m_featureExtraction == StructureMIMDSVQ::CL)
	{
		m_compLearn = new ProjectionModifierCL(nrMiddleRateUnits,0.0001,0.0005,clC);//nrMiddleRateUnits*100);
		full2->AddProjectionsEvent(m_compLearn);
	}
	else if(m_featureExtraction == StructureMIMDSVQ::BCM)
	{
		m_bcmLearn = new ProjectionModifierBCM();
		full2->AddProjectionsEvent(m_bcmLearn);
	}
	else if(m_featureExtraction == StructureMIMDSVQ::HebbAdapt)
	{
		float alpha = 0.01;
		float lambda = 0.01*0.01;

		m_bcpnn = new ProjectionModifierBcpnnOnline(alpha,lambda);
		full2->AddProjectionsEvent(m_bcpnn);

		float adampl = 1;
		float adtaudt = 0.01;
		m_adaptation = new PopulationModifierAdaptation(adampl,adtaudt);
		layer2->AddPopulationModifier(m_adaptation);
	}
	else if(m_featureExtraction == StructureMIMDSVQ::Recurr)
	{
		m_featExtrRecurrProb = 0.2;

		RandomConnectivity* randConn = new RandomConnectivity(m_featExtrRecurrProb);
		layer2->AddPre(layer2,randConn);
		randConn->SetRandomWeights(-10,-10);
		full2->SetRandomWeights(0,1);
		layer2->AddUnitsModifierToInitialize(new TransferLinear(false));//layer2->GetIncomingProjections()[1]->AddUnitsProperty(new TransferLinear(false));//layer2->AddUnitsProperty(new TransferLinear(false));

		// feedforward plasticity not utilized yet
		float alpha = 0.01;
		float lambda = 0.01*0.01;
		m_bcpnn = new ProjectionModifierBcpnnOnline(alpha,lambda);
	}
	else
	{
		m_cslLearn = new ProjectionModifierCSL();
		full2->AddProjectionsEvent(m_cslLearn);
	}

	if(m_featureExtraction != StructureMIMDSVQ::Recurr)
	{
		SoftMax* softmax = new SoftMax(1.0, SoftMax::WTA);
		m_wta = new WTA();
		//layer2->AddPopulationModifier(m_wta);
	}

	// Initializes new parts of the network (will not re-initialize old parts)
	//network->Initialize();

	m_layers.push_back(layerInput);
	m_layers.push_back(layer2);
	//m_layers.push_back(layerOutput);
}

void StructureMIMDSVQ::SetupMeters(int mpiRank, int mpiSize)//, Storage::FilePreference fileType)
{
	char* name1 = new char[30];
	sprintf(name1,"Projection_%d_n%d.csv",m_index,mpiRank);
	
	PopulationColumns* l2 = m_layers[1];

	// if run on too many nodes compared to network size, could get error here
	Meter* connMeter = new Meter(name1, Storage::CSV,Storage::Append);
	connMeter->AttachProjection(m_layers[1]->GetIncomingProjections()[0],0);//AttachUnit(layer1->GetRateUnits()[0]);
	m_network->AddMeter(connMeter);
	
	//char* name2 = "Layer2Activity_2.csv";
	char* name2 = new char[30];
	sprintf(name2,"Layer2Activity_%d.csv",m_index);

	//char* name3 = "Layer3Activity.csv";
	//char* name3 = new char[20];
	//sprintf(name3,"Layer3Activity_%d.csv",m_index);

	//char* nameTest2 = "Layer2Activity_test.csv";
	//char* nameTest3 = "Layer3Activity_test.csv";

	Meter* layerMeter = new Meter(name2, Storage::CSV);//Storage::MPI_Binary);//fileType);//Storage::CSV);
	//Meter layer3Meter(name3, Storage::CSV);

	layerMeter->AttachPopulation(m_layers[1]);
	m_network->AddMeter(layerMeter);

/*	AnalysisMDS* analysisMDS = new AnalysisMDS(m_MDS);
	m_network->AddAnalysis(analysisMDS);

	AnalysisVQ* analysisVQ = new AnalysisVQ(m_VQ);
	m_network->AddAnalysis(analysisVQ);

	AnalysisVQ* analysisVQ2 = new AnalysisVQ(this->CSLLearn());//m_VQ);
	m_network->AddAnalysis(analysisVQ2);*/

	

/*	char* name5 = new char[30];
	sprintf(name5,"mi_%d_%d.csv",m_index,mpiRank);
	Meter* miMeter = new Meter(name5,Storage::CSV, Storage::Append);
	miMeter->AttachObject((NetworkObject*)m_miHypercolumns);
	m_network->AddMeter(miMeter);*/	

	if(mpiRank == 0)
	{
		//layer3Meter.AttachPopulation(m_layers[2]);
		//m_network->AddMeter(&layer3Meter);

		char* name5 = new char[30];
		sprintf(name5,"mds_%d.csv",m_index);
		Meter* mdsMeter = new Meter(name5,Storage::CSV, Storage::Standard);
		mdsMeter->AttachObject((NetworkObject*)m_MDS, Meter::MeterPopulationModifier); // error: is not attached to layer ?
		m_network->AddMeter(mdsMeter);

		char* name4 = new char[30];
		sprintf(name4,"vq_%d.csv",m_index);
		Meter* vqMeter = new Meter(name4, Storage::CSV, Storage::Standard);
		vqMeter->AttachObject((NetworkObject*)m_VQ,Meter::MeterPopulationModifier);//(PopulationModifier*)m_VQ); // error: is not attached to layer
		m_network->AddMeter(vqMeter);
		
	}

	//}
}

void StructureMIMDSVQ::SetRecording(bool on)
{
	m_layers[1]->SetRecording(on);
	//m_VQ->SetRecording(on);
	if(m_mdsUsePearson == false)
		m_miHypercolumns->SetRecording(on);
	m_layers[1]->GetIncomingProjections()[0]->SetRecording(on,0);
}

void StructureMIMDSVQ::SetTiming(bool on)
{
/*	m_network->AddTiming(VQ());
	m_network->AddTiming(MDS());
	m_network->AddTiming(m_mdsHypercolumns);
	m_network->AddTiming(m_cslLearn);
	m_network->AddTiming(m_wta);
	m_network->AddTiming(MIHypercolumns());
	m_network->AddTiming(MIRateUnits());

	for(int i=0;i<m_layers.size();i++)
		m_network->AddTiming(m_layers[i]);
		*/
}
