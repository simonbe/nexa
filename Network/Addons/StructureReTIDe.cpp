#include "StructureReTIDe.h"
#include "NetworkUnitModifier.h"

void StructureReTIDe::Initialize(Network* network, PopulationColumns* layerInput)//,float weightStrength, float weightStrength2)
{
	m_network = network;
	//PopulationColumns* layer2 = new PopulationColumns(network,nrMiddleHypercolumns,nrMiddleRateUnits,PopulationColumns::Graded); // middle
	//layerOutput = new PopulationColumns(network,nrOutputHypercolumns,nrOutputRateUnits,PopulationColumns::Graded); // output

	RandomConnectivity* sparseConn = new RandomConnectivity(m_probRecurr, m_symmetric);//true); // true (default): symmetric
	//FullConnectivity* sparseConn = new FullConnectivity();
	sparseConn->SetRandomWeights(m_weightStrength,m_weightStrength);//-1,-1); // inhib

	layerInput->AddPre(layerInput,sparseConn); // Recurrent sparse (inhib)

/*	if(m_weightStrength3 !=0.0 || m_weightStrength3 == 0.0)
	{
		FullConnectivity* full2 = new FullConnectivity();
		full2->SetRandomWeights(m_weightStrength2,m_weightStrength2);
		layerInput->AddPre(layerInput,full2);
		FullConnectivity* full3 = new FullConnectivity();
		full3->SetRandomWeights(m_weightStrength3,m_weightStrength3);
		layerInput->AddPre(layerInput,full3);
	}
	else if(m_weightStrength2 != 0.0)
	{
		RandomConnectivity* sparseConn2 = new RandomConnectivity(m_probRecurr);
		sparseConn2->SetRandomWeights(m_weightStrength2,m_weightStrength2);
		layerInput->AddPre(layerInput,sparseConn2); // Recurrent sparse (exc)
	}
	*/
	m_layers.push_back(layerInput);
}

void StructureReTIDe::SetupStructure(Network* network, PopulationColumns* layerInput, bool useDivNormalization)
{
	TransferReTIDe* transferReTIDe = new TransferReTIDe(layerInput,m_threshold,useDivNormalization, m_maxValue);//,m_maxValue);
	//TransferReTIDe* transferReTIDe2 = new TransferReTIDe(m_threshold2);
	//TransferReTIDe* transferReTIDe3 = new TransferReTIDe(m_threshold3);
	//layerInput->GetIncomingProjections()[0]->AddUnitsProperty(transferReTIDe);//layerInput->AddUnitsProperty(transferReTIDe);
	//layerInput->GetIncomingProjections()[1]->AddUnitsProperty(transferReTIDe2);//
	//layerInput->GetIncomingProjections()[2]->AddUnitsProperty(transferReTIDe3);//
	layerInput->AddUnitsProperty(transferReTIDe);
}

void StructureReTIDe::SetupMeters(int index, Storage::SaveDataState state)
{
//	char* name1 = new char[30];
//	sprintf(name1,"Projection_%d_n%d.csv",m_index,mpiRank);
//	Meter* connMeter = new Meter(name1, Storage::CSV);

	// if run on too many nodes compared to network size, could get error here
	//connMeter->AttachProjection(m_layers[1]->GetIncomingProjections()[0]);//AttachUnit(layer1->GetRateUnits()[0]);
	//m_network->AddMeter(connMeter);

	//char* name2 = "Layer2Activity_2.csv";
	char* name2 = new char[30];
	m_index = index;
	sprintf(name2,"LayerReTideActivity_%d.csv",m_index);

	//char* name3 = "Layer3Activity.csv";
	//char* name3 = new char[20];
	//sprintf(name3,"Layer3Activity_%d.csv",m_index);

	//char* nameTest2 = "Layer2Activity_test.csv";
	//char* nameTest3 = "Layer3Activity_test.csv";

	m_layerMeter = new Meter(name2, Storage::CSV, state);
	//Meter layer3Meter(name3, Storage::CSV);

	m_layerMeter->AttachPopulation(m_layers[0]);
	m_network->AddMeter(m_layerMeter);
}

Meter* StructureReTIDe::GetMeter()
{
	return m_layerMeter;
}

void StructureReTIDe::SetRecording(bool on)
{
	m_layers[0]->SetRecording(on);
}