#pragma once


#include <iostream>
#include <vector>

#include "StructureMIMDSVQ.h"

class Network;
class PopulationColumns;

using namespace std;

class NetworkMNIST2 : public Network
{
public:

private:
	// overridden network functions
	void NetworkSetupStructure();
	void NetworkSetupMeters();
	void NetworkSetupParameters();
	void NetworkRun();

	// Train
	void TrainLayer(vector<vector<float> > trainingData, PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iterationsPatches, int iterationsFeatures);

	// other
	vector<float> toBinary(vector<float> data, int nrHc, int nrMc);
	//void ClearActivities();

	// Run
	bool m_run;
	int m_architecture;

	// Populations
	PopulationColumns* m_layerInput;
	PopulationColumns* m_layer2, *m_layer3;

	// Plasticity
	ProjectionModifierBcpnnOnline* m_bcpnn;

	// Pre-defined structures
	//vector<StructureReTIDe*> m_reTIDe;
	
	// Sizes
	int m_sizePopulation1, m_sizePopulation2, m_sizePopulation3;
	int m_nrHypercolumns2, m_nrRateUnits2, m_nrHypercolumns3, m_nrRateUnits3;
	int m_nrInputHypercolumns, m_nrInputRateUnits;

	int m_nrColumns;

	// Structures
	StructureMIMDSVQ* m_structureInput, *m_structureLayer2;

	// Meters
	Meter* m_inputMeter,*m_layer1Meter,*m_layer2Meter;

	// Probabilities
	//float m_probRecurr, m_probRecurr2, m_probForward;

	// Strengths
	float m_plasticityStrength;
};

/*
class NetworkMNIST
{
public:
	NetworkMNIST();

	void NetworkMNISTRun1(int mpiRank, int mpiSize);
	void NetworkMNISTRun2(int mpiRank, int mpiSize);
	void NetworkMNISTRun3(int mpiRank, int mpiSize);

	void NetworkMNISTRunLateralMaps(int mpiRank, int mpiSize);
	void NetworkMNISTRunLateralMaps2(int mpiRank, int mpiSize);

private:

	vector<float> toBinary(int nr, int total);
	vector<float> toBinary(vector<float> data, int nrHc, int nrMc);

	vector<PopulationColumns*> AttachPopulations(Network* network, PopulationColumns* layerInput, int nrMiddleHypercolumns, int nrMiddleRateUnits, int nrOutputHypercolumns, int nrOutputRateUnits);

	int m_architecture;
};*/