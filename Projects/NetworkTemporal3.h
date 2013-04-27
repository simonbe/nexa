#pragma once


#include <iostream>
#include <vector>

#include "StructureMIMDSVQ.h"
#include "AssignmentStrategy.h"

class Network;
class PopulationColumns;

using namespace std;


class NetworkTemporal3 : public Network
{
public:
	// setters
	void SetInputLayerSize(int size);
	void SetLayer2Size(int size);
	void SetMDSSize(int size);

	// getters

	string GetMDSInput();

private:
	// overridden network functions
	void NetworkSetupStructure();
	void NetworkInitParameters();
	void NetworkSetupMeters();
	void NetworkSetupParameters();
	void NetworkRun();

	// training process functions
	void ComputeCorrelation(const vector<vector<float> >& trainingData,  PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iter);
	void ComputeMDS(const vector<vector<float> >& trainingData,  PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iter);
	void ComputeVQ(const vector<vector<float> >& trainingData,  PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iter);
	void ExtractFeatures(const vector<vector<float> >& trainingData,  PopulationColumns* inputLayer, StructureMIMDSVQ* structure);
	
	// Train
	void TrainLayer(const vector<vector<float> >& trainingData, PopulationColumns* inputLayer, StructureMIMDSVQ* structure, int iterationsCorrs, int iterationsMDS, int iterationsVQ, int iterationsFeatures);

	// other
	vector<float> toBinary(const vector<float>& data, int nrHc, int nrMc);
	//void ClearActivities();

	


	// Run
	bool m_run;
	int m_architecture;
	bool m_verbose;
	// Populations
	PopulationColumns* m_layerInput;
	PopulationColumns* m_layer2;

	// Plasticity

	// Pre-defined structures
	//vector<StructureReTIDe*> m_reTIDe;

	// Assignment strategy
	IAssignmentStrategy* m_dataAssigner;
	
	// Sizes
	int m_sizePopulation1, m_sizePopulation2, m_sizePopulation3;
	int m_nrHypercolumns2, m_nrRateUnits2;
	int m_nrInputHypercolumns, m_nrInputRateUnits;
	int m_MDSInput;

	int m_nrColumns;

	// Structures
	StructureMIMDSVQ* m_structureInput;
	// Probabilities
	//float m_probRecurr, m_probRecurr2, m_probForward;
	Meter* m_inputMeter,*m_layer1Meter,*m_layer2Meter;
	// Strengths
	float m_plasticityStrength;
};