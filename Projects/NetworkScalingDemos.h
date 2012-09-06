#pragma once

#include <iostream>
#include <vector>
#include "StructureReTIDe.h"

class Network;
class PopulationColumns;

using namespace std;

class NetworkScalingDemos
{
public:
	NetworkScalingDemos();

	void NetworkScalingRunRecurrentSimple(bool storeData);
	void NetworkScalingRunRecurrentBCPNN(bool storeData);
	void NetworkScalingRunSpiking();
	void NetworkScalingRunMIMDSVQ();

	void RunAll();

private:

	int m_architecture;
	bool m_useBinaryFileWrite;
	bool m_allowFileWrite;
	string m_extraFilenameString;
};

class NetworkScalingStrong : public Network
{
public:

private:
	// overridden network functions
	void NetworkSetupStructure();
	void NetworkSetupMeters();
	void NetworkSetupParameters();
	void NetworkRun();

	// other
	void ClearActivities();

	// Run
	int m_run;

	// Populations
	PopulationColumns* m_layer1;
	vector<PopulationColumns*> m_layer2;
	PopulationColumns* m_layer3;

	// Plasticity
	ProjectionModifierBcpnnOnline* m_bcpnn;

	// Pre-defined structures
	vector<StructureReTIDe*> m_reTIDe;
	
	// Sizes
	int m_sizePopulation1, m_sizePopulation2, m_sizePopulation3;
	int m_nrColumns;

	// Probabilities
	float m_probRecurr, m_probRecurr2, m_probForward;

	// Strengths
	float m_plasticityStrength;
};