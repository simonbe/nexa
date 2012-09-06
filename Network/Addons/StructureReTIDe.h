#pragma once
#ifndef STRUCTRETIDE_H
#define STRUCTRETIDE_H

#include "NetworkStructure.h"
#include "Network.h"
#include "Meter.h"

class StructureReTIDe : public NetworkStructure
{
public:

	StructureReTIDe(float probRecurr, float weightStrength, float threshold = 1.0, bool symmetric = true, float maxValue = -1, float weightStrength2 = 0.0, float threshold2 = 1.0, float weightStrength3 = 0.0, float threshold3 = 1.0)
	{
		m_probRecurr = probRecurr;
		m_threshold = threshold;
		m_maxValue = maxValue;
		m_weightStrength = weightStrength;
		m_weightStrength2 = weightStrength2;
		m_threshold2 = threshold2;
		m_weightStrength3 = weightStrength3;
		m_threshold3 = threshold3;
		m_symmetric = symmetric;
	}

	void Initialize(Network* network, PopulationColumns* layerInput);
	void SetupStructure(Network* network, PopulationColumns* layerInput, bool useDivNormalization = false);
	void SetupMeters(int index, Storage::SaveDataState state);
	void SetRecording(bool on);

	PopulationColumns* GetLayer(int index)
	{
		return m_layers[index];
	}

	Meter* GetMeter();

protected:

	bool m_symmetric;

	Network* m_network;
	vector<PopulationColumns*> m_layers;
	float m_probRecurr;
	float m_threshold, m_threshold2,m_threshold3;
	float m_weightStrength,m_weightStrength2,m_weightStrength3;

	Meter* m_layerMeter;
	float m_maxValue;

	int m_index;
};

#endif