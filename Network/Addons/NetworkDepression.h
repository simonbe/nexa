#pragma once
#ifndef NETWORKEVCONNDEP_H
#define NETWORKEVCONNDEP_H

#include "NetworkConnectionModifier.h"

// Synaptic depression

class ConnectionModifierDepression : public ConnectionModifier
{
public:
	ConnectionModifierDepression();
	ConnectionModifierDepression(float strength);
	void SetConnection(Connection* c);
	void Simulate(UnitModifier* e);
	void Modify();
	void Clear();

private:
	vector<vector<long> > m_orgWeights;
	bool m_hasStoredOrgWeights;
	float m_strength;
};

#endif
