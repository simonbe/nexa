#pragma once
#ifndef NETWORKEVCONNDEP_H
#define NETWORKEVCONNDEP_H

#include "NetworkProjectionModifier.h"

/// <summary>	Synaptic depression. </summary>

class ProjectionModifierDepression : public ProjectionModifier
{
public:
	ProjectionModifierDepression();
	ProjectionModifierDepression(float strength);
	void SetProjection(Projection* c);
	void Simulate(UnitModifier* e);
	void Modify();
	void Clear();

private:
	vector<vector<long> > m_orgWeights;
	bool m_hasStoredOrgWeights;
	float m_strength;
};

#endif
