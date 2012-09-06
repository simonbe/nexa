#pragma once
#ifndef NETWORKSTDP_H
#define NETWORKSTDP_H

#include "NetworkProjectionModifier.h"

class ProjectionModifierSTDP : public ProjectionModifier
{
public:
	ProjectionModifierSTDP(bool inhibitoryWeights);

	void Initialize(Projection* Projection);
	void SetProjection(Projection* c) {}
	void Simulate(UnitModifier* e) {}
	void Modify();
	void Clear();

private:

	float synapse_stdp(float I, float w, float input);

	bool m_firstRun;
	float Wmax, Wmin, Wadj;
	float m_scaleFactor;

	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;
	
	vector<vector<float> > preexp_decay2;
	vector<vector<float> > preexp_decay1;
	//vector<float> preexp_decay2;
	vector<float> m_prevPostValues;
	
	//map<long, map<long, SynapseStandard> > m_hashSynapses;

	Projection* m_projectionFixed;
};

#endif