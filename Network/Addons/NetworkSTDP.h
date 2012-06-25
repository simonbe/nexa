#pragma once
#ifndef NETWORKSTDP_H
#define NETWORKSTDP_H

#include "NetworkConnectionModifier.h"

class ConnectionModifierSTDP : public ConnectionModifier
{
public:
	ConnectionModifierSTDP(bool inhibitoryWeights);

	void Initialize(Connection* connection);
	void SetConnection(Connection* c) {}
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

	Connection* m_connectionFixed;
};

#endif