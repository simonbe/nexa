#pragma once
#ifndef NETWORKBCM_H
#define NETWORKBCM_H

#include "NetworkConnectionModifier.h"

// BCM, Law and Cooper (1994) implementation (derived from ibcm (Intrator and Cooper 1992))

class ConnectionModifierBCM : public ConnectionModifier
{
public:

	ConnectionModifierBCM();
	ConnectionModifierBCM(float eta, float decay, float tau);

	void SetEta(float eta);
	void Initialize(Connection* connection);
	void SetConnection(Connection* c);
	void Simulate(UnitModifier* e);
	void Modify();
	void Clear();

private:

	bool m_firstRun;
	float m_eta, m_tau, m_decay;

	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;
	vector<float> m_thresholds;

	Connection* m_connectionFixed;
};

#endif