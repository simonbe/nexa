#pragma once
#ifndef NETWORKHEBBSMP_H
#define NETWORKHEBBSMP_H

#include "NetworkConnectionModifier.h"

// Hebbian with weight normalization and/or weight decay

class ConnectionModifierHebbSimple : public ConnectionModifier
{
public:

	ConnectionModifierHebbSimple();
	ConnectionModifierHebbSimple(float etaHebb, bool normalize);
	ConnectionModifierHebbSimple(float etaHebb, bool normalize, UnitModifier* transfer);

	void SetEtaHebb(float etaHebb);

	void Initialize(Connection* connection);
	void SetConnection(Connection* c);
	void Simulate(UnitModifier* e);
	void Modify();
	void Clear();

	void SetNormalizationFactor(float factor)
	{
		m_normalizationFactor = factor;
	}

	void Normalize(bool allWeights);

private:

	bool m_firstRun;

	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;

	float m_etaHebb;
	bool m_normalize;
	float m_normalizationFactor;

	Connection* m_connectionFixed;
};

#endif