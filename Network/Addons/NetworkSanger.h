#pragma once
#ifndef NETWORKSANGER_H
#define NETWORKOJA_H

#include "NetworkConnectionModifier.h"

// Sanger's rule / Generalized Hebbian Algorithm (GHA) / PCA rule

class ConnectionModifierSanger : public ConnectionModifier
{
public:

	ConnectionModifierSanger();
	ConnectionModifierSanger(float etaHebb);

	void Initialize(Connection* connection);
	void SetConnection(Connection* c);
	void Simulate(UnitModifier* e);
	void Modify();
	void Clear();

private:

	bool m_firstRun;
	float m_etaHebb;

	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;

	Connection* m_connectionFixed;
};

#endif