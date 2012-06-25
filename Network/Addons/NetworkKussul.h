#pragma once
#ifndef NETWORKKUSSUL_H
#define NETWORKKUSSUL_H

#include "NetworkConnectionModifier.h"


class ConnectionModifierKussul : public ConnectionModifier
{
public:

	ConnectionModifierKussul();

	void Initialize(Connection* connection);

	void SetConnection(Connection* c);
	void SetAlpha(float alpha);
	
	void Simulate(UnitModifier* e);
	void Modify();

private:

	float m_alpha;
	bool m_firstRun;
	
	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;

	Connection* m_connectionFixed;
};


#endif