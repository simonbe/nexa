#pragma once
#ifndef NETWORKCORR_H
#define NETWORKCORR_H

#include "NetworkConnectionModifier.h"


class ConnectionModifierPearson : public ConnectionModifier
{
public:

	ConnectionModifierPearson()
	{
		m_eventId = 11;
	}

	void SetConnection(Connection* c)
	{
		m_connection = c;
	}

	void Initialize(Connection* connection);
	void Simulate(UnitModifier* e){};
	void Modify();

	vector<vector<float> > GetRij() { return rIJ; }

private:

	vector<float> meanI,meanJ, varianceI, varianceJ, variance_nI, variance_nJ, M2I, M2J, XmeanI, XmeanJ;
	vector<vector<float> > rIJ, meanIJ;
	
	RateUnit* m_pre;
	RateUnit* m_post;

	int m_n;
};

#endif