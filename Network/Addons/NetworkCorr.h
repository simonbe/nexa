#pragma once
#ifndef NETWORKCORR_H
#define NETWORKCORR_H

#include "NetworkProjectionModifier.h"

/// <summary>	Incremental Pearson correlation. </summary>

class ProjectionModifierPearson : public ProjectionModifier
{
public:

	ProjectionModifierPearson()
	{
		m_eventId = 11;
	}

	void SetProjection(Projection* c)
	{
		m_projection = c;
	}

	void Initialize(Projection* Projection);
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