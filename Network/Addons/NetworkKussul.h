#pragma once
#ifndef NETWORKKUSSUL_H
#define NETWORKKUSSUL_H

#include "NetworkProjectionModifier.h"

/// <summary>	Kussul learning rule (delta rule). </summary>

class ProjectionModifierKussul : public ProjectionModifier
{
public:

	ProjectionModifierKussul();

	void Initialize(Projection* Projection);

	void SetProjection(Projection* c);
	void SetAlpha(float alpha);
	
	void Simulate(UnitModifier* e);
	void Modify();

private:

	float m_alpha;
	bool m_firstRun;
	
	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;

	Projection* m_projectionFixed;
};


#endif