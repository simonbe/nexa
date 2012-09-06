#pragma once
#ifndef NETWORKFOLDIAK_H
#define NETWORKFOLDIAK_H

#include "NetworkProjectionModifier.h"

/// <summary>	Foldiak learning rule. </summary>

class ProjectionModifierFoldiak : public ProjectionModifier
{
public:

	ProjectionModifierFoldiak(float eta1, float eta2, float eta3, float alpha, float beta, bool lateral);

	void Initialize(Projection* Projection);

	void SetProjection(Projection* c);
	void SetAlpha(float alpha);

	void Simulate(UnitModifier* e);
	void Modify();

private:

	bool m_firstRun;
	bool m_lateral;

	float m_eta1,m_eta2,m_eta3,m_alpha, m_beta;
	vector<float> m_s;

	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;

	Projection* m_projectionFixed;
};


#endif