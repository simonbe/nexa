#pragma once
#ifndef NETWORKHEBBSMP_H
#define NETWORKHEBBSMP_H

#include "NetworkProjectionModifier.h"

/// <summary>	Plain Hebbian learning with weight normalization and/or weight decay </summary>

class ProjectionModifierHebbSimple : public ProjectionModifier
{
public:

	ProjectionModifierHebbSimple();
	ProjectionModifierHebbSimple(float etaHebb, bool normalize);
	ProjectionModifierHebbSimple(float etaHebb, bool normalize, UnitModifier* transfer);

	void SetEtaHebb(float etaHebb);

	void Initialize(Projection* Projection);
	void SetProjection(Projection* c);
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

	Projection* m_projectionFixed;
};

#endif