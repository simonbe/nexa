#pragma once
#ifndef NETWORKTRC_H
#define NETWORKTRC_H

#include "NetworkProjectionModifier.h"

// Hebbian + intrinsic (threshold) plasticity

class ProjectionModifierTriesch : public ProjectionModifier
{
public:

	ProjectionModifierTriesch();
	ProjectionModifierTriesch(float etaHebb, float beta, float etaIP, float mu, bool thresholded);
	void Initialize(Projection* Projection);
	void SetProjection(Projection* c);
	void SetBeta(float beta)
	{
		m_beta = beta;
	}

	void SetMu(float mu);
	void Simulate(UnitModifier* e);
	void Modify();
	void Clear();

private:

	float N(float y, bool isMaxValue); // neighborhood fcn

	///
	bool m_firstRun;
	float m_etaHebb, m_beta;

	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;

	Projection* m_projectionFixed;
};

#endif