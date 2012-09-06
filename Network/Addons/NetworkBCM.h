#pragma once
#ifndef NETWORKBCM_H
#define NETWORKBCM_H

#include "NetworkProjectionModifier.h"

/// <summary>	BCM implementation, see e.g. Law and Cooper (1994) (derived from ibcm (Intrator and Cooper 1992)) </summary>

class ProjectionModifierBCM : public ProjectionModifier
{
public:

	ProjectionModifierBCM();
	ProjectionModifierBCM(float eta, float decay, float tau);

	void SetEta(float eta);
	void Initialize(Projection* Projection);
	void SetProjection(Projection* c);
	void Simulate(UnitModifier* e);
	void Modify();
	void Clear();

private:

	bool m_firstRun;
	float m_eta, m_tau, m_decay;

	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;
	vector<float> m_thresholds;

	Projection* m_projectionFixed;
};

#endif