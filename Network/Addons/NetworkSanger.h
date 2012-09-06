#pragma once
#ifndef NETWORKSANGER_H
#define NETWORKOJA_H

#include "NetworkProjectionModifier.h"

// Sanger's rule / Generalized Hebbian Algorithm (GHA) / PCA rule

class ProjectionModifierSanger : public ProjectionModifier
{
public:

	ProjectionModifierSanger();
	ProjectionModifierSanger(float etaHebb);

	void Initialize(Projection* Projection);
	void SetProjection(Projection* c);
	void Simulate(UnitModifier* e);
	void Modify();
	void Clear();

private:

	bool m_firstRun;
	float m_etaHebb;

	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;

	Projection* m_projectionFixed;
};

#endif