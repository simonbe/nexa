#pragma once
#ifndef NETWORKBCPNN_H
#define NETWORKBCPNN_H

#include "NetworkProjectionModifier.h"

/// <summary>	Implements the BCPNN learning rule. 
/// 			See e.g. 
/// 			Sandberg A., Lansner A., Petersson K.-M. and Ekeberg . (2002): Bayesian
///				attractor networks with incremental learning. Network: Computation in
///				Neural Systems: 13(2), 179-194. </summary>

class ProjectionModifierBcpnnOnline : public ProjectionModifier
{
public:

	ProjectionModifierBcpnnOnline(float alpha, float lambda, float maxValue = -1);
	ProjectionModifierBcpnnOnline();

	void Initialize(Projection* Projection);
	void SetProjection(Projection* c);
	
	void SetImpactWeights(float impact)
	{
		m_impactWeights = impact;
	}

	void SetImpactBeta(float impact)
	{
		m_impactBeta = impact;
	}

	void SetAlpha(float alpha);
	void SetLambda(float lambda);
	
	vector<float> GetBeta()
	{
		return m_beta;
	}

	void Simulate(UnitModifier* e);
	void Modify();

	void UseNoUpdateIfNoPreActivity(bool use)
	{
		m_useNoUpdateIfNoPreActivity = use;
	}

private:

	bool m_firstRun;
	bool m_useNoUpdateIfNoPreActivity;

	// map implementation
	map<long,float> m_Ai;
	vector<float> m_Aj;
	vector<float> m_beta;
	vector<float> m_inhibBeta;
	vector<vector<float> > m_Aij; // may need to restructure if synapses moved
	//vector<map<long,float> > m_Aij;

	float m_alpha;
	float m_lambda0;

	// multiplies in beta and weights
	float m_impactWeights;
	float m_impactBeta;

	vector<long> m_postIds;

	/// <summary>	Projection it is used on. </summary>
	Projection* m_projectionFixed;
};

#endif