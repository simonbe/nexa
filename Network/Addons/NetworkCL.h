#pragma once
#ifndef NETWORKCL_H
#define NETWORKCL_H

#include "NetworkProjectionModifier.h"

/// <summary>	Competitive learning implementation with adaptive threshold (see papers by DeSieno). </summary>

class ProjectionModifierCL : public ProjectionModifier
{
public:
	ProjectionModifierCL(int sizeOutLayer, float eta, float probB, float biasC);

	void Initialize(Projection* Projection);
	void SetProjection(Projection* c);
		
	void SetAlpha(float alpha);
	void SetLambda(float lambda);

	void Simulate(UnitModifier* e);
	void Modify();

	void SetC(float value);

private:

	bool m_firstRun;

	float m_eta;
	
	// bias
	vector<float> m_p;	// probabilities for post units active
	vector<float> m_b;	// biases
	float B;			// speed in est of m_p
	float C;			// bias factor (distance a losing element can reach in order to enter the solution)
	float N;			// nr elements in competitive layer
	
	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;

	Projection* m_projectionFixed;
};


#endif