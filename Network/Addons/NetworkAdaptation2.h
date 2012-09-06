#pragma once
#ifndef NETWORKADAPT_H
#define NETWORKADAPT_H

#include "NetworkProjectionModifier.h"

class PopulationModifierAdaptation2 : public PopulationModifier
{

public:
	
	PopulationModifierAdaptation2()
	{
		m_adtaudt = 0.1;
		// values of m_alpha and m_lambda0 taken from BCPNN standard initialization values: 0.01, 0.0001, resp.
		m_alpha = 0.01; // learning rate, or clipping constant for adaption of m_Aj
		m_lambda0 = 0.001;
		m_impactBeta = 1; // in Simon's BCPNN code the value is 1, however this seems to be too high 
		m_Aj = vector<float>(0);  // in order to force initialization condition
	}

	void SetParameters(float adtaudt,float alpha,float lambda)
	{		
		m_adtaudt = adtaudt;
		m_alpha=alpha;
		m_lambda0=lambda;
	}

	void SetParameters(float adtaudt)
	{		
		m_adtaudt = adtaudt;
	}

	void Simulate();
	void Initm_Aj(float aj);  // new introduced to adapt for Anders' formulas
	vector<vector<float> > GetValuesToRecord(); 

private:
	vector<float> switchdata; // this is for outputting before WTA
	vector<float> m_adaptValues,m_adampl;
	float m_adtaudt;

	// following variable definitions needed for calculation of Bj, the beta term of BCPNN, which is included in Anders' adaptation formula
	// these need to be set on initialization
	vector<float> m_Aj; 
	float m_alpha,m_lambda0;
	// comments: 
	// m_alpha can be set either by a) to standard values by ProjectionModifierBcpnnOnline(float alpha, float lambda) or b) SetAlpha(float alpha)
	// m_lambda0 can be set by a) ProjectionModifierBcpnnOnline(float alpha, float lambda) or b) SetLambda(float lambda)
	float m_impactBeta;  
};

#endif