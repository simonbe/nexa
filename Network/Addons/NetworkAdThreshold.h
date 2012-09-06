#pragma once
#ifndef NETWORKADTHR_H
#define NETWORKADTHR_H

#include "NetworkProjectionModifier.h"

// Hebbian + intrinsic (threshold) plasticity

class AdaptiveThreshold : public PopulationModifier
{
public:

	enum ThrType
	{
		Triesch,
		BCPNNBeta,
		Average,
		Variability
	};

	AdaptiveThreshold(ThrType type)
	{
		m_firstRun = true;

		m_thresholdType = type;
		if(type == this->Triesch)
		{
			m_isThresholded = false;
			eta_ip = 0.005;//etaIP;//0.005;
			mu = 0.1;//mu_;//0.1; // desired activity
		}
		else if(type == this->BCPNNBeta)
		{
			m_alpha = 0.01;
			m_lambda0 = 0.0001;
			m_impactBeta = -2;
		}
		else if(type == this->Average)
		{
			m_impactAverage = 1;
			m_averageAim = 0.1;
			m_alpha = 0.01;
			m_lambda0 = 0.01;
		}
	}
	
	void Simulate();

	void Modify();

	void SetImpactAverage(float impact)
	{
		m_impactAverage = impact;
	}

	float GetImpactAverage()
	{
		return m_impactAverage;
	}

private:
	
	ThrType m_thresholdType;
	bool m_firstRun;

	// Triesch
	map<long,float> a,b;
	bool m_isThresholded;
	float eta_ip, mu;

	// BCPNN beta + average
	vector<float> m_Aj;
	float m_alpha, m_lambda0;
	float m_impactBeta;

	// Average
	vector<float> m_average;
	vector<float> m_averageGain;
	float m_impactAverage;
	float m_averageAim;

	// Variability
};

#endif