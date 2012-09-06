#pragma once
#ifndef NETWORKADAPT_H
#define NETWORKADAPT_H

#include "NetworkProjectionModifier.h"

class PopulationModifierAdaptation : public PopulationModifier
{

public:
	
	enum AdaptationType
	{
		Standard,
		StandardMinusBeta,
		FeatureExtraction
	};

	PopulationModifierAdaptation()
	{
		m_type = this->Standard;
		m_adtaudt = 0.1;
		m_adampl = 1;
	}

	PopulationModifierAdaptation(AdaptationType type)
	{
		m_type = type;
		m_adtaudt = 0.1;
		m_adampl = 1;
	}

	PopulationModifierAdaptation(float adampl, float adtaudt)
	{
		m_adtaudt = adtaudt;
		m_adampl = adampl;
	}

	PopulationModifierAdaptation(float adampl, float adtaudt, AdaptationType type)
	{
		m_type = type;
		m_adtaudt = adtaudt;
		m_adampl = adampl;
	}

	void SetParameters(float adampl, float adtaudt)
	{
		m_adampl = adampl;
		m_adtaudt = adtaudt;
	}

	void Simulate();
	void Modify();

	void Reset()
	{
		for(int i=0;i<m_adaptValues.size();i++)
			m_adaptValues[i] = 0;
	}

private:

	AdaptationType m_type;
	vector<float> m_adaptValues;
	float m_adtaudt, m_adampl;
};


class PopulationModifierTrace : public PopulationModifier
{
public:

	PopulationModifierTrace(float traceStrength)
	{
		m_traceStrength = traceStrength;
	}

	void Simulate();
	void Modify();

private:
	float m_traceStrength;
	vector<float> m_lastValues;
};

#endif