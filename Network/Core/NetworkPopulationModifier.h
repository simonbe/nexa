#pragma once
#ifndef NETWORKEVLAY_H
#define NETWORKEVLAY_H

#include "Network.h"
#include "NetworkPopulation.h"
#include "NetworkProjectionModifier.h"

class Population;
class ProjectionModifier;

using namespace std;

/// <summary>	Modifies the activity values or any other state variables of a population of neural units.
///				Could e.g. implement abstract functions across population or columns of a population. </summary>

class PopulationModifier : public NetworkObject
{
public:

	/// <summary>	Main function to override. </summary>

	virtual void Simulate(){}// = 0;

	/// <summary>	Run after global communication step, not overridden in any example class. </summary>
	
	virtual void Modify(){}

	virtual void Initialize(Population* population)
	{
		if(m_initialized == false) // may be moved to parent
		{
			m_population = population;
		}

		m_initialized = true;
	}

	virtual ~PopulationModifier()
	{
	}

	/// <summary>	Help function for Network. </summary>

	void AddChildProjectionModifier(ProjectionModifier* e)
	{
		m_childProjectionModifier.push_back(e);
	}

	/// <summary>	Help function for Network. </summary>

	void AddParentPopulationModifier(PopulationModifier* e)
	{
		m_parentPopulationModifier.push_back(e);
	}

	/// <summary>	Help function for Network. </summary>

	void AddChildPopulationModifier(PopulationModifier* e)
	{
		e->AddParentPopulationModifier(this);
		m_childPopulationModifier.push_back(e);
	}

	virtual vector<vector<float> >* GetOutput()
	{
		return NULL;
	}

	Population* GetPopulation()
	{
		return m_population;
	}
	
protected:

	vector<ProjectionModifier*> m_childProjectionModifier;
	vector<PopulationModifier*> m_parentPopulationModifier;
	vector<PopulationModifier*> m_childPopulationModifier;
	Population* m_population;
};

/// <summary>	 SoftMax function as a population modifier, sharpness decided by m_G.
///			Observe: The exp used can lead to NaN if too high numbers are fed in.
///			Now these numbers are reduced by default but no warning is issued. </summary>

class SoftMax : public PopulationModifier
{
public:

	enum SoftMaxType
	{
		Standard,
		ProbWTA,
		WTA,
		WTAThresholded,
		KSOFT // new implementation, free after k-WTA
	};


	SoftMax()
	{
		m_G = 1.0;
		m_type = Standard;
		m_name = "SoftMax";
	}

	SoftMax(float G, SoftMaxType type)
	{
		m_G = G;
		m_type = type;
		m_name = "SoftMax";
	}

	void SetType(SoftMaxType type)
	{
		m_type = type;
	}

	SoftMax(SoftMaxType type, float threshold) // G not used
	{
		m_threshold = threshold;
		//m_G = G;
		m_type = type;
	}


	void Simulate();
	vector<double> Function(vector<double> data, float G);
	vector<double> WTAProb(vector<double> data);
	vector<double> tempWTAFunction(vector<double> data);
	vector<double> ksoftwinners(vector<double> data,int k=10);

private:

	float m_threshold;
	SoftMaxType m_type;
	float m_G;
	bool m_probabilisticWTA;
};

/// <summary>	 Implements a winner-take-all function as a population modifier.
/// Works on the columns if a population has a columnar organization. </summary>

class WTA : public PopulationModifier
{
public:

	WTA()
	{
		m_name = "WTA";
	}

	void Simulate();
	vector<float> wta(vector<float> data);

private:

};

/// <summary>	 Test function - winner-take-all with a threshold. Use normal WTA in standard cases. </summary>

class WTAThreshold : public PopulationModifier
{
public:

	WTAThreshold()
	{
		m_name = "WTAThreshold";
	}

	void Simulate();
	vector<float> Function(vector<float> data);

private:

	WTA m_wta;
};

/// <summary>	Divisive normalization. </summary>

class DivisiveNormalization : public PopulationModifier
{
public:

	DivisiveNormalization()
	{
		m_name = "DivisiveNormalization";
	}

	void Simulate();
	//vector<float> Function(vector<float> data);

private:

};

/// <summary>	Threshold as a population modifier, sets activity of unit to 0 unless over threshold. </summary>

class Threshold : public PopulationModifier
{
public:

	Threshold(float threshold)
	{
		m_name = "Threshold";
		m_threshold = threshold;
		m_numActiveUnits = 0;
	}

	int GetNumActiveUnits(){
		return m_numActiveUnits;
	}

	void Simulate();

private:

	float m_threshold;
	int m_numActiveUnits;
};

#endif