#pragma once
#ifndef NETWORKEVLAY_H
#define NETWORKEVLAY_H

#include "Network.h"
#include "NetworkPopulation.h"
#include "NetworkConnectionModifier.h"

class Population;
class ConnectionModifier;

using namespace std;

class PopulationModifier : public NetworkObject
{
public:
	virtual void Simulate(){}// = 0;
	virtual void Modify(){}

	virtual void Initialize(Population* population)
	{
		if(m_initialized == false) // move to parent
		{
			m_population = population;
		}

		m_initialized = true;
	}

	virtual ~PopulationModifier()
	{
		/*for(int i=0;i<m_childConnectionModifier.size();i++)
			delete m_childConnectionModifier[i];

		for(int i=0;i<m_parentPopulationModifier.size();i++)
			delete m_parentPopulationModifier[i];

		for(int i=0;i<m_childPopulationModifier.size();i++)
			delete m_childPopulationModifier[i];*/
	}

	void AddChildConnectionModifier(ConnectionModifier* e)
	{
		m_childConnectionModifier.push_back(e);
	}

	void AddParentPopulationModifier(PopulationModifier* e)
	{
		m_parentPopulationModifier.push_back(e);
	}

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

	/*virtual void Reset()
	{
	}*/
	
protected:

	vector<ConnectionModifier*> m_childConnectionModifier;
	vector<PopulationModifier*> m_parentPopulationModifier;
	vector<PopulationModifier*> m_childPopulationModifier;
	Population* m_population;
};

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
	std::vector<double> Function(std::vector<double> data, float G);
	std::vector<double> WTAProb(std::vector<double> data);
	std::vector<double> tempWTAFunction(std::vector<double> data);
	std::vector<double> ksoftwinners(std::vector<double> data,int k=10);

private:

	float m_threshold;
	SoftMaxType m_type;
	float m_G;
	bool m_probabilisticWTA;
};

class WTA : public PopulationModifier
{
public:

	WTA()
	{
		m_name = "WTA";
	}

	void Simulate();
	std::vector<float> Function(std::vector<float> data);

private:

};


class WTAThreshold : public PopulationModifier
{
public:

	WTAThreshold()
	{
		m_name = "WTAThreshold";
	}

	void Simulate();
	std::vector<float> Function(std::vector<float> data);

private:

	WTA m_wta;
};


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