#pragma once
#include "Network.h"
#include "NetworkProjections.h"

class FullConnectivity;

/// <summary>	Network parameters. 
/// 	
///	Used so that it is possible to specify an input vector to a network construction
/// and let the network run independent runs for all parameters.
/// Ex:
/// paramMaxWeightValues = {0.2,0.3,0.4};
/// fullConnectivity->SetRandomWeights(0, paramMaxWeightValues);
/// 
/// extend:  with e.g. boost.bind, boost.function, replace type with template,
/// 		optimize parameters wrt fitness function
/// 
/// </summary>

class NetworkParameters
{

public:

	NetworkParameters()
	{
		m_currentStep = 0;
		m_totValues = 0;
	}

	void Reset();

	void AddParameters(FullConnectivity* object, void (FullConnectivity::*function)(float), vector<float> values);
	bool SetNextParameters();
	bool ParametersLeft();

private:
	
	typedef void (FullConnectivity::*fptr)(float);

	vector<FullConnectivity*> m_objects;
	vector<fptr> m_functions;
	vector<vector<float> > m_values;
	int m_currentStep;
	int m_totValues;

	// Strategies: All combinations
	vector<float> AllCombinationsGetParameters();

	// Strategy: all combinations
	vector<int> m_allCombIndexes;
};