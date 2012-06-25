#include "NetworkParameters.h"

void NetworkParameters::AddParameters(FullConnectivity* object, void (FullConnectivity::*function)(float), vector<float> values)
{
	m_objects.push_back(object);
	m_functions.push_back(function);
	m_values.push_back(values);
	m_totValues+=values.size();
}

void NetworkParameters::Reset()
{
	m_objects.clear();
	m_functions.clear();
	m_values.clear();
	m_totValues = 0;

}

bool NetworkParameters::ParametersLeft()
{
	if(m_currentStep >= m_totValues)
		return false;
	else
		return true;
}

bool NetworkParameters::SetNextParameters()
{
	if(m_functions.size()==0)
		return false;

	if(m_currentStep >= m_totValues)
		return false;
	else
	{
		//auto func = [] () { cout << "Hello world"; };
		//func(); // now call the function

		// only all combinations possible atm
		vector<float> currentParameters = AllCombinationsGetParameters();

		if(currentParameters.size()==0)
			return false; // completed
		else
		{
			for(int i=0;i<m_functions.size();i++)
			{
				(m_objects[i]->*m_functions[i])(currentParameters[i]);
			}
		}

		m_currentStep++;
	}
}

// change to more than one variable
vector<float> NetworkParameters::AllCombinationsGetParameters()
{
	vector<float> indexes;

	if(m_currentStep<m_values[0].size())
		indexes.push_back(m_values[0][m_currentStep]);
	
	return indexes;
	/*
	vector<int> indexes;
	for(int i=0;i<m_parameterValues;i++)
	{

		for(int j=0;j<m_parameterValues[i];j++)
		{

		}
	}*/

}