#pragma once
#ifndef NETWORKBCPNN_H
#define NETWORKBCPNN_H

#include "NetworkConnectionModifier.h"

class ConnectionModifierBcpnnOnline : public ConnectionModifier
{
public:

	ConnectionModifierBcpnnOnline(float alpha, float lambda, float maxValue = -1);
	ConnectionModifierBcpnnOnline();

	void Initialize(Connection* connection);
	void SetConnection(Connection* c);
	
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
	//BCPNNOnline m_bcpnnOnline;
	/*vector<float> m_Ai;
	vector<float> m_Aj;
	vector<vector<float> > m_Aij;
	vector<float> m_beta; // in unitevent
	vector<float> m_inhibBeta; // in case of inhibitory
	*/

	// hash implementation (allows for moving any synapses)
	map<long,float> m_Ai;
	vector<float> m_Aj;
	vector<float> m_beta; // in unitevent
	vector<float> m_inhibBeta; // in case of inhibitory
	vector<vector<float> > m_Aij; // watch out, need to restructure if synapses moved
	//vector<map<long,float> > m_Aij;

	float m_alpha;
	float m_lambda0;

	// multiplies in beta and weights
	float m_impactWeights;
	float m_impactBeta;

	vector<long> m_postIds;
	//vector<vector<long> > m_idsPre;

	Connection* m_connectionFixed;
};

#endif