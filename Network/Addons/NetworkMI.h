#pragma once
#ifndef NETWORKMI_H
#define NETWORKMI_H

#include "NetworkConnectionModifier.h"



class ConnectionModifierMIHypercolumn : public ConnectionModifier
{
public:

	ConnectionModifierMIHypercolumn()
	{
		m_firstRun = true;
		sprn = prn = 1;
//		miDij = 0;
		m_eventId = 3;
		
		MAX_FLOAT = std::numeric_limits<float>::max();
		MIN_FLOAT = -MAX_FLOAT;
		EPS_FLOAT = numeric_limits<float>::min();
		m_name = "ConnectionModifierMIHypercolumn";
	}

	/*void AddChildEventConnMC(ConnectionModifierMIRateUnit* e) // corresponding to all outgoing connections from a minicolumn
	{
		m_presMImc.push_back(e);
	}*/

	void SetConnection(Connection* c);
	void Simulate(UnitModifier* e) {}
	void Modify();

	vector<vector<float> > GetDij() { return miDij; }

	vector<vector<float> > GetValuesToRecord();

private:

	bool m_firstRun;
	float MAX_FLOAT, MIN_FLOAT, EPS_FLOAT;

	float taupdt;
	float sprn,prn;

	vector<vector<float> > miDij, hmi, hji;

	//vector<ConnectionModifierMIRateUnit*> m_presMImc;
	//vector<ConnectionModifierMIRateUnit*> m_postsMImc;

	Connection* m_connectionFixedHCs;
};

class ConnectionModifierMIRateUnit : public ConnectionModifier
{
public:

	ConnectionModifierMIRateUnit(ConnectionModifier* eventHc);

	~ConnectionModifierMIRateUnit();

	void SetConnection(Connection* c);

	void Initialize(Connection* connection);
	void Simulate(UnitModifier* e){};
	void Modify();

	vector<float>* GetCi(){ return &miCi; }
	vector<float>* GetCj(){ return &miCj; }
	vector<vector<float> >* GetCij(){ return &miCij; }

	float GetCi(int i)
	{ 
		return miCi[i]; 
	}
	
	float GetCj(int j)
	{ 
		return miCj[j]; 
	}
	
	float GetCij(int i,int j)
	{ 
		return miCij[i][j]; 
	}


	Connection* GetConnection()
	{
		return (Connection*)m_connectionFixedMCs;
	}

	void Reset()
	{
		int i = miCi.size();
		int j = miCj.size();
		miCi = vector<float>(i,0.01f);
		miCj = vector<float>(j,0.01f);

		miCij = vector<vector<float> >(i,vector<float>(j, 0.01f));
	}

private:

	float taupdt;
	float sprn,prn;

	//float miCi,miCj,miCij;
	vector<float> miCi,miCj;
	vector<vector<float> > miCij;

	Connection* m_connectionFixedMCs;
};


#endif