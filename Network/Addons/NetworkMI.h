#pragma once
#ifndef NETWORKMI_H
#define NETWORKMI_H

#include "NetworkProjectionModifier.h"



class ProjectionModifierMIHypercolumn : public ProjectionModifier
{
public:

	ProjectionModifierMIHypercolumn()
	{
		m_firstRun = true;
		sprn = prn = 1;
//		miDij = 0;
		m_eventId = 3;
		
		MAX_FLOAT = std::numeric_limits<float>::max();
		MIN_FLOAT = -MAX_FLOAT;
		EPS_FLOAT = numeric_limits<float>::min();
		m_name = "ProjectionModifierMIHypercolumn";
	}

	/*void AddChildEventConnMC(ProjectionModifierMIRateUnit* e) // corresponding to all outgoing Projections from a minicolumn
	{
		m_presMImc.push_back(e);
	}*/

	void SetProjection(Projection* c);
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

	//vector<ProjectionModifierMIRateUnit*> m_presMImc;
	//vector<ProjectionModifierMIRateUnit*> m_postsMImc;

	Projection* m_projectionFixedHCs;
};

class ProjectionModifierMIRateUnit : public ProjectionModifier
{
public:

	ProjectionModifierMIRateUnit(ProjectionModifier* eventHc);

	~ProjectionModifierMIRateUnit();

	void SetProjection(Projection* c);

	void Initialize(Projection* Projection);
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


	Projection* GetProjection()
	{
		return (Projection*)m_projectionFixedMCs;
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

	Projection* m_projectionFixedMCs;
};


#endif