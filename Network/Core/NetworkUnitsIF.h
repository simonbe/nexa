#pragma once
#ifndef NETWORKUNITSIF_H
#define NETWORKUNITSIF_H


//#include "Network.h"
#include "NetworkUnits.h"

#if GSL_AVAILABLE==1

#include <gsl/gsl_odeiv.h>

#endif

class Unit;

using namespace std;


class UnitIF : public Unit
{
public:
	UnitIF(float km = 0.0041, float vth = 60000, float vahp = 0, float vk = 0, float ik = 0);

	void SimulateEventQueue();
	UnitModifier* CreateEvent() { return CreateEvent(0.0); }
	UnitModifier* CreateEvent(float time);
	
	vector<vector<float> > GetValuesToRecord()
	{
		return m_recordedValues;
	}

	void SetUnitIdLocalHypercolumn(int id)
	{
		m_localIdHypercolumn = id;
	}

	void SetHypercolumnId(int id)
	{
		m_hypercolumnId = id;
	}

	void AddHypercolumn(Hypercolumn* h)
	{
		m_hypercolumns.push_back(h);
	}

private:

	vector<float> m_storedValues;
	int m_hypercolumnId; // parent hypercolumn id (used?)
	int m_localIdHypercolumn;
	vector<Hypercolumn*> m_hypercolumns;

	// naming from matlab script
	float Km,Vth,Vahp,Vk,Ik;
	
};
#if GSL_AVAILABLE==1

// Model brought in with small modifications from NEST
// Equation solving exactly as in NEST

class UnitadEIF : public Unit
{
public:
private:
};
#endif

#endif