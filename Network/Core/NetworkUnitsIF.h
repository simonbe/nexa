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
	UnitIF();

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


class UnitadEIF : public Unit
{
 // not included atm in repository.
};




#endif

#endif