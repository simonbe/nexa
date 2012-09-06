#pragma once
#ifndef NETWORKEVCON_H
#define NETWORKEVCON_H

#include <string>
#include <limits>
#include "Network.h"
#include "NetworkUnits.h"
#include "NetworkPopulationModifier.h"
//#include "NetworkUnitModifier.h"
#include "NetworkPopulation.h"
#include "NetworkProjections.h"
#include "NetworkObject.h"

class PopulationModifier;
class Projection;
class RateUnit;
class Hypercolumn;
class UnitModifier;

using namespace std;

/// <summary>	Projection modifier: Changes state variables of a projection.
/// 			Examples include synaptic plasticity and synaptic depression.  </summary>

class ProjectionModifier : public NetworkObject
{
public:

	ProjectionModifier()
	{
		m_value = 0;
		m_eventId = -1;
		m_transferFunction = NULL;
	}

	virtual ~ProjectionModifier();

	virtual void Initialize(Projection* Projection);
	virtual void Simulate(UnitModifier* e) = 0;//{ }
	virtual void Modify() = 0;

	double GetValue()
	{
		return m_value;
	}

	virtual void SetProjection(Projection* c);

	virtual Projection* GetProjection()
	{
		return m_projection;
	}

	int GetEventId()
	{
		return m_eventId;
	}

	void AddParentProjectionModifier(ProjectionModifier* e)
	{
		m_parentProjectionModifier.push_back(e);
	}

	void AddChildProjectionModifier(ProjectionModifier* e)
	{
		m_childProjectionModifier.push_back(e);
	}

	virtual void AddParentPopulationModifier(PopulationModifier* e)
	{
		m_parentPopulationModifier.push_back(e);
	}

	virtual void SaveState()
	{
	}

	virtual void LoadState()
	{
	}

	// specifically implemented for sensor drift application
	virtual void DriftAdaptStep(float eta)
	{
	}

	UnitModifier *GetTransferFunction();

	void SetTransferFunction(UnitModifier* transfer)
	{
		m_transferFunction = transfer;
	}

	/*virtual void Reset()
	{
	}*/

protected:

	double m_value;
	int m_eventId;

	//Unit* m_pre;
	//Unit* m_post;
	
	Projection* m_projection;
	UnitModifier* m_transferFunction;
	
	vector<ProjectionModifier*> m_parentProjectionModifier;
	vector<ProjectionModifier*> m_childProjectionModifier;
	vector<PopulationModifier*> m_parentPopulationModifier; // allows event layer to have direct access to this event Projection
};



#endif