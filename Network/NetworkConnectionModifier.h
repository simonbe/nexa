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
#include "NetworkConnections.h"
#include "NetworkObject.h"

class PopulationModifier;
class Connection;
class RateUnit;
class Hypercolumn;
class UnitModifier;

using namespace std;

class ConnectionModifier : public NetworkObject
{
public:

	ConnectionModifier()
	{
		m_value = 0;
		m_eventId = -1;
		m_transferFunction = NULL;
	}

	virtual ~ConnectionModifier();

	virtual void Initialize(Connection* connection);
	virtual void Simulate(UnitModifier* e) = 0;//{ }
	virtual void Modify() = 0;

	double GetValue()
	{
		return m_value;
	}

	virtual void SetConnection(Connection* c);

	virtual Connection* GetConnection()
	{
		return m_connection;
	}

	int GetEventId()
	{
		return m_eventId;
	}

	void AddParentConnectionModifier(ConnectionModifier* e)
	{
		m_parentConnectionModifier.push_back(e);
	}

	void AddChildConnectionModifier(ConnectionModifier* e)
	{
		m_childConnectionModifier.push_back(e);
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
	
	Connection* m_connection;
	UnitModifier* m_transferFunction;
	
	vector<ConnectionModifier*> m_parentConnectionModifier;
	vector<ConnectionModifier*> m_childConnectionModifier;
	vector<PopulationModifier*> m_parentPopulationModifier; // allows event layer to have direct access to this event connection
};



#endif