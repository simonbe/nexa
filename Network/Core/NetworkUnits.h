#pragma once
#ifndef NETWORKUNITS_H
#define NETWORKUNITS_H

#include <math.h>
#include <algorithm>
#include <deque>
#include "NetworkProjections.h"
#include "NetworkUnitModifier.h"
#include "NetworkPopulation.h"
#include "NetworkProjectionModifier.h"


class UnitModifier;
class Projection;
class Population;
class ProjectionModifier;
class Hypercolumn;

using namespace std;

/// <summary>	A neural unit. </summary>

class Unit : public NetworkObject
{
public:

	Unit()
	{
		m_processId = 0;
		m_nrEventsSent = 0;
		m_nrEventsReceived = 0;

		m_unitTypeInt = 0;
		m_unitType = "minicolumn"; // default type
	}

	virtual ~Unit()
	{

	}

	virtual void ClearMemory()
	{
		m_recordedValues.clear();
	}

	bool IsLocal()
	{
		return m_isLocal;
	}

	void SetLocal(bool isLocal)
	{
		m_isLocal = isLocal;
	}

	string GetType()
	{
		return m_unitType;
	}

	void SetNodeId(int node)
	{
		m_processId = node;
	}

	int GetNodeId()
	{
		return m_processId;
	}

	virtual long GetUnitIdLocal()
	{
		return m_unitIdLocal;
	}

	long GetUnitId()
	{
		return m_unitId;
	}

	void SetUnitId(long id)
	{
		m_unitId = id;
	}

	string GetUnitType()
	{
		return m_unitType;
	}

	void SetUnitIdLocal(long localId)
	{
		m_unitIdLocal = localId;
	}

	void ClearEventsIncoming()
	{
	}

	vector<UnitModifier*> GetEvents()
	{
		return m_eventsOutgoing;
	}

	void SetEventOutgoing(UnitModifier* e)
	{
		m_eventOutgoing = e;
	}

	void AddEventOutgoing();

	void AddEventIncoming(long eventIndex);

	virtual void SimulateEventQueue() = 0;
	virtual void Simulate();

	virtual void SimulateMisc()
	{
	}

	virtual UnitModifier* CreateEvent() = 0;
	virtual bool IsNewEvent()
	{
		return m_isNewEvent;
	}

	void IsNewEvent(bool isNewEvent)
	{
		m_isNewEvent = isNewEvent;
		
		for(unsigned int i=0;i<m_eventsOutgoing.size();i++)
		{
			delete m_eventsOutgoing[i];
		}

		m_eventsOutgoing.clear();
	}

	virtual UnitModifier* CreateEvent(float value) = 0;

	void SetPopulation(Population* net)
	{
		m_population = net;
	}

	Population* GetPopulation()
	{
		return m_population;
	}

	void AddUnitModifier(UnitModifier* p);

	void AddPostProcess(int process)
	{
		m_postProcesses.push_back(process);
		/*vector<int>::iterator result;
		result = find( m_postProcesses.begin(), m_postProcesses.end(), node );

		if( result == m_postProcesses.end() )
			m_postProcesses.push_back(node);*/
	}

	vector<int>* GetPostProcesses()
	{
		return &m_postProcesses;
	}

	bool IsPostNode(int node)
	{
		vector<int>::iterator result;
		result = find( m_postProcesses.begin(), m_postProcesses.end(), node );

		if( result == m_postProcesses.end() )
			return false;
		else return true;
	}

	virtual float GetValue()
	{
		return m_value;
	}

	virtual void SetValue(float value)
	{
		m_value = value;
	}

protected:

	vector<UnitModifier*> m_eventsOutgoing;

	float m_value;
	float m_subThreshValue, m_subThreshDrive; // TODO: used for gradedThresholded unit - move it

	bool m_isNewEvent;
	vector<vector<float> > m_recordedValues; // TODO: may get moved to network object class
	string m_unitType;
	int m_unitTypeInt;

	bool m_isLocal;
	vector<Unit*> m_pres;
	vector<Unit*> m_posts;
	vector<Projection*> m_preProjections;
	vector<Projection*> m_postProjections;
	vector<int> m_postProcesses; // mpi post processes

	long m_unitId;
	int m_processId;
	long m_unitIdLocal;

	vector<ProjectionModifier*> m_eventsProjections;

	UnitModifier* m_eventOutgoing;

	vector<int> m_sentToNode;
	int m_maxNode;

	Population* m_population;

	vector<UnitModifier*> m_unitProperties; // global level - for all incoming Projections
	int m_nrEventsSent, m_nrEventsReceived;

	/// <summary>	Stores incoming data if a delay is used. </summary>
	vector<vector<long> > m_eventsIncoming;
};

/// <summary>	Rate unit. 
/// 			+ (optionally) adaptation and trace</summary>

class RateUnit : public Unit
{
public:

	RateUnit(bool useThreshold = false)
	{
		m_firstRun = true;
		m_isTransferCSL = false; // TODO: move to CSL (special check for CSL for optimization)
		m_value = 0;
		m_beta = 0;
		m_inhibBeta = 0;
		m_unitType = "minicolumn";
		m_nrEventsSent = 0;
		m_nrEventsReceived = 0;
		m_useTrace = false; // TODO: make global for population
		m_useAdaptation = false; // "
		m_isNewEvent = true;
		m_useThreshold = useThreshold; // TODO: make globalcould be global
		m_subThreshValue = 0; // TODO: replace to also do allocation here

		m_noUpdatingCurrentTimeStep = false;
		m_name="RateUnit";
	}

	~RateUnit()
	{
	}

	void SimulateEventQueue();

	float GetValue()
	{
		return m_value;
	}

	float GetValueFromGroup(vector<int> nodeIndexes);

	void SetValue(float value)
	{
		m_value = value;
	}

	void SetSubThresholdDrive(float value)
	{
		m_subThreshDrive = value;
	}

	void SetSubThresholdValue(float value) // valid if unit GradedThresholded is used
	{
		m_subThreshValue = value;
	}

	float GetSubThresholdValue() // valid if unit GradedThresholded is used
	{
		return m_subThreshValue;
	}
	
	float GetSubThresholdDrive() // valid if unit GradedThresholded is used
	{
		return m_subThreshDrive;
	}

	UnitModifier* CreateEvent(float value);
	UnitModifier* CreateEvent();
	
	void AddHypercolumn(Hypercolumn* h)
	{
		m_hypercolumns.push_back(h);
	}

	vector<Hypercolumn*> GetHypercolumns()
	{
		return m_hypercolumns;
	}

	int GetUnitIdLocalInHypercolumn()
	{
		return m_localIdHypercolumn;
	}

	void SetUnitIdLocalHypercolumn(int id)
	{
		m_localIdHypercolumn = id;
	}

	void SetHypercolumnId(int id)
	{
		m_hypercolumnId = id;
	}

	int GetHypercolumnId()
	{
		return m_hypercolumnId;
	}

	// TODO: move (specific to BCPNN)
	void SetBeta(float value)
	{
		m_beta = value;
	}

	// TODO: move (specific to BCPNN)
	void AddBeta(float value)
	{
		m_beta += value;
	}

	// TODO: move (specific to BCPNN)
	void AddInhibBeta(float value)
	{
		m_inhibBeta += value;
	}

	// TODO: move (specific to BCPNN)
	float GetBeta()
	{
		return m_beta;
	}

	vector<vector<float> > GetValuesToRecord();

	void AddGains();

	void UseAdaptation(bool useAdaptation) { m_useAdaptation = useAdaptation; }

	void SetTrace(bool useTrace, float traceStrength)
	{
		m_useTrace = useTrace;
		m_traceStrength = traceStrength;
	}

	void SetAdaptation(float ampl, float tau, bool on) { m_adaptationAmpl = ampl; m_adaptationTau = tau; m_useAdaptation = on; }

	void SimulateMisc();

	void Dispose() { }

	/// <summary>	Forces no update of unit in the coming time step. </summary>
	///
	/// <param name="noUpdating">	if true no updating. </param>

	void SetNoUpdatingCurrentTimeStep(bool noUpdating)
	{
		m_noUpdatingCurrentTimeStep = noUpdating;
	}

private:

	// Functions
	void SimulateUnitProperties(vector<UnitModifier*> unitProperties, vector<UnitModifier*> eus, vector<float> weights);
	void SimulateUnitPropertiesV2(vector<UnitModifier*>* unitProperties, vector<float>* values, vector<float>* weights, vector<long>* hypercolumnIds); // v2, optimized

	bool m_firstRun;
	bool m_isTransferCSL;

	bool m_useAdaptation;
	bool m_useTrace;
	float m_adaptationValue;
	float m_traceStrength;
	float m_traceValueLast;
	
	float m_adaptationAmpl;
	float m_adaptationTau;
	bool m_useThreshold;

	bool m_noUpdatingCurrentTimeStep;
	
	float m_beta;		// TODO: move (used in bcpnn)
	float m_inhibBeta;	// TODO: move

	int m_localIdHypercolumn; // local parent hypercolumn id
	int m_hypercolumnId; // parent hypercolumn id

	vector<Hypercolumn*> m_hypercolumns;
	vector<vector<int> > m_incomingHypercolumnIndexes;
};

/// <summary>	Hypercolumn. A column containing a number of (typically) RateUnits (minicolumns). 
/// 			If a silent threshold set in constructor, activity of all minicolumns forced to zero if none over threshold value. </summary>

class Hypercolumn : public Unit
{
public:

	Hypercolumn()
	{
		m_unitType = "hypercolumn";
		m_isSilent = false;
		m_useSilent = false;
	}

	Hypercolumn(float silentThreshold)
	{
		m_unitType = "hypercolumn";
		m_useSilent = true;
		m_isSilent = false;
		m_silentThreshold = silentThreshold;
	}

	~Hypercolumn()
	{
	}

	void SimulateEventQueue();

	UnitModifier* CreateEvent(float value);
	UnitModifier* CreateEvent();

	void AddRateUnit(Unit* m)
	{
		m_minicolumns.push_back((RateUnit*)m);
	}

	vector<RateUnit*> GetRateUnits()
	{
		return m_minicolumns;
	}
	
	void SetUseSilent(bool useSilent)
	{
		m_useSilent = useSilent;
	}

	void SetSilentThreshold(float threshold)
	{
		m_silentThreshold = threshold;
	}

	bool IsSilent()
	{
		return m_isSilent;
	}
	vector<vector<float> > GetValuesToRecord();

	// works as a buffer to include also non-local values from minicolumns that are not residing on this process
	vector<float> GetValues()
	{
		if(m_values.size()>0)
			return m_values;
		else // all local values
		{
			vector<float> data(this->GetRateUnits().size());

			if(this->m_useSilent == false) // standard
			{
				for(int i=0;i<this->GetRateUnits().size();i++)
				{
					data[i] = this->GetRateUnits()[i]->GetValue();
				}
				return data;
			}
			else
			{
				float maxVal = -1e8;
				for(int i=0;i<this->GetRateUnits().size();i++)
				{
					data[i] = this->GetRateUnits()[i]->GetValue();
					if(data[i]>maxVal)
						maxVal = data[i];
				}

				if(maxVal > this->m_silentThreshold)
					return data;
				else
				{
					return vector<float>(0);
				}
			}
			
		}
	}

	void SetValues(vector<float> values)
	{
		m_values = values; // not necessary

		vector<RateUnit*> units = this->GetRateUnits();
		for(int j=0;j<units.size();j++)
		{
			units[j]->SetValue(values[units[j]->GetUnitIdLocalInHypercolumn()]);//data[j]);
		}
	}

	int GetTotalNrRateUnits()
	{
		return m_totalNrRateUnits;
	}

	void SetTotalNrRateUnits(int tot)
	{
		m_totalNrRateUnits = tot;
	}

private:

	int m_totalNrRateUnits; // global (nr including non-local minicolumns)
	vector<RateUnit*> m_minicolumns; // local minicolumns
	vector<float> m_values;

	bool m_useSilent; // if we should be able to use silent
	bool m_isSilent; // current state (TODO: check if needed)
	float m_silentThreshold;
};

#endif