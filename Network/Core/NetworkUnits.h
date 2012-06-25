#pragma once
#ifndef NETWORKUNITS_H
#define NETWORKUNITS_H

#include <math.h>
#include <algorithm>
#include <deque>



//#include "Network.h"
//#include "Logger.h"
#include "NetworkConnections.h"
#include "NetworkUnitModifier.h"
#include "NetworkPopulation.h"
#include "NetworkConnectionModifier.h"


class UnitModifier;
class Connection;
class Population;
class ConnectionModifier;
class Hypercolumn;

using namespace std;

class Unit : public NetworkObject
{
public:

	Unit()
	{
		m_nodeId = 0;
		m_nrEventsSent = 0;
		m_nrEventsReceived = 0;

		m_unitTypeInt = 0;
		m_unitType = "minicolumn"; // default
	}

	virtual ~Unit()
	{

	}

	virtual void ClearMemory()
	{
		m_recordedValues.clear();
	}

	virtual void AddPre(Unit* unit, Connection* c)
	{
		/*m_pres.push_back(unit);
		m_preConnections.push_back(c);
		m_hashIdConnection[unit->GetUnitId()] = c;*/
	}

	virtual void AddPre(Unit* unit)
	{
		m_pres.push_back(unit);
	}

	virtual void AddPost(Unit* unit, Connection* c)
	{
		m_posts.push_back(unit);
		m_postConnections.push_back(c);
		//m_hashIdConnection[unit->GetUnitId()] = c;
	}

	virtual void AddPostConnection(long long connectionId)
	{
	}

	vector<Unit*> GetPres()
	{
		return m_pres;
	}

/*	Connection* GetConnectionFromUnitId(long unitId)
	{
		map<long,Connection*>::iterator itr;

		if ( (itr = m_hashIdConnection.find(unitId)) == m_hashIdConnection.end()) 
			return NULL;
		else
			return m_hashIdConnection[unitId];
	}
	*/

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
		m_nodeId = node;
	}

	int GetNodeId()
	{
		return m_nodeId;
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

	/*vector<UnitModifier*> GetEventsIncoming()
	{
		return m_eventsIncoming;
	}*/

	void ClearEventsIncoming()
	{
		/*for(int j=0;j<m_eventsIncoming.size();j++)
		{
			//if(m_eventsIncoming[j]!=0)
			//{
				delete m_eventsIncoming[j];
			//	m_eventsIncoming[j] = 0;
			//}
		}

		m_eventsIncoming.clear();
		m_nrEventsReceived = 0;*/
		
//		if(m_eventsIncoming.size()>0)
//			m_eventsIncoming.erase(m_eventsIncoming.begin(),m_eventsIncoming.begin()+1);
		
		//m_eventsIncoming.clear();
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

	//void AddEventIncoming(UnitModifier* e)
	//{
		/*if(m_eventsIncoming.size()>m_nrEventsReceived)
		{
			m_eventsIncoming[m_nrEventsReceived] = e;
		}
		else*/
	//		m_eventsIncoming.push_back(e);


	//	m_nrEventsReceived++;
	//}

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

	virtual void ReceiveEvents();
	virtual void SendReceive();
	virtual void SendReceiveNext();

	bool AlreadySentEventToNode(int nodeId);

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

//	UnitModifier* GetUnitModifier(int id);
//	UnitModifier* GetUnitModifier(string name);

protected:

	vector<UnitModifier*> m_eventsOutgoing;

	float m_value;
	bool m_isNewEvent;
	vector<vector<float> > m_recordedValues; // be in network object instead?
	string m_unitType;
	int m_unitTypeInt;

	bool m_isLocal;
	vector<Unit*> m_pres;
	vector<Unit*> m_posts;
	vector<Connection*> m_preConnections;
	vector<Connection*> m_postConnections;
	vector<int> m_postProcesses; // mpi post processes

	long m_unitId;
	int m_nodeId;
	long m_unitIdLocal;

	//vector<UnitModifier*> m_eventsIncoming;
	//vector<vector<long> > m_eventsIncoming;
//	vector<vector<long> > m_eventsIncoming; //

	vector<ConnectionModifier*> m_eventsConnections;

	UnitModifier* m_eventOutgoing;

	vector<int> m_sentToNode;
	int m_maxNode;

	Population* m_population;
	// Coordinates m_coordinates;

	vector<UnitModifier*> m_unitProperties; // global level - for all incoming connections
	int m_nrEventsSent, m_nrEventsReceived;
};

// RateUnit + (optionally) adaptation and trace
class RateUnit : public Unit
{
public:

	RateUnit(bool useThreshold)
	{
		m_firstRun = true;
		m_isTransferCSL = false; // special check for CSL (optimization)
		m_value = 0;
		m_beta = 0;
		m_inhibBeta = 0;
		m_unitType = "minicolumn";
		m_nrEventsSent = 0;
		m_nrEventsReceived = 0;
		m_useTrace = false; // should be global and not state variable
		m_useAdaptation = false; // "
		m_isNewEvent = true;
		m_useThreshold = useThreshold; // should be global
		if(m_useThreshold)
			m_subThreshValue = 0; // replace to also do allocation here

		m_noUpdatingCurrentTimeStep = false;
		m_name="RateUnit";
	}

	~RateUnit()
	{
		// deallocate subthreshold val etc
	}

	void SimulateEventQueue();

	// override
	float GetValue()
	{
		/*if(m_useTrace)
			return m_trace[m_traceTimeSteps];
		else*/
			return m_value;
	}

/*	virtual float GetValue()
	{
		return m_value;
	}*/

	float GetValueFromGroup(vector<int> nodeIndexes);

	void SetValue(float value)
	{
		m_value = value;
	}

	void SetSubThresholdValue(float value) // valid if unit GradedThresholded is used
	{
		m_subThreshValue = value;
	}

	float GetSubThresholdValue() // valid if unit GradedThresholded is used
	{
		return m_subThreshValue;
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

	void SetBeta(float value)
	{
		m_beta = value;
	}

	void AddBeta(float value)
	{
		m_beta += value;
	}

	void AddInhibBeta(float value)
	{
		m_inhibBeta += value;
	}

	float GetBeta()
	{
		return m_beta;
	}

	vector<vector<float> > GetValuesToRecord();

	void AddGains();

	void UseAdaptation(bool useAdaptation) { m_useAdaptation = useAdaptation; }
	
	//void UseTrace(bool useTrace) { m_useTrace = useTrace; }
	//void SetTraceTimeSteps(int traceTimeSteps);

	void SetTrace(bool useTrace, float traceStrength)
	{
		m_useTrace = useTrace;
		m_traceStrength = traceStrength;
	}

	void SetAdaptation(float ampl, float tau, bool on) { m_adaptationAmpl = ampl; m_adaptationTau = tau; m_useAdaptation = on; }

	// override
	void SimulateMisc();

	void Dispose()
	{
	//	m_trace.clear();
	}

	void SetNoUpdatingCurrentTimeStep(bool noUpdating)
	{
		m_noUpdatingCurrentTimeStep = noUpdating;
	}

private:

	// Functions

	void SimulateUnitProperties(vector<UnitModifier*> unitProperties, vector<UnitModifier*> eus, vector<float> weights);
	void SimulateUnitPropertiesV2(vector<UnitModifier*>* unitProperties, vector<float>* values, vector<float>* weights, vector<long>* hypercolumnIds); // v2, optimized
	// v3, with delays

	// Variables
	bool m_firstRun;
	bool m_isTransferCSL;

	bool m_useAdaptation;
	bool m_useTrace;
	float m_adaptationValue;
	//map<int,float> m_trace; // map<timestep constant, trace value>
	float m_traceStrength;
	float m_traceValueLast;
	
	float m_adaptationAmpl;
	float m_adaptationTau;
	bool m_useThreshold;

	bool m_noUpdatingCurrentTimeStep;

	float m_subThreshValue; // used for gradedThresholded unit - can avoid allocation

	
	float m_beta;		// used in bcpnn, (change location)
	float m_inhibBeta;	// " (change location)

	int m_localIdHypercolumn;
	int m_hypercolumnId; // parent hypercolumn id

	vector<Hypercolumn*> m_hypercolumns;
	vector<vector<int> > m_incomingHypercolumnIndexes;

	// adaptation + trace

	//int m_traceTimeSteps;

};

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

	// overloaded
	/*long GetUnitIdLocal()
	{
		return m_minicolumns[0]->GetUnitIdLocal();
	}*/

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
		m_values = values;
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

	int m_totalNrRateUnits; // global
	vector<RateUnit*> m_minicolumns; // local
	vector<float> m_values;

	bool m_useSilent; // if we should be able to use silent
	bool m_isSilent; // current state
	float m_silentThreshold;
};


/*class UnitIF : public Unit
{
public:
	UnitIF();

	// functions for columnar structure - should merge this and minicolumn
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

	bool IsNewEvent() { return true; }

	vector<vector<float> > GetValuesToRecord();

private:

	int m_localIdHypercolumn;
	int m_hypercolumnId; // parent hypercolumn id
	vector<Hypercolumn*> m_hypercolumns;

	void update(const long from, const long to);

	void SimulateEventQueue();
	
	UnitModifier* CreateEvent(float value) { return CreateEvent(); };
	UnitModifier* CreateEvent();

	 struct Parameters_ {
      // Membrane capacitance in pF.
      double C_;
    
      // Membrane time constant in ms. 
      double Tau_; 

      // Time constant of synaptic current in ms.
      double tau_syn_;
      
      // Refractory period in ms. 
      double TauR_;

      // Resting potential in mV. 
      double U0_;

      // Reset value of the membrane potential, in mV.
      //    @note Value is relative to resting potential U0_.
      double V_reset_;

      // Threshold in mV. 
      //    @note Value is relative to resting potential U0_.
      double Theta_;

      // External current in pA
      double I_e_;

      Parameters_();  //!< Sets default parameter values

      //void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      //void set(const DictionaryDatum&);  //!< Set values from dicitonary
    };

    // ----------------------------------------------------------------

    // State variables

    struct State_ {
      double y0_; //!< Constant current
      double y1_;  
      double y2_;
      double y3_; //!< This is the membrane potential RELATIVE TO RESTING POTENTIAL.

      int    r_;  //!< number of refractory steps remaining

      State_();  //!< Default initialization
      
      //void get(DictionaryDatum&, const Parameters_&) const;
      //void set(const DictionaryDatum&, const Parameters_&);
    };
    
    // ---------------------------------------------------------------- 

    // Buffers

    struct Buffers_ {
      // buffers and summs up incoming spikes/currents
      //RingBuffer spikes_;
      //RingBuffer currents_;


      // Buffer for membrane potential.

      //AnalogDataLogger<PotentialRequest> potentials_;
    };
    
    // ---------------------------------------------------------------- 

    
    // Internal variables of the model.

    struct Variables_ { 
      // Amplitude of the synaptic current.
	  //      This value is chosen such that a post-synaptic potential with
	  //      weight one has an amplitude of 1 mV.
      
     double PSCInitialValue_;
     int    RefractoryCounts_;  //!< refractory time in steps
    
     double P11_;   
     double P21_;
     double P22_;
     double P31_;
     double P32_;
     double P30_;
     double P33_;
   };

   // ---------------------------------------------------------------- 

   //
   //  @defgroup iaf_neuron_data
   //  Instances of private data structures for the different types
   //  of data pertaining to the model.
   //  @note The order of definitions is crucial: Moving Variables_
   //        to the very end increases simulation time for brunel-2.sli
   //        from 72s to 81s on a Mac, Intel Core 2 Duo 2.2GHz, g++ 4.0.1 -O3
   //  @{
       
   Parameters_ P_;
   State_      S_;
   Variables_  V_;
   Buffers_    B_;
   // @}
};

*/



#endif

/*void Pop::updsup() {
	if (noiseampl<=0)
		for (int j=0; j<nunit; j++) dsup[j] += (gain*sup[j] - dsup[j])* taumdt;
	else {
		if (simstep==0) cerr << "#" ;
		for (int j=0; j<nunit; j++)
			dsup[j] += (gain*(sup[j] +
			noiseampl*(poisson(noiseintens) - noiseintens)) -
			dsup[j])* taumdt;
	}
}*/

