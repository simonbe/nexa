#pragma once

#include "Network.h"
#include "Analysis.h"

class Network;
class AnalysisTiming;

using namespace std;

/// <summary>	Generic class all network objects inherit from.
/// 			Contains flags for enabling/disabling the object and its timing, recording of data. </summary>

class NetworkObject
{
public:

	NetworkObject()
	{
		m_on = true;
		m_recording = false;
		m_initialized = false;
		m_timing = false;
		m_name = "";

		m_mpiCommunicator = NULL;//NETWORK_COMM_WORLD;
	}

	bool IsOn()
	{
		return m_on;
	}

	bool IsRecording()
	{
		return m_recording;
	}

	bool IsRecording(string recordParameter)
	{
		/*if(recordParameter.compare(m_recordParameter) == 0)
			return true;
		else return false;*/
		return m_recordParameters[recordParameter];
	}

	bool IsTiming()
	{
		return m_timing;
	}

	void SwitchOnOff(bool on)
	{
		m_on = on;
	}

	Network* network()
	{
		return m_network;
	}

	void network(Network* net)
	{
		m_network = net;
	}

	virtual vector<vector<float> > GetValuesToRecord() 
	{ 
		return vector<vector<float> >(0); 
	}

	void SetRecording(bool on)
	{
		m_recording = on;
	}

	bool IsInitialized()
	{
		return m_initialized;
	}

	void SetInitialized(bool initialized)
	{
		m_initialized = initialized;
	}

	virtual void Dispose() // manual destructor used for reset
	{
		m_initialized = false;
	}

	void SetRecordingParameter(string recordParameter)
	{
		m_recordParameters[recordParameter] = true;
	}
	
	// Timing functions - will be moved to timing (possibly parent) class
	virtual void SetTiming(bool on, Network* network)
	{
		m_network = network;

		if(on == true && m_timing == false)
		{
			m_analysisTiming = new AnalysisTiming(m_network);
			m_analysisTiming->SetTiming(on);
		}
		else if(on == false && m_timing == true)
		{
			delete m_analysisTiming;
		}

		m_timing = on;
	}

	void SetName(string name)
	{
		m_name = name;
	}

	string GetName()
	{
		return m_name;
	}

	AnalysisTiming* Timing()
	{
		return m_analysisTiming;
	}

	void TimingStart(string name)
	{
		if(m_timing)
		{
			m_analysisTiming->TimingStart(name);
		}
	}

	void TimingStop(string name)
	{
		if(m_timing)
		{
			m_analysisTiming->TimingStop(name);
		}
	}

	AnalysisTiming* GetAnalysisTiming()
	{
		return m_analysisTiming;
	}

	void SetMPICommunicator(MPI_Comm* comm)
	{
		m_mpiCommunicator = comm;
	}

	MPI_Comm GetMPICommunicator()
	{
		return *m_mpiCommunicator;
	}

protected:

	MPI_Comm* m_mpiCommunicator; // may be moved away from here - not by default initialized (NULL) but set from outside.

	Network* m_network;
	bool m_on, m_recording, m_initialized; // flags for enabling/disabling object and recording
	map<string,bool> m_recordParameters; // which parameters that should have their values recorded (only used in MDS atm).

	// Timing/Timer of process time for this object 
	AnalysisTiming* m_analysisTiming;
	bool m_timing;

	string m_name; // not efficient, change
};