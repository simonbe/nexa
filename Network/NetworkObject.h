#pragma once

#include "Network.h"
#include "Analysis.h"

class Network;
class AnalysisTiming;

using namespace std;

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

		//m_timingTot = vector<double>(1,0.0); // move
		//m_timingStart = vector<double>(1,0.0);
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

	virtual std::vector<std::vector<float> > GetValuesToRecord() 
	{ 
		return std::vector<std::vector<float> >(0); 
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

	MPI_Comm* m_mpiCommunicator; // not efficient, change

	Network* m_network;
	bool m_on, m_recording, m_initialized;
	map<string,bool> m_recordParameters;

	// timer variables, create new parent timing class? since not interesting for all network objects
	// or move into timing class, initialize if turned on

	AnalysisTiming* m_analysisTiming;
	bool m_timing;

	string m_name; // not efficient, change
};