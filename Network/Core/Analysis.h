#pragma once

#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <mpi.h>

using namespace std;

class Network;
class Population;

class Analysis
{
	// derive from NetworkObject instead ?
public:

	Analysis()
	{
		m_name = "Analysis";
		m_on = true;
	}

	Analysis(Network* network)
	{
		m_on = true;
		m_name = "Analysis";
		m_network = network;
	}

	virtual ~Analysis()
	{
	}

	virtual void Initialize(){}
	virtual void Simulate() = 0;
	virtual void Finalize(){}

	void SwitchOnOff(bool on)
	{
		m_on = on;
	}

	bool IsOn()
	{
		return m_on;
	}

	virtual void Dispose()
	{
	}

	string GetName()
	{
		return m_name;
	}

	vector<vector<float> > GetResults()
	{
		return m_results;
	}

	vector<vector<float> > GetResults2()
	{
		return m_results2;
	}

	void SetName(string name)
	{
		m_name = name;
	}

	void network(Network* network)
	{
		m_network = network;
	}

	void SetOwnFileWriting(bool use, string filename)
	{
		m_useOwnFileWriting = use;
		m_ownFilename = filename;
	}

	void DoOwnFileWriting();


protected:

	string m_ownFilename;
	bool m_useOwnFileWriting;
	bool m_useSparseStorage;
	float m_sparseThreshold;

	Network* m_network;
	bool m_on;
	string m_name;
	vector<vector<float> > m_results;
	vector<vector<float> > m_results2;
};

class AnalysisTiming : public Analysis
{
public:

	AnalysisTiming(Network* network)
	{
		m_name = "Timing";
		m_network = network;
	}

	void Initialize() 
	{
		m_name = "Timings";		
	}

	void Simulate() {}
	void Finalize();

	void SetTiming(bool on)
	{
		m_timing = on;
	}

	bool IsTiming()
	{
		return m_timing;
	}

	// could also reduce to one fcn instead of start/stop
	virtual void TimingStart(string name)
	{
		if(m_timing)
		{
			m_timingStart[name] = MPI_Wtime();
		}
	}

	virtual void TimingStop(string name)
	{
		if(m_timing)
		{
			double time = MPI_Wtime() - m_timingStart[name];
			m_timingTot[name] += time;

			if(time>m_timingMax[name])
				m_timingMax[name] = time;
			
			if(fabs(m_timingMin[name])<EPS)
				m_timingMin[name] = 1e8;
			if(time<m_timingMin[name])
				m_timingMin[name] = time;
		}
	}

	double GetTimingTotal(string name)
	{
		return m_timingTot[name];
	}

	virtual string GetTimingString()
	{
		stringstream ss;
		map<string,double>::iterator p;
		
		int index = 0;
		for(p = m_timingTot.begin(); p != m_timingTot.end(); p++)
		{
			ss<<p->first<<", ";
			ss<<p->second<<", ";
			if(m_meanValues.size()>index)
				ss<<m_meanValues[index]<<", ";
			if(m_stdValues.size()>index)
				ss<<m_stdValues[index]<<", ";
			if(m_timingMax.size()>index)
				ss<<m_timingMax[p->first]<<", ";
			if(m_meanMaxValues.size()>index)
				ss<<m_meanMaxValues[index]<<", ";
			if(m_stdMaxValues.size()>index)
				ss<<m_stdMaxValues[index]<<", ";
			if(m_min.size()>index && m_locMin.size()>index && m_max.size()>index && m_locMax.size()>index)
				ss<<m_min[index]<<", "<<m_locMin[index]<<", "<<m_max[index]<<", "<<m_locMax[index];
			ss<<"\n";
			index++;
		}

		return ss.str();
	}

private:

	bool m_timing;
	map<string,double> m_timingStart;
	map<string,double> m_timingTot; // could be changed to vector
	map<string,double> m_timingMax;
	map<string,double> m_timingMin;
	vector<float> m_meanValues;
	vector<float> m_stdValues;
	vector<float> m_meanMaxValues;
	vector<float> m_stdMaxValues;
	vector<float> m_meanMinValues;
	vector<float> m_stdMinValues;
	vector<int> m_locMin;
	vector<int> m_locMax;
	vector<float> m_min;
	vector<float> m_max;
};

class AnalysisLayer : public Analysis
{
public:

	AnalysisLayer() // don't use
	{
		m_name = "AnalysisLayer";
	}

	AnalysisLayer(Population* Population); // use

	AnalysisLayer(vector<Population*> Populations); // use


	void SetSamplingRate(int startSamplingTimeStep,int samplingTimeSteps)
	{
		m_samplingTimeSteps = samplingTimeSteps;
		m_samplingStartTimeStep = startSamplingTimeStep;
	}

	virtual void Simulate();
	// to implement: help virtual functions for analysis on 
	// - all data
	// - local data

protected:

	vector<Population*> m_populations;
	
	bool m_continue;
	int m_samplingTimeSteps, m_currentSamplingTimeStep, m_samplingStartTimeStep;

};

class AnalysisProjection : public Analysis
{
public:

private:

};

class FisherRatio : public AnalysisLayer
	// from Luis' OB evaluation
{
public:

	FisherRatio(Population* Population);

	FisherRatio(Population* Population, bool useIncrementalVersion);

	FisherRatio(Population* Population1, Population* Population2, bool useIncrementalVersion);

	void SetCurrentPatternClass(int index)
	{
		m_currentPatternClass = index;
	}
	
	void Simulate();
	void Finalize();
	
private:

	void CalcFisherRatio();
	bool m_useIncrementalVersion;

	vector<vector<vector<float> > > m_data;
	int m_currentPatternClass;
};

class AnalysisDistance : public AnalysisLayer
	// using abs distance atm, not euclidean
{

public:

	enum Measure
	{
		Abs,
		Euclidean,
		Hamming
	};

	AnalysisDistance(Population* Population, Measure measure)
	{
		m_populations.push_back(Population);
		m_useOwnFileWriting = false;
		m_useSparseStorage = false;
		m_name = "Distance";
		m_measure = measure;
		m_nrConcatenatedItems = -1;
		m_nrPreviousTimeSteps = -1;
		m_samplingTimeSteps = -1;
		m_currentSamplingTimeStep = -1;
		m_samplingStartTimeStep = -1;
		m_collectRepresentations = false;
	}

	AnalysisDistance(Population* Population, bool useOwnFileWriting, bool useSparseStorage, float sparseThreshold, Measure measure)
	{
		m_populations.push_back(Population);
		m_useOwnFileWriting = useOwnFileWriting;
		m_useSparseStorage = useSparseStorage;
		m_sparseThreshold = sparseThreshold;
		m_name = "Distance";
		m_measure = measure;
		m_nrConcatenatedItems = -1;
		m_nrPreviousTimeSteps = -1;
		m_samplingTimeSteps = -1;
		m_currentSamplingTimeStep = -1;
		m_samplingStartTimeStep = -1;
		m_collectRepresentations = false;
	}

	void SetNrConcatenatedItems(int nrItems)
	{
		m_nrConcatenatedItems = nrItems;
	}
	
	void SetNrPreviousTimeSteps(int nrSteps)
	{
		m_nrPreviousTimeSteps = nrSteps;
	}

	void SetOwnFileWriting(bool use, string filename)
	{
		m_useOwnFileWriting = use;
		m_ownFilename = filename;
	}

	void AddSpecificRepresentation(vector<float> values)
	{
		//if(m_specificRepresentations.size()==0) // fix (!), correct why last ends up first
			m_specificRepresentations.push_back(values);
		//else
//			m_specificRepresentations.insert(m_specificRepresentations.end()-1,values);
	}

	void AddSpecificRepresentation(vector<float> values, int correspondingClassIndex)
	{
		m_specificRepresentations.push_back(values);
		m_specificCorrespondingClassIndexes.push_back(correspondingClassIndex);
	}

	vector<vector<float> > GetCollectedRepresentations()
	{
		return m_specificRepresentations;
	}

	void SetCollectRepresentations(bool on)
	{
		m_collectRepresentations = on;
	}

	bool GetCollectRepresentations()
	{
		return m_collectRepresentations;
	}

	void Simulate();
	void Finalize();

	// Segmentation specific functions
	vector<float> GetSegmOutputResult(vector<vector<float> > targets, int classStep, int dynStep);
	void FinalizeSegmResults(vector<vector<float> > targets, int classStep, int dynStep);
	//

private:

	int m_nrConcatenatedItems;
	int m_nrPreviousTimeSteps;

	bool m_collectRepresentations; // collects the incoming values in m_specificRepresentations if on

	Measure m_measure;
	vector<vector<float> > m_buffer;
	vector<float> m_buffer2;

	vector<vector<float> > m_specificRepresentations;
	vector<int> m_specificCorrespondingClassIndexes;
};