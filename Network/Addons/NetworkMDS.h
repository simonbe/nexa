#pragma once
#ifndef NETWORKMDS_H
#define NETWORKMDS_H
#define MDS_EPS 1e-3
#define MDS_NUM 10
#include "NetworkProjectionModifier.h"

class ProjectionModifier;

/// <summary>	Multi-dimensional scaling. </summary>

class LayerMDS : public PopulationModifier
{
public:

	LayerMDS(int MIDimension, int MDSDimension, Network* net)
	{
		this->network(net);
		m_mdsDimension = MDSDimension;
		m_miDimension = MIDimension;
		m_prevDist=0;
		m_distUnchanged=0;
		m_Xi = vector<vector<float> >(MIDimension,vector<float>(MDSDimension));
		m_diffXi = vector<vector<float> >(MIDimension,vector<float>(MDSDimension));
		srand(0); // make sure all get same start matrix

		for(int i=0;i<MIDimension;i++)
			for(int j=0;j<MDSDimension;j++)
			{
				m_Xi[i][j] = (rand()/(float)RAND_MAX) - 0.5; // nodes have same seed
			}

		this->network()->SetSeed(); // go back to default
		m_name = "LayerMDS";
	}
	
	void Simulate();

	vector<vector<float> >* GetCurrentXi()
	{
		return &m_Xi;
	}

	void AddDiffXi(vector<vector<float> > xi)
	{
		for(int i=0;i<xi.size();i++)
			for(int j=0;j<xi[0].size();j++)
				m_diffXi[i][j] += xi[i][j];
	}

	void SetLocalStress(float localStress)
	{
		m_localStress = localStress;
	}

	void AddDiffXi(vector<float> xi, int location)
	{
		for(int i=0;i<xi.size();i++)
		{
			m_diffXi[location][i] += xi[i];
		}
	}

	int GetDimension()	{ return m_mdsDimension; }

	vector<vector<float> >* GetOutput()
	{
		return &m_Xi;
	}

	vector<vector<float> > GetValuesToRecord() 
	{
		//vector<float> time(1);
		//time[0] = this->network()->GetCurrentTimeStep();
		if(IsRecording("stress") == true)
		{
			vector<vector<float> > data = m_recordedValues; // pop back instead with iterator

			for(int i=0;i<m_recordedValues.size();i++)
			{
				m_recordedValues[i].clear();
			}
			m_recordedValues.clear();

			return data;
		}
		else if(IsRecording())
		{
			vector<vector<float> > data;
			//data.push_back(time);

			for(int i=0;i<m_Xi.size();i++)
				data.push_back(m_Xi[i]);

			return data;//m_Xi;
		}
	}

	void Reset() // Used to be Dispose: will work as a reset here (!)
	{
		for(int i=0;i<m_miDimension;i++)
			for(int j=0;j<m_mdsDimension;j++)
			{
				m_Xi[i][j] = (rand()/(float)RAND_MAX) - 0.5; // nodes have same seed
			}

		//m_diffXi = vector<vector<float> >(m_miDimension,vector<float>(m_mdsDimension));

		m_initialized = false;
	}

	float GetCurrentStress()
	{
		return m_totStress; // total change and not stress ?
	}

	float GetCurrentStressChange()
	{
		return m_totChange;
	}

private:

	vector<float> m_Xm;
	vector<vector<float> > m_Xi;
	vector<vector<float> > m_diffXi;
	
	float m_totChange;
	float m_totStress;
	float m_localStress;
	float m_prevDist;
	int m_distUnchanged;
	int m_mdsDimension, m_miDimension;
	vector<vector<float> > m_recordedValues;
};

/// <summary>	Analysis of a running Multi-dimensional scaling. Retrieves specific values from class for easy storage on disk. </summary>

class AnalysisMDS : public AnalysisLayer
{
public:

	AnalysisMDS(LayerMDS* mds)
	{
		m_layerMDS = mds;
		m_name = "MDS-StressChange-and-Stress";
	}

	void Simulate()
	{
		vector<float> f(3);
		f[0] = m_layerMDS->GetPopulation()->network()->GetCurrentTimeStep();
		f[1] = m_layerMDS->GetCurrentStressChange();
		f[2] = m_layerMDS->GetCurrentStress();
		
		m_results.push_back(f);
	}

	void SetMDSLayer(LayerMDS* mds)
	{
		m_layerMDS = mds;
	}

private:

	LayerMDS* m_layerMDS;
};

/// <summary>	Multi-dimensional scaling. Used in combination with LayerMDS. </summary>

class ProjectionModifierMDS : public ProjectionModifier
{
public:

	ProjectionModifierMDS()
	{
		ProjectionModifierMDS(false);
		m_name = "ProjectionModifierMDS";
	}

	ProjectionModifierMDS(bool usePearson)
	{
		sprn = prn = 1;
		m_MDSK = 2e-2;//2e-5;//2e-2; // Change in SetUpdateSize
		m_eventId = 4;
		m_firstRun = true;
		m_usePearson = usePearson;
		if(usePearson)
			m_name = "ProjectionModifierMDS_Pearson"; // Uses Pearson correlation as input
		else
			m_name = "ProjectionModifierMDS_MI"; // Uses mutual information as input

		m_useThreshold = false;
		m_threshold = 0.9;
	}
		
	~ProjectionModifierMDS()
	{
 		//delete m_Xi;
	}

	/// <summary>	If threshold is used, points will only affect each other if they are within threshold distance. </summary>
	///
	/// <remarks>	Post Lazarus, 9/4/2012. </remarks>
	///
	/// <param name="useThreshold">	true to use threshold. </param>
	/// <param name="threshold">   	Threshold value. </param>

	void SetUseThreshold(bool useThreshold, float threshold)
	{
		m_useThreshold = useThreshold;
		m_threshold = threshold;
	}

	void AddParentPopulationModifier(PopulationModifier* e);
	void SetProjection(Projection* c);
	void Simulate(UnitModifier* e) {}
	void Modify();

	vector<float> GetXm(){ return Xm; }
	void TranslateToOrigin(vector<float> XmTot);
	float GetMeands() { return meands; }

	void Dispose();
	//~ProjectionModifierMDS();

	/// <summary>	Magnitude of update in one iteration of MDS. </summary>
	///
	/// <param name="mdsk">	Update value. </param>

	void SetUpdateSize(float mdsk)
	{
		m_MDSK = mdsk;
	}

private:

	bool m_usePearson; // uses pearson's corr instead of mi
	bool m_useThreshold;
	float m_threshold;

	bool m_firstRun;
	float taupdt;
	float sprn,prn;
	float m_MDSK;

	int m_mdsDim;

	vector<float> Xi_i, Xi_j, oldXi_i, Xm;
	vector<vector<float> >* m_Xi;
	int srcnhyp;
	float meands;

	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;

	Projection* m_projectionFixed;
};

#endif