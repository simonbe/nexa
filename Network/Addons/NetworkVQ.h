#pragma once
#ifndef NETWORKVQ_H
#define NETWORKVQ_H

#include "NetworkPopulationModifier.h"
//#include "Parameter.h"

using namespace std;

class LayerVQ;
class CSL;


class CSL
{
public:

	CSL()
	{
		m_epsilon = 0.01;//0.001;
		m_eta = 0.01f;//0.0001;//0.5;//0.001;
	}

	void Initialize(vector<vector<float> >* x, int nrCodeVectors, int mpiRank, int mpiSize, bool printOutResult, NetworkObject* parent, MPI_Comm comm);

	/*~CSL()
	{

	}*/

	vector<float> NewCodeVector(vector<float> oldCodeVector, vector<float> inputVector);
	vector<float> NewCodeVectorVector(vector<float> oldCodeVector, vector<float> inputVector, vector<float> eta);

	void Selection(vector<int> indexes, vector<float> D);
	vector<int> GetHighestIndexes(bool highestAndLowest, int nrValues, vector<float> data);
	unsigned int GetSeed();
	void ParallelDistortionCalculation(float *distortion, std::vector<float> *D);
	void ParallelCompetitiveLearning();
	void Step(vector<vector<float> >* x);
	float RRValue(vector<float> x1, vector<float> x2);
	int GetClosestCodeVectorIndex(std::vector<float> data, std::vector<std::vector<float> > codeVectors);
	std::vector<std::vector<float> > GetClosestCodeVectorMultipleIndexesAndValues(vector<float> data, vector<vector<float> > codeVectors, int maxSize);
	vector<int> GetClosestCodeVectorIndexes(std::vector<std::vector<float> > data);
	vector<float> GetClosestCodeVectorIndexAndValue(vector<float> data, vector<vector<float> > codeVectors);

	void SaveState();
	void LoadState();

	void SetMPIParameters(int mpiRank, int mpiSize)
	{
		m_mpiRank = mpiRank;
		m_mpiSize = mpiSize;
	}

	/*void SetInputData(vector<vector<float> >* x)
	{
		m_x = x;
	}*/

	vector<vector<int> > GetTotalWinnersUnits()
	{
		return m_totalWinnersUnits;
	}

	vector<vector<int> > GetTotalWinnersUnits(int nrOverlaps);

	vector<vector<float> >* GetCodeVectors() // make m_c pointer
	{
		return &m_c;
	}

	void SetEta(float eta, bool useSelection)
	{
		m_eta = eta;
		m_useSelection = useSelection;
	}

	void SetEtaVector(vector<float> eta, bool useSelection)
	{
		m_etaVector = eta;
		m_useSelection = useSelection;
	}

	// application specific
	void AddToCodeVector(int cvIndex, vector<float> values)
	{
		for(int i=0;i<values.size();i++)
		{
			m_c[cvIndex][i] += values[i];
		}
	}

	//void Reset();
	void DriftAdaptStep(vector<float> x, float eta);
	void DriftAdaptStepVector(vector<float> x, vector<float> eta);

	float GetDistortion() { return m_currentDistortion; }
	int GetSelection()
	{
		if(m_S.size()>0)
			return m_S[m_M];
		else return -1;
	}

private:

	MPI_Comm m_comm; // the nodes involved in this instance of csl
	NetworkObject* m_parent;

	bool m_printOutResult;

	int m_nrCodeVectors;
	vector<vector<int> > m_totalWinnersUnits;

	// CSL

	bool m_useSelection;
	float m_epsilon;
	float m_eta;
	vector<float> m_etaVector;
	float m_distortion;

	int m_iterations;
	int m_M;
	int m_N; // current nr of codewords
	int m_k; // dimension of input and codevectors
	float m_currentDistortion;
	int m_mpiSize,m_mpiRank;

	std::vector<std::vector<float> >* m_x; // input data
	std::vector<std::vector<float> > m_c; // codebook
	std::vector<std::vector<float> > m_cSaved; // saved state codebook

	vector<int> m_S; // Selection size per time step (declining)

	vector<int> m_winners;
	vector<int> m_total_winners;
};


class ConnectionModifierCSL : public ConnectionModifier
{
public:

	ConnectionModifierCSL();

	~ConnectionModifierCSL()
	{
		for(int i=0;i<m_csl.size();i++)
			delete m_csl[i];

		for(int i=0;i<m_comms.size();i++)
			delete m_comms[i];
	}

	void Initialize(Connection* connection);
	void SetConnection(Connection* c);

	void SetMaxPatterns(int maxPatterns)
	{
		m_maxPatterns = maxPatterns;
	}

	void Simulate(UnitModifier* e);
	void Modify();
	void UpdateWeights(int i, int currentIndex);

	void Reset()
	{
		for(int i=0;i<m_csl.size();i++)
		{
			delete m_csl[i]; // ok ? (!)
		}

		m_csl.clear();

		for(int i=0;i<m_comms.size();i++)
			delete m_comms[i];

		m_comms.clear();

		m_hcIndexes.clear();
		m_hcAllIndexes.clear();
		m_mcReprHc.clear();
		m_mpiLocSize.clear();
		m_mpiLocRank.clear();

		for(int i=0;i<m_x.size();i++)
			m_x[i].clear();
		m_x.clear();

		m_nrCodeVectors.clear();
		m_firstRun.clear();
		m_idsPost.clear();
		for(int i=0;i<m_idsPre.size();i++)
			m_idsPre[i].clear();

		m_idsPre.clear();
		m_idsIndexInHc.clear();
	}

	void DriftAdaptStep(float eta);
	void DriftAdaptStepVector(vector<float> eta);

	// SetEta/SetEtaVector are currently application specific (sensor drift)
	void SetEta(float eta, bool useSelection)
	{
		m_csl[0]->SetEta(eta, useSelection);
	}

	void SetEta(float eta)
	{
		m_csl[0]->SetEta(eta, true);
	}

	void SetEtaVector(vector<float> eta, bool useSelection)
	{
		m_csl[0]->SetEtaVector(eta, useSelection);
	}

	vector<vector<float> >* GetCodeVectors()
	{
		return m_csl[0]->GetCodeVectors();
	}

	// application specific
	void AddToCodeVector(int cvIndex, vector<float> values)
	{
		m_csl[0]->AddToCodeVector(cvIndex, values);
	}

	void SaveState();
	void LoadState();

	// only gives info about one of the connection sets atm.
	float GetDistortion()
	{
		return m_csl[0]->GetDistortion();
	}

	int GetSelection()
	{
		return m_csl[0]->GetSelection();
	}

private:

	vector<CSL*> m_csl;
	vector<int> m_hcIndexes; // indexes of the output hcs
	vector<int> m_hcAllIndexes; // index of each mc in hc - which it belongs to (in m_csl)
	vector<long> m_mcReprHc; // First minicolumn id representing hypercolumn (used to retrieve input values - all other mcs in hc will have same input)
	vector<MPI_Comm*> m_comms; // all comms used in the output hcs
	vector<int> m_mpiLocSize;
	vector<int> m_mpiLocRank;

	vector<vector<vector<float> > > m_x; // semi-online
	int m_maxPatterns;
	vector<int> m_nrCodeVectors;

	vector<bool> m_firstRun;
	
	vector<long> m_idsPost;
	vector<vector<long> > m_idsPre;
	map<long,int> m_idsIndexInHc;

	Connection* m_connectionFixed;
};

// shell class to let LayerVQ modify connections (currently only allowing ConnectionModifier to do it)
class ConnectionModifierVQ : public ConnectionModifier
{
public:

	ConnectionModifierVQ(LayerVQ* layerVQ);

	void Simulate(UnitModifier* e) {};
	void Modify(); // sends forward to LayerVQ::ModifyConnection

private:
	LayerVQ* m_layerVQ;
};



class LayerVQ : public PopulationModifier
{
public:

	enum VQType
	{
		VQStandard,
		VQCSL
	};

	LayerVQ(int nrCodeVectors, LayerVQ::VQType type)
	{
		m_nrCodeVectors = nrCodeVectors;
		m_firstRun = true;
		m_type = type;
		m_isRateUnitPreConnections = false;
		m_hasCheckedPreConnections = false;

		m_csl = NULL;

		if(type == LayerVQ::VQStandard)
		{
			trgnhyp = nrCodeVectors;
			m_MDSK = 0.5;//2e-2f;
			m_NWIN = 2;
			prn = sprn = 1.0;
			m_RLIM = 0.65f;
		}
		else if(type == LayerVQ::VQCSL)
		{
			m_csl = new CSL();
		}

		m_eventConnectionVQ = new ConnectionModifierVQ(this);

		m_name = "LayerVQ";
		m_nrOverlaps = 1;
	}

	void SetNrOverlaps(int nrOverlaps)
	{
		m_nrOverlaps = nrOverlaps;
	}

	~LayerVQ()
	{
		//delete m_eventConnectionVQ; // added remotely
		delete m_csl;
	}

	void Simulate();
	void ModifyConnections(Connection* connections);
	float GetError(){
		return m_csl->GetDistortion();
	}

	ConnectionModifierVQ* GetConnectionModifier();
	vector<vector<float> > GetValuesToRecord();


	float GetDistortion()
	{
		return m_distortion;
	}

	int GetSelection()
	{
		if(m_S.size()>0)
			return m_S[m_M];
		else
			return -1;
	}

	CSL* GetCSL()
	{
		return m_csl;
	}

//	Parameter *lambda;
private:

	void InitiateVQUnits();
	ConnectionModifierVQ* m_eventConnectionVQ;

	int m_nrCodeVectors;
	LayerVQ::VQType m_type;
	vector<vector<float> >* m_data;
	bool m_firstRun;
	vector<vector<int> > m_totalWinnersUnits;
	bool m_isRateUnitPreConnections;
	bool m_hasCheckedPreConnections;

	// standard
	void updvqwin(int u);
	void updvq();
	int trgnhyp,srcnhyp;
	float m_MDSK;
	int m_NWIN;
	int m_MDSDIM;
	vector<float> vqd, vqn, vqdst;
	vector<vector<float> > vqu, vquDiff;
	vector<vector<int> > oldhuse, huse;
	long simstep;
	float prn, sprn;
	vector<int> vqwin;
	float m_RLIM;
	float euclid(vector<float> v1,vector<float> v2,int n);
	vector<long> m_mpiParts;
	int m_sizeThisNode; // (absolute) size of input data node handles
	int m_mpiNrNodes;
	bool m_distFracChanged;
	float m_distFrac;

	vector<vector<int> > m_winnersUnits; // inputs belonging to which winning units (this node takes care of), index-based

	// CSL

	CSL* m_csl;

	float m_epsilon;
	float m_eta;
	int m_iterations;
	int m_M;
	int m_N; // current nr of codewords
	int m_k; // dimension of input and codevectors
	float m_currentDistortion;
	int m_mpiSize,m_mpiRank;

	float m_distortion;

	std::vector<std::vector<float> > m_x; // input data
	std::vector<std::vector<float> > m_c; // codebook
	vector<int> m_S; // Selection size per time step (declining)

	vector<int> m_winners;
	vector<int> m_total_winners;

	void CSL_Initiate();
	vector<float> CSL_NewCodeVector(vector<float> oldCodeVector, vector<float> inputVector);
	void CSL_Selection(vector<int> indexes, vector<float> D);
	vector<int> CSL_GetHighestIndexes(bool highestAndLowest, int nrValues, vector<float> data);
	unsigned int GetSeed();
	void CSL_ParallelDistortionCalculation(float *distortion, std::vector<float> *D);
	void CSL_ParallelCompetitiveLearning();
	void CSL_Step();
	float CSL_RRValue(vector<float> x1, vector<float> x2);
	int CSL_GetClosestCodeVectorIndex(std::vector<float> data, std::vector<std::vector<float> > codeVectors);
	vector<float> CSL_GetClosestCodeVectorIndexAndValue(vector<float> data, vector<vector<float> > codeVectors);

	// overlapping winners
	int m_nrOverlaps;

	// could place in parent
	vector<vector<float> > m_recordedValues;
};

class AnalysisVQ: public AnalysisLayer
{
public:

	AnalysisVQ(LayerVQ* vq)
	{
		m_layerVQ = vq;
		m_name = "VQ-Layer-Distortion/Selection";
		m_eventConnectionCSL = NULL;
	}

	AnalysisVQ(ConnectionModifierCSL* vq)
	{
		m_eventConnectionCSL = vq;
		m_name = "VQ-Connection-Distortion/Selection";
		m_layerVQ = NULL;
	}

	void Simulate()
	{
		if(m_layerVQ == NULL)
		{
			vector<float> f(2);
			if(m_eventConnectionCSL == NULL)
			{
				f[0] = -1;
				f[1] = -1;
			}
			else
			{
				f[0] = m_eventConnectionCSL->GetDistortion();
				f[1] = (float)m_eventConnectionCSL->GetSelection();
			}


			m_results.push_back(f);
		}
		else
		{
			vector<float> f(2);
			//f[0] = m_layerVQ->Layer()->network()->GetCurrentTimeStep();
			f[0] = m_layerVQ->GetDistortion();
			f[1] = (float)m_layerVQ->GetSelection();

			m_results.push_back(f);
		}
	}

	void SetVQLayer(LayerVQ* vq)
	{
		m_layerVQ = vq;
	}

private:

	LayerVQ* m_layerVQ;
	ConnectionModifierCSL* m_eventConnectionCSL;
};

#endif
