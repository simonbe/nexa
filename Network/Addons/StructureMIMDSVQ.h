#pragma once
#ifndef STRUCTMIMDSVQ_H
#define STRUCTMIMDSVQ_H

#include "NetworkStructure.h"

#include "NetworkBCPNN.h"
#include "NetworkKussul.h"
#include "NetworkCL.h"
#include "NetworkMDS.h"
#include "NetworkMI.h"
#include "NetworkVQ.h"
#include "NetworkBCM.h"
#include "NetworkCorr.h"
#include "NetworkAdaptation.h"
#include "Network.h"
#include "Meter.h"
//#include "KernelConnectivity.h"

class PopulationModifierAdaptation;

class StructureMIMDSVQ : public NetworkStructure
{
public:

	enum FeatureExtraction
	{
		CSL,
		CL,
		BCM,
		HebbAdapt,
		Recurr
	};

	StructureMIMDSVQ()
	{
		m_index = 0;
		m_featureExtraction = StructureMIMDSVQ::CSL;
		m_mdsUsePearson = false;
		m_mdsDimension = 10;
	}

	StructureMIMDSVQ(FeatureExtraction featExtractor, bool usePearson)
	{
		m_index = 0;
		m_featureExtraction = featExtractor;
		m_mdsUsePearson = usePearson;
		m_mdsDimension = 10;
	}

	void SetFeatureExtraction(FeatureExtraction method)
	{
		m_featureExtraction = method;
	}

	FeatureExtraction GetFeatureExtraction()
	{
		return m_featureExtraction;
	}

	void SetMDSMeasure(bool usePearson)
	{
		m_mdsUsePearson = usePearson;
	}

	bool UsePearson()
	{
		return m_mdsUsePearson;
	}

	void SetMDSDimension(int dimension)
	{
		m_mdsDimension = dimension;
	}
	
	void SetupStructure(Network* network, PopulationColumns* layerInput, int nrMiddleHypercolumns, int nrMiddleRateUnits, bool addInputLayerToNetwork, bool useSilentHypercolumns = false, float silentHypercolumnsThreshold = 0.0);//, int nrOutputHypercolumns, int nrOutputRateUnits);
	void SetupStructure(Network* network, int nrInputHypercolumns, int nrInputRateUnits, int nrMiddleHypercolumns, int nrMiddleRateUnits, bool addInputLayerToNetwork);//, int nrOutputHypercolumns, int nrOutputRateUnits);
	void SetupMeters(int mpiRank, int mpiSize);//, Storage::FilePreference fileType);

	vector<PopulationColumns*> Layers()
	{
		return m_layers;
	}
	
	LayerMDS* MDS()
	{
		return m_MDS;
	}

	LayerVQ* VQ()
	{
		return m_VQ;
	}

	ProjectionModifierMDS* MDSHypercolumns()
	{
		return m_mdsHypercolumns;
	}

	// this way of calling will be changed
	ProjectionModifierCL* CompLearn()
	{
		return m_compLearn;
	}

	ProjectionModifierCSL* CSLLearn()
	{
		return m_cslLearn;
	}

	ProjectionModifierBCM* BCMLearn()
	{
		return m_bcmLearn;
	}

	PopulationModifierAdaptation* Adaptation()
	{
		return m_adaptation;
	}

	ProjectionModifier* FeatureExtractor()
	{
		if(m_featureExtraction == StructureMIMDSVQ::CL)
			return m_compLearn;
		else if(m_featureExtraction == StructureMIMDSVQ::BCM)
			return m_bcmLearn;
		else if(m_featureExtraction == StructureMIMDSVQ::HebbAdapt)
			return m_bcpnn;
		else if(m_featureExtraction == StructureMIMDSVQ::Recurr)
			return m_bcpnn;
		else
			return m_cslLearn;
	}

	ProjectionModifierMIHypercolumn* MIHypercolumns()
	{
		return m_miHypercolumns;
	}

	ProjectionModifierMIRateUnit* MIRateUnits()
	{
		return m_miRateUnits;
	}

	ProjectionModifierPearson* Pearson()
	{
		return m_pearson;
	}

	void SetIndex(int index)
	{
		m_index = index;
	}

	void SetRecording(bool on);

	PopulationColumns* GetLayer(int index)
	{
		return m_layers[index];
	}

	void SetTiming(bool on);

protected:

	Network* m_network;
	vector<PopulationColumns*> m_layers;
	LayerMDS* m_MDS;
	LayerVQ* m_VQ;
	ProjectionModifierMDS* m_mdsHypercolumns;
	ProjectionModifierMIHypercolumn* m_miHypercolumns;
	ProjectionModifierMIRateUnit* m_miRateUnits;
	ProjectionModifierCL* m_compLearn;
	ProjectionModifierCSL* m_cslLearn;
	ProjectionModifierBCM* m_bcmLearn;
	ProjectionModifierBcpnnOnline* m_bcpnn;
	ProjectionModifierPearson* m_pearson; // between minicolumns
	PopulationModifierAdaptation* m_adaptation;

	float m_featExtrRecurrProb;

	WTA* m_wta;
	SoftMax* m_softmax;

	FeatureExtraction m_featureExtraction;
	bool m_mdsUsePearson;
	int m_mdsDimension;

	int m_index;
};


#endif