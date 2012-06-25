#pragma once


#include "NetworkBCPNN.h"
#include "NetworkKussul.h"
#include "NetworkCL.h"
#include "NetworkMDS.h"
#include "NetworkMI.h"
#include "NetworkVQ.h"
#include "Network.h"
#include "Meter.h"
#include "Storage.h"

class NetworkTests
{
public:
	NetworkTests();
	void NetworkTestBCPNNRecurrent();
	void NetworkTestOR2ORN2MT(int mpiRank, int mpiSize);
	void NetworkTestOlfCortex(int mpiRank, int mpiSize);
	void NetworkTestInclTimingBCPNNRecurrent();
	void NetworkTestInclTimingMIMDSVQVisual();
	void NetworkTestMDSVQ(int mpiRank, int mpiSize);
	void NetworkTestHierarchyLayers(int mpiRank, int mpiSize);
	void NetworkTestMNISTRecurrent(int mpiRank, int mpiSize);
	void NetworkTestMNISTClassification(int mpiRank, int mpiSize);
	void NetworkTestIF(int mpiRank,int mpiSize);
	void NetworkTestPearson(int mpiRank, int mpiSize);
	void NetworkTestTrieschAndFoldiak(int mpiRank, int mpiSize);
	void NetworkTestSanger(int mpiRank, int mpiSize);

	void NetworkTestSwitching(int mpiRank, int mpiSize);
	//void NetworkTestSensorModelInput(int mpiRank, int mpiSize);
private:
	vector<float> toBinary(int nr, int total);
	vector<float> toBinary(vector<float> data, int nrHc, int nrMc);
	int m_architecture;
};