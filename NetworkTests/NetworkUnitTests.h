#pragma once

#include "Network.h"
#include "NetworkBCPNN.h"

using namespace std;

class UnitTests
{
public:
	void Run();

private:

};


class NetworkProjectionTests : public Network
{
public:

/*	 NetworkProjectionTests();
    ~NetworkProjectionTests();*/
private:

	int m_nrHypercolumns;
	int m_nrRateUnits;

	// overridden
	void NetworkSetupStructure();
	void NetworkSetupMeters();
	void NetworkSetupParameters();
	void NetworkRun();
	
};


class NetworkSynapsesTests : public Network
{
public:
private:
	void NetworkSetupStructure();
	void NetworkRun();

	ProjectionModifierBcpnnOnline* m_bcpnn;
};
