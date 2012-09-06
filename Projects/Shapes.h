#include <mpi.h>
#include "Network.h"
#include "NetworkBCPNN.h"
#include "NetworkAdaptation.h"
#include "NetworkDepression.h"
#include "StructureReTIDe.h"

class NetworkEventHandlerLocal; //see def below
class StructureReTIDe;

class ShapesNetwork : public Network
{
public:

	vector<RateUnit*> minicolumns;
	int numEvents;
	int rank;
	int numEventsTotal;

	vector< vector<float> > ReadPatternSequence(const char*);

private:

	ProjectionModifierBcpnnOnline* recConnBcpnn;
	PopulationModifierAdaptation* recConnAdaptation;
	ProjectionModifierDepression* recConnDepression;
	FullConnectivity* recConn;
	Threshold* layer1Threshold;

	StructureReTIDe* reTIDe;

	double timestep;

	void Train();
	void Test();
	void NormalizeWeights();
	void RemoveUnusedSynapses();

	void NetworkSetupStructure();
	void NetworkSetupMeters();
	void NetworkSetupParameters();
	void NetworkRun();
};
