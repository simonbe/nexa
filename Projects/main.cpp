#undef SEEK_SET
#undef SEEK_END
#undef SEEK_CUR

//#define WIN32

//#include "NetworkUnitTests.h"
//#include "NetworkDemoVis1.h"
#include "Network.h"
//#include "NetworkTests.h"
#include "NetworkMNIST.h"

#include "stubs.h"

//#include "D:/Databases/Olfaction/Bernhard/ob_output/ob_output/NetworkOlfaction_BK.h"
//#include "OB_OCTX_Connectivity/NetworkOlfaction_BK.h"
//#include "OlfactionSystem.h"

//#include "NetworkScalingDemos.h"
#include <mpi.h>
//#include <VisItControlInterface_V2.h>
//#include <VisItDataInterface_V2.h>
//#pragma comment(lib, "Ws2_32.lib")
#if VISIT_AVAILABLE == 1
#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>
#pragma comment(lib, "Ws2_32.lib")
#endif

// memory leaks detection
//#define _CRTDBG_MAP_ALLOC
//#include <crtdbg.h>
//#include <stdlib.h>

// tail -n 10 core.0 | addr2line -e project1

using namespace std;

//class NetworkDemoVis1;
//class Network;

/*void SimulationArguments(int argc, char **argv)
{
    int i;
    for (i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-dir") == 0 &&
                (i + 1) < argc)
        {
            VisItSetDirectory(argv[i + 1]);
            ++i;
        }
        else if (strcmp(argv[i], "-options") == 0 &&
                (i + 1) < argc)
        {
            VisItSetOptions(argv[i + 1]);
            ++i;
        }
#ifdef VISIT_CONTROL_INTERFACE_V2_H
        else if (strcmp(argv[i], "-trace") == 0 &&
                (i + 1) < argc)
        {
#ifdef PARALLEL
            int rank;
            char *tmpfile = NULL;
            tmpfile = (char*) malloc(strlen(argv[i + 1]) + 10);
            MPI_Comm_rank(NETWORK_COMM_WORLD, &rank);
            sprintf(tmpfile, "%s.%04d", argv[i + 1], rank);
            VisItOpenTraceFile(tmpfile);
            free(tmpfile);
#else
            VisItOpenTraceFile(argv[i + 1]);
#endif
            ++i;
        }
#endif
    }
}
*/

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
	int mpiSize = -1, mpiRank = -1; // not set here anymore

#if MUSIC_AVAILABLE == 1
	//This has to be done first of all. Calls MPI_Init.
	//MUSIC::Setup* setup = new MUSIC::Setup(argc, argv);
	//Tell the simulator to use the MUSIC communicator instead of //NETWORK_COMM_WORLD
	//Network::MPIComm = setup->communicator();
#endif

//	NetworkIFTests ifTests;
//	ifTests.NetworkTestSTDP();

	//ifTests.NetworkTestSingleIF(mpiRank,mpiSize);
	//ifTests.NetworkTestNetworkIF(mpiRank,mpiSize);
	//ifTests.NetworkTestNetworkIF2(mpiRank,mpiSize);
	//ifTests.NetworkTestRetinaIF(mpiRank,mpiSize);
	//ifTests.NetworkTestNetworkIFInterneuron(mpiRank,mpiSize);
	//ifTests.NetworkTestNetworkIFDecorrelation();

//	NetworkTests tests;
//	tests.NetworkTestInclTimingBCPNNRecurrent();
	//tests.NetworkTestSanger(mpiRank,mpiSize);
	//tests.NetworkTestMNISTClassification(mpiRank,mpiSize);//NetworkTestIF(mpiRank,mpiSize);
	//tests.NetworkTestTrieschAndFoldiak(mpiRank,mpiSize);
	//tests.NetworkTestOlfCortex(mpiRank,mpiSize);
	//tests.NetworkTestOR2ORN2MT(mpiRank,mpiSize);

//	NetworkMNIST2 mnist2;
//	mnist2.Run();
	//NetworkMNIST mnist;
	//mnist.NetworkMNISTRun3(mpiRank,mpiSize);
	//mnist.NetworkMNISTRun1(mpiRank,mpiSize);
	//mnist.NetworkMNISTRunLateralMaps(mpiRank,mpiSize);
	//mnist.NetworkMNISTRunLateralMaps2(mpiRank,mpiSize);

//NetworkScalingDemos scaling;
//	scaling.RunAll();
//	NetworkScalingStrong scalingStrong;
//	scalingStrong.Run();

//	NetworkOlfactionScenariosClassification_PH olf;
//	olf.Start();

//		NetworkOlfaction olf;
//		olf.NetworkOlfactionRunWithOBModel2();
//	olf.NetworkOlfactionRunWithOBModel1();

	//olf.NetworkOlfactionRun3(mpiRank,mpiSize);
	//olf.NetworkOlfactionRunWithOBModel1();
//	olf.NetworkOlfactionRun3(mpiRank,mpiSize);

	//NetworkSensorDrift sensors;
	//sensors.NetworkSensorDriftRun1(mpiRank,mpiSize);
	//sensors.NetworkSensorDriftExamples(mpiRank,mpiSize);
	//sensors.NetworkSensorDriftRunCCCompare(mpiRank,mpiSize);
	//sensors.NetworkSensorDriftCollectiveRuns(mpiRank,mpiSize);

	//NetworkSetupConns setupConns;
	//setupConns.Example1();

//	NetworkTopDownDecorr decorr;
//	decorr.RunRule1();//RunHistoryDependencePlot();

//	decorr.RunSequence1();
//	decorr.RunTestClassification();
//	decorr.RunTestExcRecurr();
	//decorr.Run3();
//	decorr.Run3DataRandToCorr();
//	decorr.Run3DataRandToCorr_2ndRun();
//	NetworkTopDownDecorr_RandToCorr decorrRandToCorr;
//	decorrRandToCorr.Run();
//	NetworkTopDownDecorr_Invariance decorrInv;
//	decorrInv.Run();

//	NetworkTopDownDecorr_SequencePrediction seq;
//	seq.Run();
//	NetworkFridaySemantics semantics;
//	semantics.Run();

//	ShapesNetwork shapes;
//	shapes.Run();
//	NetworkProjectionTests testsConns;
//	testsConns.Run();

//	NetworkSynapsesTests tests;
//	tests.Run();

//	decorr.Run3DataInvariance_2ndRun();
//	decorr.RunTestHS();
//	decorr.RunTestHSInhib();
//	NetworkArithmetic arithmetic;
//	arithmetic.Run();
//	NetworkTES tes;
//	tes.Run();

//	NetworkTemporal temporal;
//	temporal.Run();

//	VisItSetupEnvironment();
//	NetworkDemoVis1 vis1;
//	vis1.Run();

	
//	NetworkDemoVis1 vis1;
//	vis1.Run();

	
//	NetworkDemoVis0 vis0;
//	vis0.Run();
	
//	UnitTests test;
//	test.Run();
	
	//NetworkDemoVis1 vis1;
	//vis1.Run();

//	NetworkDemoVis1 demoVis1;
//	demoVis1.Run();
//	OlfactionSystem olfsys;
//	olfsys.test();	

	
	/*if (argc < 7) {
		cout << "Error: Not enough arguments given!" << endl;
		cout << "Example usage: mpirun -np 8 ./TestNetworkOlfaction 1 80 4 10 FileName.dat" << endl;
	}*/
/*	
	int nrPatterns = 25;
	int nrInputHypercolumns = 25;
	int nrInputRateUnits = 8;
	int nrMiddleHypercolumns = 5;
	int nrMiddleRateUnits = 25;
	std::string filename = "C:/CurrentProjects/Network/Projects/OB_OCTX_Connectivity/mit_response_normalized_25_25_8_3_25";//"C:/CurrentProjects/Network/Projects/OB_OCTX_Connectivity/mit_response_normalized";//mit_response_normalized_np3";//mit_pattern_response_10x48";
	NetworkOlfaction_BK test;
	test.TrainWithSpikingOBData(nrPatterns, nrInputHypercolumns, nrInputRateUnits, nrMiddleHypercolumns, nrMiddleRateUnits, filename);
	*/

/*	NetworkOlfaction_BK test;
	test.TrainWithSpikingOBData();*/
	MPI_Finalize();

//	_CrtDumpMemoryLeaks();

	return 0;
}