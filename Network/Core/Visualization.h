#pragma once

#include "Network.h"
#define PARALLEL

#if VISIT_AVAILABLE == 1
#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>
#pragma comment(lib, "Ws2_32.lib")
#endif

//#include "stubs.h"
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <mpi.h>

using namespace std;

#define SIM_RUNNING 1
#define SIM_STOPPED 0


struct simulation_data
{
    int     cycle;
    double  time;
    int     runMode;
    int     done;
    
    int nrHypercolumns;
    int nrRateUnits;
    float   *values;
    int     rMeshDims[3];
    float   *rmeshX, *rmeshY;
//    NetworkDemoVis1 *network;
    
#ifdef PARALLEL
    int     par_rank;
    int     par_size;
#endif
    simulation_data();
    ~simulation_data();
};

void simulation_data_ctor(simulation_data *sim);

void simulation_data_dtor(simulation_data *sim);


// Handles setting up and managing Projection with VisIt in-situ visualization (the libsim library)

class Visualization
{
	// derive from NetworkObject instead ?
public:

#if VISIT_AVAILABLE == 1

	Visualization()
	{
		m_name = "Visualization";
		m_on = true;
	}

#endif

private:

#if VISIT_AVAILABLE == 1

	bool m_on;
	string m_name;

#endif
};

