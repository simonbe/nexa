#include "Network.h"

#if VISIT_AVAILABLE == 1

#include <VisItControlInterface_V2.h>
#include <VisItDataInterface_V2.h>
#pragma comment(lib, "Ws2_32.lib")

#include "stubs.h"

#define SIM_RUNNING 1
#define SIM_STOPPED 0
#define VISIT_COMMAND_PROCESS 0
#define VISIT_COMMAND_SUCCESS 1
#define VISIT_COMMAND_FAILURE 2

class NetworkDemoVis1 : public Network
{
public:
    NetworkDemoVis1();
    ~NetworkDemoVis1();
    
private:  
    simulation_data sim;  
    vector<RateUnit*> minicolumns;
    int step;
    
    static void ControlCommandCallback(const char *cmd, const char *args, void *cbdata);
    void SimulationArguments(int argc, char **argv);
    void NetworkSetupStructure();
    void NetworkSetupMeters();
    void NetworkSetupParameters();
    void NetworkRun();
    void simulate_one_timestep();      
    
    static const char *cmd_names[];
    void MeshConfigure();
#ifdef PARALLEL
    int ProcessVisItCommand();
    static int visit_broadcast_int_callback(int *value, int sender);
    static int visit_broadcast_string_callback(char *str, int len, int sender);
#endif
    static void BroadcastSlaveCommand(int *command);
    static void SlaveProcessCallback();
    
    visit_handle static SimGetDomainList(const char *name, void *cbdata);
    visit_handle static SimGetMesh(int domain, const char *name, void *cbdata);
    visit_handle static SimGetVariable(int domain, const char *name, void *cbdata);
};
#endif