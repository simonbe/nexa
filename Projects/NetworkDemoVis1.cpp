#include "NetworkDemoVis1.h"

#if VISIT_AVAILABLE == 1
const char *cmd_names[] = {"Run", "Stop", "Step", "Update", "Reset"};

#define PARALLEL

void NetworkDemoVis1::ControlCommandCallback(const char *cmd, const char *args, void *cbdata)
{
    simulation_data *sim = (simulation_data *) cbdata;

    if (strcmp(cmd, "Stop") == 0)
        sim->runMode = SIM_STOPPED;
    else if (strcmp(cmd, "Step") == 0)
        sim->network->simulate_one_timestep();
    else if (strcmp(cmd, "Run") == 0)
        sim->runMode = SIM_RUNNING;
    else if (strcmp(cmd, "Update") == 0)
    {
        VisItTimeStepChanged();
        VisItUpdatePlots();
    }
    else if (strcmp(cmd, "Reset") == 0)
    {
        sim->network->Reset();
        sim->cycle = 0; sim->time = 0;      
        VisItTimeStepChanged();
        VisItUpdatePlots();
    }
}

NetworkDemoVis1::NetworkDemoVis1()
{
    sim.network = this;
    /* Initialize environment variables. */
#ifdef PARALLEL
	
    MPI_Comm_rank(NETWORK_COMM_WORLD, &sim.par_rank);
    MPI_Comm_size(NETWORK_COMM_WORLD, &sim.par_size);

    /* Install callback functions for global communication. */
    VisItSetBroadcastIntFunction(visit_broadcast_int_callback);
    VisItSetBroadcastStringFunction(visit_broadcast_string_callback);
    /* Tell VSIL whether the simulation is parallel. */
    VisItSetParallel(sim.par_size > 1);
    VisItSetParallelRank(sim.par_rank);
	cout<<"Parallel VisIt."<<sim.par_rank; cout.flush();
#endif

    /* Write out .sim2 file that VisIt uses to connect. Only do it
     * on processor 0.
     */
#ifdef PARALLEL
    if (sim.par_rank == 0)
#endif
    {
        VisItInitializeSocketAndDumpSimFile("Network",
                                            "Parallel simulation visualization",
                                            "D:\\Network\\", NULL, NULL, NULL);
    }
}

// Defines network structure

void NetworkDemoVis1::NetworkSetupStructure()
{
    // total nrHypercolumns*nrRateUnits neural units
    sim.nrHypercolumns = 8;
    sim.nrRateUnits = 50;

    PopulationColumns* layer1 = new PopulationColumns(this, sim.nrHypercolumns, sim.nrRateUnits, PopulationColumns::Graded);
    this->AddPopulation(layer1); // this population/layer will have index 0

    // Population fully connected to itself by random synaptic weights
    FullConnectivity* full1 = new FullConnectivity();
    full1->SetRandomWeights(0, 10);
    layer1->AddPre(layer1, full1);

    // Winner-take-all operation
    WTA* wta = new WTA();
    layer1->AddPopulationModifier(wta);
    this->SetSeed(10);
}

// Setup outputs

void NetworkDemoVis1::NetworkSetupMeters()
{
    // example of file output
    Meter* layerMeter = new Meter("activity.csv", Storage::CSV, Storage::Standard);
    //layerMeter->AttachPopulation(this->GetLayer(0));
    //this->AddMeter(layerMeter);

    // also check simulation times
    this->SetTiming(true, this);
}

void NetworkDemoVis1::NetworkSetupParameters()
{
    // not used
}

visit_handle
SimGetMetaData(void *cbdata);

void NetworkDemoVis1::MeshConfigure()
{
    int hstart = minicolumns[0]->GetHypercolumnId();
    int hend = minicolumns[minicolumns.size() - 1]->GetHypercolumnId();
    int mstart = minicolumns[0]->GetUnitIdLocalInHypercolumn();
    int mend = minicolumns[minicolumns.size() - 1]->GetUnitIdLocalInHypercolumn();
    int count = minicolumns.size();

    sim.rMeshDims[1] = hend - hstart + 2;
    sim.rMeshDims[0] = mend - mstart + 2;
    sim.rMeshDims[2] = 1;
    sim.rmeshX = new float[sim.rMeshDims[0]];
    sim.rmeshY = new float[sim.rMeshDims[1]];

    for (int i = hstart; i <= hend + 1; i++)
    {
        sim.rmeshY[i - hstart] = float(i);
    }

    for (int i = mstart; i <= mend + 1; i++)
    {
        sim.rmeshX[i - mstart] = float(i);
    }

	cout<<"\n"<<sim.par_rank<<": "<<sim.rMeshDims[0]<<","<<sim.rMeshDims[1];cout.flush();
}

void NetworkDemoVis1::NetworkRun()
{
    int blocking, visitstate, err = 0;
    int nrTimeSteps = 100;

    // local minicolumns/neural units on this process for population 0
    minicolumns = ((PopulationColumns*)this->GetLayer(0))->GetRateUnits();
    
    MeshConfigure();
    
    sim.values = new float[minicolumns.size()];

    // simulate
    for (int i = 0; !sim.done && err == 0; i++)
    {
		cout.flush();
		cerr.flush();
        blocking = (sim.runMode == VISIT_SIMMODE_RUNNING) ? 0 : 1;
#ifdef PARALLEL
        if (sim.par_rank == 0)
#endif
            visitstate = VisItDetectInput(blocking, -1);
#ifdef PARALLEL
        MPI_Bcast(&visitstate, 1, MPI_INT, 0, NETWORK_COMM_WORLD);
#endif
        /* Do different things depending on the output from VisItDetectInput. */
        if (visitstate >= -5 && visitstate <= -1)
        {
            fprintf(stderr, "Can't recover from error!\n");
            err = 1;
			cerr.flush();
        }
        else if (visitstate == 0)
        {
		//	cout<<"'";cout.flush();
            simulate_one_timestep();
        }
        else if (visitstate == 1)
        {
			VisItSetupEnvironment();
            /* VisIt is trying to connect to sim. */
            if (VisItAttemptToCompleteProjection() == VISIT_OKAY)
            {
                fprintf(stderr, "VisIt connected\n");
                VisItSetCommandCallback(ControlCommandCallback, (void*) &sim);
                VisItSetSlaveProcessCallback(SlaveProcessCallback);
                VisItSetGetMetaData(SimGetMetaData, (void*) &sim);
                VisItSetGetMesh(SimGetMesh, (void*) &sim);
                VisItSetGetVariable(SimGetVariable, (void*) &sim);
                VisItSetGetDomainList(SimGetDomainList, (void*) &sim);
            }
            else
                fprintf(stderr, "VisIt did not connect\n");
        }
        else if (visitstate == 2)
        {
            /* VisIt wants to tell the engine something. */
            sim.runMode = VISIT_SIMMODE_STOPPED;
#ifdef PARALLEL
            if (!ProcessVisItCommand())
#else
            if (!VisItProcessEngineCommand())
#endif
            {
				cout<<"VisItDisconnect.";cout.flush();
                /* Disconnect on an error or closed Projection. */
                VisItDisconnect();
                /* Start running again if VisIt closes. */
                sim.runMode = VISIT_SIMMODE_RUNNING;
            }
			//else sim.runMode = VISIT_SIMMODE_RUNNING;
        }
    }
    sim.done = true;

	cout<<"\nError: "<<err;cout.flush();
    // store data afterwards for this case
    this->RecordAll();
    this->StoreTimings();
}

visit_handle SimGetMetaData(void *cbdata)
{
    visit_handle md = VISIT_INVALID_HANDLE;
    simulation_data *sim = (simulation_data *) cbdata;

    /* Create metadata. */
    if (VisIt_SimulationMetaData_alloc(&md) == VISIT_OKAY)
    {
        int i;
        visit_handle m1 = VISIT_INVALID_HANDLE;
        visit_handle vmd = VISIT_INVALID_HANDLE;
        visit_handle cmd = VISIT_INVALID_HANDLE;

        /* Set the simulation state. */
        VisIt_SimulationMetaData_setMode(md, (sim->runMode == SIM_STOPPED) ?
                                         VISIT_SIMMODE_STOPPED : VISIT_SIMMODE_RUNNING);
        VisIt_SimulationMetaData_setCycleTime(md, sim->cycle, sim->time);

        /* Set the first mesh's properties.*/
        if (VisIt_MeshMetaData_alloc(&m1) == VISIT_OKAY)
        {
            /* Set the mesh's properties.*/
            VisIt_MeshMetaData_setName(m1, "Mesh2d Layer 1");
            VisIt_MeshMetaData_setMeshType(m1, VISIT_MESHTYPE_RECTILINEAR);
            VisIt_MeshMetaData_setTopologicalDimension(m1, 2);
            VisIt_MeshMetaData_setSpatialDimension(m1, 2);
            VisIt_MeshMetaData_setNumDomains(m1, sim->par_size);
            VisIt_MeshMetaData_setXLabel(m1, "RateUnits");
            VisIt_MeshMetaData_setYLabel(m1, "Hypercolumns");

            VisIt_SimulationMetaData_addMesh(md, m1);
        }

        /* Add a zonal scalar variable on mesh2d. */
        if (VisIt_VariableMetaData_alloc(&vmd) == VISIT_OKAY)
        {
            VisIt_VariableMetaData_setName(vmd, "value");
            VisIt_VariableMetaData_setMeshName(vmd, "Mesh2d Layer 1");
            VisIt_VariableMetaData_setType(vmd, VISIT_VARTYPE_SCALAR);
            VisIt_VariableMetaData_setCentering(vmd, VISIT_VARCENTERING_ZONE);
            VisIt_SimulationMetaData_addVariable(md, vmd);
        }

        /* Add some custom commands. */
        for (int i = 0; i < sizeof (cmd_names) / sizeof (const char *); ++i)
        {
            visit_handle cmd = VISIT_INVALID_HANDLE;
            if (VisIt_CommandMetaData_alloc(&cmd) == VISIT_OKAY)
            {
                VisIt_CommandMetaData_setName(cmd, cmd_names[i]);
                VisIt_SimulationMetaData_addGenericCommand(md, cmd);
            }
        }
    }

    return md;
}

visit_handle NetworkDemoVis1::SimGetMesh(int domain, const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    simulation_data *sim = (simulation_data *) cbdata;

    if (strcmp(name, "Mesh2d Layer 1") == 0)
    {
        if (VisIt_RectilinearMesh_alloc(&h) != VISIT_ERROR)
        {
            int minRealIndex[3], maxRealIndex[3];
            minRealIndex[0] = minRealIndex[1] = minRealIndex[2] = 0;

            maxRealIndex[0] = sim->rMeshDims[0] - 1;
            maxRealIndex[1] = sim->rMeshDims[1] - 1;
            maxRealIndex[2] = 0;
            visit_handle hxc, hyc;
            VisIt_VariableData_alloc(&hxc);
            VisIt_VariableData_alloc(&hyc);
            VisIt_VariableData_setDataF(hxc, VISIT_OWNER_SIM, 1, sim->rMeshDims[0], sim->rmeshX);
            VisIt_VariableData_setDataF(hyc, VISIT_OWNER_SIM, 1, sim->rMeshDims[1], sim->rmeshY);
            VisIt_RectilinearMesh_setCoordsXY(h, hxc, hyc);

            VisIt_RectilinearMesh_setRealIndices(h, minRealIndex, maxRealIndex);
        }
    }

    return h;
}

void NetworkDemoVis1::simulate_one_timestep()
{
    ++sim.cycle;
    sim.time += 1;
    // inject some random data every nth iteration
    if (step % 10 == 0)
    {
        vector<float> randData(minicolumns.size());
        for (int m = 0; m < randData.size(); m++) randData[m] = rand() / (float(RAND_MAX) + 1);
        this->GetLayer(0)->SetValuesAll(randData, true);
    }

    this->Simulate(); // 1 timestep

    // example of how to retrieve current activity (local for this process)
    //vector<float> values(minicolumns.size());
    //cout << sim.par_rank << ": ";cout << i << ": ";
    for (int j = 0; j < minicolumns.size(); j++)
        sim.values[j] = (minicolumns[j]->GetValue());
//    VisItTimeStepChanged();
//    VisItUpdatePlots();
    step++;
}

visit_handle NetworkDemoVis1::SimGetVariable(int domain, const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    int nComponents = 1, nTuples = 0;
    simulation_data *sim = (simulation_data *) cbdata;

    if (VisIt_VariableData_alloc(&h) == VISIT_OKAY)
    {
        if (strcmp(name, "value") == 0)
        {
            nTuples = (sim->rMeshDims[0] - 1) * (sim->rMeshDims[1] - 1);
            VisIt_VariableData_setDataF(h, VISIT_OWNER_SIM, nComponents, nTuples, sim->values);
        }
        else
        {
            VisIt_VariableData_free(h);
            h = VISIT_INVALID_HANDLE;
        }
    }
	cout<<"SimGetVariable ("<<sim->par_rank<<")";cout.flush();
    return h;
}

/*
visit_handle NetworkDemoVis1::SimGetDomainList(const char *name, void *cbdata)
{
    visit_handle h = VISIT_INVALID_HANDLE;
    if (VisIt_DomainList_alloc(&h) != VISIT_ERROR)
    {
        visit_handle hdl;
        int i, *iptr = NULL;
        simulation_data *sim = (simulation_data *) cbdata;

        iptr = (int *) malloc(sizeof (int));
        *iptr = sim->par_rank;

        //if (VisIt_VariableData_alloc(&hdl) == VISIT_OKAY)
        //{
			VisIt_VariableData_alloc(&hdl);
            VisIt_VariableData_setDataI(hdl, VISIT_OWNER_VISIT, 1, 1, iptr);
            VisIt_DomainList_setDomains(h, sim->par_size, hdl);
        //}
			cout<<"SimGetDomainList ("<<sim->par_rank<<")";cout.flush();
    }
    return h;
}
*/

visit_handle NetworkDemoVis1::SimGetDomainList(const char *name, void *cbdata)
{
	int par_size, par_rank;
	visit_handle h = VISIT_INVALID_HANDLE;
	MPI_Comm_rank(NETWORK_COMM_WORLD, &par_rank);
	MPI_Comm_size(NETWORK_COMM_WORLD, &par_size);
	if(VisIt_DomainList_alloc(&h) != VISIT_ERROR)
	{
		visit_handle hdl;
		VisIt_VariableData_alloc(&hdl);
		VisIt_VariableData_setDataI(hdl, VISIT_OWNER_COPY, 1, 1, &par_rank);
		VisIt_DomainList_setDomains(h, par_size, hdl);
	}
	return h;
}

NetworkDemoVis1::~NetworkDemoVis1()
{
    if (sim.rmeshX) delete [] sim.rmeshX;
    if (sim.rmeshY) delete [] sim.rmeshY;
    if (sim.values) delete [] sim.values;
}

int NetworkDemoVis1::visit_broadcast_int_callback(int *value, int sender)
{
    return MPI_Bcast(value, 1, MPI_INT, sender, NETWORK_COMM_WORLD);
}

int NetworkDemoVis1::visit_broadcast_string_callback(char *str, int len, int sender)
{
    return MPI_Bcast(str, len, MPI_CHAR, sender, NETWORK_COMM_WORLD);
}

void NetworkDemoVis1::BroadcastSlaveCommand(int *command)
{
#ifdef PARALLEL
    MPI_Bcast(command, 1, MPI_INT, 0, NETWORK_COMM_WORLD);
#endif
}

/* Callback involved in command communication. */
void NetworkDemoVis1::SlaveProcessCallback()
{
    int command = VISIT_COMMAND_PROCESS;
    BroadcastSlaveCommand(&command);
}

#ifdef PARALLEL
int NetworkDemoVis1::ProcessVisItCommand()
{
    int command;
    if (sim.par_rank == 0)
    {
		cout<<"Command = ";cout.flush();
        int success = VisItProcessEngineCommand();
		cout<<success<<" ";cout.flush();
        if (success)
        {
            command = VISIT_COMMAND_SUCCESS;
            BroadcastSlaveCommand(&command);
            return 1;
        }
        else
        {
            command = VISIT_COMMAND_FAILURE;
            BroadcastSlaveCommand(&command);
            return 0;
        }
    }
    else
    {
        /* Note: only through the SlaveProcessCallback callback
         * above can the rank 0 process send a VISIT_COMMAND_PROCESS
         * instruction to the non-rank 0 processes. */
        while (1)
        {
			
            BroadcastSlaveCommand(&command);
			cout<<"Command (slave) = "<<command;cout.flush();
			switch (command)
            {
            case VISIT_COMMAND_PROCESS:
                VisItProcessEngineCommand();
                break;
            case VISIT_COMMAND_SUCCESS:
                return 1;
            case VISIT_COMMAND_FAILURE:
				cout<<"fail.";cout.flush();
                return 0;
            }
        }
    }
    return 1;
}
#endif
#endif