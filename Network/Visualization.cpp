#include "Visualization.h"

#if VISIT_AVAILABLE == 1

const char *cmd_names[] = {"Run", "Stop", "Step", "Update", "Reset"};

simulation_data::simulation_data()
{
    cycle = 0;
    time = 0.;
    runMode = 1; /* VISIT_SIMMODE_RUNNING */
    done = 0;
#ifdef PARALLEL
    par_rank = 0;
    par_size = 1;
#endif
}

simulation_data::~simulation_data()
{
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



#endif