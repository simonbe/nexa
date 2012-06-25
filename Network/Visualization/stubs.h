/* 
 * File:   stubs.h
 * Author: markosp
 *
 * Created on December 20, 2011, 2:16 PM
 */

#define PARALLEL

#ifndef STUBS_H
#define	STUBS_H

class NetworkDemoVis1;//_prace;

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
    //NetworkDemoVis1_prace *network;
	NetworkDemoVis1 *network;
    
#ifdef PARALLEL
    int     par_rank;
    int     par_size;
#endif
    simulation_data();
    ~simulation_data();
};

void simulation_data_ctor(simulation_data *sim);

void simulation_data_dtor(simulation_data *sim);

#endif	/* STUBS_H */

