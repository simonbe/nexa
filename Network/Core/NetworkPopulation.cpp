#include <math.h>
#include "Network.h"
#include "NetworkPopulation.h"
#include "MPIDistribution.h"
#include "NetworkUnitsIF.h"


Population::Population() 
{
	m_mpiDistribution = new MPIDistribution(this);
	m_keepValues = false;//false;
	m_parallelizationScheme = MPIDistribution::ParallelizationDefault;
}

void Population::ResetLocalities()
{
	TimingStart("ResetLocalities");

	for(int i=0;i<m_units.size();i++)
	{
		m_units[i]->SetLocal(true);	
		m_units[i]->IsNewEvent(true);
	}

	TimingStop("ResetLocalities");
}

Population::~Population()
{
	for(int i=0;i<m_units.size();i++)
		delete m_units[i];

	for(int i=0;i<m_projectionIncoming.size();i++)
	{
		if(m_projectionIncoming[i] != 0)
		{
			delete m_projectionIncoming[i];
			m_projectionIncoming[i] = NULL;
		}
	}

	for(int i=0;i<m_preConnectivitys.size();i++)
	{
		if(m_preConnectivitys[i] != NULL)
		{
			delete m_preConnectivitys[i];
			m_preConnectivitys[i] = NULL;
		}
	}
}

/// <summary>	Calls all projection modifiers for the incoming projections of this population. 
/// 			Also call population modifiers' Modify function but this is rarely overridden for any population modifier
/// 			as they update in their Simulate functions instead.</summary>

void Population::Modify()
{
	// projection modifiers, e.g. synaptic plasticity

	for(int i=0;i<m_projectionIncoming.size();i++)
	{
		m_projectionIncoming[i]->ModifyProjection();
	}

	// typically seldom utilized, currently adaptation and trace classes use it
	// but that should be possible to put in their Simulate functions instead so consider removing this
	for(int i=0;i<m_populationModifiers.size();i++)
		m_populationModifiers[i]->Modify(); 
}

/// <summary>	Resets the localities. Not used.
/// 			TODO: Check dependencies and remove. </summary>

void PopulationColumns::ResetLocalities()
{
	TimingStart("ResetLocalities");
	// minicolumns

	for(int i=0;i<m_minicolumns.size();i++)
	{
		m_minicolumns[i]->SetLocal(true);
		m_minicolumns[i]->IsNewEvent(true); // needed?
	}

	// hypercolumns - all by default local now
	for(int i=0;i<m_hypercolumns.size();i++)
	{
		m_hypercolumns[i]->SetLocal(true);//false);
		m_hypercolumns[i]->IsNewEvent(true); // ok?
	}

	TimingStop("ResetLocalities");
}

/// <summary>	Adds a units modifier for all incoming Projections for layer </summary>
///
/// <param name="p">	UnitModifier to add to population. </param>

void Population::AddUnitsModifierToInitialize(UnitModifier* p)
{
	if(m_units.size()>0)	// called after initialize
	{
		for(int i=0;i<m_units.size();i++)
		{
			m_units[i]->AddUnitModifier(p);
		}
	}
	else					// called before initialize
	{
		m_unitPropertiesToInitialize.push_back(p);
	}
}

/// <summary>	Constructor for columnar population (parts may get moved into Population::Population) </summary>
///
/// <param name="net">						  	network (this). </param>
/// <param name="nrHypercolumns">			  	Number of hypercolumns/macrocolumns - set to 1 to get a non-columnar population. </param>
/// <param name="nrRateUnits">				  	Nr of (by default rate) units in each hypercolumn. </param>
/// <param name="unitType">					  	Neural unit type. Default RateUnit. </param>
/// <param name="parallelizationScheme">	  	The parallelization scheme. </param>
/// <param name="useSilentHypercolumns">	  	True to allow a hypercolumn to have no activity (default false, TODO: move away from constructor and into a population modifier) </param>
/// <param name="silentHypercolumnsThreshold">	The silent hypercolumns threshold (default not used). </param>

PopulationColumns::PopulationColumns(Network* net, unsigned long nrHypercolumns, unsigned long nrRateUnits, UnitType unitType, MPIDistribution::ParallelizationSchemeLayer parallelizationScheme, bool useSilentHypercolumns, float silentHypercolumnsThreshold)
{
	m_network = net;
	m_keepValues = false;//false;

	m_unitType = unitType;
	m_parallelizationScheme = parallelizationScheme;
	m_useSilentHypercolumns = useSilentHypercolumns;
	m_silentHypercolumnsThreshold = silentHypercolumnsThreshold;
	m_firstRun = true;
	m_initialized = false;
	m_nrUnits = 0;
	m_nrRateUnits.clear();
	m_localHypercolumnIndexes.clear();
	m_localRateUnits.clear();
	m_mcHcIndexes.clear();
	m_nodeHcIndexes.clear();
	m_nodeHcIds.clear();

	if(parallelizationScheme == MPIDistribution::ParallelizationDefault || parallelizationScheme == MPIDistribution::ParallelizationSplitPopulations)
	{
		m_network->MPIAddNrUnitsParallelizationDefault(nrHypercolumns*nrRateUnits);
	}

	m_nrHypercolumns = nrHypercolumns;
	m_nrRateUnitsInHypercolumn = nrRateUnits;	
}

/// <summary>	Main initialization for columnar population (may get moved into Population::Initialize 
/// 			</summary>
///
///  TODO: Guarantee to always scale properly by changing loop around m_nrRateUnits.size()
///

void PopulationColumns::Initialize()
{
	///////// Prevously put in constructor

	double timeStart;
	if(this->network()->MPIGetNodeId() == 0)
	{
		timeStart = MPI_Wtime();
		cout<<"Setting up object PopulationColumns...";
		cout.flush();
	}

	vector<int> mcsInHc;

	for(int i=0;i<m_nrHypercolumns;i++)
	{
		m_nrRateUnits.push_back(m_nrRateUnitsInHypercolumn);
		m_nrUnits+=m_nrRateUnitsInHypercolumn;
		mcsInHc.push_back(m_nrRateUnitsInHypercolumn);
	}

	m_populationType = Population::FixedWeightsAndRateUnits;
	TimingStart("LayerInitDiv");
	vector<vector<vector<long> > > div;

	if(this->m_parallelizationScheme == MPIDistribution::ParallelizationDefault || this->m_parallelizationScheme == MPIDistribution::ParallelizationSplitPopulations)
	{
		// if not more processes than total nr populations, divide by total units
		//if(network()->GetNrPopulations()<=network()->MPIGetNrProcs())
		//	div = m_mpiDistribution->DivideEqualByTotalUnits(mcsInHc,network()->MPIGetCurrentUnitParallelizationDefault(),network()->MPIGetNrUnitsParallelizationDefault(),network()->MPIGetNodeId(),network()->MPIGetNrProcs());
		//else
			div = m_mpiDistribution->DivideEqualByPopulations(mcsInHc,network()->MPIGetCurrentUnitParallelizationDefault(),network()->MPIGetNrUnitsParallelizationDefault(),network()->MPIGetNodeId(),network()->MPIGetNrProcs());

		// if more processes than total nr populations, divide by populations
		network()->MPIAddCurrentUnitParallelizationDefault(m_nrUnits);

		vector<int> allUseProcess(network()->MPIGetNrProcs());
		int useProcess = 1;
		if( div[1][0][0] == 0 &&  div[1][0][1] == 0)
			useProcess = 0;

		MPI_Allgather(&useProcess,1,MPI_INT,&allUseProcess[0],1,MPI_INT,NETWORK_COMM_WORLD);

		for(int i=0;i<allUseProcess.size();i++)
		{
			if(allUseProcess[i] == 1)
				m_mpiProcessesUsed.push_back(i);
		}
	}
	else // DivideByHypercolumns
		div = m_mpiDistribution->DivideEqualByHypercolumn(mcsInHc,m_parallelizationScheme,network()->MPIGetNodeId(),network()->MPIGetNrProcs());//m_nrUnits,network()->MPIGetNodeId(),network()->MPIGetNrProcs());//(m_nrRateUnits,m_network->MPIGetNodeId(),m_network->MPIGetNrProcs());
	

	TimingStop("LayerInitDiv");
	
	m_nodeHcIndexes = div[0];
	m_localUnitIndexesInterval = div[1][0];
	div.clear();
	mcsInHc.clear();

	//m_hashMcHcIndexes

	int currentTotUnit = 0;
	vector<int> currentNodeIndexesLayer;
	int lastMaxUsedInterval = 0;

	currentTotUnit = m_localUnitIndexesInterval[0];

	TimingStart("LayerInitHcIndexes");
	for(int i=0;i<m_nodeHcIndexes.size();i++)
	{
		for(int j=0;j<m_nodeHcIndexes[i].size();j++)
		{
			bool exist = false;
			for(int k=0;k<m_nodeLayerIndexes.size();k++)
			{
				if(m_nodeLayerIndexes[k] == m_nodeHcIndexes[i][j])
				{
					exist = true;
					break;
				}
			}

			if(exist == false)
				m_nodeLayerIndexes.push_back(m_nodeHcIndexes[i][j]);
		}
	}

	TimingStop("LayerInitHcIndexes");

	TimingStart("LayerInitLocalUnitIndexes");
	for(int i=m_localUnitIndexesInterval[0];i<m_localUnitIndexesInterval[1];i++)
	{
		int cHc = i/m_nrRateUnitsInHypercolumn;

		if(find(m_localHypercolumnIndexes.begin(),m_localHypercolumnIndexes.end(),cHc) == m_localHypercolumnIndexes.end())
			m_localHypercolumnIndexes.push_back(cHc);
	}
	
	sort(m_localHypercolumnIndexes.begin(), m_localHypercolumnIndexes.end()); // not needed with current default division
	TimingStop("LayerInitLocalUnitIndexes");
	
	m_unitIndexesInterval.clear(); // fix, may change logic

	if(m_network->MPIGetNodeId() == 0)
	{
		cout<<"\n(process 0) LocalHcIndex: "<<m_localHypercolumnIndexes.size()<<", localHcIndexes: "<<m_localHypercolumnIndexes.size() <<", nodeLayerIndexes: "<<m_nodeLayerIndexes.size()<<".\n";
		cout.flush();
	}
	

	if(this->network()->MPIGetNodeId() == 0)
	{
		cout<<"ok. (Time (process 0): "<<MPI_Wtime()-timeStart<<")\n";
		cout.flush();
	}

	/// end previously in constructor
	////////////////////////////////////////////////////

	//// Initialize the individual units

	if(this->network()->MPIGetNodeId() == 0)
	{
		cout<<"Initializing local hypercolumns, ";
		cout.flush();
	}

	if(m_initialized == false || m_initialized == true) // re-initializing (as a result of added populations etc) should only initialize certain parts, so may need to change
	{
		long localId = 0;

		TimingStart("HypercolumnsInit");

		for(long l=0;l<m_nrRateUnits.size();l++)
		{
			bool makeLocal = false;

			for(int i=0;i<m_nrRateUnits[l];i++)
			{
				int lId = localId+i;
				if(lId >= m_localUnitIndexesInterval[0] && lId < m_localUnitIndexesInterval[1])
				{
					makeLocal = true;
					break;
				}
			}

			long unitId = network()->GetNextUnitId();
			if(l==0)
				m_unitsStartId = unitId + m_nrRateUnits.size(); // assumes the startid starts first with the minicolumns and there is a strict ordering

			int hcId = network()->GetNextHypercolumnId(); // not used anymore?

			if(makeLocal == true)
			{
				Hypercolumn* h;
				
				if(this->m_useSilentHypercolumns == true)
					h = new Hypercolumn(this->m_silentHypercolumnsThreshold);
				else
					h = new Hypercolumn();

				
				h->SetUnitId(unitId);//hcId);
				h->SetTotalNrRateUnits(m_nrRateUnits[l]);

				//h->SetNodeId(this->MPIBelongsToNode(localId)); // will be set to same node as the first minicolumn belongs to
				h->SetUnitIdLocal(l); // will be ordered first - may be problematic in vector (local) index based lookups
				h->network(this->network());
				
				this->AddUnit(h);

				localId+=m_nrRateUnits[l];
			}
			else
			{	
				localId+=m_nrRateUnits[l];
			}
		}

		TimingStop("HypercolumnsInit");
		
		localId = 0;

		if(this->network()->MPIGetNodeId() == 0)
		{
			cout<<"minicolumns...";
			cout.flush();
		}

		int currentHc = 0;
			
		TimingStart("RateUnitsInit");

		if(this->network()->MPIGetNodeId() == 0)
		{
			cout<<"Interval process 0 = ["<<m_localUnitIndexesInterval[0]<<", "<<m_localUnitIndexesInterval[1]<<"] ";
			cout.flush();
		}
		
		for(long l=0;l<m_nrRateUnits.size();l++)
		{
			bool hcUsed = false;

			for(int i=0;i<m_nrRateUnits[l];i++)
			{
				int lId = localId+i;
				if(lId >= m_localUnitIndexesInterval[0] && lId < m_localUnitIndexesInterval[1])
				{
					hcUsed = true;
					break;
				}
			}

			for(int i=0;i<m_nrRateUnits[l];i++)
			{
				bool makeLocal = false;

				if(localId >= m_localUnitIndexesInterval[0] && localId < m_localUnitIndexesInterval[1])
					makeLocal = true;

				long id = network()->GetNextUnitId();

				if(l==0 && i==0)
					m_startRateUnitId = id;

				if(makeLocal == true)
				{
					Hypercolumn* h = this->GetHypercolumns()[currentHc];

					if(l==0 && i==0)
						m_startRateUnitId = id;

					// TODO: reduce to T_ etc

					if(m_unitType == PopulationColumns::Graded || m_unitType == PopulationColumns::GradedThresholded)
					{
						RateUnit* m;
						if(m_unitType == PopulationColumns::Graded)
							m = new RateUnit(false);
						else if(m_unitType == PopulationColumns::GradedThresholded)
							m = new RateUnit(true);
						
						m->SetUnitId(id);
						m->SetUnitIdLocal(localId);
						m->SetUnitIdLocalHypercolumn(i);
						m->SetHypercolumnId(h->GetUnitId());//hcId);
						this->SetUnitIdLocalId(id,localId);
						
						m->SetPopulation(this);
						m->network(this->network());

						if(m_useTrace == true)
						{
							m->SetTrace(true,m_traceStrength);
						}

						this->AddUnit(m);
						h->AddRateUnit(m);
						m->AddHypercolumn(h);
					}
					else if(m_unitType == PopulationColumns::adEIF)
					{
#if GSL_AVAILABLE==1
						UnitadEIF* m = new UnitadEIF();
						m->SetUnitId(id);//m_network->GetNextUnitId());
						m->SetUnitIdLocal(localId);
						m->SetUnitIdLocalHypercolumn(i);
						m->SetHypercolumnId(h->GetUnitId());
						//m->AddUnitModifier(t);
						this->SetUnitIdLocalId(id,localId);

						m->SetNodeId(this->MPIBelongsToNode(localId));
						m->Population(this);
						m->network(this->network());

						this->AddUnit(m);
						h->AddRateUnit(m);
						m->AddHypercolumn(h);
#endif
					}
					else if(m_unitType == PopulationColumns::IF)
					{
						UnitIF* m = new UnitIF();
						m->SetUnitId(id);//network()->GetNextUnitId());
						m->SetUnitIdLocal(localId);
						m->SetUnitIdLocalHypercolumn(i);
						m->SetHypercolumnId(h->GetUnitId());//hcId);
						this->SetUnitIdLocalId(id,localId);
						//m->AddUnitModifier(t);

						m->SetNodeId(this->MPIBelongsToNode(localId));
						m->SetPopulation(this);
						m->network(this->network());

						this->AddUnit(m);
						h->AddRateUnit(m);
						m->AddHypercolumn(h);
					}
				}

				localId++;
			}

			if(hcUsed == true)
				currentHc++;
		}

		TimingStop("RateUnitsInit");
		
		// could put in unit properties here

		ResetLocalities();

		for(int i=0;i<m_minicolumns.size();i++)
		{
			// global/non-Projection-specific unit properties/transfer functions
			for(int j=0;j<m_unitPropertiesToInitialize.size();j++)
				m_minicolumns[i]->AddUnitModifier(m_unitPropertiesToInitialize[j]);
			
			if(m_minicolumns[i]->IsLocal())
				m_localRateUnits.push_back(m_minicolumns[i]);
		}
	}

	if(this->network()->MPIGetNodeId() == 0)
	{
		cout<<"ok ("<<this->GetRateUnits().size()<<" nr local mcs, "<<this->GetHypercolumns().size()<<" nr local hcs, "<<this->GetUnitsAll()->size()<<" tot units). ";
		cout.flush();
	}
}

/// <summary>	Additional initializations not included in standard Initialize. </summary>

void PopulationColumns::InitializeProjectionsEventsAndParameters()
{
	
	for(int i=0;i<m_pres.size();i++)
	{
		m_preConnectivitys[i]->network(this->network());
		m_preConnectivitys[i]->Initialize();
	}

	for(int i=0;i<m_populationModifiers.size();i++)
	{
		m_populationModifiers[i]->network(this->network());
		m_populationModifiers[i]->Initialize(this);
	}

	m_initialized = true;
}

/// <summary>	Main population simulate (may get moved into Population::Simulate) </summary>

void PopulationColumns::Simulate()
{
	TimingStart(m_name);

	if(IsOn())
	{
		TimingStart("SimulateUnits");

		int count = 0;

		for(int i=0;i<m_minicolumns.size();i++)
		{
			m_minicolumns[i]->Simulate();
			count++;
		}

		for(int i=0;i<m_hypercolumns.size();i++)
		{
			m_hypercolumns[i]->Simulate();
			count++;
		}

		TimingStop("SimulateUnits");
	}
	else
	{
		for(int i=0;i<m_minicolumns.size();i++)
			m_minicolumns[i]->ClearEventsIncoming();
	}

	// Note: PopulationModifiers are run even if population is switched off (is put on top of population, not "inside" a population)
	TimingStart("PopulationModifiers");
	for(int i=0;i<m_populationModifiers.size();i++)
	{
		if(m_populationModifiers[i]->IsOn())
			m_populationModifiers[i]->Simulate();
	}
	TimingStop("PopulationModifiers");

	for(int j=0;j<m_minicolumns.size();j++)
	{
		if(!(fabs(m_minicolumns[j]->GetValue()) < EPS))
		{
			m_minicolumns[j]->IsNewEvent(true); // currently removes prev existing events (TODO: check if always ok.)
			m_minicolumns[j]->AddEventOutgoing();
		}
		else
			m_minicolumns[j]->IsNewEvent(false);
	}

	if(IsRecording())
	{
		TimingStart("RecordingLayer");

/*		if(m_meterLayer->GetSaveType() == Storage::MPI_Binary)//m_network->GetFilePreference() == Network::MPI_Binary)
		{
			vector<float> values = vector<float>(m_minicolumns.size());
			for(int j=0;j<m_minicolumns.size();j++)
			{
				values[j] = m_minicolumns[j]->GetValue();
			}

			m_recordedValues.push_back(values);
		}
		else
		{*/
			vector<float> values = vector<float>(m_minicolumns.size());
			for(int j=0;j<m_minicolumns.size();j++)
			{
				values[j] = m_minicolumns[j]->GetValue();
			}

			m_recordedValues.push_back(values);
		//}
		
		TimingStop("RecordingLayer");
	}

	// for graded, go through all incoming Projections and retrieve events
	// for spiking, could allow any event to come during this timestep

	TimingStop(m_name);
}

void Population::ClearEventsIncoming()
{
	// assumes even distribution
	for(int i=m_localUnitIndexesInterval[0];i< m_localUnitIndexesInterval[1]; i++)
	{
		m_units[i]->ClearEventsIncoming();
	}
}

void PopulationColumns::ClearEventsIncoming()
{
	for(int i=0;i<m_minicolumns.size();i++)
	{
		m_minicolumns[i]->ClearEventsIncoming();
	}
}

void PopulationColumns::SetValue(int localIndex, float value)
{
	for(int i=0;i<m_minicolumns.size();i++)
	{
		if(m_minicolumns[i]->GetUnitIdLocal() == localIndex)
		{
			m_minicolumns[i]->SetValue(value);
			break;
		}
	}
}

void PopulationColumns::SetValuesAll(vector<float> values, bool alsoSimulate)
{
	if(m_minicolumns.size() == values.size())
	{
		for(int i=0;i<m_minicolumns.size();i++)
		{
			m_minicolumns[i]->SetValue(values[i]);
			if(alsoSimulate == false)
				((RateUnit*)m_minicolumns[i])->SetNoUpdatingCurrentTimeStep(true); // test this
		}
	}
	else
	{
		int index = 0;
		for(int i=0;i<m_minicolumns.size();i++)//values.size();i++)
		{
			//		if(index==values.size())
			//			index-=values.size(); // repeat pattern

			if(m_minicolumns[index]->GetUnitIdLocal()<values.size())
			{
				((RateUnit*)m_minicolumns[i])->SetValue(values[m_minicolumns[index]->GetUnitIdLocal()]);//index]);
			}
			
			if(alsoSimulate == false)
				((RateUnit*)m_minicolumns[i])->SetNoUpdatingCurrentTimeStep(true); // test this

			index++;
		}
	}

	if(alsoSimulate == true)
		this->m_keepValues = true;
	else
		this->m_keepValues = false;

	// clear also all incoming values
}

vector<float> Population::GetValuesLocal()
{
	vector<float> f(m_units.size());

	for(int i=0;i<m_units.size();i++)
	{
		f[i] = m_units[i]->GetValue();
	}

	return f;
}

vector<float> PopulationColumns::GetValuesLocal()
{
	vector<float> f(m_minicolumns.size());

	for(int i=0;i<m_minicolumns.size();i++)
	{
		f[i] = ((RateUnit*)m_minicolumns[i])->GetValue();
	}

	return f;
}


void PopulationColumns::Reset()
{
	for(int i=0;i<m_minicolumns.size();i++)
	{
		((RateUnit*)m_minicolumns[i])->SetValue(0.0);
	}
}

void PopulationColumns::SetValuesLocal(vector<float> values)
{
	//	if(values.size()!=m_units.size()) // should always be the same size as values
	//		Logger::print("Warning: values not same size as units.");

	long start = m_localUnitIndexesInterval[0]; // relative
	long end = m_localUnitIndexesInterval[1];

	if(values.size()>end-start) // absolute
	{
		start = 0;
		end = values.size();
	}

	for(long i=start;i<end;i++) 
	{
		((RateUnit*)m_minicolumns[i])->SetValue(values[i-start]);
	}
}


/// <summary>	Connects this population to another population. </summary>
///
/// <param name="layer">				Other (pre-)population. </param>
/// <param name="ProjectionType">   	Connectivity used. </param>
/// <param name="calledFromAddPost">	True if called from AddPre within class (default false). </param>

void Population::AddPre(Population* layer, Connectivity* ProjectionType,  bool calledFromAddPost)
{
	m_pres.push_back(layer);
	m_preConnectivitys.push_back(ProjectionType);
	ProjectionType->SetPreAndPost(layer,this);
	ProjectionType->network(network());

	if(calledFromAddPost == false)
		layer->AddPost(this, ProjectionType, true);	
}

/// <summary>	Connects this population to another population. </summary>
///
/// <param name="layer">		   	Other (post-)population. </param>
/// <param name="ProjectionType">  	Connectivity used. </param>
/// <param name="calledFromAddPre">	True if called from AddPost within class (default false). </param>

void Population::AddPost(Population* layer, Connectivity* ProjectionType,  bool calledFromAddPre)
{
	m_posts.push_back(layer);
	m_postConnectivitys.push_back(ProjectionType);
	ProjectionType->SetPreAndPost(this,layer);
	ProjectionType->network(network());

	if(calledFromAddPre == false)
		layer->AddPre(this,ProjectionType,true);
}

vector<Unit*> Population::GetLocalUnits()
{
	vector<Unit*> localUnits;

	for(int i=0;i<m_units.size();i++)
		if(m_units[i]->IsLocal())
			localUnits.push_back(m_units[i]);
	return localUnits;
}

/// <summary>	Gets local indexes. 
/// 			TODO: Change name.</summary>
///
/// <returns>	Local indexes of units on this process. (now should be all m_units) </returns>
///

vector<int> PopulationColumns::GetMPIDistribution()
{
	vector<int> unitIndexes(m_minicolumns.size());

	for(int i=0;i<m_minicolumns.size();i++)
	{
		unitIndexes[i] = m_minicolumns[i]->GetUnitIdLocal();
	}

	return unitIndexes;
}

/// <summary>	Get ids of hypercolumns that are local on process. </summary>
///
/// <param name="processId">	Process id to get indexes for. </param>
///
/// <returns>	Hypercolumn ids/processes. </returns>

vector<int> PopulationColumns::GetMPIDistributionHypercolumns(int processId)
{
	vector<int> unitIndexes;

	for(int i=0;i<m_hypercolumns.size();i++)
	{
		if(m_hypercolumns[i]->GetNodeId() == processId)
			unitIndexes.push_back(i);
	}

	return unitIndexes;
}

/// <summary>	Adds unit to population. </summary>
///
/// <param name="unit">	[in,out] If non-null, the unit. </param>

void PopulationColumns::AddUnit(Unit* unit)
{
	m_units.push_back(unit);
	network()->AddUnit(unit);

	if(unit->GetType().compare("minicolumn")==0)
	{
		m_minicolumns.push_back((RateUnit*)unit);
	}
	else if(unit->GetType().compare("hypercolumn")==0)
	{
		m_hypercolumns.push_back((Hypercolumn*)unit);
	}
}

/// <summary>	Adds a population modifier (e.g. WTA) to population. </summary>
///
/// <param name="eventLayer">	The population modifier (e.g. WTA). </param>

void Population::AddPopulationModifier(PopulationModifier* eventLayer)
{
	eventLayer->network(this->network()); // move?
	eventLayer->Initialize(this);
	m_populationModifiers.push_back(eventLayer);
	m_network->AddNetworkObjectToDelete(eventLayer); // to avoid duplicate deletion - network responsible for the deletion of this object
}

void Population::Dispose()
{
	DisposeParent();
}

void Population::DisposeParent()
{
	for(int i=0;i<m_units.size();i++)
	{
		m_units[i]->Dispose(); // can put in destructor
		delete m_units[i];
	}

	for(int i=0;i<m_projectionOutgoing.size();i++) // contains hash weights, will be rebuilt
	{
		m_projectionOutgoing[i]->Dispose();
		delete m_projectionOutgoing[i];
	}
	
	m_projectionOutgoing.clear();
	m_projectionIncoming.clear();
	m_units.clear();

	for(int i=0;i<m_preConnectivitys.size();i++)
	{
		m_preConnectivitys[i]->SetInitialized(false);
	}

	m_nrUnits = 0;
}

void PopulationColumns::Dispose()
{
	DisposeParent();
	m_minicolumns.clear();
	m_hypercolumns.clear();
	m_localRateUnits.clear();
}

vector<int> PopulationColumns::GetUnitIdLocals()
{
	vector<int> ids(m_minicolumns.size());
	for(int i=0;i<m_minicolumns.size();i++)
	{
		ids[i] = ((RateUnit*)m_minicolumns[i])->GetUnitIdLocal();
	}

	return ids;
}

/// <summary>	Sets hypercolumns in population to use silent hypercolumns,
/// 			i.e. forces a hypercolumn to have no activity in its units if none over a threshold. </summary>
///
/// <param name="useSilent">	true to use silent hypercolumns. </param>
/// <param name="threshold">	The threshold determining if a hypercolumn should be silent (if no activity above). </param>

void PopulationColumns::SetSilentHypercolumns(bool useSilent, float threshold)
{
	m_useSilentHypercolumns = useSilent;
	m_silentHypercolumnsThreshold = threshold;
	for(int i=0;i<m_hypercolumns.size();i++)
	{
		((Hypercolumn*)m_hypercolumns[i])->SetUseSilent(useSilent);
		((Hypercolumn*)m_hypercolumns[i])->SetSilentThreshold(threshold);
	}
}