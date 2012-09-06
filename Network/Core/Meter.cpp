#include <sstream>
#include <stdio.h>
#include <string.h>
#include "Meter.h"
#include "Storage.h"


Meter::Meter(char* filename, Storage::FilePreference fileType)
{
	m_filename = filename;
	m_firstRun = true;
	m_fileType = fileType;
	m_useExtraFilenameString = false;
}

Meter::Meter(char* filename, Storage::FilePreference fileType, Storage::SaveDataState state)
{
	m_filename = filename;
	m_firstRun = true;
	m_fileType = fileType;
	m_useExtraFilenameString = false;

	m_saveState = state;
}

void Meter::AttachUnit(Unit* unit)
{
	m_attachedUnits.push_back(unit);
	unit->SetRecording(true);
	m_meterTypes.push_back(Meter::MeterUnit);
}

void Meter::AttachProjection(Projection* Projection,int samplingRate)
{
	m_attachedProjections.push_back(Projection);
	Projection->SetRecording(true,samplingRate);
	m_meterTypes.push_back(Meter::MeterProjection);
	m_additionalInfo.push_back(Projection->PreLayer()->GetLayerId());
	m_additionalInfo.push_back(Projection->PostLayer()->GetLayerId());
}

void Meter::AttachPopulation(Population* layer)
{
	m_attachedLayers.push_back(layer);
	layer->SetRecording(true);
	layer->SetAttachedMeter(this);
	m_meterTypes.push_back(Meter::MeterLayer);
	m_additionalInfo.push_back(layer->GetLayerId());
}

void Meter::AttachPopulationModifier(PopulationModifier* layer)
{
	m_attachedPopulationModifiers.push_back(layer);
	layer->SetRecording(true);
	m_meterTypes.push_back(Meter::MeterPopulationModifier);
	m_additionalInfo.push_back(layer->GetPopulation()->GetLayerId());
}

void Meter::AttachObject(NetworkObject* object, MeterType type)
{
	m_attachedObjects.push_back(object);
	object->SetRecording(true);
	m_meterTypes.push_back(type);

	if(type == Meter::MeterPopulationModifier)
	{
		m_additionalInfo.push_back(((PopulationModifier*)object)->GetPopulation()->GetLayerId());
	}
	else if(type == Meter::MeterLayer)
	{
		m_additionalInfo.push_back(((Population*)object)->GetLayerId());
		((Population*)object)->SetAttachedMeter(this);
	}
	else if(type == Meter::MeterProjection)
	{
		m_additionalInfo.push_back(((Projection*)object)->PreLayer()->GetLayerId());
		m_additionalInfo.push_back(((Projection*)object)->PostLayer()->GetLayerId());
	}
}

/// <summary>	Stores all recorded data (from added meters) to disk.</summary>
///
/// <param name="timestep">	The current timestep (may write as well). </param>

void Meter::RecordAll(float timestep)
{
	if(m_firstRun)
		m_storage = new Storage();
	
	int tot = m_attachedUnits.size() + m_attachedLayers.size() + m_attachedProjections.size() + m_attachedPopulationModifiers.size() + m_attachedObjects.size();
	
	if(tot>0)
	{
		vector<vector<float> > dataToWrite;//((tot);
		
		int index = 0;

		Storage::SaveDataState state;
		state = Storage::Standard;

		if(m_saveState == Storage::Append)
			if(m_firstRun == false)
				state = m_saveState;

		// currently all data fetched from same function name - this could be specified instead (or function could be parameterized)

		for(int i=0;i<m_attachedUnits.size();i++)
		{
			vector<vector<float> > val = m_attachedUnits[i]->GetValuesToRecord();
			for(int j=0;j<val.size();j++)
				dataToWrite.push_back(val[j]);

			m_attachedUnits[i]->ClearMemory();
		}

		for(int i=0;i<m_attachedLayers.size();i++)
		{
			vector<vector<float> > val = m_attachedLayers[i]->GetValuesToRecord();

			if(val.size()>0)
			{
				// special case - all processes stored separately during runtime, put together
				int sizeItem = 0; 
				sizeItem = val[0].size();
				int totSize = 0;
				MPI_Allreduce(&sizeItem,&totSize,1,MPI_INT,MPI_SUM,NETWORK_COMM_WORLD); // assumes all processes part of layer

				vector<int> rcounts(this->network()->MPIGetNrProcs());
				vector<int> displs(this->network()->MPIGetNrProcs());

				// counts from each node
				MPI_Allgather( &sizeItem, 1, MPI_INT, &rcounts[0], 1, MPI_INT, NETWORK_COMM_WORLD);

				for(int i=1;i<rcounts.size();i++)
				{
					displs[i] = rcounts[i-1]+displs[i-1];
				}

				for(int j=0;j<val.size();j++)
				{
					vector<float> data(totSize);
					MPI_Allgatherv( &val[j][0],val[j].size(),MPI_FLOAT,&data[0],&rcounts[0],&displs[0],MPI_FLOAT,NETWORK_COMM_WORLD);

					if(this->network()->MPIGetNodeId() == 0)
					{
						dataToWrite.push_back(data);
					}
				}

				m_attachedLayers[i]->ClearMemory();
			}
		}

		for(int i=0;i<m_attachedPopulationModifiers.size();i++)
		{
			vector<vector<float> > val = m_attachedPopulationModifiers[i]->GetValuesToRecord();
			for(int j=0;j<val.size();j++)
				dataToWrite.push_back(val[j]);
		}

		for(int i=0;i<m_attachedProjections.size();i++)
		{
			vector<vector<float> > val = m_attachedProjections[i]->GetValuesToRecord();
			for(int j=0;j<val.size();j++)
				dataToWrite.push_back(val[j]);
		}

		for(int i=0;i<m_attachedObjects.size();i++)
		{
			vector<vector<float> > val = m_attachedObjects[i]->GetValuesToRecord();
			for(int j=0;j<val.size();j++)
				dataToWrite.push_back(val[j]);
		}

		char* filename = new char[100];
		stringstream ss2;
		if(this->m_useExtraFilenameString == true)
		{
			ss2<<"E"<<this->m_filenameExtraString<<"_"<<this->m_filename;
			strcpy(filename,ss2.str().c_str());
		}
		else
		{
			strcpy(filename,m_filename);
		}


		if(m_fileType == Storage::CSV)
		{
			TimingStart("MeterSaveCSV");
			stringstream ss;
			ss<<"t="<<timestep<<"\n";

			if(dataToWrite.size()>0)
			{
				m_storage->SaveDataFloatCSV(filename,dataToWrite,state, "");
			}

			TimingStop("MeterSaveCSV");
		}
		else if(m_fileType == Storage::MPI_Binary)
		{
			if(dataToWrite.size()>0)
			{
				if(this->m_mpiCommunicator == NULL)
					m_storage->SaveDataFloatMPIBin(filename,dataToWrite,this->network()->MPIGetNodeId(),this->network()->MPIGetNrProcs(),NETWORK_COMM_WORLD);//this->GetMPICommunicator()); // typically assumes a NETWORK_COMM_WORLD distribution
				else
					m_storage->SaveDataFloatMPIBin(filename,dataToWrite,this->network()->MPIGetNodeId(),this->network()->MPIGetNrProcs(),this->GetMPICommunicator()); // typically assumes a NETWORK_COMM_WORLD distribution
			}
		}
		else if(m_fileType == Storage::MPI_HDF5)
		{
			// TODO: Implement
		}

		delete[] filename;
	}

	if(m_firstRun)
		m_firstRun = false;
}