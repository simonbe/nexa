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

void Meter::AttachConnection(Connection* connection,int samplingRate)
{
	m_attachedConnections.push_back(connection);
	connection->SetRecording(true,samplingRate);
	m_meterTypes.push_back(Meter::MeterConnection);
	m_additionalInfo.push_back(connection->PreLayer()->GetLayerId());
	m_additionalInfo.push_back(connection->PostLayer()->GetLayerId());
}

void Meter::AttachLayer(Population* layer)
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
	else if(type == Meter::MeterConnection)
	{
		m_additionalInfo.push_back(((Connection*)object)->PreLayer()->GetLayerId());
		m_additionalInfo.push_back(((Connection*)object)->PostLayer()->GetLayerId());
	}
}

void Meter::RecordAll(float timestep)
{
	if(m_firstRun)
		m_storage = new Storage();
	
	int tot = m_attachedUnits.size() + m_attachedLayers.size() + m_attachedConnections.size() + m_attachedPopulationModifiers.size() + m_attachedObjects.size();
	
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
				//for(int j=0;j<val.size();j++)
				//	dataToWrite.push_back(val[j]);

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

		for(int i=0;i<m_attachedConnections.size();i++)
		{
			vector<vector<float> > val = m_attachedConnections[i]->GetValuesToRecord();
			for(int j=0;j<val.size();j++)
				dataToWrite.push_back(val[j]);
		}

		for(int i=0;i<m_attachedObjects.size();i++)
		{
			vector<vector<float> > val = m_attachedObjects[i]->GetValuesToRecord();
			for(int j=0;j<val.size();j++)
				dataToWrite.push_back(val[j]);
		}

		/*if(m_firstRun)
		state = Storage::Open;
		else
		state = Storage::Append;*/

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
			//filename = m_filename;
		}


		if(m_fileType == Storage::CSV)
		{
			TimingStart("MeterSaveCSV");
			stringstream ss;
			ss<<"t="<<timestep<<"\n";

			if(dataToWrite.size()>0)// || m_firstRun == true)
			{
				m_storage->SaveDataFloatCSV(filename,dataToWrite,state, "");//ss.str());	
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
			RecordUnitsHDF5();
		}

		delete[] filename;
	}

	if(m_firstRun)
		m_firstRun = false;
}

void Meter::RecordUnitsHDF5()
{
	vector<vector<float> > values(m_attachedUnits.size());

	for(int i=0;i<m_attachedUnits.size();i++)
	{
		vector<vector<float> > val = m_attachedUnits[i]->GetValuesToRecord();
		//values[i] = val;

		/*if(m_firstRun)
		m_Storage->SaveDataFloatHDF5(m_filename,"temp",val,Storage::Open);
		else
		m_Storage->SaveDataFloatHDF5(m_filename,"temp",val,Storage::Append);*/

	}
}