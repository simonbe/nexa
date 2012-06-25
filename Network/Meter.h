#pragma once
#ifndef METER_H
#define METER_H

#include "Network.h"
#include "NetworkPopulationModifier.h"
#include "Storage.h"

class Connection;
class Unit;
class Population;
class PopulationModifier;
//class Network;

using namespace std;

class Meter : public NetworkObject
{
public:

/*	enum FileType
	{
		CSV,
		MPI_BIN,
		HDF5 // not used atm
	};*/

	enum MeterType // could be determined automatically instead of specified
	{
		MeterUnit,
		MeterLayer,
		MeterPopulationModifier,
		MeterConnection,
		MeterObject
	};

	Meter(char* filename, Storage::FilePreference fileType);//FileType fileType);
	Meter(char* filename, Storage::FilePreference fileType, Storage::SaveDataState state);

	~Meter()
	{
		delete m_storage;
		//delete m_filename;
	}

	void SetFilename(char* filename)
	{
		m_filename = filename;
	}

	void AttachUnit(Unit* unit);
	void AttachLayer(Population* layer);
	void AttachPopulationModifier(PopulationModifier* layer);
	void AttachObject(NetworkObject* object, MeterType type=MeterObject);
	void AttachConnection(Connection* conn, int samplingRate); // 0 samplingrate if only last result to be stored, otherwise weights evolution
	void RecordAll(float timestep);

	vector<int> GetMeterStructure()
	{
		vector<int> structure;
		return structure;
	}

	char* GetFilename()
	{
		return m_filename;
	}


	vector<MeterType> GetMeterTypes()
	{
		return m_meterTypes;
	}

	vector<int> GetAdditionalInfo()
	{
		return m_additionalInfo;
	}

	/*void Network(Network* network)
	{
		m_network = network;
	}*/

	Storage::FilePreference GetSaveType()
	{
		return m_fileType;
	}

	void SetExtraFilenameString(char* extraString)
	{
		m_useExtraFilenameString = true;
		m_filenameExtraString = extraString;
	}

private:

	//Network* m_network;

	vector<int> m_additionalInfo;
	bool m_firstRun;
	void RecordUnitsHDF5();
	char* m_filename;
	char* m_filenameExtraString;
	bool m_useExtraFilenameString;
	Storage* m_storage;
	Storage::FilePreference m_fileType;//FileType m_fileType;
	Storage::SaveDataState m_saveState;
	vector<MeterType> m_meterTypes;

	vector<Unit*> m_attachedUnits;
	vector<Population*> m_attachedLayers;
	vector<PopulationModifier*> m_attachedPopulationModifiers;
	vector<Connection*> m_attachedConnections;
	vector<NetworkObject*> m_attachedObjects;
};

#endif