#pragma once
#ifndef STORAGE5_H
#define STORAGE5_H

#define HDF5_AVAILABLE 0

#define ARCHITECTURE_BGL 1

#include <string>
#include "csv.h"
//#include "Logger.h"

#if HDF5_AVAILABLE==1
	#include "hdf5.h" // hdf5d.lib or hdf5.lib
	//#include "hdf5_hl.h" // hdf5_hld.lib or hdf5_hl.lib
#endif

using namespace std;

/* 
/////////////////  Comments on setting up hdf5 ////////////////////////////

downloads
http://www.hdfgroup.org/HDF5/release/obtain5.html
http://www.hdfgroup.org/ftp/HDF5/current/src/

for debug: 
hdf5d.lib;hdf5_cppd.lib;hdf5_hld.lib;hdf5_hl_cppd.lib;libszipd.lib;zlibd.lib
necessary now only: hdf5d.lib;hdf5_hld.lib;libszipd.lib;zlibd.lib

for release: 
hdf5.lib;hdf5_cpp.lib;hdf5_hl.lib;hdf5_hl_cpp.lib;libszip.lib;zlib.lib

libs need to use the same version of libraries: debug/release, dynamic/static, single/multi-threaded
http://msdn.microsoft.com/en-us/library/abx4dbyh%28VS.80%29.aspx

windows build instructions
http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL_Windows.txt

download zlib and szip and include the directories with header files in the solution properties when compiling
first build zlib, szip, then use debug lib binaries for linking hdf5

remember when compiling, always to do a complete re-build both in network and in libs
take sub-projects hdf5 and hdf5_hl for compilation (rest not needed)

problem description and solutions: http://support.microsoft.com/kb/148652
workaround: ignore libraries in linker settings, e.g. MSVCR100D.dll, msvcrtd.lib;
solution: [[/verbose:lib]] in project settings

finally working when ignoring some VS libraries! 

for hebb, hdf5 not installed already: http://www.pdc.kth.se/resources/software/installed-software
configuration instructions with cygwin: http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL_Cygwin.txt
szip binaries for aix: ftp://ftp.hdfgroup.org/lib-external/szip/2.0/bin/aix/

o  SZIP Library is not available for Crays SV1 and T3E
http://www.hdfgroup.org/ftp/HDF/prev-releases/ReleaseFiles/RELEASE-4.2r0.txt

*/


class Sample{
public:
	Sample(){
		data.resize(0);
		set=Undefined;
		compoundnr=-1;
	}
	int getCompound(){
		return this->compoundnr;
	}
	int getSampleNr(){
		return this->samplenr;
	}
	vector<float> getData(){
		return data;
	}	
	float getSensor(int index){
		return data[index];
	}
	enum SetType{
		Training,
		Test,		
		Validation,
		Interim,
		Undefined
	};
	float concentrations[3];
	int samplenr;	
	SetType set;
	int compoundnr;
/*	hvl_t data_hvl_t;
	void savedata(){
		data_hvl_t.p = &data[0];
		float dummy;
		data_hvl_t.len = sizeof( dummy )*data.size();	
	}
	*/
	vector<float> data;	
};

class UPCdata{
	// UPC data treatment
	public:
	UPCdata(char* filename,int filetype=0){
		//if(filetype==0)
			readCSVfile(filename);
		//else
			//readHDF5file(filename); // just for test
			// readHDFfile(filename); 
	}
	//UPCdata(){ // empty constructor, do nothing (for now)
	//}
	void readCSVfile(char* filename);
	void writeHDFfile(char* filename);
	void readHDF5file(char* filename);
	void writeHDF5file(char* filename);
	void readHDFfile(char* filename,int count=14640);

	vector<Sample> getTraining(){
		vector<Sample> training;
		for(int i=0;i<samples.size();i++)
			if(samples[i].set==Sample::Training)
				training.push_back(samples[i]);

		return training;
	}
	vector<Sample> getValidation(){
		vector<Sample> validation;
		for(int i=0;i<samples.size();i++)
			if(samples[i].set==Sample::Validation)
				validation.push_back(samples[i]);

		return validation;
	}
	vector<Sample> getTest(){
		vector<Sample> test;
		for(int i=0;i<samples.size();i++)
			if(samples[i].set==Sample::Test)
				test.push_back(samples[i]);

		return test;
	}
	vector<string> getCompounds(){
		return compounds;
	}
	vector<float> getSensor(int index){
		// gets data for one column, i.e. one sensor
		vector<float> sensor;
		for(int i=0;i<samples.size();i++)			
				sensor.push_back(samples[i].getSensor(index));	
	}
	vector<vector< float> > getData(){
		vector<vector<float> > data;
		for(int i=0;i<samples.size();i++)			
				data.push_back(samples[i].getData());	
	}
	vector<vector< float> > getData(Sample::SetType set){
		// return data for either training, testing, or validation
		vector<vector<float> > data;
		for(int i=0;i<samples.size();i++)
			    if(samples[i].set==set)
				data.push_back(samples[i].getData());
	}
private:	
	vector<string> compounds;
	vector<Sample> samples;	
	struct csv_parser m_csvParser;
};


class Storage
{
public:	

	enum FilePreference
	{
		CSV,
		MPI_Binary,
		MPI_HDF5
	};

	enum SaveDataState
	{
		Standard,
		Open,
		Append,
		Close
	};

	vector<vector<float> > LoadDataFloatCSV(char* filename, int nrItems, bool keepOpen);
	vector<vector<float> > LoadDataFloatCSVNextItems(int nrItems, bool close);
	vector<vector<vector<float> > > LoadDataFloatBinSparse(char* filename, int nrItems, bool keepOpen);

	vector<vector<float> > LoadDataFloatMPIBin(char* filename, int nrItems, int startColumn, int endColumn, MPI_Comm comm);

	void SeparateTrainingTest(vector<vector<float> > data,vector<vector<float> > &training,vector<vector<float> > &test,float trainingper=.9);

	//vector<vector<float> > LoadDataFloatHDF5(char* filename, char* datasetname, int startItem, int endItem);
	//vector<vector<float> > LoadDataFloatHDF5(char* filename, char* datasetname, vector<int> columnsLoad, int startItem, int endItem);

	void SaveDataFloatCSV(char* filename, vector<vector<float> > data, SaveDataState state, string first);
	void SaveDataFloatMPIBin(char* filename, vector<vector<float> > data, int mpiRank, int mpiSize, MPI_Comm comm);//, SaveDataState state, string first);
	//void SaveDataFloatHDF5(char* filename, char* datasetname, vector<float> data, SaveDataState state);
	//void SaveDataFloatHDF5(char* filename, long datasetStartId, vector<vector<float> > data, SaveDataState state);
	//void CreateFileHDF5(char* filename);

	void SetMPIParameters(int mpiRank,int mpiSize)
	{
		m_mpiRank = mpiRank;
		m_mpiSize = mpiSize;
	}

	// Alternative load class into strings
	vector<string> GetFileData(string filename);

private:
	int m_mpiRank,m_mpiSize;
	struct csv_parser m_csvParser;

	FILE *m_fp;

	// kept for faster appending in SaveDataFloatHDF5
/*	hid_t       m_file, m_dataset;         // file and dataset handles
    hid_t       m_datatype, m_dataspace;   // handles
    hsize_t     m_dimsf[1];              // dataset dimensions*/
};

#endif
