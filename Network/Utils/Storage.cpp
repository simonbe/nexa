#include <iostream>
#include <string>
#include <mpi.h>
//#include <mpio.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sstream>


#include "Storage.h"


void UPCdata::readHDF5file(char* filename){
#if HDF5_AVAILABLE==1
	int     matrix[5][5];
	hid_t   file, dataset;           // File and dataset identifiers 
	herr_t  status;

	file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset = H5Dopen1(file,"dataset1");
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                  H5P_DEFAULT, matrix);


    for (int i=0; i < 5; i++) {
        for(int j=0; j < 5; j++) 
			printf("%3d  \n", matrix[i][j]);
        printf("\n");
    }
	
#endif
}

/*
void UPCdata::writeHDF5file(char* filename){
#define dim1 5
#define dim2 5
#define ndims 2

	int     matrix[dim1][dim2];
	hid_t   file, dataset;           // File and dataset identifiers 
	herr_t  status;
	hid_t   fid; 
   for (int i = 0; i < dim1; i++) {
       for (int j = 0; j < dim2; j++)
		matrix[i][j] = i+j;
    }

<<<<<<< .mine
vector<vector<float> > Storage::LoadDataFloatHDF5(char* filename, char* datasetname, int startItem, int endItem)
{
	vector<int> allColumns;
	return LoadDataFloatHDF5(filename,datasetname,allColumns,startItem,endItem);
}

vector<vector<float> > Storage::LoadDataFloatHDF5(char* filename, char* datasetname, vector<int> columnsLoad, int startItem, int endItem)
{
	hid_t       file;                        // handles 
    hid_t       dataset;
    hid_t       filespace;
    hid_t       memspace;
    hid_t       cparms;
    hsize_t     dims[2];                     // dataset and chunk dimensions
    hsize_t     chunk_dims[2];
    hsize_t     col_dims[1];
    hsize_t     count[2];
    hsize_t     offset[2];
=======
    hsize_t fdim[] = {dim1, dim2}; 
		
	file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);  // Create a file.	
	fid = H5Screate_simple(ndims, fdim, NULL);// Create dataspace for the dataset in the file. Important first parameter: dimensions of the dataset
	dataset = H5Dcreate1(file, "/dataset2", H5T_NATIVE_INT, fid, H5P_DEFAULT); // Create dataset 
	status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix); // write it into the file.
>>>>>>> .r628

	status = H5Sclose(fid); 
    status = H5Dclose(dataset);
	status = H5Fclose(file);

	char key;
	cin >> key;

	return;
}  
*/


void UPCdata::writeHDF5file(char* filename){
#define LENGTH 5
#if HDF5_AVAILABLE==1

    typedef struct s1_t {
		int    a;
		float  b;
		double c; 
    } s1_t;

    s1_t       s1[LENGTH];
    hid_t      s1_tid;     /* File datatype identifier */

    hid_t      file, dataset, space; /* Handles */
    herr_t     status;
    hsize_t    dim[] = {LENGTH};   /* Dataspace dimensions */


    /*
     * Initialize the data
     */
    for (int i = 0; i< LENGTH; i++) {
        s1[i].a = i;
        s1[i].b = i*i;
        s1[i].c = 1./(i+1);
    }

    /*
     * Create the data space.
     */
    space = H5Screate_simple(1, dim, NULL);  // one dimensional in this case

    /*
     * Create the file.
     */
    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create the memory datatype. 
     */
    s1_tid = H5Tcreate (H5T_COMPOUND, sizeof(s1_t));
    H5Tinsert(s1_tid, "a_name", HOFFSET(s1_t, a), H5T_NATIVE_INT);
    H5Tinsert(s1_tid, "c_name", HOFFSET(s1_t, c), H5T_NATIVE_DOUBLE);
    H5Tinsert(s1_tid, "b_name", HOFFSET(s1_t, b), H5T_NATIVE_FLOAT);

    /* 
     * Create the dataset.
     */
    dataset = H5Dcreate1(file, "/dataset3", s1_tid, space, H5P_DEFAULT);

    /*
     * Write data to the dataset; 
     */
    status = H5Dwrite(dataset, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, s1);

    /*
     * Release resources
     */
    H5Tclose(s1_tid);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);
#endif
}



void UPCdata::writeHDFfile(char* filename){
#if HDF5_AVAILABLE == 1
#define NRECORDS 50
	//const int NRECORDS=samples.size();
	const int NFIELDS = 5;
	hsize_t    chunk_size = 10;  // ??	

	typedef struct _samplestruct{
		float data[1020];  // having difficulties here; dynamic allocation? 
		float concentrations[3];
		int samplenr;
		int set;
		int compoundnr;		
	}samplestruct_t;
		
	cout << "converting data format" << endl;
	samplestruct_t* sampledata=NULL;
	sampledata=new samplestruct_t[NRECORDS];
	for(int i=0;i<NRECORDS;i++){
		for(int j=0;j<1020;j++)
			sampledata[i].data[j]=samples[i].data[j];
		sampledata[i].concentrations[0]=samples[i].concentrations[0];
		sampledata[i].concentrations[1]=samples[i].concentrations[1];
		sampledata[i].concentrations[2]=samples[i].concentrations[2];
		sampledata[i].samplenr=samples[i].samplenr;
		sampledata[i].set=samples[i].set;
		sampledata[i].compoundnr=samples[i].compoundnr;
	}
	
	size_t dst_offset[NFIELDS] = {
		HOFFSET(samplestruct_t, data), 
		HOFFSET(samplestruct_t, concentrations), 
		HOFFSET(samplestruct_t, samplenr), 
		HOFFSET(samplestruct_t, set),
		HOFFSET(samplestruct_t, compoundnr), 
	};
	size_t dst_sizes[NFIELDS] = { 
		sizeof( sampledata[0].data),
		sizeof(sampledata[0].concentrations),
		sizeof(sampledata[0].samplenr),
		sizeof(sampledata[0].set),
		sizeof(sampledata[0].compoundnr)
	};
	
	hsize_t  dim[] = {NRECORDS};   /* Dataspace dimensions */
	hid_t space = H5Screate_simple(1, dim, NULL);  // one dimensional in this case
	hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); // Create a new file using default properties. 
	
	const hsize_t datadim[] = {1020}; 
	const hsize_t concdim[] = {3};    
	hid_t dataarray = H5Tarray_create(H5T_NATIVE_FLOAT, 1, datadim);
	hid_t concarray = H5Tarray_create(H5T_NATIVE_FLOAT, 1, concdim);   
    hid_t s_tid = H5Tcreate (H5T_COMPOUND, sizeof(samplestruct_t));
	H5Tinsert(s_tid, "Data", HOFFSET(samplestruct_t, data), dataarray);
    H5Tinsert(s_tid, "Concentrations", HOFFSET(samplestruct_t, concentrations), concarray);
    H5Tinsert(s_tid, "samplenr", HOFFSET(samplestruct_t, samplenr), H5T_NATIVE_INT);
    H5Tinsert(s_tid, "set", HOFFSET(samplestruct_t, set), H5T_NATIVE_INT);
	H5Tinsert(s_tid, "compoundnr", HOFFSET(samplestruct_t, compoundnr), H5T_NATIVE_INT);

	hid_t dataset = H5Dcreate1(file, "/dataset3", s_tid, space, H5P_DEFAULT);
	herr_t status = H5Dwrite(dataset, s_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, sampledata);

	delete [] sampledata;
	H5Tclose(dataarray);
	H5Tclose(concarray);
    H5Tclose(s_tid);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);

#endif
}


void UPCdata::readHDFfile(char* filename, int count){ // count: 11040 for segmentation and 14640 for classification
#if HDF5_AVAILABLE==2
    // based on http://www.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html

	// maybe first deleting all data? does it do that automatically or is there some possible memory leak? 
	// I need to know first the data set sizes... so I create one dummy sample... this is only a hack 
	samples.resize(0); 
	/*
	Sample dummy; 
	dummy.data.resize(1020,0); 
	dummy.concentrations[0]=1;dummy.concentrations[1]=1;dummy.concentrations[2]=1; 
	dummy.samplenr=1; 
	dummy.set=Sample::Training; 
	dummy.compoundnr=1; 
	*/

	const int NFIELDS = 5;
	hsize_t    chunk_size = 10;  // ??	
	
	typedef struct _samplestruct{
		int concentrations[3];
		int samplenr;
		int set;
		int compoundnr;
		float data[1020];  // having difficulties here; dynamic allocation? 
	}samplestruct;
		
	samplestruct *sampledata;
	void *sampleptr = malloc(NRECORDS*sizeof(samplestruct));
	sampledata = (samplestruct*)sampleptr;
	if (sampleptr == NULL) {
       // Memory could not be allocated, the program should handle the error here as appropriate. 
	   cout << "memory could not be allocated!";
	   return;
	} 

	// ... 
	// free(samplestruct);
	
/*	size_t dst_size =  sizeof(Sample);
	size_t dst_offset[NFIELDS] = { 
		HOFFSET(Sample, concentrations), 
		HOFFSET(Sample, samplenr), 
		HOFFSET(Sample, set), 
		HOFFSET(Sample, compoundnr), 
		HOFFSET(Sample, data)
		};
	size_t dst_sizes[NFIELDS] = { sizeof( samples[0].data),
		sizeof(samples[0].concentrations),
		sizeof(samples[0].samplenr),
		sizeof(samples[0].set),
		sizeof(samples[0].compoundnr)};
		*/

	size_t dst_size =  sizeof(samplestruct);
	size_t dst_offset[NFIELDS] = {HOFFSET(samplestruct, concentrations), 
		HOFFSET(samplestruct, samplenr), 
		HOFFSET(samplestruct, set), 
		HOFFSET(samplestruct, compoundnr),
		HOFFSET(samplestruct, data), 
	};
	size_t dst_sizes[NFIELDS] = { sizeof( sampledata[0].data),
		sizeof(sampledata[0].concentrations),
		sizeof(sampledata[0].samplenr),
		sizeof(sampledata[0].set),
		sizeof(sampledata[0].compoundnr)
	};

	hid_t      file_id;  
	int        *fill_data = NULL;
	int        compress  = 1;
	herr_t     status;
	
	//const hsize_t datadim[] = {samples[0].data.size()}; 
	const hsize_t datadim[] = {1020*sizeof(float)}; 
	const hsize_t concdim[] = {3};    

	hid_t floatarray = H5Tarray_create(H5T_NATIVE_FLOAT, 1, datadim);
	hid_t concarray = H5Tarray_create(H5T_NATIVE_FLOAT, 1, concdim);    

	// Initialize field_type 
	hid_t      field_type[NFIELDS]; 
	field_type[0] = floatarray;
	field_type[1] = concarray;
	field_type[2] = H5T_NATIVE_INT;
	field_type[3] = H5T_NATIVE_INT;
	field_type[4] = H5T_NATIVE_INT;
	
	// opening file
	file_id=H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( file_id < 0 ){
		std::cout << "UPCdata::readHDFfile(): Error opening " << filename <<": " << file_id << std::endl;
		return;
	}
	
	//const int NRECORDS=samples.size(); // how to get the data set size? 
	// is this the way? 
	hid_t dataset_id=H5Dopen1(file_id, "Dataset");
	hsize_t size = H5Dget_storage_size(dataset_id);
	
	//vector<Sample> buf(NRECORDS+1);
	//vector<Sample> buf(static_cast<int>(size+1), 0x00);

	status=H5TBread_table(file_id, "Dataset", dst_size, dst_offset, dst_sizes, sampledata);

	cout << "converting data format" << endl;
	for(int i=0;i<NRECORDS;i++){
		Sample dummysample;
		dummysample.concentrations[0]=sampledata[i].concentrations[0];
		dummysample.concentrations[1]=sampledata[i].concentrations[1];
		dummysample.concentrations[2]=sampledata[i].concentrations[2];
		dummysample.samplenr=sampledata[i].samplenr;
		dummysample.set=(Sample::SetType)sampledata[i].set;
		dummysample.compoundnr=sampledata[i].compoundnr;
		dummysample.data.resize(1020);
		for(int j=0;j<1020;j++)
			dummysample.data[j]=sampledata[i].data[j];
		samples.push_back(dummysample);
	}

	// close type(s)??
	// close the file /
	H5Fclose( file_id );	
#endif
}


/*
int writehdf5(char* filename,char* datasetname, vector<vector<float> > data){	
	// open file
	herr_t  status;
	hid_t fh = H5Fcreate(filename, 1, 1, 1);  // H5ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT

	// create and initialize the essential components of a dataset for writing to a file
	hid_t    dataset, datatype, dataspace;   // declare identifiers 
	hsize_t dimsf[2];
	dimsf[0] = data[0].size();
	dimsf[1] = data.size();
	dataspace = H5Screate_simple(1, dimsf, NULL);  // RANK??
	datatype = H5Tcopy(H5T_NATIVE_FLOAT);
	status = H5Tset_order(datatype, H5T_ORDER_LE);
	dataset = H5Dcreate(fh, datasetname, datatype, dataspace, H5P_DEFAULT);

	// would this instead work? 
	// const void* data2=&data.begin; 
	const void* data2;
	copy(data.begin(), data.end(), data2);

	// write a dataset to a new file
	status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data2);

	status = H5Fclose(fh); 
	return 0;
}
*/

struct counts {
	vector<float> fData;
  long unsigned fields;
  long unsigned rows;
};

void cbColumn (void *s, size_t len, void *data) 
{ 
	float f = atof((char*)s);
	((struct counts *)data)->fData.push_back(f);
	s = const_cast<char*>("                                                                                      ");
	//((struct counts *)data)->fields++; 
}
void cbRow (int c, void *data) 
{ 
	((struct counts *)data)->rows++; 
}

/*
  FILE * pFile;
  long lSize;
  char * buffer;
  size_t result;

  pFile = fopen ( "myfile.bin" , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (1);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (2);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); exit (3);}

  // terminate
  fclose (pFile);
  free (buffer);
  
*/

void UPCdata::readCSVfile(char* filename){
	/* UPC dataset is of an un-orthodox format
	 each row, starting from second, contain 1020 data columns and 6 description columns
	 data columns correspond to 60 ORNs per 17 OR types. 
	 description columns denote compound, set (TrainingSet, TestSet, ValidationSet), conc1, conc2, conc3, sample (minute)
	 I use c-style parameter handling in that vector parameters are handed over by pointers (c++ references to be exact). 
	 */

	cout << "reading file " << filename << endl;

	// empty output vectors
	compounds.resize(0);
	samples.resize(0);
	
	/*
	std::ifstream fh(filename, ifstream::in);
	if (fh.is_open()){
		if(fh.good()){
			string line;
			getline(fh,line); // ignore first line
		}
	*/

	char buf[10000];
	FILE *fp = fopen( filename, "rb" );
	fgets(buf, 10000, fp); // ignore first line
	while(fgets(buf, 10000, fp)){
	//	while(fh.good() ){
			// initalizations for file reading
			string line(buf);
			Sample sample;	
			
			struct counts c = {vector<float>(0),0, 0};

			if (csv_init(&m_csvParser, CSV_STRICT) != 0){
				fprintf(stderr, "Failed to initialize csv parser\n");
				return;
				//exit(EXIT_FAILURE);
			}

			//getline(fh,line);
			int bytes_read=line.size();
			if (csv_parse(&m_csvParser, buf, bytes_read, cbColumn, cbRow, &c) != bytes_read){
				fprintf(stderr, "Error while parsing file: %s\n", csv_strerror(csv_error(&m_csvParser)));
			}			
			csv_fini(&m_csvParser, cbColumn, cbRow, &c);			

			if(c.fData.size()>0){
				for(int i=0;i<c.fData.size()-6;i++)
					sample.data.push_back(c.fData[i]);

				sample.samplenr=c.fData[c.fData.size()-1]; 
				sample.concentrations[0]=c.fData[c.fData.size()-4]; 
				sample.concentrations[1]=c.fData[c.fData.size()-3]; 
				sample.concentrations[2]=c.fData[c.fData.size()-2]; 

				// now get the compound and the set without using regular expressions
				// in order to avoid problems in code compatibility if neither boost nor tr1 are available
				string str=line;
				string setstr="";
				for(int i=0;i<5;i++){
					int found=str.find_last_of(",");
					str=str.substr(0,found);
					if(i==3){
						setstr=str.substr(str.find_last_of(",")+2,str.size()-1);
					}
				}						
				string compound=str.substr(str.find_last_of(",")+2,str.size()-1);						

				if(setstr.compare("TrainingSet\"")==0)
					sample.set=Sample::Training;
				else if(setstr.compare("TestSet\"")==0)			
					sample.set=Sample::Test;
				else if(setstr.compare("ValidationSet\"")==0)
					sample.set=Sample::Validation;
				else if(setstr.compare("InterimSet\"")==0)
					sample.set=Sample::Interim;

				if(sample.set==Sample::Undefined)
					fprintf(stderr, "Storage::ReadUPCfile(): Could not match set\n");
						
				for(int i=0;i<compounds.size();i++){
					if(compounds[i]==compound){
						sample.compoundnr=i;
						break;
					}
				}
				if(sample.compoundnr<0){
					compounds.push_back(compound);
					sample.compoundnr=compounds.size()-1;
				}
				samples.push_back(sample);			
			}
		}  // line treatment
    //fh.close();
	fclose(fp);

	//} else 
	//	cerr<<"Error: Could not open file. ("<<filename<<")\n";

	cout<<"finished reading " << filename << endl;
	return;
}


vector<vector<float> > Storage::LoadDataFloatCSV(char* filename, int nrItems, bool keepOpen)
{
	vector<vector<float> > out;
	size_t bytes_read;
	char buf[1024];
	struct counts c = {vector<float>(0),0, 0};

	if (csv_init(&m_csvParser, CSV_STRICT) != 0) 
	{
		fprintf(stderr, "Failed to initialize csv parser\n");
		return out;
		//exit(EXIT_FAILURE);
	}

	//cout<<"Loading "<<filename<<"\n";
	m_fp = fopen(filename, "r");

	if(m_fp==NULL)
	{
		cerr<<"Error: Could not open file. ("<<filename<<")\n";
	}

	cout.flush();

	long oldRows = 0;

	while ((bytes_read=fread(buf, 1, 1, m_fp)) > 0) 
	{//((bytes_read=fread(buf, 1, 1024, m_fp)) > 0) {
		if (csv_parse(&m_csvParser, buf, bytes_read, cbColumn, cbRow, &c) != bytes_read) {
			fprintf(stderr, "Error while parsing file: %s\n", csv_strerror(csv_error(&m_csvParser)));
		}
		else{
//			cout<<".";cout.flush();
		}

		if(c.rows>oldRows){
			out.push_back(c.fData);
			c.fData.clear();
			oldRows = c.rows;
		}

		if(c.rows>=nrItems)
			break;
    }

	if(keepOpen == false)
	{
		csv_fini(&m_csvParser, cbColumn, cbRow, &c);
		fclose(m_fp);
	}

	if(m_mpiRank == 0)
	{
		if(out.size()>0)
			cout<<"Loaded "<<out.size()<<" nr items of dimension "<<out[0].size()<<".\n";
		else
			cout<<"Loaded no items.\n";
		cout.flush();
	}

	return out;
}

void Storage::SeparateTrainingTest(vector<vector<float> > data,vector<vector<float> > &training,vector<vector<float> > &test,float trainingper){
	// splits data matrix into two parts, training and test, with parts corresponding to trainingper,1-trainingper
	// out(1) is training data, out(2) is test data, sampling is random without repetition
	// if((trainingper<0) || (trainingper>1)) error	
	int trainingex=data.size()*trainingper;
	int testex=data.size()-trainingex;
	for(int i=0;i<data.size();i++){
		vector<float> out;
		for(int j=0;j<data[i].size();j++)
			out.push_back(data[i][j]);
		training.push_back(out);
	}
	for(int i=0;i<testex;i++){
		int ind=rand() % training.size();
		test.push_back(training[ind]);
		training.erase(training.begin()+ind);
	}
	vector<vector<vector<float> > > out;	
}

vector<vector<float> > Storage::LoadDataFloatCSVNextItems(int nrItems, bool close)
{
	vector<vector<float> > out;

	size_t bytes_read;
	char buf[1024];
	struct counts c = {vector<float>(0),0, 0};

	long oldRows = 0;

	while ((bytes_read=fread(buf, 1, 1, m_fp)) > 0) {//((bytes_read=fread(buf, 1, 1024, m_fp)) > 0) {
		if (csv_parse(&m_csvParser, buf, bytes_read, cbColumn, cbRow, &c) != bytes_read) {
			
			fprintf(stderr, "Error while parsing file: %s\n", csv_strerror(csv_error(&m_csvParser)));
		}
		
		if(c.rows>oldRows)
		{
			out.push_back(c.fData);
			c.fData.clear();
			oldRows = c.rows;
		}

		if(c.rows>=nrItems)
			break;
    }

	if(close == true)
		csv_fini(&m_csvParser, cbColumn, cbRow, &c);

	
	return out;
}

void cbWriteColumn (void *s, size_t i, void *outfile) {
  csv_fwrite((FILE *)outfile, s, i);
  fputc(',',(FILE *)outfile);
}

void cbWriteRow (int c, void *outfile) {
  fseek((FILE *)outfile, -1, SEEK_CUR);
  fputc('\n', (FILE *)outfile);
}

// nrItems need to be exactly specified
vector<vector<float> > Storage::LoadDataFloatMPIBin(char* filename, int nrItems, int startColumn, int endColumn, MPI_Comm comm)//vector<int> indexes, MPI_Comm comm)
{
	double timeStart;
	if(m_mpiRank == 0)
	{
		cout<<"Loading "<<filename<<"...";cout.flush();
		timeStart = MPI_Wtime();
	}

	vector<vector<float> > data(nrItems);

	MPI_File fh;
	MPI_File_open(comm,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&fh);

	vector<float> tempData(nrItems*(endColumn-startColumn));

	if(endColumn-startColumn == 0)
	{
		cout<<"(E) endColumn-startColumn == 0\n";cout.flush();
	}

	MPI_Status status;
	MPI_File_read_at(fh,startColumn*nrItems*sizeof(MPI_REAL4),&tempData[0],(endColumn-startColumn)*nrItems,MPI_REAL4,&status);//MPI_FLOAT,&status);
	//MPI_File_read_at(fh,startColumn*nrItems*sizeof(MPI_FLOAT),&tempData[0],(endColumn-startColumn)*nrItems,MPI_FLOAT,&status);

	
	//for(int i=0;i<(endColumn-startColumn);i++)
	for(int i=0;i<nrItems;i++)
	{
		vector<float> f(endColumn-startColumn);
		data[i] = f;
	}

	int index = 0;
	for(int j=0;j<(endColumn-startColumn);j++)
	{
		for(int i=0;i<nrItems;i++)
		{
			data[i][j] = tempData[index];
			index++;
		}
	}

	MPI_File_close(&fh);

	if(m_mpiRank == 0)
	{
		if(data.size()>0)
			cout<<"Loaded "<<data.size()<<" items of size "<<data[0].size()<<". (Time (process 0): "<<MPI_Wtime()-timeStart<<")\n";
		else
			cout<<"Warning: Loaded no items from filename: "<<filename<<"\n";
		
		cout.flush();
	}

	return data;
}

// this is currently slower than it could be because we write in chunks (data.size()..) with barriers in between - instead we could take care of correct format of matrix on the read-in side and write everything at once, or assume same format of data and only do one barrier
void Storage::SaveDataFloatMPIBin(char* filename, vector<vector<float> > data, int mpiRank, int mpiSize, MPI_Comm comm)
{
	MPI_File fh;
	// assumes specific data distribution for processes

	MPI_File_delete(filename,MPI_INFO_NULL);
	MPI_File_open(comm,filename,MPI_MODE_RDWR|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);
	int globalOffset = 0;

	// (!) Currently assumes same size of all items in data
	int size = 0;
	for(int i=0;i<data.size();i++) // assumes same nr items each process
	{
		size+= data[i].size();
	}

	vector<float> tempData(size);
	int index=0;
	for(int i=0;i<data.size();i++) // assumes same nr items each process
	{
		for(int j=0;j<data[i].size();j++)
		{
			tempData[index] = data[i][j];
			index++;
		}
	}

	vector<int> allSizes(mpiSize);
	MPI_Allgather(&size,1,MPI_INT,&allSizes[0],1,MPI_INT,comm);
	int startPos = 0;
	for(int j=0;j<mpiRank;j++)
		startPos+=allSizes[j];

	//for(int i=0;i<data.size();i++) // assumes same nr items each process
	//{
	//	/*int size = data[i].size(); // each item can be of different size
	//	vector<int> allSizes(mpiSize);
	//	MPI_Allgather(&size,1,MPI_INT,&allSizes[0],1,MPI_INT,comm);
	//	int startPos = 0;
	//	for(int j=0;j<mpiRank;j++)
	//		startPos+=allSizes[j];*/

	//	MPI_Status status;
	//	MPI_File_write_at(fh,(startPos+globalOffset)*sizeof(MPI_FLOAT),&data[i][0],data[i].size(),MPI_FLOAT,&status);
	//	//if(i==0)
	//	//MPI_File_write_at(fh,(startPos+globalStartPos),&data[i][0],data[i].size(),MPI_FLOAT,&status);

	//	for(int j=0;j<mpiSize;j++)
	//		globalOffset+=allSizes[j];
	//	//globalOffset = 0;
	//	//for(int j=mpiRank;j<mpiSize;j++)
	//	//	globalOffset+=allSizes[j];
	//}

	MPI_Status status;
	if(size>0)
		MPI_File_write_at(fh,(startPos+globalOffset)*sizeof(MPI_FLOAT),&tempData[0],size,MPI_FLOAT,&status);//&data[0][0],size,MPI_FLOAT,&status);
		//MPI_File_write_at(fh,(startPos+globalOffset)*sizeof(MPI_FLOAT),&data[0][0],size,MPI_FLOAT,&status);
	
	MPI_File_close(&fh);
}

void Storage::SaveDataFloatCSV(char* filename, vector<vector<float> > data, SaveDataState state, string first)
{
	if(state == Storage::Standard || state == Storage::Open)
	{
		// open
		/*if (csv_init(&m_csvParser, CSV_STRICT) != 0) 
		{
			fprintf(stderr, "Failed to initialize csv parser\n");
			return;
			//exit(EXIT_FAILURE);
		}
		*/
		m_fp = fopen(filename, "wb");
	}
	else if(state == Storage::Append)
	{
		m_fp = fopen(filename, "ab");
	}

	if(m_fp)
	{
		if(first.size()>0)
			fputs(first.c_str(),m_fp);

		for(int i=0;i<data.size();i++)
		{
			for(int j=0;j<data[i].size();j++)
			{
				char sz[20];
				int nr = sprintf(sz, "%G", data[i][j]);
				size_t bytes = nr*sizeof(char);
				fputs(sz,m_fp);
				//csv_fwrite(m_fp, sz, bytes);
				if(j!=data[i].size()-1)
					fputc(',',m_fp);
				/*
				csv_parse(&m_csvParser, sz, bytes, cbWriteColumn, cbWriteRow, m_fp);
				*/
			}
			//if(i!=data.size()-1)
			fputc('\n', m_fp);
		}

		//if(state == Storage::Standard || state == Storage::Close)
		//{
			// close
			//csv_fini(&m_csvParser, cbWriteColumn, cbWriteRow, m_fp);
			//csv_free(&m_csvParser);
			fclose(m_fp);
		//}
	}
	else cout<<"Error saving file "<<filename;
}

vector<vector<vector<float> > > Storage::LoadDataFloatBinSparse(char* filename, int nrItems, bool keepOpen)
{
	// slow version
	unsigned long long position = -1;
	ifstream* in = new ifstream;

	int extraPosition = 0;

	in->open(filename,ifstream::binary);

	if (!in)
	{
		cout<<"Error: Could not open data location file (LoadDataFloatBin): "<<filename<<".\n";
		//Logger::print(ss.str().c_str());
		cout.flush();
		return vector<vector<vector<float> > >(0);// out;
    }

	/*int minRow = rowsToLoad[0];
	int maxRow = rowsToLoad[0];

	for(int i=0;i<rowsToLoad.size();i++)
	{
		if(rowsToLoad[i]<minRow)
			minRow = rowsToLoad[i];
		if(rowsToLoad[i]>maxRow)
			maxRow = rowsToLoad[i];
	}*/

	// 1. first number specify nr of items/patterns in file
	float nrPatterns;
	//in->seekg( 4 );
	in->read( (char *) &nrPatterns, sizeof( float) ); // ? needs to be flipped in case of different file system architecture ?

	// 2. second number specify nr hypercolumns N
	float nrHypercolumns;
	in->read( (char *) &nrHypercolumns, sizeof( float ) ); // = 1 for semantic

	vector<int> nrMinicolumns(nrHypercolumns);
	// 3. N next numbers specify minicolumns in each hypercolumn

	int totalMinicolumns = 0;
	float nrMc = 1;
	for(int i=0;i<nrHypercolumns;i++)
	{
		in->read( (char *) &nrMc, sizeof(float));
		/*if(m_specialLimitNrValues>0)
			nrMinicolumns[i] = (int)nrMc*m_specialLimitNrValues;
		else*/
		nrMinicolumns[i] = (int)nrMc;

		totalMinicolumns += (int)nrMc;
	}

	if(nrItems>-1 && nrItems<=nrPatterns)
		nrPatterns = nrItems;

	vector<vector<vector<float> > >out(nrPatterns);//,vector<float>(0));

	//out = ublas::compressed_matrix<float>(nrPatterns,totalMinicolumns);//totalMinicolumns,nrPatterns);//rowsToLoad.size(),nrHypercolumns);

	bool finished = false;


	int index = 0;
	while(finished == false)
	{
		vector<float> vals(3000);

		if(!in->read((char*) &vals[0], sizeof(float)*3000))
		{
			finished = true;
		}
		else
		{
			vector<float> val(3);
			vector<float> dataTwo(2);

			for(int i=0;i<vals.size();i+=3)
			{
				val[0] = vals[i]-1;
				dataTwo[0] = vals[i+1]-1;
				dataTwo[1] = vals[i+2];

				if(val[0]+1>nrPatterns)
				{
					finished = true;
					i = vals.size();
				}
				else
				{
					/*if(vals[1]>maxRow+1)  // assuming 1-based indexes, starting at (1,1), otherwise change = to minRow place for zero-based
					{//	finished = true;
					}
					else
					{
					if(vals[1]>=minRow) // assuming linear (non-excluding) distribution of indexes
					{
					out.push_back(vals[0]-minRow-1,vals[1]-1,vals[2]);//out(vals[1]-1,vals[0]-minRow-1) = vals[2];
					index++;

					if(index%10000 == 0)
					{
					cout<<index<<" ";
					}
					}
					}
					*/

					out[val[0]].push_back(dataTwo);//[vals[1]-1] = vals[2];
				}
			}
		}
	}

//	ublas::compressed_matrix<float> a = ublas::compressed_matrix<float>(nrPatterns,totalMinicolumns);
//	a = ublas::trans(out);
	//out = ublas::trans(out);

	return out;
}


vector<string> Storage::GetFileData(string filename)
{
	string line;
	vector<string> allStrings;

	ifstream myfile(filename);
	if (myfile.is_open())
	{
		while ( myfile.good() )
		{
			getline (myfile,line);
			if(line.length()>0)
			{
				if(line[0]!='#')
				{
					allStrings.push_back(line);
				}
			}
		}
		myfile.close();
	}

	return allStrings;
}