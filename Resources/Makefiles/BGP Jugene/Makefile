# This is a simple makefile that compiles multiple C++ source files

# set the names here to be the names of your source files with the
# .cxx or .cpp replaced by .o
# Be *** SURE *** to put the .o files here rather than the source files

ProjectObjects = main.o Analysis.o NetworkScalingDemos.o NetworkOlfactionScenariosClassification_PH.o NetworkOlfactionScenariosClassification.o NetworkOlfactionScenariosTemplate.o ParameterSet_OlfactionPH.o ParameterSet.o NetworkAdThreshold.o NetworkIFTests.o StructureReTIDe.o NetworkTopDownDecorr.o KernelConnectivity.o NetworkSanger.o NetworkAdaptation2.o EarlyOlfSys.o NetworkSensorDrift.o NetworkAdaptation.o NetworkBCM.o NetworkHebbSimple.o NetworkOlfaction.o NetworkOlfactionConcInv.o NetworkOlfactionMixtures.o NetworkTriesch.o NetworkCorr.o DataSources.o Meter.o NetworkMNIST.o StructureMIMDSVQ.o NetworkCL.o NetworkFoldiak.o TransferFunctions.o NetworkTests.o MPIDistribution.o NetworkBCPNN.o NetworkKussul.o NetworkConnections.o NetworkEventConnection.o NetworkEventUnit.o NetworkEventLayer.o NetworkLayer.o Network.o NetworkMDS.o NetworkMI.o NetworkSTDP.o NetworkUnits.o NetworkUnitsIF.o NetworkVQ.o ParamReader.o Storage.o csv.o 

#------------ no need to change between these lines -------------------
CFLAGS = -O3 -qstrict -qarch=450d -qtune=450 -DMPICH_IGNORE_CXX_SEEK 
BGP_SYS = /bgsys/drivers/ppcfloor/comm
CXX_INCLUDES = -I$(BGP_SYS)/include \
		-I/bgsys/local/gsl/include

CXX_LIBS = -L$(BGP_SYS)/lib -L/bgsys/local/gsl/lib \
	-lcxxmpich.cnk -lmpich.cnk -ldcmf.cnk -ldcmfcoll.cnk -lpthread -lrt -L/bgsys/drivers/ppcfloor/runtime/SPI -lSPI.cna -lgsl -lgslcblas -lm


.SUFFIXES: .cxx .cpp

.cxx.o:
	bgxlC $(CFLAGS) $(CXX_INCLUDES) -c $<

.cpp.o:
	bgxlC $(CFLAGS) $(CXX_INCLUDES) -c $<

.c.o:
	bgxlC $(CFLAGS) $(CXX_INCLUDES) -c $<


#------------ no need to change between these lines -------------------


#------------ targets --------------------------------------------
# describe how to create the targets - often there will be only one target

all: project1

project1: $(ProjectObjects)
	bgxlC $(CFLAGS) $(CXX_INCLUDES) $(ProjectObjects) -o project1 $(CXX_LIBS)

clean:
	rm -f $(ProjectObjects) project1

#------------ dependencies --------------------------------------------
# put the .o that depends on a .h, then colon, then TAB, then the .h

MPIDistribution.o: MPIDistribution.h
Network.o: Network.h
NetworkBCPNN.o: NetworkBCPNN.h
NetworkKussul.o: NetworkKussul.h
NetworkConnections.o: NetworkConnections.h
NetworkEventConnection.o: NetworkEventConnection.h
NetworkEventLayer.o: NetworkEventLayer.h
NetworkEventUnit.o: NetworkEventUnit.h
NetworkLayer.o: NetworkLayer.h
NetworkMDS.o: NetworkMDS.h
NetworkMI.o: NetworkMI.h
NetworkSTDP.o: NetworkSTDP.h
NetworkUnits.o: NetworkUnits.h
NetworkUnitsIF.o: NetworkUnitsIF.h
NetworkVQ.o: NetworkVQ.h
Storage.o: Storage.h
Meter.o: Meter.h
NetworkUnitsIF.o: NetworkUnitsIF.h
csv.o: csv.h
KernelConnectivity.o: KernelConnectivity.h
NetworkCL.o: NetworkCL.h
TransferFunctions.o: TransferFunctions.h
NetworkFoldiak.o: NetworkFoldiak.h
NetworkMNIST.o: NetworkMNIST.h
StructureMIMDSVQ.o: StructureMIMDSVQ.h
NetworkTriesch.o: NetworkTriesch.h
NetworkCorr.o: NetworkCorr.h
DataSources.o: DataSources.h
NetworkOlfaction.o: NetworkOlfaction.h
NetworkHebbSimple.o: NetworkHebbSimple.h
NetworkBCM.o: NetworkBCM.h
NetworkAdaptation.o: NetworkAdaptation.h
NetworkAdaptation2.o: NetworkAdaptation2.h
NetworkSensorDrift.o: NetworkSensorDrift.h
NetworkIFTests.o: NetworkIFTests.h
NetworkSanger.o: NetworkSanger.h
EarlyOlfSys.o: EarlyOlfSys.h
NetworkTopDownDecorr.o: NetworkTopDownDecorr.h
StructureReTIDe.o: StructureReTIDe.h
ParamReader.o: ParamReader.h
Analysis.o: Analysis.h 
NetworkAdThreshold.o: NetworkAdThreshold.h
ParameterSet.o: ParameterSet.h
ParameterSet_OlfactionPH.o: ParameterSet_OlfactionPH.h 
NetworkOlfactionScenariosTemplate.o: NetworkOlfactionScenariosTemplate.h
NetworkOlfactionScenariosClassification.o: NetworkOlfactionScenariosClassification.h
NetworkOlfactionScenariosClassification_PH.o: NetworkOlfactionScenariosClassification_PH.h
NetworkScalingDemos.o: NetworkScalingDemos.h
main.o:
