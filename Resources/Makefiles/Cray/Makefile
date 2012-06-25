# This is a simple makefile that compiles multiple C++ source files

# set the names here to be the names of your source files with the
# .cxx or .cpp replaced by .o
# Be *** SURE *** to put the .o files here rather than the source files

ProjectObjects = main.o NetworkScalingDemos.o NetworkAdThreshold.o Analysis.o OlfactionSystem.o NetworkOlfactionScenariosClassification.o NetworkOlfactionScenariosTemplate.o NetworkOlfactionScenariosClassification_PH.o NetworkIFTests.o StructureReTIDe.o NetworkTopDownDecorr.o KernelConnectivity.o NetworkSanger.o NetworkAdaptation2.o EarlyOlfSys.o NetworkSensorDrift.o NetworkAdaptation.o NetworkBCM.o NetworkHebbSimple.o NetworkOlfaction.o NetworkOlfactionConcInv.o NetworkOlfactionMixtures.o NetworkTriesch.o NetworkCorr.o DataSources.o Meter.o NetworkMNIST.o StructureMIMDSVQ.o NetworkCL.o NetworkFoldiak.o TransferFunctions.o NetworkTests.o MPIDistribution.o NetworkBCPNN.o NetworkKussul.o NetworkConnections.o NetworkEventConnection.o NetworkEventUnit.o NetworkEventLayer.o NetworkLayer.o Network.o NetworkMDS.o NetworkMI.o NetworkSTDP.o NetworkUnits.o NetworkUnitsIF.o NetworkVQ.o ParamReader.o Storage.o csv.o ParameterSet.o ParameterSet_OlfactionPH.o 

#------------ no need to change between these lines -------------------
CFLAGS = -O3
CXX_INCLUDES = 
CXX_LIBS = 

.SUFFIXES: .cxx .cpp

.cxx.o:
	CC $(CFLAGS) $(CXX_INCLUDES) -c $<

.cpp.o:
	CC $(CFLAGS) $(CXX_INCLUDES) -c $<


#------------ no need to change between these lines -------------------


#------------ targets --------------------------------------------
# describe how to create the targets - often there will be only one target

all: projectCray

projectCray: $(ProjectObjects)
	CC $(CFLAGS) $(CXX_INCLUDES) $(ProjectObjects) -o projectCray $(CXX_LIBS)

clean:
	rm -f $(ProjectObjects) projectCray

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
OlfactionSystem.o: OlfactionSystem.h
NetworkOlfactionScenariosClassification.o: NetworkOlfactionScenariosClassification.h
NetworkOlfactionScenariosTemplate.o: NetworkOlfactionScenariosTemplate.h
NetworkOlfactionScenariosClassification_PH.o: NetworkOlfactionScenariosClassification_PH.h
ParameterSet.o: ParameterSet.h  
ParameterSet_OlfactionPH.o: ParameterSet_OlfactionPH.h 
NetworkOlfactionScenariosClassification_PH.o: NetworkOlfactionScenariosClassification_PH.h
Analysis.o: Analysis.h
NetworkAdThreshold.o: NetworkAdThreshold.h
NetworkScalingDemos.o: NetworkScalingDemos.h
main.o: