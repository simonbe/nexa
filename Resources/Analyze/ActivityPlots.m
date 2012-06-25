layer2 = 'C:\\CurrentProjects\\Network\\debug\\Layer2Activity_triesch.csv';
layer3 = 'C:\\CurrentProjects\\Network\\debug\\Layer3Activity.csv';
layer2test = 'C:\\CurrentProjects\\Network\\debug\\Layer2Activity_test.csv';
layer3test = 'C:\\CurrentProjects\\Network\\debug\\Layer3Activity_test.csv';
outValues = 'C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_trainingLabels.csv';
outValuesTest = 'C:\\CurrentProjects\\Network\\Databases\\MNIST\\MNIST_testingLabels.csv';

figure(23);
ActivityPlotLayer(layer2, 'none',1)%outValues,1);
%figure(22);
%ActivityPlotLayer(layer3, 'none',1);
%figure(23);
%ActivityPlotLayerWTA(layer3, outValues,2);
%figure(24);
%ActivityPlotLayer(layer2test, outValuesTest,1);
%figure(25);
%ActivityPlotLayerWTA(layer3test, outValuesTest,2);