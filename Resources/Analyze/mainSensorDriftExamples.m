clear all;
clf;

run = 5;


%%%%%%%%%% Input data

if run == 1
    dataFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example1In.csv';
	dataOutFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example1Out.csv';
    dataTimeFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example1Time.csv';
elseif run == 2
    dataFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example2In.csv';
	dataOutFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example2Out.csv';
    dataTimeFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example2Time.csv';   
elseif run == 3
    dataFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example3In.csv';
	dataOutFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example3Out.csv';
    dataTimeFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example3Time.csv';   
elseif run == 5
    dataFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example5In.csv';
	dataOutFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example5Out.csv';
    dataTimeFilename = 'C:\\CurrentProjects\\Network\\Databases\\Sensor_drift\\artificial\\example5Time.csv';   
end


dataIn = load(dataFilename);
dataOut = load(dataOutFilename);
dataTime = load(dataTimeFilename);

%%%%%%%%%% cc results

pcDataIn = dataIn(1:2:end);
pcDataIn = pcDataIn;

trainSetMean = mean(pcDataIn);
trainSetFrac = 1/(max(pcDataIn)-trainSetMean);
pcDataIn = (pcDataIn - trainSetMean) * trainSetFrac;

x = 1:length(pcDataIn);
timeMean = mean(x);
timeFrac = 1/(max(x)-timeMean);
x=(x-timeMean) * timeFrac;

pcdata = [x;(pcDataIn)'];
%figure
%plot(pcdata(1,:),pcdata(2,:))
cov(pcdata')
[v d] = eig(cov(pcdata'))
[coeff, score, latent, tsquare] = princomp(pcdata');
figure(30);
subplot(2,1,2);
hold on;
plot([0 v(1,1)],[0 v(2,1)],'g-'); % or coeff
plot([0 v(1,2)],[0 v(2,2)],'g-');

subplot(2,1,1);
plot(x,pcDataIn);

%%%%%%%%%% Simulation results
indir = 'C:\CurrentProjects\Network\debug\';
resInFilename1 = 'Connection_0_n0.cvs';
resInFilename2 = 'Connection_0_n1.cvs';
[time resultCV1 sizeXY] = LoadWithTime('C:\CurrentProjects\Network\debug\Connection_0_n0.csv');%[indir resInFilename1]);
[time resultCV2 sizeXY] = LoadWithTime('C:\CurrentProjects\Network\debug\Connection_0_n1.csv');%[indir resInFilename2]);

resultCV1 = cell2mat(resultCV1);
resultCV2 = cell2mat(resultCV2);

startItem = 105*2+1;
resultCV1 = resultCV1(startItem:end);
resultCV2 = resultCV2(startItem:end);

figure(run);
clf;
subplot(2,1,1);
hold on;
plot(dataTime(1:2:end),dataIn(1:2:end),'r','MarkerSize',1);
plot(dataTime(2:2:end),dataIn(2:2:end),'b','MarkerSize',1);

plot(1:length(resultCV1)/4,resultCV1(1:4:end)','.-r','MarkerSize',10);
plot(2:length(resultCV2)/4+1,resultCV2(1:4:end)','.-b','MarkerSize',10);

subplot(2,1,2);
hold on;
plot([0 v(1,1)],[0 v(2,1)],'g-'); % or coeff
plot([0 v(1,2)],[0 v(2,2)],'g-');
