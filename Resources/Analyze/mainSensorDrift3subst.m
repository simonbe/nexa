run = 1;
directory = 'c:\CurrentProjects\Network\debug\';
trainResults = load([directory 'trainResultsRun' num2str(run) '.csv']);
testResults = load([directory 'testResultsRun' num2str(run) '.csv']);
testResultsSets = load([directory 'testResultsSetsRun' num2str(run) '.csv']);

figure(31);
hold off;
plot(trainResults(trainResults(:,2)==0,1),trainResults(trainResults(:,2)==0,3))
hold on;
plot(testResults(testResults(:,2)==0,1),testResults(testResults(:,2)==0,3),'r')
title('Training + testing results NO Train adaptation');

figure(30);
hold off;
plot(trainResults(:,3));
hold on;
plot(testResults(:,3),'r');
title('Training + testing results ALL');

A=[1 0.960938 0.953488 0.875 0.953125 0.829457 0.953125 0.945736 0.84375 0.953125];
x=1:10;

figure(32);
plot(x,A);
xlabel('Time');
ylabel('Classification ratio');
title('V 2');

r1 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun1.csv');
r2 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun2.csv');
r2Sets = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsSetsRun2.csv');
r3 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun3.csv');
r4 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun4.csv');
r5 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun5.csv');
%r6 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun6.csv');
r7 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun7.csv');
r8 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun8.csv');
r9 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun9.csv');
r10 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun10.csv');
r11 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun11.csv');
r112 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun112.csv');
r12 = load('c:\CurrentProjects\Network\Output\SensorDrift\new\testResultsRun12.csv');

figure(33);
clf;

subplot(5,1,1);
plot(r1(:,2),r1(:,4))
legend(num2str(r1(1,1)))

subplot(5,1,2);
hold on;
D=r2(:,1)==3;
plot(r2(D,2),r2(D,4),'y')
A=r2(:,1)==6;
plot(r2(A,2),r2(A,4))
B=r2(:,1)==9;
plot(r2(B,2),r2(B,4),'magenta')
C=r2(:,1)==15;
plot(r2(C,2),r2(C,4),'black')
plot(r3(:,2),r3(:,4),'r')
plot(r4(:,2),r4(:,4),'g')
legend('Code vectors = 3','6','9','15','30','60');

subplot(5,1,3);
plot(r3(:,2),r3(:,4))

subplot(5,1,4);
plot(r4(:,2),r4(:,4))

subplot(5,1,5);
plot(r5(:,2),r5(:,4))

% continued

figure(34)
hold on;
D=r2(:,1)==3;
plot(r2(D,2),r2(D,4),'y')
A=r2(:,1)==6;
plot(r2(A,2),r2(A,4))
B=r2(:,1)==9;
plot(r2(B,2),r2(B,4),'magenta')
C=r2(:,1)==15;
plot(r2(C,2),r2(C,4),'black')
plot(r3(:,2),r3(:,4),'r')
plot(r4(:,2),r4(:,4),'g')
legend('# codevectors = 3','6','9','15','30','60');
xlabel('adaptation');
ylabel('correct');

figure(341);
clf;
hold on;
B=r2(:,1)==9;
plot(r2(B,2),100*r2(B,4),'blue')
C=r2(:,1)==15;
plot(r2(C,2),100*r2(C,4),'black')
%plot(r3(:,2),r3(:,4),'r')
plot(r4(:,2),100*r4(:,4),'red')
legend('Number of codevectors = 9','15','60');
title('Total classification rates for subset of data');
xlabel('Drift correction parameter');
ylabel('Classification rates (%)');
axis([0,0.35,60,100]);

figure(35);
clf;
hold on;
sets15=r2Sets(C,:);
sets9=r2Sets(B,:);
plot(sets9(1,:),'-.');
plot(sets9(5,:));
plot(sets15(1,:),'r-.');
plot(sets15(5,:),'r');
legend('# codevectors = 9, adaptation = 0', '# codevectors = 9, adaptation = 0.05','# codevectors = 15, adaptation = 0','# codevectors = 15, adaptation = 0.05');
ylabel('correct');
xlabel('time');
axis([0,10,0.6,1]);
figure(36);

subplot(2,1,1);
clf;
hold on;
plot(r11(1:61,3),r11(1:61,4));
%plot(r11(62:122,3),r11(62:122,4),'r-');
plot(r11(123:183,3),r11(123:183,4),'r-.');
%plot(r11(184:244,3),r11(184:244,4),'r');
plot(r11(245:305,3),r11(245:305,4),'g-');
%plot(r11(306:366,3),r11(306:366,4),'g-.');
plot(r11(367:427,3),r11(367:427,4),'g');
plot(r11(428:488,3),r11(428:488,4),'y');
plot(r11(489:549,3),r11(489:549,4),'y-.');
%legend('adaptation = 0','0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08')
legend('Drift correction = 0','0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08')
title('Prediction removal on input layer');


figure(361);

subplot(2,1,1);
clf;
hold on;
plot(r11(1:61,3),r11(1:61,4));
%plot(r11(62:122,3),r11(62:122,4),'r-');
plot(r11(123:183,3),r11(123:183,4),'r-');
%plot(r11(184:244,3),r11(184:244,4),'r');
%plot(r11(245:305,3),r11(245:305,4),'g-');
%plot(r11(306:366,3),r11(306:366,4),'g-.');
plot(r11(367:427,3),r11(367:427,4),'black');
%plot(r11(428:488,3),r11(428:488,4),'y');
plot(r11(489:549,3),r11(489:549,4),'g');
%legend('adaptation = 0','0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08')
legend('Drift correction = 0','0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08')
title('Prediction removal on input layer');
xlabel('Removal strength');
ylabel('Classification rates (%)');

figure(37);
clf;
hold on;
A=0:max(r2(C,2))/(length(r10(:,4))-1):max(r2(C,2));
plot(r2(C,2),100*r2(C,4));
plot(A,100*r10(:,4),'r');
legend('Number of codevectors = 15, standard drift correction','Sensor-dependent correction');
title('Sensor-dependent drift correction');
xlabel('Drift correction strength');
ylabel('Classification rates (%)');

% different strengths of input gains
%A = [];
%for j=1:length(r8)/20
%    for i=0:19
%        x = r8((j-1)*20+i,2);
%        y = r8((j-1)*20+i,4);
%        A(i,:,:) = [A(i,:,:),x,y];
%    end
%end

% different strengths of etavector