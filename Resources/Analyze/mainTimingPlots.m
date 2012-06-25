% scaling plots - will re

A = [16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144];
% if exist include

% WS11 scaling plots
ws11 = 1;

directory = 'C:\CurrentProjects\Network\Output\JUGENE\';%'C:\CurrentProjects\Network\Output\Olfaction_Pawel\';%'C:\Projects\Network\NetworkResults5\';%'C:\Projects\Network\Outputs\Timings\recurr\small\';

X = [];
Y = [];

for i=1:length(A)
    filename = [directory 'NetworkTimings' num2str(A(i)) '_rec21.txt'];
    if(exist(filename))
        X = [X A(i)];
        [labels data] = LoadTimingInformation(filename);
        Y = [Y data];
        
    end
end


% plot all of the timings
figure(75);
clf;
hold on;
subplot(2,1,1)
plot(X,Y);
%for j=1:size(Y,1)
%    plot(X,Y(j,:));
%end

legend(labels)
subplot(2,1,2)
hold on;
Yscaled = 1./(Y./repmat(Y(:,1),1,size(Y,2)))
plot(X,Yscaled)
plot(X,X./X(1),'-.b')
legend(labels)

n = 9;
% plot the n largest timings
figure(76);
clf;
hold on;
[values,indexes] = sort(Y,1,'descend')
YLarge = Y(indexes(1:n,1),:);
labelsLarge = labels(indexes(1:n,1),:);
subplot(2,1,1);
plot(X,YLarge);
legend(labelsLarge);
subplot(2,1,2)
hold on;
YscaledLarge = 1./(YLarge./repmat(YLarge(:,1),1,size(YLarge,2)))
plot(X,YscaledLarge')
plot(X,X./X(1),'-.b')
legend(labelsLarge)

%WS11 plots
if(ws11 == 1)
    
    YWS11 = Y(4,:) + Y(2,:);
    YscaledWS11 = 1./(YWS11./repmat(YWS11(:,1),1,size(YWS11,2)))
    figure(77);
    clf;
    subplot(1,2,1);
    hold on;
    plot(X,X./X(1),'-.b');
    plot(X,YscaledWS11,'r+-');
    ylabel('Speedup');
    xlabel('Processes');
    title('Memory storage and retrieval simulation');%'Medium-sized network - 65536, 131072, 262144 processes');
    legend('linear','6.6x10^7 minicolumns, 3.9x10^{11} synapses');
    
    subplot(1,2,2);
    hold on;
    InitTime = Y(2,:);
    AnalysisTime = Y(1,:);
    PlastTime = Y(3,:);
    CommTime = Y(6,:);
    UnitsTime = Y(14,:);
    EventsTime = Y(8,:);
    TotTime = InitTime+PlastTime+CommTime+UnitsTime+EventsTime+AnalysisTime;
    plot(X,100*InitTime./TotTime,'.-');
    plot(X,100*PlastTime./TotTime,'-o');
    plot(X,100*CommTime./TotTime,'-x');
    plot(X,100*UnitsTime./TotTime,'-*');
    plot(X,100*EventsTime./TotTime,'-.');
    plot(X,100*AnalysisTime./TotTime,'-+');
    ylabel('% Time spent');
    xlabel('Processes');
    legend('Building network','Plasticity','MPI Communication','Simulating neural units','Population operations (WTA)','Analysis');
    title('Fraction time spent in parts of simulation');
    
end