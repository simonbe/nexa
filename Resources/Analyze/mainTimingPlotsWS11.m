mainTimingPlots

% add larger network run

directory2 = 'C:\Projects\Network\NetworkResults6\';%'C:\Projects\Network\Outputs\Timings\recurr\small\';

X2 = [];
Y2 = [];

for i=1:length(A)
    filename = [directory2 'NetworkTimings' num2str(A(i)) '.txt'];
    if(exist(filename))
        X2 = [X2 A(i)];
        [labels2 data2] = LoadTimingInformation(filename);
        Y2 = [Y2 data2];
        
    end
end

YWS11_2 = Y2(4,:) + Y2(2,:);
YscaledWS11_2 = 1./(YWS11_2./repmat(YWS11_2(:,1),1,size(YWS11_2,2)))
figure(78);
clf;
subplot(1,2,1);
hold on;
plot(X2,X2./X2(1),'-.b');
plot(X2,YscaledWS11_2,'rx-');
ylabel('Speedup');
xlabel('Processes');
title('Medium-sized network - 8192,16384, 32768 and 65536 processes');
legend('linear','6.6x10^6 minicolumns, 1.3x10^{10} synapses');

figure(78);
subplot(1,2,2);
hold on;
InitTime = Y2(2,:);
PlastTime = Y2(3,:);
CommTime = Y2(5,:);
UnitsTime = Y2(9,:);
EventsTime = Y2(7,:);
TotTime = InitTime+PlastTime+CommTime+UnitsTime+EventsTime;
plot(X2,100*InitTime./TotTime,'.-');
plot(X2,100*PlastTime./TotTime,'-o');
plot(X2,100*CommTime./TotTime,'-x');
plot(X2,100*UnitsTime./TotTime,'-*');
plot(X2,100*EventsTime./TotTime,'-.');
ylabel('% Time spent');
xlabel('Processes');
legend('Building network','Plasticity','MPI Communication','Simulating neural units','Population operations (WTA)');
title('Fraction time spent in parts of simulation');