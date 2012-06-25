function MDSTimePlots(mdsData, sizeData)
%mdsData = load('C:\\CurrentProjects\\Network\\debug\\mds_0_2.csv');
%mdsData = load(filename);

%sizeData = 25;
dims = size(mdsData,2);
timeSteps = size(mdsData,1)/sizeData;

TimeData = zeros(timeSteps,sizeData,dims);

for i=1:timeSteps
    TimeData(i,:,:) = mdsData((i-1)*sizeData+1:i*sizeData,:);
end

figure;

for i=1:timeSteps
    
    A = squeeze(TimeData(i,:,:));
    
    if(mod(i,1) == 0)
        scatter(A(:,1),A(:,2));
        
        %refresh plot
        %axis([-0.1 0.1 -0.1 0.2])
        title(['T = ' num2str(i)])
        %refresh(3)
        %pause(0.0005)
        drawnow
    end
end
