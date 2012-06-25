X = repmat(spikingLayer1(:,1),size(spikingLayer1,2));

figure(70);
hold on;
for i=1:size(spikingLayer1,2)-1
    scatter(spikingLayer1(:,1),repmat(i,length(spikingLayer1(:,1)),1), spikingLayer1(:,1+i));
end