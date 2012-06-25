function GroupPlots(filename)
%groupIndexes = LoadGrouping('C:\\CurrentProjects\\Network\\debug\\vqGroups_0.csv');
groupIndexes = LoadGrouping(filename);

sizeX=28;
sizeY=28;

img = zeros(sizeX,sizeY);

for i=1:length(groupIndexes)
    
    for j=1:length(groupIndexes{i})
        val = groupIndexes{i}(j);
        y = floor(val/sizeY)+1;
        x = mod(val,sizeX);
        if(x==0)
            x=sizeX;
        end
        
        img(x,y) = i;
    end
end

%figure(14);
imagesc(img);