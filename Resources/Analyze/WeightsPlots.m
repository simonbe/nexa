weightValues = LoadGrouping('C:\\CurrentProjects\\Network\\debug\\Connection_triesch_n0.csv');

sizeX=10;
sizeY=10;
img = zeros(sizeX,sizeY);

for i=1:length(weightValues)
    figure(i);
    colormap(jet);

    index = 1;
    for j=1:1:length(weightValues{i})
        val = weightValues{i}(j);
        y = floor((index-1)/sizeY)+1;
        x = mod(index,sizeX);
        if(x==0)
            x=sizeX;
        end

        img(x,y) = val;
        index = index+1;
    end
    
    imagesc(img');
    colorbar
end

