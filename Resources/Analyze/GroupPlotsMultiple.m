clear all;

groupIndexes0 = []; groupIndexes1 = []; groupIndexes2 = []; 
groupIndexes0 = LoadGrouping('C:\\CurrentProjects\\Network\\debug\\vqGroups_0.csv');
groupIndexes1 = LoadGrouping('C:\\CurrentProjects\\Network\\debug\\vqGroups_1.csv');
groupIndexes2 = LoadGrouping('C:\\CurrentProjects\\Network\\debug\\vqGroups_2.csv');

sizeX=28;
sizeY=28;

img = zeros(sizeX,sizeY);
img2 = zeros(sizeX,sizeY);

%% layer 1
for i=1:length(groupIndexes0)
    
    for j=1:length(groupIndexes0{i})
        val = groupIndexes0{i}(j);
        y = floor(val/sizeY)+1;
        x = mod(val,sizeX);
        if(x==0)
            x=sizeX;
        end
        
        img(x,y) = i;
    end
end

figure(14);
imagesc(img);

%% layer 2

for i=1:length(groupIndexes1)
    
    indexes = [];
    for j=1:length(groupIndexes1{i})
        
        n = groupIndexes0{groupIndexes1{i}(j)+1};
        for k=1:length(n)
            indexes = [indexes n(k)];
        end
        
    end
    
    groupIndexesNew{i} = indexes;
end

for i=1:length(groupIndexesNew)
    
    for j=1:length(groupIndexesNew{i})
        val = groupIndexesNew{i}(j);
        y = floor(val/sizeY)+1;
        x = mod(val,sizeX);
        if(x==0)
            x=sizeX;
        end
        
        img2(x,y) = i;
    end
end

figure(15);
imagesc(img2);


%% layer 3

if(length(groupIndexes2)>0)
    
    for i=1:length(groupIndexes2)
        
        indexes = [];
        for j=1:length(groupIndexes2{i})
            
            n = groupIndexesNew{groupIndexes2{i}(j)+1};
            for k=1:length(n)
                indexes = [indexes n(k)];
            end
            
        end
        
        groupIndexesNew2{i} = indexes;
    end
    
    for i=1:length(groupIndexesNew2)
        
        for j=1:length(groupIndexesNew2{i})
            val = groupIndexesNew2{i}(j);
            y = floor(val/sizeY)+1;
            x = mod(val,sizeX);
            if(x==0)
                x=sizeX;
            end
            
            img3(x,y) = i;
        end
    end
    
    figure(16);
    imagesc(img3);
end