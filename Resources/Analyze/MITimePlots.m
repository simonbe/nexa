indir = 'C:\\CurrentProjects\\Network\\debug\\mds';

files = dir(indir);

sizeData = 25;

n=1;
TotData = [];

firstRun = 1;

for i=1:length(files)
    [path name ext] = fileparts(files(i).name);

    if(strcmp(ext,'.csv') == 1)
        
        A = load([indir '\\' name ext]);
        
        timesteps = length(A)/sizeData;
        
        if(firstRun == 1)
            TotData = zeros(timesteps,sizeData,sizeData);
            firstRun = 0;
        end
        
        for j=1:timesteps
             TotData(j,n,:) = A((j-1)*sizeData+1:j*sizeData);
        end
        
        n=n+1;
    end
end