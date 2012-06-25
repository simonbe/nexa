indir = 'C:\CurrentProjects\Network\Analyze\Activities_set1_1\';
files = dir(indir);

for i=1:length(files)
    [path name ext] = fileparts(files(i).name);

    if(strcmp(ext,'.csv') == 1)
        disp(name)
        [time data] = LoadWithTime([indir name '.csv']);
    end
end