indir = 'C:\CurrentProjects\Network\Databases\Olfaction\set2\';
files = dir(indir);

for i=1:length(files)
    [path name ext] = fileparts(files(i).name);

    if(strcmp(ext,'.dat') == 1)
        disp(name);
        A = load([indir name '.dat']);
        dlmwrite([name '.csv'] ,A,'delimiter',',','precision','%.6f');
    end
end