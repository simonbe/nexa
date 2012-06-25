function [A count] = LoadBinary(filename, unixFile)

% nr values assumed ordered as - unit1: nr of time steps, unit 2: nr of
% time steps ...

if(unixFile == 1)
    fid = fopen(filename,'r','b');
else
    fid = fopen(filename,'r');
end

[A count] = fread(fid,1e9,'real*4');


fclose(fid);

