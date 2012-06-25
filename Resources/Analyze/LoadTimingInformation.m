function [labels data] = LoadTimingInformation(filename)

fid = fopen(filename);
C = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',200, 'delimiter',',');
fclose(fid);

labels = char(C{1});
data = str2num(char(C{3}));