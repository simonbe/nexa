function [nrNodes nrLayers NrHypercolumns NrMinicolumns] = LoadNetworkDetails(filename)

fid = fopen(filename)
C = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s',100, 'delimiter',' ');
fclose(fid);

% computer details
nrNodes = str2num(char(C{2}(2)));

% layers details
nrLayers = str2num(char(C{2}(3)));
pos = 4;

NrHypercolumns = [];
NrMinicolumns = [];

for i=1:nrLayers
    id = str2num(char(C{2}(pos)));
    type = str2num(char(C{4}(pos)));
    nrUnits = str2num(char(C{6}(pos)));
    nrHcs = str2num(char(C{8}(pos)));
    nrMcs = [];
    for j=1:nrHcs
        nrMcs=[nrMcs str2num(char(C{9+j}(pos)))];
    end
    
    NrHypercolumns = [NrHypercolumns nrHcs];
    NrMinicolumns = [NrMinicolumns nrMcs];
    pos = pos+1;
end

t=1