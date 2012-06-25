% Loads csv data/time separately
function [time data sizeXY] = LoadWithTimeV2(filename, dataSizeX, dataSizeY)

time = [];
data = [];

fid = fopen(filename);

tline = fgetl(fid);
index = 1;
while ischar(tline)
    
    f = strread(tline,'%f','delimiter',',');
    data{index} = f(2:end);
    time = [time f(1)];

    tline = fgetl(fid);
    index = index+1;
end

fclose(fid);
sizeXY = 0;
%data = dlmread(filename);
%raw = data;
% 
% if nargin == 1
%     dataSizeX = size(data,2);
% end
% 
% nrZero = sum(data'==0);
% indexesTime = find(nrZero == dataSizeX-1); % only one element neq 0
% 
% time = data(indexesTime,1);
% data(indexesTime,:) = [];
% 
% if nargin == 1
%     if(length(time)>1)
%         dataSizeY = indexesTime(2) - indexesTime(1) - 1;
%     else
%         dataSizeY = size(data,1) - 1;
%     end
% end
% 
% sizeXY = [dataSizeX dataSizeY]