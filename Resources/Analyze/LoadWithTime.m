% Loads csv data/time separately
function [time data sizeXY] = LoadWithTime(filename, dataSizeX, dataSizeY)

time = [];
data = [];

fid = fopen(filename);

tline = fgetl(fid);
index = 1;
while ischar(tline)
    
    isTimeLine = 1;
    if(length(tline) == 0)
        isTimeLine = 0;
    elseif(strcmp(tline(1:1),'t') == 0)
        isTimeLine = 0;
    end
    
    if(isTimeLine == 1)%strcmp(tline(1:1),'t'))
        time = [time str2num(tline(3:end))];
    else
        f = strread(tline,'%f','delimiter',',');
        data{index} = f;
        index = index+1;
    end
    tline = fgetl(fid);
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