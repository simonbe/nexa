function [labels data] = LoadTimingInformation2(searchstring,xvalues)
% Given a searchstring, takes all files that match the string, sorts them 
% in ascending order, and gives out label and data information about 
% profiling (output from Timings functions in Network). This assumes that
% profiling information uses the same profiling parameters (it could run
% different sizes, e.g.). 
% The optional second parameter indicates if you want a plot to be generated. 
% (off by default). It should give a vector with values for the changed
% parameter for each profile. 
% Example: 
% [labels,data]=LoadTimingInformation2('..\Debug\benjamin\NetworkProjects\p
% rofiling\NetworkTimings*',[5,50,200]);

% some variables
linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':',...
'-.','--','-',':','-.','--','-',':','-.','-',':','-.','--','-',':','-.','--','-',':','-',':'));

MarkerEdgeColors=['r','c','b','k','g','m','y',...
'b','m','k','c','g','r','b','m','k','c','g',...
'r','c','b','k','g','m','r','b','m','k','c','g','r','r','c','b','k','g','m','r','b','m','k','c'];
Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.',...
'+','*','o','x','^','<','h','.','>','p','s','d','v',...
'o','x','+','*','s','d','v','^','<','>','p','h','.','+','*','o','x','^','<','h','.','>'];

[pathstr, name, ext] = fileparts(searchstring);

if length(pathstr)>0 
% the path is relative, but misses a ., at least for cases I looked at. 
% Might need tweaking
    concstr=[pathstr,'/',name,ext];
    files=dir(concstr);
else
    files=dir([name,ext]);
end

if numel(files)==0 
    disp('no files found')
    return
end

% you could convert files just doing filenames={files.name};
% however, I prefer making sure, there are no dirs in the structure
for i=1:numel(files) 
     if files(i).isdir==0
         if exist('filenames','var')
            filenames{numel(filenames)+1}=files(i).name;
         else
            filenames{1}=files(i).name;
         end         
     end
end

[filenames] = sort_nat(filenames);
filenames % to see if the order is right

data=[];
for i=1:numel(filenames)
    filename=filenames{i};
    [labels data2] = LoadTimingInformation([pathstr,'\',filename]);
    data=[data,data2];
end

if (nargin==2)
    figure;
    hold on
    for i=1:size(data,1)
        plot(xvalues,data(i,:),[linestyles{i} Markers(i) MarkerEdgeColors(i)]);
    end    
    plot(xvalues,data');    
    legend(labels,'Location','NorthWest');
    yLabel('Time');
end


