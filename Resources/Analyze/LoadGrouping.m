function data = LoadGrouping(filename)

fid = fopen(filename);

tline = fgets(fid);

data = []
i=1;

while ischar(tline) 
    parts = str2double(regexp(tline,',','split'));
    data{i} = parts;
    
    tline = fgets(fid);
    i=i+1;
end

fclose(fid);