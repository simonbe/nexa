function data = LoadGroupingParallel(filename)

fid = fopen(filename, 'r');
a = fscanf(fid, '%c', inf);

str = '';

Collection = [];
indexRow = 1;
indexColumn = 1;
for i=1:length(a)
    if(strcmp(a(i),',') || strcmp(a(i),';'))
        num = str2num(str);
        if(~isempty(num))
            Collection{indexRow}(indexColumn) = num;
            str = '';
            indexColumn = indexColumn+1;

            if(strcmp(a(i),';'))
                indexRow=indexRow+1;
                indexColumn = 1;
            end
        end
    else
        str = [str a(i)];
    end
end

data=Collection;