function A = LoadTriangularMatrix(filename)

raw = textread(filename,'%s,%s,%s,%s,%s','bufsize', 1024*500-1);%textread(filename,'%s,%s,%s,%s,%s','bufsize', 1024*500-1);

l = length(str2num(raw{length(raw)}));
A = zeros(length(raw),l);

for i=1:length(raw)
    data = str2num(raw{i});
    k = length(data);
    A(i,1+(i-k):length(data)+(i-k)) = data;
    A(1+(i-k):length(data)+(i-k),i) = data;
end