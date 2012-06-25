function A = MakeSparseMatrix(data)

A = zeros(5880,5880);%sparse(0);

for i=1:3:length(data)
    A(data(i)+1,data(i+1)+1) = data(i+2);
end