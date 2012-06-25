%This m-file is a short example of how to read and write a binary file 
%in MATLAB using low level routines

%write file 
%fid = fopen('square_mat.bin','wb')   %open file in binary write mode
%fwrite(fid, 'This is a square matrix', 'char');   %insert a header
%fwrite(fid, [1 2 3; 4 5 6; 7 8 9]', 'int32');   %write data to file
%fclose(fid);   %close file


%read in the same file
fid = fopen('DistanceMatrixfMRI.dat','rb');   %open file
[fname, mode, mformat] = fopen(fid);
mformat
x = fread(fid, 100, 'char', 'ieee-le')
dataFirst = fread(fid,3,'float64')
%nrItems = fread(fid, 20, 'short')   %read in the header
%nrHypercolumns = fread(fid, 1, 'int32')
%fread(fid, 1, 'int32')
%fread(fid, 1, 'int32')
%fread(fid, 1, 'int32')
%title = char(bintitle')   
%data = fread(fid, [3 inf], 'int32')  %read in the data
%data_tranpose = data'   %must transpose data after reading in
fclose(fid);   %close file

