A = load('C:\CurrentProjects\Network\Databases\Olfaction\mtF_sdbl.dat');

dlmwrite('mtF_sdbl.csv',A,'delimiter',',','precision','%.6f');