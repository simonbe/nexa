indir = 'C:\CurrentProjects\Network\Release\PostLazarus\NetworkProjects\';
%close all
figure
subplot(4,2,1);
[timeData connectionData rawData] = LoadWithTimeV2([indir 'ligands.csv']);
%connectionData{1}=zeros(size(connectionData{2}));
C=log(cell2mat(connectionData));
ligands=C;
plot(C');title('ligands');xlabel('time');ylabel('log ligand concentration');
subplot(4,2,2);
imagesc(C);colorbar;
subplot(4,2,3);
[timeData connectionData rawData] = LoadWithTimeV2([indir 'receptors.csv']);
C=cell2mat(connectionData);ors=C;
plot(ors');title('ORs');xlabel('time');ylabel('OR activation');
subplot(4,2,4);
imagesc(C);colorbar;
subplot(4,2,5);
[timeData connectionData rawData] = LoadWithTimeV2([indir 'orn_activations.csv']);
C=cell2mat(connectionData);orns=C;
plot(orns');title('ORNs');xlabel('time');ylabel('ORN frequency');
subplot(4,2,6);
imagesc(C);colorbar;
subplot(4,2,7);
[timeData connectionData rawData] = LoadWithTimeV2([indir 'mt_activations.csv']);
C=cell2mat(connectionData);mts=C;
plot(mts');title('MTs');xlabel('time');ylabel('MT frequency');
subplot(4,2,8);
imagesc(C);colorbar;
