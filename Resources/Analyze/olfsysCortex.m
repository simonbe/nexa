% analyze olfactory system to cortex
indir = '..\debug\benjamin\NetworkProjects\';
iterations=20;
figure

%% Ligands
subplot(4,2,1);
[timeData connectionData rawData] = LoadWithTimeV2([indir 'ligands.txt']);
C=log(cell2mat(connectionData));
ligands=C;
plot(C');title('ligands');xlabel('time');ylabel('log ligand concentration');
subplot(4,2,2);
imagesc(C);colorbar;
%% ORs
subplot(4,2,3);
[timeData connectionData rawData] = LoadWithTimeV2([indir 'receptors.txt']);
C=cell2mat(connectionData);ors=C;
plot(ors');title('ORs');xlabel('time');ylabel('OR activation');
subplot(4,2,4);
imagesc(C);colorbar;

%% ORNs
subplot(4,2,5);
[timeData connectionData rawData] = LoadWithTimeV2([indir 'orn_activations.txt']);
C=cell2mat(connectionData);
orns=C;
plot(orns');title('ORNs');xlabel('time');ylabel('ORN activation');
subplot(4,2,6);
imagesc(C);colorbar;

%% MTs
% subplot(4,2,5);
% [timeData connectionData rawData] = LoadWithTimeV2([indir 'mt_activations.txt']);
% C=cell2mat({connectionData{1:end/2}});mts=C;
% plot(mts');title('MTs');xlabel('time');ylabel('MT activation');title('MTs');xlabel('time');
% subplot(4,2,6);
% imagesc(C);colorbar;
% subplot(4,2,7);
% [timeData connectionData rawData] = LoadWithTimeV2([indir 'Layer2Activity_0.csv']);
% C=cell2mat({connectionData{1:end}});mts=C;
% plot(mts');title('MTs');xlabel('time');ylabel('layer1 activation');title('MTs');xlabel('time');
% subplot(4,2,8);
% imagesc(C);colorbar;

%% cortex
subplot(4,2,7);
[timeData connectionData rawData] = LoadWithTimeV2([indir 'mimdsvq_output.txt']);
C=cell2mat(connectionData);cortex=C;
plot(cortex');title('MIMDSVQ');xlabel('time');ylabel('MIMDSVQ activation');title('MIMDSVQ');xlabel('time');
subplot(4,2,8);
imagesc(C);colorbar;