indir = '..\debug\benjamin\NetworkTests\';
close all
figure
subplot(2,1,1);
[timeData connectionData rawData] = LoadWithTime([indir 'Layer1ActivityWTA.csv']);
C=cell2mat(connectionData);
imagesc(C)
ylabel('units')
colorbar
title('layer activity before switching')
subplot(2,1,2);
[timeData connectionData rawData] = LoadWithTime([indir 'Layer1Activity.csv']);
C=cell2mat(connectionData);
imagesc(C)
xlabel('time')
ylabel('units')
colorbar
title('layer activity after WTA')
subplot(2,1,2);
%hgsave(1,'switching_normexp')
%hgsave(1,'switching_normWTA')