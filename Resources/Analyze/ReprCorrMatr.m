function ReprCorrMatr(activityValues)


nrItems = size(activityValues,1);
A = corr(activityValues',activityValues');%zeros(nrItems,nrItems);

imagesc(A);
