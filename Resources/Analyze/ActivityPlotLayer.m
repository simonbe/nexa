function ActivityPlotLayer(timeValues, activityValues, filenameCorrectActivities, startIndex)%filenameActivities, filenameCorrectActivities, startIndex)

%if(exist(filenameActivities))
    
    %activityValues = load(filenameActivities);
    activityValues = activityValues(startIndex:end,:);
    rawActivityValues = activityValues;
    
    correctExists = 0;
    if(exist(filenameCorrectActivities))
        outValues = load(filenameCorrectActivities);
        correctExists = 1;
    end
    
    if(correctExists == 1)
        for j=1:size(activityValues,1)
            index = find(activityValues(j,:)==1);
            activityValues(j,index) = outValues(j)+1;
        end
    end
    
    size(activityValues)
    
    activityValue(activityValues==0) = -1;
    
    subplot(2,1,1);
    imagesc(activityValues);
    colorbar;
    
    xlabel('Unit nr');
    ylabel('Time');
    %title(filenameActivities);
    %colorbar;

    % hist
    subplot(2,1,2);
    if(correctExists == 1)
        maxVal = max(outValues)
        stacked = [];
        for i=1:maxVal+1
            stacked = [stacked; sum(activityValues(:,:)==i)];
        end
        
        D = sum(stacked);
        M = max(stacked);
        
        s = 0;
        
        for i=1:length(D)
            if(D(i)>0)
                s = s + M(i)/D(i);
            end
        end
        
        tot = sum(M)/sum(D);
        MdivD = M/D
        meanSame = s/sum(D~=0);
        
        size(stacked);
        bar(stacked','stack');
        
        title(['Mean same activity = ' num2str(meanSame) ', Tot max correct = ' num2str(tot)]);
    else
        bar(sum(rawActivityValues(:,:))/size(rawActivityValues,1));%bar(sum(rawActivityValues(:,:)==1));
    end
    
    xlabel('Unit nr');
    ylabel('Nr active');
%end