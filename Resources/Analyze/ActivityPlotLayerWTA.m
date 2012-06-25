function histCorrect = ActivityPlotLayerWTA(filenameActivities, filenameCorrectActivities, startIndex)

if(exist(filenameActivities))
    
    activityValues = load(filenameActivities);
    activityValues = activityValues(startIndex:end,:);
    histCorrect = zeros(size(activityValues,1),1);
    resultValues = zeros(size(activityValues));
    
    correctExists = 0;
    if(exist(filenameCorrectActivities))
        outValues = load(filenameCorrectActivities);
        correctExists = 1;
    end

    totCorrect = 0;
    if(correctExists == 1)
        for j=1:size(activityValues,1)
            [val index] = max(activityValues(j,:));
            resultValues(j,index) = outValues(j)+1;
            
            if(index == outValues(j)+1)
                totCorrect=totCorrect+1;
            end
            
            histCorrect(j) = totCorrect/j;
        end
    end
    
    subplot(2,1,1);
    imagesc(resultValues');
    xlabel('Time');
    ylabel('Unit nr');
    colorbar;
    
    subplot(2,1,2);
    plot(histCorrect);
    
    xlabel('Time');
    ylabel('Correct');
    
    
    
    
end