function GroupTimePlot2d(vqData, sizeData, background)

dims = size(vqData,2);
timeSteps = dims/sizeData;

figure;

for i=1:timeSteps
    
    if(mod(i,1) == 0)
        for j=1:sizeData
            index = (i-1)*sizeData+j;
            
            if(size(vqData{index})>0)
                plot(vqData{index},j,'.');
            end
            
            hold on;
        end
        
        if(length(background)>0)
            plot(background);
        end
        hold off;
        
        title(['T = ' num2str(i)])
        drawnow
    end
end
