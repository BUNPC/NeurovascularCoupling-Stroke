function avg = trialAverage(data, numTrials, numFrame)

datareshape = zeros([numTrials,size(data,1),size(data,2),numFrame/numTrials]);
for t = 1:numTrials
    datareshape(t,:,:,:) = data(:,:,(numFrame/numTrials*t)-(numFrame/numTrials-1):numFrame/numTrials*t);
end
avg = zeros([size(data,1),size(data,2),numFrame/numTrials]);
for i = 1:numFrame/numTrials 
    avg(:,:,i) = mean(squeeze(datareshape(:,:,:,i)),1);  
end

end