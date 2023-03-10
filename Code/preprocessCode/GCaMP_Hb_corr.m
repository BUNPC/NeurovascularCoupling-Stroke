%% Correlation between response magnitude of GCaMP and Hb

SS75 = '../MouseData/SS75/';
SS76 = '../MouseData/SS76/';
SS77 = '../MouseData/SS77/';
SS78 = '../MouseData/SS78/';
SS79 = '../MouseData/SS79/';
SS80 = '../MouseData/SS80/';
SS81 = '../MouseData/SS81/';
SS82 = '../MouseData/SS82/';
SS83 = '../MouseData/SS83/';
SS84 = '../MouseData/SS84/';
SS85 = '../MouseData/SS85/';
SS93 = '../MouseData/SS93/';
mouse = {SS75, SS76, SS77, SS78, SS79, SS80, SS81, SS82, SS83, SS84, SS85, SS93};
mousename = {'SS75', 'SS76', 'SS77', 'SS78', 'SS79', 'SS80', 'SS81', 'SS82', 'SS83', 'SS84', 'SS85', 'SS93'};
thresh = [0.65 0.75 0.85];
thGCaMP{12,5,20,3} = zeros(128);
thHbT{12,5,20,3} = zeros(128);
thHbO{12,5,20,3} = zeros(128);
thHbR{12,5,20,3} = zeros(128);

for m = 1:length(mouse)
    m
[timepoints, img, mask] = animals_contra(mousename{m});
load(mask)
for time = 1:5
    time
respmask = maskSFDI(time).unaff_mask;       
load([timepoints{time},'/','act_IOSI_ipsi.mat'])
load([timepoints{time},'/','act_GCaMPcorr_ipsi.mat'])

deltaGCaMPcorr = respmask.*deltaGCaMPcorr;
HbO = respmask.*HbO;
HbR = respmask.*HbR;
HbT = respmask.*HbT;
deltaGCaMPcorr(deltaGCaMPcorr==0) = NaN;
HbT(HbT==0) = NaN;
HbO(HbO==0) = NaN;
HbR(HbR==0) = NaN;


for t = 1:numTrials
    Greshape(t,:,:,:) = deltaGCaMPcorr(:,:,(numFrame/numTrials*t)-(numFrame/numTrials-1):numFrame/numTrials*t);
    HbTreshape(t,:,:,:) = HbT(:,:,(numFrame/numTrials*t)-(numFrame/numTrials-1):numFrame/numTrials*t);
    HbOreshape(t,:,:,:) = HbO(:,:,(numFrame/numTrials*t)-(numFrame/numTrials-1):numFrame/numTrials*t);
    HbRreshape(t,:,:,:) = HbR(:,:,(numFrame/numTrials*t)-(numFrame/numTrials-1):numFrame/numTrials*t);
end

rStart = baseFrames;
rEnd = 2*baseFrames;
Gresp = squeeze(mean(Greshape(:,:,:,rStart:rEnd),4));
HbTresp = squeeze(mean(HbTreshape(:,:,:,rStart:rEnd),4));
HbOresp = squeeze(mean(HbOreshape(:,:,:,rStart:rEnd),4));
HbRresp = squeeze(mean(HbRreshape(:,:,:,rStart:rEnd),4));

for trial = 1:20
    
    for th = 1:length(thresh)
    
    cutG = squeeze(Gresp(trial,:,:));
    peakG = prctile(cutG(:),99);
    cutG(cutG<=thresh(th)*peakG) = NaN;
    cutG(isnan(cutG)) = 0;
    mGCaMP(m,time,trial,th) = mean(nonzeros(cutG(:)));    

    cutHbT = squeeze(HbTresp(trial,:,:));
    peakHbT = prctile(cutHbT(:),99);
    cutHbT(cutHbT<=thresh(th)*peakHbT) = NaN;
    cutHbT(isnan(cutHbT)) = 0;
    mHbT(m,time,trial,th) = mean(nonzeros(cutHbT(:)));
    
    cutHbO = squeeze(HbOresp(trial,:,:));
    peakHbO = prctile(cutHbO(:),99);
    cutHbO(cutHbO<=thresh(th)*peakHbO) = NaN;
    cutHbO(isnan(cutHbO)) = 0;
    mHbO(m,time,trial,th) = mean(nonzeros(cutHbO(:)));
    
    cutHbR = squeeze(HbRresp(trial,:,:));
    peakHbR = prctile(cutHbR(:),1);
    cutHbR(cutHbR>=thresh(th)*peakHbR) = NaN;
    cutHbR(isnan(cutHbR)) = 0;
    mHbR(m,time,trial,th) = mean(nonzeros(cutHbR(:)));
    
    BWG = cutG;
    BWG(BWG>0)=1;
    BWG(BWG<0)=0;
    
    BWT = cutHbT;
    BWT(BWT>0)=1;
    BWT(BWT<0)=0;
    
    BWO = cutHbO;
    BWO(BWO>0)=1;
    BWO(BWO<0)=0;
    
    BWR = cutHbR;
    BWR(BWR>0)=0;
    BWR(BWR<0)=1;
    
    if isempty(dice(BWG,BWT)) == 1
        simT(m,time,trial,th) = 0;
    else
        simT(m,time,trial,th) = dice(BWG,BWT);
    end
    if isempty(dice(BWG,BWO)) == 1
        simO(m,time,trial,th) = 0;
    else
        simO(m,time,trial,th) = dice(BWG,BWO);
    end
    if isempty(dice(BWG,BWR)) == 1
        simR(m,time,trial,th) = 0;
    else
        simR(m,time,trial,th) = dice(BWG,BWR);
    end

    
    areaHbT(m,time,trial,th) = length(find(cutHbT~=0));
    areaHbO(m,time,trial,th) = length(find(cutHbO~=0));
    areaHbR(m,time,trial,th) = length(find(cutHbR~=0));
    areaG(m,time,trial,th) = length(find(cutG~=0));
        
    thGCaMP{m,time,trial,th} = cutG;
    thHbT{m,time,trial,th} = cutHbT;
    thHbO{m,time,trial,th} = cutHbO;
    thHbR{m,time,trial,th} = cutHbR;
    
    end
end
end
end
mGCaMP(isnan(mGCaMP)) = 0;
mGCaMP = 100.*mGCaMP;
mHbT = 1e6.*mHbT;
mHbO = 1e6.*mHbO;
mHbR = 1e6.*mHbR;
areaG = 0.0027.*areaG; % convert to pixel area in mm^2
areaHbO = 0.0027.*areaHbO;
areaHbR = 0.0027.*areaHbR;
areaHbT = 0.0027.*areaHbT;

save('GCaMP_Hb_corr_unimpUnaff.mat','mGCaMP','mHbT','mHbO','mHbR','BWG','BWT','BWO','BWR','simT','simO','simR',...
    'areaHbT','areaHbO','areaHbR','areaG','thGCaMP','thHbT','thHbO','thHbR','-v7.3')

