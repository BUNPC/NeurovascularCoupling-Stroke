%% Rewrite the RSFC code as a script to run all animals together

%% Forelimb connectivity select forelimb position

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
mouse = {SS75, SS76, SS77, SS78, SS79, SS80, SS81, SS82, SS83, SS84, SS85 SS93};
mousename = {'SS75', 'SS76', 'SS77', 'SS78', 'SS79', 'SS80', 'SS81', 'SS82', 'SS83', 'SS84', 'SS85', 'SS93'};

% Find forelimb area from evoked GCaMP response map pre-stroke, find center
% and then take 0.25mm around that pixel (so +-5 pixels on each side from
% center pixel.
% figure
for m = 1:length(mouse)
    [timepoints, img, mask] = animals(mousename{m});
    load(mask)
    for t = 1:5
        load([timepoints{t},'/','act_IOSI_ipsi.mat'])
        load([timepoints{t},'/','act_GCaMPcorr_ipsi.mat'])
        % Trial average hemodynamics and GCaMP (corrected and uncorrected)
        HbOavg = trialAverage(HbO, numTrials, numFrame);
        HbRavg = trialAverage(HbR, numTrials, numFrame);
        HbTavg = trialAverage(HbT, numTrials, numFrame);
        GCaMPcorr = trialAverage(deltaGCaMPcorr, numTrials, numFrame);
        rStart = baseFrames;
        rEnd = 2*baseFrames;
        respmask = maskSFDI(t).aff_mask;
        GCaMPcorrresp = respmask.*squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
        GCaMPcorrresp(isnan(GCaMPcorrresp)) = 0;
        cutGCaMPcorrresp = GCaMPcorrresp; 
        peakG = prctile(cutGCaMPcorrresp(:),99);
        cutGCaMPcorrresp(cutGCaMPcorrresp<=0.75*peakG) = 0;
        cutGCaMPcorrresp(cutGCaMPcorrresp>0.75*peakG) = 1;
        stats_ipsi(m,t) = regionprops(cutGCaMPcorrresp,'basic');
        
        % Select bregma and lambda
%         figure
%         imagesc(img(t).img)
%         axis image
%         colormap gray
%         [breg_y{m,t},breg_x{m,t}] = ginput(2);
        
%         subplot(3,4,m)
%         imagesc(cutGCaMPcorrresp)
%         axis image
%         axis off
%         hold on
%         plot(stats(m).Centroid(1),stats(m).Centroid(2),'o')
    end
end
     
for m = 1:length(mouse)
    [timepoints, img, mask] = animals_contra(mousename{m});
    load(mask)
    for t = 1:5
        load([timepoints{t},'/','act_IOSI_ipsi.mat'])
        load([timepoints{t},'/','act_GCaMPcorr_ipsi.mat'])
        % Trial average hemodynamics and GCaMP (corrected and uncorrected)
        HbOavg = trialAverage(HbO, numTrials, numFrame);
        HbRavg = trialAverage(HbR, numTrials, numFrame);
        HbTavg = trialAverage(HbT, numTrials, numFrame);
        GCaMPcorr = trialAverage(deltaGCaMPcorr, numTrials, numFrame);
        rStart = baseFrames;
        rEnd = 2*baseFrames;
        respmask = maskSFDI(t).unaff_mask;
        GCaMPcorrresp = respmask.*squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
        GCaMPcorrresp(isnan(GCaMPcorrresp)) = 0;
        cutGCaMPcorrresp = GCaMPcorrresp; 
        peakG = prctile(cutGCaMPcorrresp(:),99);
        cutGCaMPcorrresp(cutGCaMPcorrresp<=0.75*peakG) = 0;
        cutGCaMPcorrresp(cutGCaMPcorrresp>0.75*peakG) = 1;
        stats_contra(m,t) = regionprops(cutGCaMPcorrresp,'basic');
    end
end

save('bregma_lambda.mat','breg_x','breg_y','stats_ipsi','stats_contra','mouse')

%% Calculate connectivity metrics for all mice

load('bregma_lambda.mat')
for m = 1:length(mouse)
m
baseline = [mouse{m},'RestingState/Baseline'];
day2 = [mouse{m},'RestingState/Day2'];
week1 = [mouse{m},'RestingState/Week1'];
week2 = [mouse{m},'RestingState/Week2'];
week4 = [mouse{m},'RestingState/Week4'];
mask = [mouse{m},'brainmaskSFDI.mat'];
pathname = {baseline, day2, week1, week2, week4};
for t = 1:length(pathname)
    t
    dataDir = pathname{t};   
    load(mask)
    load([dataDir,'/','gcampcorr_v2_highrsfc_I.mat'])
    [nY,nX,nLam,nT] = size(I);
    brain_mask = maskSFDI(t).aff_mask + maskSFDI(t).unaff_mask;
    brainIdx = find(brain_mask > 0);
    if t == 1 || t == 2
        % Global signal regression
        yAll = reshape(I(:,:,1,:),[nY*nX,nT]);
        y1 = yAll(brainIdx,:);
        yMean = mean(y1,1);
        a = pinv(yMean*yMean')*yMean*y1';
        ynew = y1'-yMean'*a;
    else
        % Multi region signal regression
        infarctIdx = find(maskSFDI(t).stroke_mask == 1);
        yAll = reshape(I(:,:,1,:),[nY*nX,nT]);
        y1 = yAll(brainIdx,:);
        yinfarct = yAll(infarctIdx,:);
        yinfarctMean = mean(yinfarct,1);
        noninfarctIdx = setdiff(brainIdx,infarctIdx);
        ynoninfarct = yAll(noninfarctIdx,:);
        ynoninfarctMean = mean(ynoninfarct,1);
        yMean = [ynoninfarctMean; yinfarctMean];
        a = pinv(yMean*yMean')*yMean*y1';
        ynew = y1'-yMean'*a;
    end
    Chb1 = zeros(nY*nX,nT);
    Chb1(brainIdx,:) = ynew';
    Chb1 = reshape(Chb1, [nY, nX, nT]);
    deltaGCaMPcorr = Chb1;

    
    load([dataDir,'\','1_v2_highrsfc_I.mat'])   
    if t == 1 || t == 2
        % Global signal regression
        yAll = reshape(I(:,:,1,:),[nY*nX,nT]);
        y1 = yAll(brainIdx,:);
        yMean = mean(y1,1);
        a = pinv(yMean*yMean')*yMean*y1';
        ynew = y1'-yMean'*a;
    else
        % Multi region signal regression
        infarctIdx = find(maskSFDI(t).stroke_mask == 1);
        yAll = reshape(I(:,:,1,:),[nY*nX,nT]);
        y1 = yAll(brainIdx,:);
        yinfarct = yAll(infarctIdx,:);
        yinfarctMean = mean(yinfarct,1);
        noninfarctIdx = setdiff(brainIdx,infarctIdx);
        ynoninfarct = yAll(noninfarctIdx,:);
        ynoninfarctMean = mean(ynoninfarct,1);
        yMean = [ynoninfarctMean; yinfarctMean];
        a = pinv(yMean*yMean')*yMean*y1';
        ynew = y1'-yMean'*a;
    end
    Chb1 = zeros(nY*nX,nT);
    Chb1(brainIdx,:) = ynew';
    Chb1 = reshape(Chb1, [nY, nX, nT]);
    HbO = Chb1;
    
    
    % Calculate seed-based analysis for forelimb seed
    pSeed = round(stats_ipsi(m,t).Centroid);
    yG = deltaGCaMPcorr;
    h = 10; %Seed dimension 
    ySeed = squeeze(mean(mean(yG(pSeed(2)+[-floor(h/2):floor(h/2)],pSeed(1)+[-floor(h/2):floor(h/2)],:),1),2));
    yG = reshape(yG,[nY*nX nT]);
    yG = yG(brainIdx,:);
    yCCG = corr(yG.',ySeed);
    Chb1CC = zeros(nY,nX);
    Chb1CC(brainIdx) = yCCG;
    GCaMP_forelimb_new{m,t} = Chb1CC;        
    pSeed = round(stats_ipsi(m,t).Centroid);
    yH = HbO;
    h = 10; %Seed dimension 
    ySeed = squeeze(mean(mean(yH(pSeed(2)+[-floor(h/2):floor(h/2)],pSeed(1)+[-floor(h/2):floor(h/2)],:),1),2));
    yH = reshape(yH,[nY*nX nT]);
    yH = yH(brainIdx,:);
    yCCH = corr(yH.',ySeed);
    Chb1CC = zeros(nY,nX);
    Chb1CC(brainIdx) = yCCH;
    HbO_forelimb_new{m,t} = Chb1CC;    
    
    pSeed = round(stats_ipsi(m,1).Centroid);
    yG = deltaGCaMPcorr;
    h = 10; %Seed dimension 
    ySeed = squeeze(mean(mean(yG(pSeed(2)+[-floor(h/2):floor(h/2)],pSeed(1)+[-floor(h/2):floor(h/2)],:),1),2));
    yG = reshape(yG,[nY*nX nT]);
    yG = yG(brainIdx,:);
    yCCG = corr(yG.',ySeed);
    Chb1CC = zeros(nY,nX);
    Chb1CC(brainIdx) = yCCG;
    GCaMP_forelimb_conn{m,t} = Chb1CC;        
    pSeed = round(stats_ipsi(m,1).Centroid);
    yH = HbO;
    h = 10; %Seed dimension 
    ySeed = squeeze(mean(mean(yH(pSeed(2)+[-floor(h/2):floor(h/2)],pSeed(1)+[-floor(h/2):floor(h/2)],:),1),2));
    yH = reshape(yH,[nY*nX nT]);
    yH = yH(brainIdx,:);
    yCCH = corr(yH.',ySeed);
    Chb1CC = zeros(nY,nX);
    Chb1CC(brainIdx) = yCCH;
    HbO_forelimb_conn{m,t} = Chb1CC;    
    
    % Calculate seed-based analysis for contralateral forelimb
    pSeed = round(stats_contra(m,1).Centroid);
    yG = deltaGCaMPcorr;
    h = 10; %Seed dimension 
    ySeed = squeeze(mean(mean(yG(pSeed(2)+[-floor(h/2):floor(h/2)],pSeed(1)+[-floor(h/2):floor(h/2)],:),1),2));
    yG = reshape(yG,[nY*nX nT]);
    yG = yG(brainIdx,:);
    yCCG = corr(yG.',ySeed);
    Chb1CC = zeros(nY,nX);
    Chb1CC(brainIdx) = yCCG;
    GCaMP_forelimb_conn_contra{m,t} = Chb1CC;        
    pSeed = round(stats_contra(m,1).Centroid);
    yH = HbO;
    h = 10; %Seed dimension 
    ySeed = squeeze(mean(mean(yH(pSeed(2)+[-floor(h/2):floor(h/2)],pSeed(1)+[-floor(h/2):floor(h/2)],:),1),2));
    yH = reshape(yH,[nY*nX nT]);
    yH = yH(brainIdx,:);
    yCCH = corr(yH.',ySeed);
    Chb1CC = zeros(nY,nX);
    Chb1CC(brainIdx) = yCCH;
    HbO_forelimb_conn_contra{m,t} = Chb1CC;   
    
    % Calculate global connectivity matrix
    yG = deltaGCaMPcorr;
    yG = reshape(yG,[nY*nX nT]);
    yG = yG(brainIdx,:);
    yCC = corr(yG.',yG.');
    yCCp = yCC;
    yCCp(yCC<=0) = nan;    
    yCCavg = mean(abs(yCCp),1,'omitnan');
    BOC = zeros(nY*nX,1);    
    BOC(brainIdx) = yCCavg';
    BOC_GCaMP{m,t} = reshape(BOC, [nY, nX]);

    yH = HbO;
    yH = reshape(yH,[nY*nX nT]);
    yH = yH(brainIdx,:);
    yCC = corr(yH.',yH.');
    yCCp = yCC;
    yCCp(yCC<=0) = nan;    
    yCCavg = mean(abs(yCCp),1,'omitnan');
    BOC = zeros(nY*nX,1);    
    BOC(brainIdx) = yCCavg';
    BOC_HbO{m,t} = reshape(BOC, [nY, nX]);
    
    
    % Calculate interhemispheric connectivity
    h = 1; % seed size
    a = breg_y{m,t}(2)-breg_y{m,t}(1);
    b = breg_x{m,t}(1)-breg_x{m,t}(2);
    c = (breg_y{m,t}(1)-breg_y{m,t}(2))*breg_x{m,t}(1)+(breg_x{m,t}(2)-breg_x{m,t}(1))*breg_y{m,t}(1);
    yLRCG = zeros(nY,nX);
    yLRCH = zeros(nY,nX);
    CRG = zeros(nY,nX);
    CRH = zeros(nY,nX);
    for u = 1:nY
        for v = 1:nX
            if brain_mask(u,v) == 1
                xs = round(2*(b^2*u-a*b*v-a*c)/(a^2+b^2)-u);
                ys = round(2*(a^2*v-a*b*u-b*c)/(a^2+b^2)-v);
                if xs <1 
                    xs = 1;
                end
                if xs > nY
                    xs = nY;
                end
                if ys <1 
                    ys = 1;
                end
                if ys > nX
                    ys = nX;
                end
                CRG = corr(squeeze(deltaGCaMPcorr(xs,ys,:)),squeeze(deltaGCaMPcorr(u,v,:)));
                yLRCG(xs,ys) = CRG;
                yLRCG(u,v) = CRG;
                CRH = corr(squeeze(HbO(xs,ys,:)),squeeze(HbO(u,v,:)));
                yLRCH(xs,ys) = CRH;
                yLRCH(u,v) = CRH;
            end
        end 
    end
    LRC_GCaMP{m,t} = yLRCG;
    LRC_HbO{m,t} = yLRCH;
    
end
end

save('RSFC_analysis_highHz.mat','deltaGCaMPcorr','HbO','GCaMP_forelimb_new','HbO_forelimb_new',...
    'GCaMP_forelimb_conn','HbO_forelimb_conn','GCaMP_forelimb_conn_contra','HbO_forelimb_conn_contra',...
    'BOC_GCaMP','BOC_HbO','LRC_GCaMP','LRC_HbO','-v7.3')


%% Correlations between GCaMP and HbO

load('RSFC_analysis_lowHz.mat')

for m = 1:12
    for t = 1:5
        forelimb_new_corr{m,t} = corr(GCaMP_forelimb_new{m,t},HbO_forelimb_new{m,t});
    end
end







