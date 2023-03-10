%% Combine animals for neurovascular coupling analysis

% Load data
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

for m = 1:12
    m
baseline = [mouse{m},'FunctionalActivation/Baseline'];
day2 = [mouse{m},'FunctionalActivation/Day2'];
week1 = [mouse{m},'FunctionalActivation/Week1'];
week2 = [mouse{m},'FunctionalActivation/Week2'];
week4 = [mouse{m},'FunctionalActivation/Week4'];
pathname = {baseline, day2, week1, week2, week4};
time = 0:0.2:15-0.2;
mask1 = [mouse{m},'brainmaskSFDI.mat'];
load(mask1)
for t = 1:5
    t
    dataDir = pathname{t};
    load([dataDir,'/','act_HRF.mat'])    
    if t == 1
        HRF1 = reshape(HRF,[size(HRF,1)*size(HRF,2) size(HRF,3)]);
        mHRF1 = mean(HRF1,2);
        idx = find(mHRF1~=0 & ~isnan(mHRF1));
        HRF1 = HRF1(idx,:);
        refHRF = mean(HRF1,1);
    end
    l2normirf = zeros(size(HRF,1),size(HRF,2));
    for i = 1:size(HRF,1)
        for j = 1:size(HRF,2)
            sub = squeeze(refHRF'-squeeze(HRF(i,j,:)));
            l2normirf(i,j) = norm(sub);
        end
    end
    
    if t == 1 || t == 2
        newMask = imdilate(maskSFDI(3).stroke_mask, true(30));
        peri = abs(newMask - maskSFDI(3).stroke_mask);       
        stroke = maskSFDI(3).stroke_mask;
    else        
        newMask = imdilate(maskSFDI(t).stroke_mask, true(30));
        peri = abs(newMask - maskSFDI(t).stroke_mask);
        stroke = maskSFDI(t).stroke_mask;
    end
    forelimb = maskSFDI(1).ipsiOutline;
    contra = maskSFDI(1).contraOutline;
    
    for p = 1:3
        clear peakTime amplitudePos amplitudeNeg fwhm

        if p == 1
            newmask = stroke;
        elseif p == 2
            newmask = peri;
        elseif p == 3
            newmask = contra;
        elseif p == 4
            newmask = forelimb;
        end
        HRF_roi = HRF.*newmask;
        
        HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
        mHRF_roi = mean(HRF_roi,2);
        idx = find(mHRF_roi~=0 & ~isnan(mHRF_roi));
        HRF_roi = HRF_roi(idx,:);
        
    
        for h = 1:size(HRF_roi,1)
            idx = find(HRF_roi(h,3:27) == max(HRF_roi(h,3:27)));
            peakTime(h) = time(idx);
            amplitudePos(h) = max(HRF_roi(h,3:27));
            amplitudeNeg(h) = min(HRF_roi(h,8:27));
            halfMax = amplitudePos(h)/2;     
            if isempty(find(HRF_roi(h,3:27) >= halfMax, 1, 'first'))
                fwhm(h) = nan;
            else  
                index1 = find(HRF_roi(h,3:27) >= halfMax, 1, 'first');
                index2 = find(HRF_roi(h,3:27) >= halfMax, 1, 'last');
                fwhm(h) = time(index2) - time(index1);  
            end
        end
            
        newcorr = correlation.*newmask;        
        newcorr(newcorr==0) = NaN;
        corr_roi(p,t,m) = mean(newcorr(:),'omitnan');        
%         newL2norm = L2norm.*newmask;
%         newL2norm(newL2norm==0) = NaN;
%         L2norm_roi(p,t,m) = mean(newL2norm(:),'omitnan');        
        newL2normIRF = l2normirf.*newmask;
        newL2normIRF(newL2normIRF==0) = NaN;
        L2normIRF_roi(p,t,m) = mean(newL2normIRF(:),'omitnan');        
        peakTime_roi(p,t,m) = mean(peakTime,'omitnan');
        amplitudePos_roi(p,t,m) = mean(amplitudePos,'omitnan');
        amplitudeNeg_roi(p,t,m) = mean(amplitudeNeg,'omitnan');
        fwhm_roi(p,t,m) = mean(fwhm,'omitnan'); 
        HRF_HbT{p,t,m} = HRF_roi;
%         L2normHbT{m,t} = L2norm;
        L2normIRF{m,t} = l2normirf;
        Corr{m,t} = correlation;
    end
end

end

%% 
save('NVC_HbT_stats.mat','Corr','corr_roi','peakTime_roi','amplitudePos_roi','amplitudeNeg_roi','fwhm_roi',...
    'HRF_HbT','L2normIRF','L2normIRF_roi') %,'corr_comb_stim','corr_comb_rest')

%% Deviation from pre-stroke

mcorr_roi = mean(corr_roi,3);
stdcorr_roi = std(corr_roi,[],3); %./sqrt(size(corr_roi,3));
a = 0.05;
[h_fore(1),p(1)] = ttest2(corr_roi(1,1,:),corr_roi(1,2,:),'Alpha',a);
[h_fore(2),p(2)] = ttest2(corr_roi(1,1,:),corr_roi(1,3,:),'Alpha',a);
[h_fore(3),p(3)] = ttest2(corr_roi(1,1,:),corr_roi(1,4,:),'Alpha',a);
[h_fore(4),p(4)] = ttest2(corr_roi(1,1,:),corr_roi(1,5,:),'Alpha',a);
[h_peri(1),p(1)] = ttest2(corr_roi(2,1,:),corr_roi(2,2,:),'Alpha',a);
[h_peri(2),p(2)] = ttest2(corr_roi(2,1,:),corr_roi(2,3,:),'Alpha',a);
[h_peri(3),p(3)] = ttest2(corr_roi(2,1,:),corr_roi(2,4,:),'Alpha',a);
[h_peri(4),p(4)] = ttest2(corr_roi(2,1,:),corr_roi(2,5,:),'Alpha',a);
[h_contra(1),p(1)] = ttest2(corr_roi(3,1,:),corr_roi(3,2,:),'Alpha',a);
[h_contra(2),p(2)] = ttest2(corr_roi(3,1,:),corr_roi(3,3,:),'Alpha',a);
[h_contra(3),p(3)] = ttest2(corr_roi(3,1,:),corr_roi(3,4,:),'Alpha',a);
[h_contra(4),p(4)] = ttest2(corr_roi(3,1,:),corr_roi(3,5,:),'Alpha',a);

grp=3;
figure
roi = {'Stroke','Peri-infarct','Contralesional'};
x=1:grp;
b = bar(x,mcorr_roi);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = grp;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mcorr_roi(:,e),stdcorr_roi(:,e),stdcorr_roi(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
hold off
set(gca, 'XTickLabel', roi);
% legend('Pre-Stroke','Day2','Week1','Week2','Week4')
ylim([0 1])
title('Correlation coefficient')
set(gca,'FontSize',16)


%%
% load('NVC_stats')
grp = 4;
for m = 1:12
    corr_norm(:,:,m) = corr_roi(:,:,m)./corr_roi(:,1,m);
    peakTime_norm(:,:,m) = peakTime_roi(:,:,m)./peakTime_roi(:,1,m);
    amplitudePos_norm(:,:,m) = amplitudePos_roi(:,:,m)./amplitudePos_roi(:,1,m);
    amplitudeNeg_norm(:,:,m) = amplitudeNeg_roi(:,:,m)./amplitudeNeg_roi(:,1,m);
    fwhm_norm(:,:,m) = fwhm_roi(:,:,m)./fwhm_roi(:,1,m); 
end

mcorr_roi = mean(corr_roi,3);
stdcorr_roi = std(corr_roi,[],3)./sqrt(size(corr_roi,3));
mpeakTime_roi = mean(peakTime_roi,3);
stdpeakTime_roi = std(peakTime_roi,[],3)./sqrt(size(corr_roi,3));
mamplitudePos_roi = mean(amplitudePos_roi,3);
stdamplitudePos_roi = std(amplitudePos_roi,[],3)./sqrt(size(corr_roi,3));
mamplitudeNeg_roi = mean(amplitudeNeg_roi,3);
stdamplitudeNeg_roi = std(amplitudeNeg_roi,[],3)./sqrt(size(corr_roi,3));
mfwhm_roi = mean(fwhm_roi,3);
stdfwhm_roi = std(fwhm_roi,[],3)./sqrt(size(corr_roi,3));

mcorr_norm = mean(corr_norm,3);
stdcorr_norm = std(corr_norm,[],3)./sqrt(size(corr_norm,3));
mpeakTime_norm = mean(peakTime_norm,3);
stdpeakTime_norm = std(peakTime_norm,[],3)./sqrt(size(corr_norm,3));
mamplitudePos_norm = mean(amplitudePos_norm,3);
stdamplitudePos_norm = std(amplitudePos_norm,[],3)./sqrt(size(corr_norm,3));
mamplitudeNeg_norm = mean(amplitudeNeg_norm,3);
stdamplitudeNeg_norm = std(amplitudeNeg_norm,[],3)./sqrt(size(corr_norm,3));
mfwhm_norm = mean(fwhm_norm,3);
stdfwhm_norm = std(fwhm_norm,[],3)./sqrt(size(corr_norm,3));

figure
roi = {'Forelimb','Peri-infarct','Contralesional','Stroke'};
x=1:grp;
b = bar(x,mcorr_roi);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = grp;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mcorr_roi(:,e),stdcorr_roi(:,e),stdcorr_roi(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
hold off
set(gca, 'XTickLabel', roi);
% legend('Pre-Stroke','Day2','Week1','Week2','Week4')
ylim([0 1])
title('Correlation coefficient')
set(gca,'FontSize',16)

figure
x=1:grp;
b = bar(x,mpeakTime_roi);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = grp;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mpeakTime_roi(:,e),stdpeakTime_roi(:,e),stdpeakTime_roi(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
hold off
ylim([0 2])
set(gca, 'XTickLabel', roi);
% legend('Pre-Stroke','Day2','Week1','Week2','Week4')
title('Time to peak (sec)')
set(gca,'FontSize',16)

figure
x=1:grp;
b = bar(x,mamplitudePos_roi);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = grp;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mamplitudePos_roi(:,e),stdamplitudePos_roi(:,e),stdamplitudePos_roi(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
hold off
ylim([0 3e-5])
set(gca, 'XTickLabel', roi);
% legend('Pre-Stroke','Day2','Week1','Week2','Week4')
title('Peak amplitude')
set(gca,'FontSize',16)

figure
x=1:grp;
b = bar(x,mamplitudeNeg_roi);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = grp;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mamplitudeNeg_roi(:,e),stdamplitudeNeg_roi(:,e),stdamplitudeNeg_roi(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
hold off
ylim([-1.5e-5 0])
set(gca, 'XTickLabel', roi);
% legend('Pre-Stroke','Day2','Week1','Week2','Week4')
title('Post-stimulus undershoot')
set(gca,'FontSize',16)

figure
x=1:grp;
b = bar(x,mfwhm_roi);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = grp;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mfwhm_roi(:,e),stdfwhm_roi(:,e),stdfwhm_roi(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
hold off
ylim([0 1.5])
set(gca, 'XTickLabel', roi);
% legend('Pre-Stroke','Day2','Week1','Week2','Week4')
title('Width at half-max (sec)')
set(gca,'FontSize',16)

%% Comparing evoked and resting state correlations

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
corr_comb_rest = [];
for m = 1:12
    m
baseline = [mouse{m},'RestingState/Baseline'];
day2 = [mouse{m},'RestingState/Day2'];
week1 = [mouse{m},'RestingState/Week1'];
week2 = [mouse{m},'RestingState/Week2'];
week4 = [mouse{m},'RestingState/Week4'];
pathname = {baseline, day2, week1, week2, week4};
time = -5+0.2:0.2:10;
mask1 = [mouse{m},'brainmaskSFDI.mat'];
load(mask1)
for t = 1
    dataDir = pathname{t};
    load([dataDir,'/','rest_HRF.mat'])
    newmask = maskSFDI(t).aff_mask + maskSFDI(t).unaff_mask;
    newcorr = correlation.*newmask;
    newcorr = reshape(newcorr,[128*128 1]);
    newcorr(newcorr==0) = [];
    newcorr(isnan(newcorr)) = [];
    corr_comb_rest = cat(1,corr_comb_rest,newcorr);
end
end

corr_comb_stim = [];
for m = 1:12
    m
baseline = [mouse{m},'FunctionalActivation/Baseline'];
day2 = [mouse{m},'FunctionalActivation/Day2'];
week1 = [mouse{m},'FunctionalActivation/Week1'];
week2 = [mouse{m},'FunctionalActivation/Week2'];
week4 = [mouse{m},'FunctionalActivation/Week4'];
pathname = {baseline, day2, week1, week2, week4};
time = -5+0.2:0.2:10;
mask1 = [mouse{m},'brainmaskSFDI.mat'];
load(mask1)
for t = 1
    dataDir = pathname{t};
    load([dataDir,'/','act_HRF.mat'])
    newmask = maskSFDI(t).aff_mask + maskSFDI(t).unaff_mask;
    newcorr = correlation.*newmask;
    newcorr = reshape(newcorr,[128*128 1]);
    newcorr(newcorr==0) = [];
    newcorr(isnan(newcorr)) = [];
    corr_comb_stim = cat(1,corr_comb_stim,newcorr);
end
end

%% Rest vs Stim IRF

mouse = {SS75, SS76, SS77, SS78, SS79};

for m = 1:5
    baseline_rest = [mouse{m},'RestingState/Baseline/rest_HRF'];
    baseline_act = [mouse{m},'FunctionalActivation/Baseline/act_HRF'];
    pathname = {baseline_rest, baseline_act};
    for p = 1:2
        load(pathname{p})
        clear peakTime amplitudePos amplitudeNeg fwhm
        figure(1)
        imagesc(correlation)
        axis image
        colormap jet
        caxis([0 1])
        axis off
        h = imrect;
        fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
        setPositionConstraintFcn(h,fcn);
        pos = wait(h);
        delete(h);
        x1(1) = int32(pos(1));
        x1(2) = x1(1)+int32(pos(3));
        y1(1) = int32(pos(2));
        y1(2) = y1(1)+int32(pos(4));

        HRF_roi = HRF(y1(1):y1(2),x1(1):x1(2),:);
        HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
        mHRF(p,:,m) = mean(HRF_roi,1);
        stdHRF(p,:,m) = std(HRF_roi,1);
    end
end

%
rest = mean(mHRF(1,28:53,:),3);
act = mean(mHRF(2,28:53,:),3);
time = 0:0.2:5;
figure
subplot(1,3,1)
for m = 1:6
    plot(time,mHRF(1,28:53,m))
    hold on
end
ylim([-1e-5 3e-5])
xlim([0 5])
xlabel('Time (sec)')
title('Rest IRFs')
set(gca,'FontSize',16)
legend('SS75', 'SS76', 'SS77', 'SS78', 'SS79', 'SS80')
subplot(1,3,2)
for m = 1:6
    plot(time,mHRF(2,28:53,m))
    hold on
end
ylim([-1e-5 3e-5])
xlim([0 5])
xlabel('Time (sec)')
title('Evoked IRFs')
set(gca,'FontSize',16)
legend('SS75', 'SS76', 'SS77', 'SS78', 'SS79', 'SS80')

subplot(1,3,3)
plot(time,rest,'LineWidth',2)
hold on
plot(time,act,'LineWidth',2)
ylim([-1e-5 3e-5])
xlim([0 5])
xlabel('Time (sec)')
title('Mean IRFs')
set(gca,'FontSize',16)
legend('Mean rest IRF','Mean stimulation IRF')

%% Stim IRF before and after stroke
% Load data
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

for m = 1:7
    baseline = [mouse{m},'FunctionalActivation\Baseline\Ipsi_ForepawStim_IOSI\act_HRF2'];
    day2 = [mouse{m},'FunctionalActivation\Day2\Ipsi_ForepawStim_IOSI\act_HRF2'];
    week1 = [mouse{m},'FunctionalActivation\Week1\Ipsi_ForepawStim_IOSI\act_HRF2'];
    week2 = [mouse{m},'FunctionalActivation\Week2\Ipsi_ForepawStim_IOSI\act_HRF2'];
    week4 = [mouse{m},'FunctionalActivation\Week4\Ipsi_ForepawStim_IOSI\act_HRF2'];

    pathname = {baseline, day2}; %, week1, week2, week4};
    time = -5+0.2:0.2:10;
    for t = 1:2
        load(pathname{t})
        clear peakTime amplitudePos amplitudeNeg fwhm
        figure(1)
        imagesc(correlation)
        axis image
        colormap jet
        caxis([0 1])
        axis off
        h = imrect;
        fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
        setPositionConstraintFcn(h,fcn);
        pos = wait(h);
        delete(h);
        x1(1) = int32(pos(1));
        x1(2) = x1(1)+int32(pos(3));
        y1(1) = int32(pos(2));
        y1(2) = y1(1)+int32(pos(4));

        HRF_roi = HRF(y1(1):y1(2),x1(1):x1(2),:);
        HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
        mHRF(t,:,m) = mean(HRF_roi,1);
        stdHRF(t,:,m) = std(HRF_roi,1);
    end
end

for t = 1:2
    figure
    for m = 1:7
        plot(time(5:end-10),mHRF(t,8:end-10,m),'LineWidth',2)
        hold on
    end
    ylim([-2e-5 3e-5])
    xlim([-4 8])
    xlabel('Time (sec)')
    set(gca,'FontSize',16)
    legend('Mouse 1', 'Mouse 2', 'Mouse 3', 'Mouse 4', 'Mouse 5', 'Mouse 6', 'Mouse 7')
end

%% Stim IRF across rois

mouse = {SS75, SS76, SS77, SS78, SS79, SS80};
color = ['b','y','g'];
for m = 1
    baseline = [mouse{m},'FunctionalActivation\Baseline\Ipsi_ForepawStim_IOSI\act_HRF2'];
    day2 = [mouse{m},'FunctionalActivation\Day2\Ipsi_ForepawStim_IOSI\act_HRF2'];
%     week1 = [mouse{m},'FunctionalActivation\Week1\Ipsi_ForepawStim_IOSI\act_HRF2'];
%     week2 = [mouse{m},'FunctionalActivation\Week2\Ipsi_ForepawStim_IOSI\act_HRF2'];
%     week4 = [mouse{m},'FunctionalActivation\Week4\Ipsi_ForepawStim_IOSI\act_HRF2'];

%     pathname = {baseline, day2, week1, week2, week4};
    pathname = {baseline};
    time = -5+0.2:0.2:10;
    for t = 1
        load(pathname{t})
        clear peakTime amplitudePos amplitudeNeg fwhm
        figure(1)
        subplot(1,2,1)
        imagesc(correlation)
        axis image
        colormap jet
        caxis([-1 1])
        axis off
        hold on
        for r = 1:3
            figure(1)
            subplot(1,2,1)
            h = imrect;
            fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
            setPositionConstraintFcn(h,fcn);
            pos = wait(h);
            delete(h);
            plot([pos(1),pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)],[pos(2),pos(2)+pos(4),pos(2)+pos(4),pos(2),pos(2)],color(r),'LineWidth',3);
            x1(1) = int32(pos(1));
            x1(2) = x1(1)+int32(pos(3));
            y1(1) = int32(pos(2));
            y1(2) = y1(1)+int32(pos(4));

            HRF_roi = HRF(y1(1):y1(2),x1(1):x1(2),:);
            HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
            figure(1)
            subplot(1,2,2)
            mHRF(t,:,m) = mean(HRF_roi,1);
            stdHRF(t,:,m) = std(HRF_roi,1);
            plot(time(5:end-10),mHRF(t,8:end-10,m),'LineWidth',2)
            hold on
            plot(time, HRF_roi(:,4:end),color(r))
            
        end
        xlabel('Time (sec)')
        title('IRF')
        ylim([-1.5*10^-5 3.5*10^-5])
        hold on
        plot([0 0],[-1.5*10^-5 3.5*10^-5],'k--')
    end
end

%%
figure
subplot(1,5,1)
imagesc(correlation)
axis image
colormap jet
caxis([-1 1])
axis off
title('Correlation')
colorbar

subplot(1,5,2)
imagesc(peakTime)
axis image
colormap jet
caxis([0 2])
axis off
title('Time to peak (sec)')
colorbar

subplot(1,5,3)
imagesc(amplitudePos)
axis image
colormap jet
caxis([0 3e-5])
axis off
title('Peak positive amplitude')
colorbar

subplot(1,5,4)
imagesc(amplitudeNeg)
axis image
colormap jet
caxis([-1e-5 0])
axis off
title('Peak negative amplitude')
colorbar

subplot(1,5,5)
imagesc(fwhm)
axis image
colormap jet
caxis([0 3])
axis off
title('Width at half-max (sec)')
colorbar
