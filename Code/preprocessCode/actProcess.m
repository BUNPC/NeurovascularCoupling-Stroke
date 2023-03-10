%% Activation area metrics

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

% Magnitude of response (without threshold)
%
for m = 1:length(mouse)
    m
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
        if t == 1 || t == 2
        forelimb = maskSFDI(1).ipsiOutline;
        newMask = imdilate(maskSFDI(1).ipsiOutline, true(30));
        peri = abs(newMask - maskSFDI(1).ipsiOutline);
        stroke = forelimb;
        else
        forelimb = maskSFDI(1).ipsiOutline;
        newMask = imdilate(maskSFDI(t).stroke_mask, true(30));
        peri = abs(newMask - maskSFDI(t).stroke_mask);
        stroke = maskSFDI(t).stroke_mask;
        end
        respmask = peri;
        GCaMPcorrresp = respmask.*squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
        GCaMPcorrresp(GCaMPcorrresp==0) = NaN;
        HbOresp = respmask.*squeeze(mean(HbOavg(:,:,rStart:rEnd),3));
        HbOresp(HbOresp==0) = NaN;
        HbRresp = respmask.*squeeze(mean(HbRavg(:,:,rStart:rEnd),3));
        HbRresp(HbRresp==0) = NaN;
        HbTresp = respmask.*squeeze(mean(HbTavg(:,:,rStart:rEnd),3));
        HbTresp(HbTresp==0) = NaN;
        
        magHbT(m,t) = mean(HbTresp(:),'omitnan');
        magHbO(m,t) = mean(HbOresp(:),'omitnan');
        magHbR(m,t) = mean(HbRresp(:),'omitnan');
        magG(m,t) = mean(GCaMPcorrresp(:),'omitnan');
    end
end
magG = 100.*magG; % convert percent change
magHbO = 1e6.*magHbO; % convert to micro molar
magHbR = 1e6.*magHbR;
magHbT = 1e6.*magHbT;

for t = 1:5
    mmagG(t) = mean(magG(:,t));
    smagG(t) = std(magG(:,t));
    mmagHbO(t) = mean(magHbO(:,t));
    smagHbO(t) = std(magHbO(:,t));
    mmagHbR(t) = mean(magHbR(:,t));
    smagHbR(t) = std(magHbR(:,t));
    mmagHbT(t) = mean(magHbT(:,t));
    smagHbT(t) = std(magHbT(:,t));
end
mmag = [mmagG; mmagHbT; mmagHbO; mmagHbR];
smag = [smagG; smagHbT; smagHbO; smagHbR];
a = 0.05;
[hG(1),pG(1)] = ttest2(magG(:,1),magG(:,2),'Alpha',a);
[hG(2),pG(2)] = ttest2(magG(:,1),magG(:,3),'Alpha',a);
[hG(3),pG(3)] = ttest2(magG(:,1),magG(:,4),'Alpha',a);
[hG(4),pG(4)] = ttest2(magG(:,1),magG(:,5),'Alpha',a);
[hT(1),pG(1)] = ttest2(magHbT(:,1),magHbT(:,2),'Alpha',a);
[hT(2),pG(2)] = ttest2(magHbT(:,1),magHbT(:,3),'Alpha',a);
[hT(3),pG(3)] = ttest2(magHbT(:,1),magHbT(:,4),'Alpha',a);
[hT(4),pG(4)] = ttest2(magHbT(:,1),magHbT(:,5),'Alpha',a);
[hO(1),pG(1)] = ttest2(magHbO(:,1),magHbO(:,2),'Alpha',a);
[hO(2),pG(2)] = ttest2(magHbO(:,1),magHbO(:,3),'Alpha',a);
[hO(3),pG(3)] = ttest2(magHbO(:,1),magHbO(:,4),'Alpha',a);
[hO(4),pG(4)] = ttest2(magHbO(:,1),magHbO(:,5),'Alpha',a);
[hR(1),pG(1)] = ttest2(magHbR(:,1),magHbR(:,2),'Alpha',a);
[hR(2),pG(2)] = ttest2(magHbR(:,1),magHbR(:,3),'Alpha',a);
[hR(3),pG(3)] = ttest2(magHbR(:,1),magHbR(:,4),'Alpha',a);
[hR(4),pG(4)] = ttest2(magHbR(:,1),magHbR(:,5),'Alpha',a);

figure
b = bar(1:4,mmag);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = 4;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mmag(:,e),smag(:,e),smag(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'GCaMP','HbT','HbO','HbR'};
set(gca, 'XTickLabel', group, 'FontSize', 16);
ylabel('\DeltaF/F (%)   \DeltaHb (\muM)')
set(b, 'FaceAlpha', 1)
% legend('Pre-stroke','Day2','Week1','Week2','Week4')
ylim([-3 5])


%%
load('respMag.mat')
respMag(1).impPeri = magG;
respMag(2).impPeri = magHbO;
respMag(3).impPeri = magHbR;
respMag(4).impPeri = magHbT;
save('respMag.mat','respMag')

%% Area and magnitude of response (with thresholding)

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
mouse = {SS75, SS76, SS77, SS78, SS79, SS80, SS81, SS82, SS83, SS84, SS93};
mousename = {'SS75', 'SS76', 'SS77', 'SS78', 'SS79', 'SS80', 'SS81', 'SS82', 'SS83', 'SS84', 'SS93'};

thresh = [0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95]; 
for m = 1:length(mouse)
    m
    [timepoints, img, mask] = animals(mousename{m});
    load(mask)
    for t = 1:5
        load(timepoints{t})
        rStart = baseFrames;
        rEnd = 2*baseFrames;
        
        for th = 1:length(thresh)
        
        respmask = maskSFDI(t).aff_mask;
        
        GCaMPcorrresp = respmask.*squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
        GCaMPcorrresp(GCaMPcorrresp==0) = NaN;
        HbOresp = respmask.*squeeze(mean(HbOavg(:,:,rStart:rEnd),3));
        HbOresp(HbOresp==0) = NaN;
        HbRresp = respmask.*squeeze(mean(HbRavg(:,:,rStart:rEnd),3));
        HbRresp(HbRresp==0) = NaN;
        HbTresp = respmask.*squeeze(mean(HbTavg(:,:,rStart:rEnd),3));
        HbTresp(HbTresp==0) = NaN;
        
        cutGCaMPcorrresp = GCaMPcorrresp; %imgaussfilt(GCaMPcorrresp,1); 
        cutHbO = HbOresp; %imgaussfilt(HbTresp,2);
        cutHbR = HbRresp; %imgaussfilt(HbTresp,2);
        cutHbT = HbTresp; %imgaussfilt(HbTresp,2);
        
        if t == 1
            peakG = prctile(cutGCaMPcorrresp(:),99);
            peakHbO = prctile(cutHbO(:),99);
            peakHbR = prctile(cutHbR(:),1);
            peakHbT = prctile(cutHbT(:),99);
        end
        cutGCaMPcorrresp(cutGCaMPcorrresp<=thresh(th)*peakG) = NaN;
        cutHbO(cutHbO<=thresh(th)*peakHbO) = NaN;
        cutHbR(cutHbR>=thresh(th)*peakHbR) = NaN;       
        cutHbT(cutHbT<=thresh(th)*peakHbT) = NaN;
        
        GCaMP{m,t,th} = cutGCaMPcorrresp;
        HbT{m,t,th} = cutHbT;
        HbO{m,t,th} = cutHbO;
        HbR{m,t,th} = cutHbR;
        areaHbT(m,t,th) = length(find(~isnan(cutHbT)));
        areaHbO(m,t,th) = length(find(~isnan(cutHbO)));
        areaHbR(m,t,th) = length(find(~isnan(cutHbR)));
        areaG(m,t,th) = length(find(~isnan(cutGCaMPcorrresp)));
        magHbT(m,t,th) = sum(cutHbT(:),'omitnan')./areaHbT(m,t,th);
        magHbO(m,t,th) = sum(cutHbO(:),'omitnan')./areaHbO(m,t,th);
        magHbR(m,t,th) = sum(cutHbR(:),'omitnan')./areaHbR(m,t,th);
        magG(m,t,th) = sum(cutGCaMPcorrresp(:),'omitnan')./areaG(m,t,th);
        
        end
    end
end

magG = 100.*magG; % convert percent change
magG(isnan(magG)) = 0;
magHbO = 1e6.*magHbO; % convert to micro molar
magHbR = 1e6.*magHbR;
magHbT = 1e6.*magHbT;
magHbO(isnan(magHbO)) = 0;
magHbR(isnan(magHbR)) = 0;
magHbT(isnan(magHbT)) = 0;


areaG = 0.0027.*areaG; % convert to pixel area in mm^2
areaHbO = 0.0027.*areaHbO;
areaHbR = 0.0027.*areaHbR;
areaHbT = 0.0027.*areaHbT;



%% Magnitude of the response

% Absolute magnitude
th = 6;
for t = 1:5
    mmagG(t) = mean(magG(:,t,th));
    smagG(t) = std(magG(:,t,th));
    mmagHbO(t) = mean(magHbO(:,t,th));
    smagHbO(t) = std(magHbO(:,t,th));
    mmagHbR(t) = mean(magHbR(:,t,th));
    smagHbR(t) = std(magHbR(:,t,th));
    mmagHbT(t) = mean(magHbT(:,t,th));
    smagHbT(t) = std(magHbT(:,t,th));
end
mmag = [mmagG; mmagHbT; mmagHbO; mmagHbR];
smag = [smagG; smagHbT; smagHbO; smagHbR];
a = 0.05;
[hG(1),pG(1)] = ttest2(magG(:,1,th),magG(:,2,th),'Alpha',a);
[hG(2),pG(2)] = ttest2(magG(:,1,th),magG(:,3,th),'Alpha',a);
[hG(3),pG(3)] = ttest2(magG(:,1,th),magG(:,4,th),'Alpha',a);
[hG(4),pG(4)] = ttest2(magG(:,1,th),magG(:,5,th),'Alpha',a);
[hT(1),pG(1)] = ttest2(magHbT(:,1,th),magHbT(:,2,th),'Alpha',a);
[hT(2),pG(2)] = ttest2(magHbT(:,1,th),magHbT(:,3,th),'Alpha',a);
[hT(3),pG(3)] = ttest2(magHbT(:,1,th),magHbT(:,4,th),'Alpha',a);
[hT(4),pG(4)] = ttest2(magHbT(:,1,th),magHbT(:,5,th),'Alpha',a);
[hO(1),pG(1)] = ttest2(magHbO(:,1,th),magHbO(:,2,th),'Alpha',a);
[hO(2),pG(2)] = ttest2(magHbO(:,1,th),magHbO(:,3,th),'Alpha',a);
[hO(3),pG(3)] = ttest2(magHbO(:,1,th),magHbO(:,4,th),'Alpha',a);
[hO(4),pG(4)] = ttest2(magHbO(:,1,th),magHbO(:,5,th),'Alpha',a);
[hR(1),pG(1)] = ttest2(magHbR(:,1,th),magHbR(:,2,th),'Alpha',a);
[hR(2),pG(2)] = ttest2(magHbR(:,1,th),magHbR(:,3,th),'Alpha',a);
[hR(3),pG(3)] = ttest2(magHbR(:,1,th),magHbR(:,4,th),'Alpha',a);
[hR(4),pG(4)] = ttest2(magHbR(:,1,th),magHbR(:,5,th),'Alpha',a);

figure
b = bar(1:4,mmag);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = 4;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mmag(:,e),smag(:,e),smag(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'GCaMP','HbT','HbO','HbR'};
set(gca, 'XTickLabel', group, 'FontSize', 16);
ylabel('\DeltaF/F (%)   \DeltaHb (\muM)')
set(b, 'FaceAlpha', 1)
legend('Pre-stroke','Day2','Week1','Week2','Week4')
ylim([-3 5])

% % Magnitude across threshold
figure
fig = gcf;
options.x_axis = thresh;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = [0 0 0]; 
options.color_line = [0 0 0]; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(magG(:,1,:)),options)
hold on
options.color_area = [0.5 0 0]; 
options.color_line = [0.5 0 0]; 
hold on
plot_areaerrorbar(squeeze(magG(:,2,:)),options)
options.color_area = [1 0 0]; 
options.color_line = [1 0 0]; 
hold on
plot_areaerrorbar(squeeze(magG(:,3,:)),options)
options.color_area = [1 0.5 0]; 
options.color_line = [1 0.5 0]; 
hold on
plot_areaerrorbar(squeeze(magG(:,4,:)),options)
options.color_area = [1 1 0]; 
options.color_line = [1 1 0]; 
hold on
plot_areaerrorbar(squeeze(magG(:,5,:)),options)
xlabel('Threshold')
ylabel('\DeltaF/F (%)')
set(gca, 'FontSize', 16);
set(b, 'FaceAlpha', 1)
ylim([-5 8])
xlim([0.5 0.95])

figure
fig = gcf;
options.x_axis = thresh;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = [0 0 0]; 
options.color_line = [0 0 0]; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(magHbT(:,1,:)),options)
hold on
options.color_area = [0.5 0 0]; 
options.color_line = [0.5 0 0]; 
hold on
plot_areaerrorbar(squeeze(magHbT(:,2,:)),options)
options.color_area = [1 0 0]; 
options.color_line = [1 0 0]; 
hold on
plot_areaerrorbar(squeeze(magHbT(:,3,:)),options)
options.color_area = [1 0.5 0]; 
options.color_line = [1 0.5 0]; 
hold on
plot_areaerrorbar(squeeze(magHbT(:,4,:)),options)
options.color_area = [1 1 0]; 
options.color_line = [1 1 0]; 
hold on
plot_areaerrorbar(squeeze(magHbT(:,5,:)),options)
xlabel('Threshold')
ylabel('\DeltaHbT (\muM)')
set(gca, 'FontSize', 16);
set(b, 'FaceAlpha', 1)
ylim([-5 8])
xlim([0.5 0.95])

figure
fig = gcf;
options.x_axis = thresh;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = [0 0 0]; 
options.color_line = [0 0 0]; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(magHbO(:,1,:)),options)
hold on
options.color_area = [0.5 0 0]; 
options.color_line = [0.5 0 0]; 
hold on
plot_areaerrorbar(squeeze(magHbO(:,2,:)),options)
options.color_area = [1 0 0]; 
options.color_line = [1 0 0]; 
hold on
plot_areaerrorbar(squeeze(magHbO(:,3,:)),options)
options.color_area = [1 0.5 0]; 
options.color_line = [1 0.5 0]; 
hold on
plot_areaerrorbar(squeeze(magHbO(:,4,:)),options)
options.color_area = [1 1 0]; 
options.color_line = [1 1 0]; 
hold on
plot_areaerrorbar(squeeze(magHbO(:,5,:)),options)
xlabel('Threshold')
ylabel('\DeltaHbO (\muM)')
set(gca, 'FontSize', 16);
set(b, 'FaceAlpha', 1)
ylim([-5 8])
xlim([0.5 0.95])

figure
fig = gcf;
options.x_axis = thresh;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = [0 0 0]; 
options.color_line = [0 0 0]; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(magHbR(:,1,:)),options)
hold on
options.color_area = [0.5 0 0]; 
options.color_line = [0.5 0 0]; 
hold on
plot_areaerrorbar(squeeze(magHbR(:,2,:)),options)
options.color_area = [1 0 0]; 
options.color_line = [1 0 0]; 
hold on
plot_areaerrorbar(squeeze(magHbR(:,3,:)),options)
options.color_area = [1 0.5 0]; 
options.color_line = [1 0.5 0]; 
hold on
plot_areaerrorbar(squeeze(magHbR(:,4,:)),options)
options.color_area = [1 1 0]; 
options.color_line = [1 1 0]; 
hold on
plot_areaerrorbar(squeeze(magHbR(:,5,:)),options)
xlabel('Threshold')
ylabel('\DeltaHbR (\muM)')
set(gca, 'FontSize', 16);
set(b, 'FaceAlpha', 1)
ylim([-5 8])
xlim([0.5 0.95])

%% Area of the response

% Absolute area
th = 6;
for t = 1:5
    mareaG(t) = mean(areaG(:,t,th));
    sareaG(t) = std(areaG(:,t,th));
    mareaHbO(t) = mean(areaHbO(:,t,th));
    sareaHbO(t) = std(areaHbO(:,t,th));
    mareaHbR(t) = mean(areaHbR(:,t,th));
    sareaHbR(t) = std(areaHbR(:,t,th));
    mareaHbT(t) = mean(areaHbT(:,t,th));
    sareaHbT(t) = std(areaHbT(:,t,th));
end
marea = [mareaG; mareaHbT; mareaHbO; mareaHbR];
sarea = [sareaG; sareaHbT; sareaHbO; sareaHbR];
a = 0.05;
[hG(1),pG(1)] = ttest2(areaG(:,1,th),areaG(:,2,th),'Alpha',a);
[hG(2),pG(2)] = ttest2(areaG(:,1,th),areaG(:,3,th),'Alpha',a);
[hG(3),pG(3)] = ttest2(areaG(:,1,th),areaG(:,4,th),'Alpha',a);
[hG(4),pG(4)] = ttest2(areaG(:,1,th),areaG(:,5,th),'Alpha',a);
[hT(1),pG(1)] = ttest2(areaHbT(:,1,th),areaHbT(:,2,th),'Alpha',a);
[hT(2),pG(2)] = ttest2(areaHbT(:,1,th),areaHbT(:,3,th),'Alpha',a);
[hT(3),pG(3)] = ttest2(areaHbT(:,1,th),areaHbT(:,4,th),'Alpha',a);
[hT(4),pG(4)] = ttest2(areaHbT(:,1,th),areaHbT(:,5,th),'Alpha',a);
[hO(1),pG(1)] = ttest2(areaHbO(:,1,th),areaHbO(:,2,th),'Alpha',a);
[hO(2),pG(2)] = ttest2(areaHbO(:,1,th),areaHbO(:,3,th),'Alpha',a);
[hO(3),pG(3)] = ttest2(areaHbO(:,1,th),areaHbO(:,4,th),'Alpha',a);
[hO(4),pG(4)] = ttest2(areaHbO(:,1,th),areaHbO(:,5,th),'Alpha',a);
[hR(1),pG(1)] = ttest2(areaHbR(:,1,th),areaHbR(:,2,th),'Alpha',a);
[hR(2),pG(2)] = ttest2(areaHbR(:,1,th),areaHbR(:,3,th),'Alpha',a);
[hR(3),pG(3)] = ttest2(areaHbR(:,1,th),areaHbR(:,4,th),'Alpha',a);
[hR(4),pG(4)] = ttest2(areaHbR(:,1,th),areaHbR(:,5,th),'Alpha',a);

figure
b = bar(1:4,marea);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = 4;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,marea(:,e),sarea(:,e),sarea(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'GCaMP','HbT','HbO','HbR'};
set(gca, 'XTickLabel', group, 'FontSize', 16);
ylabel('Response Area (mm^2)')
set(b, 'FaceAlpha', 1)
legend('Pre-stroke','Day2','Week1','Week2','Week4')
% ylim([-5 10])
% set(b(1), 'FaceAlpha', 1)
% set(b(2), 'FaceAlpha', 0.8)
% set(b(3), 'FaceAlpha', 0.6)
% set(b(4), 'FaceAlpha', 0.4)
% set(b(5), 'FaceAlpha', 0.2)


% Normalized area
% figure
% b = bar(1:4,marea./marea(:,1));
% b(1).FaceColor = 'k';
% b(2).FaceColor = 'r';
% b(3).FaceColor = 'g';
% b(4).FaceColor = 'b';
% b(5).FaceColor = 'm';
% group = {'GCaMP','HbO','HbR','HbT'};
% set(gca, 'XTickLabel', group, 'FontSize', 16, 'FaceAlpha', 0.5);
% ylim([0 1.2])
% set(b, 'FaceAlpha', 0.8)

% Area across threshold

figure
fig = gcf;
options.x_axis = thresh;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = [0 0 0]; 
options.color_line = [0 0 0]; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(areaG(:,1,:)),options)
hold on
options.color_area = [0.5 0 0]; 
options.color_line = [0.5 0 0]; 
hold on
plot_areaerrorbar(squeeze(areaG(:,2,:)),options)
options.color_area = [1 0 0]; 
options.color_line = [1 0 0]; 
hold on
plot_areaerrorbar(squeeze(areaG(:,3,:)),options)
options.color_area = [1 0.5 0]; 
options.color_line = [1 0.5 0]; 
hold on
plot_areaerrorbar(squeeze(areaG(:,4,:)),options)
options.color_area = [1 1 0]; 
options.color_line = [1 1 0]; 
hold on
plot_areaerrorbar(squeeze(areaG(:,5,:)),options)
xlabel('Threshold')
ylabel('\DeltaF/F Response Area (mm^2)')
set(gca, 'FontSize', 16);
set(b, 'FaceAlpha', 1)
xlim([0.5 0.95])

figure
fig = gcf;
options.x_axis = thresh;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = [0 0 0]; 
options.color_line = [0 0 0]; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(areaHbT(:,1,:)),options)
hold on
options.color_area = [0.5 0 0]; 
options.color_line = [0.5 0 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbT(:,2,:)),options)
options.color_area = [1 0 0]; 
options.color_line = [1 0 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbT(:,3,:)),options)
options.color_area = [1 0.5 0]; 
options.color_line = [1 0.5 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbT(:,4,:)),options)
options.color_area = [1 1 0]; 
options.color_line = [1 1 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbT(:,5,:)),options)
xlabel('Threshold')
ylabel('\DeltaHbT Response Area (mm^2)')
set(gca, 'FontSize', 16);
set(b, 'FaceAlpha', 1)
xlim([0.5 0.95])

figure
fig = gcf;
options.x_axis = thresh;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = [0 0 0]; 
options.color_line = [0 0 0]; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(areaHbO(:,1,:)),options)
hold on
options.color_area = [0.5 0 0]; 
options.color_line = [0.5 0 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbO(:,2,:)),options)
options.color_area = [1 0 0]; 
options.color_line = [1 0 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbO(:,3,:)),options)
options.color_area = [1 0.5 0]; 
options.color_line = [1 0.5 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbO(:,4,:)),options)
options.color_area = [1 1 0]; 
options.color_line = [1 1 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbO(:,5,:)),options)
xlabel('Threshold')
ylabel('\DeltaHbO Response Area (mm^2)')
set(gca, 'FontSize', 16);
set(b, 'FaceAlpha', 1)
xlim([0.5 0.95])

figure
fig = gcf;
options.x_axis = thresh;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = [0 0 0]; 
options.color_line = [0 0 0]; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(areaHbR(:,1,:)),options)
hold on
options.color_area = [0.5 0 0]; 
options.color_line = [0.5 0 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbR(:,2,:)),options)
options.color_area = [1 0 0]; 
options.color_line = [1 0 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbR(:,3,:)),options)
options.color_area = [1 0.5 0]; 
options.color_line = [1 0.5 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbR(:,4,:)),options)
options.color_area = [1 1 0]; 
options.color_line = [1 1 0]; 
hold on
plot_areaerrorbar(squeeze(areaHbR(:,5,:)),options)
xlabel('Threshold')
ylabel('\DeltaHbR Response Area (mm^2)')
set(gca, 'FontSize', 16);
set(b, 'FaceAlpha', 1)
xlim([0.5 0.95])



%% Plotting


th = 6;
time = 3;
Gcomb = reshape(mGCaMP(:,time,:,th),[size(mGCaMP,1)*size(mGCaMP,3) 1]);
HbTcomb = reshape(mHbT(:,time,:,th),[size(mHbT,1)*size(mHbT,3) 1]);
HbOcomb = reshape(mHbO(:,time,:,th),[size(mHbO,1)*size(mHbO,3) 1]);
HbRcomb = reshape(mHbR(:,time,:,th),[size(mHbR,1)*size(mHbR,3) 1]);
figure
subplot(2,3,1)
[p,S] = polyfit(Gcomb,HbTcomb,1);
f = polyval(p,Gcomb);
plot(Gcomb,HbTcomb,'o',Gcomb,f,'-')
[r,p] = corrcoef(Gcomb,HbTcomb);
z = 0.5*(log(1+r) - log(1-r));
title(['HbT: r = ',num2str(r(1,2)),', z = ',num2str(z(1,2)),', p = ',num2str(p(1,2))])

subplot(2,3,2)
[p,S] = polyfit(Gcomb,HbOcomb,1);
f = polyval(p,Gcomb);
plot(Gcomb,HbOcomb,'o',Gcomb,f,'-')
[r,p] = corrcoef(Gcomb,HbOcomb);
z = 0.5*(log(1+r) - log(1-r));
title(['HbO: r = ',num2str(r(1,2)),', z = ',num2str(z(1,2)),', p = ',num2str(p(1,2))])

subplot(2,3,3)
[p,S] = polyfit(Gcomb,HbRcomb,1);
f = polyval(p,Gcomb);
plot(Gcomb,HbRcomb,'o',Gcomb,f,'-')
[r,p] = corrcoef(Gcomb,HbRcomb);
z = 0.5*(log(1+r) - log(1-r));
title(['HbR: r = ',num2str(r(1,2)),', z = ',num2str(z(1,2)),', p = ',num2str(p(1,2))])

subplot(2,3,4)
bar(mean(simT(:)))
ylim([0 1])
subplot(2,3,5)
bar(mean(simO(:)))
ylim([0 1])
subplot(2,3,6)
bar(mean(simR(:)))
ylim([0 1])


for m = 1:length(mouse)
figure(101)
subplot(3,4,m)
[p,S] = polyfit(mGCaMP(m,:),mHbT(m,:),1);
f = polyval(p,mGCaMP(m,:));
plot(mGCaMP(m,:),mHbT(m,:),'o',mGCaMP(m,:),f,'-')
[r,p] = corrcoef(mGCaMP(m,:),mHbT(m,:));
z = 0.5*(log(1+r) - log(1-r));
title(['r = ',num2str(r(1,2)),' z = ',num2str(z(1,2))])

figure(102)
subplot(3,4,m)
[p,S] = polyfit(mGCaMP(m,:),mHbO(m,:),1);
f = polyval(p,mGCaMP(m,:));
plot(mGCaMP(m,:),mHbO(m,:),'o',mGCaMP(m,:),f,'-')
[r,p] = corrcoef(mGCaMP(m,:),mHbO(m,:));
z = 0.5*(log(1+r) - log(1-r));
title(['r = ',num2str(r(1,2)),' z = ',num2str(z(1,2))])

figure(103)
subplot(3,4,m)
[p,S] = polyfit(mGCaMP(m,:),mHbR(m,:),1);
f = polyval(p,mGCaMP(m,:));
plot(mGCaMP(m,:),mHbR(m,:),'o',mGCaMP(m,:),f,'-')
[r,p] = corrcoef(mGCaMP(m,:),mHbR(m,:));
z = 0.5*(log(1+r) - log(1-r));
title(['r = ',num2str(r(1,2)),' z = ',num2str(z(1,2))])


end

