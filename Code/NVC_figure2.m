%% Stroke NVC figure 2
% Smrithi Sunil
% BOAS Lab

%% Figure 2b

clear
mousename = {'SS75', 'SS76', 'SS77', 'SS78', 'SS79', 'SS80', 'SS81', 'SS82', 'SS83', 'SS84', 'SS85', 'SS93'};
% Magnitude of response (without threshold) for contralateral stimulation
for m = 1:length(mousename)
    m
    [timepoints, mask] = animals(mousename{m});
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
 
        arespmask = maskSFDI(t).aff_mask;
        aGCaMPcorrresp = arespmask.*squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
        aGCaMPcorrresp(aGCaMPcorrresp==0) = NaN;
        aHbOresp = arespmask.*squeeze(mean(HbOavg(:,:,rStart:rEnd),3));
        aHbOresp(aHbOresp==0) = NaN;
        aHbRresp = arespmask.*squeeze(mean(HbRavg(:,:,rStart:rEnd),3));
        aHbRresp(aHbRresp==0) = NaN;
        aHbTresp = arespmask.*squeeze(mean(HbTavg(:,:,rStart:rEnd),3));
        aHbTresp(aHbTresp==0) = NaN;

        urespmask = maskSFDI(t).unaff_mask;
        uGCaMPcorrresp = urespmask.*squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
        uGCaMPcorrresp(uGCaMPcorrresp==0) = NaN;
        uHbOresp = urespmask.*squeeze(mean(HbOavg(:,:,rStart:rEnd),3));
        uHbOresp(uHbOresp==0) = NaN;
        uHbRresp = urespmask.*squeeze(mean(HbRavg(:,:,rStart:rEnd),3));
        uHbRresp(uHbRresp==0) = NaN;
        uHbTresp = urespmask.*squeeze(mean(HbTavg(:,:,rStart:rEnd),3));
        uHbTresp(uHbTresp==0) = NaN;
        
        amagHbT(m,t) = mean(aHbTresp(:),'omitnan');
        amagHbO(m,t) = mean(aHbOresp(:),'omitnan');
        amagHbR(m,t) = mean(aHbRresp(:),'omitnan');
        amagG(m,t) = mean(aGCaMPcorrresp(:),'omitnan');
        
        umagHbT(m,t) = mean(uHbTresp(:),'omitnan');
        umagHbO(m,t) = mean(uHbOresp(:),'omitnan');
        umagHbR(m,t) = mean(uHbRresp(:),'omitnan');
        umagG(m,t) = mean(uGCaMPcorrresp(:),'omitnan');
    end
end

magG = 100.*amagG; % convert percent change
magHbO = 1e6.*amagHbO; % convert to micro molar
magHbR = 1e6.*amagHbR;
magHbT = 1e6.*amagHbT;

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
a = 0.01;
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
set(b, 'FaceAlpha', 1)

% For contralateral affected 
sigline2([0.69 0.85],[],2.6)
sigline2([0.69 1],[],2.9)
sigline2([0.69 1.15],[],3.2)
sigline2([0.69 1.3],[],3.5)

sigline2([1.69 1.85],[],1.9)
sigline2([1.69 2],[],2.2)
sigline2([1.69 2.15],[],2.5)
sigline2([1.69 2.3],[],2.8)

sigline2([2.69 2.85],[],3.9)
sigline2([2.69 3],[],4.2)
sigline2([2.69 3.15],[],4.5)
sigline2([2.69 3.3],[],4.8)
    
group = {'GCaMP','HbT','HbO','HbR'};
ylabel('\DeltaF/F (%)   \DeltaHb (\muM)')
ylim([-4 6])
legend({'Pre-stroke','Day2','Week1','Week2','Week4'})
set(gca, 'XTickLabel', group, 'FontSize', 24);
saveas(gcf,'figures/figure2b_top.png')

magG = 100.*umagG; % convert percent change
magHbO = 1e6.*umagHbO; % convert to micro molar
magHbR = 1e6.*umagHbR;
magHbT = 1e6.*umagHbT;

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
a = 0.01;
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
set(b, 'FaceAlpha', 1)
group = {'GCaMP','HbT','HbO','HbR'};
ylabel('\DeltaF/F (%)   \DeltaHb (\muM)')
ylim([-4 6])
set(gca, 'XTickLabel', group, 'FontSize', 24);
saveas(gcf,'figures/figure2b_bottom.png')

%% Supplementary figure 4a

for t = 1:5
    mmagG(t) = mean(amagG(:,t));
    mmagHbO(t) = mean(amagHbO(:,t));
    mmagHbR(t) = mean(amagHbR(:,t));
    mmagHbT(t) = mean(amagHbT(:,t));
end
mmag = [mmagG; mmagHbT; mmagHbO; mmagHbR];
mmag = 100.*mmag./abs(mmag(:,1));

figure
b = bar(1:4,mmag);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
group = {'GCaMP','HbT','HbO','HbR'};
set(gca, 'XTickLabel', group, 'FontSize', 24);
ylabel({'\DeltaF/F (%)   \DeltaHb (\muM)';'Normalized to pre-stroke'})
set(b, 'FaceAlpha', 1)
legend('Pre-stroke','Day2','Week1','Week2','Week4')
ylim([-120 120])
saveas(gcf,'figures/suppfig4a.png')

%% Figure 2d
clear
% Magnitude of response (without threshold) for ipsilateral stimulation
mousename = {'SS75', 'SS76', 'SS77', 'SS78', 'SS79', 'SS80', 'SS81', 'SS82', 'SS83', 'SS84', 'SS85', 'SS93'};

for m = 1:length(mousename)
    m
    [timepoints, mask] = animals(mousename{m});
    load(mask)
    for t = 1:5        
        load([timepoints{t},'/','act_IOSI_contra.mat'])
        load([timepoints{t},'/','act_GCaMPcorr_contra.mat'])
        % Trial average hemodynamics and GCaMP (corrected and uncorrected)
        HbOavg = trialAverage(HbO, numTrials, numFrame);
        HbRavg = trialAverage(HbR, numTrials, numFrame);
        HbTavg = trialAverage(HbT, numTrials, numFrame);
        GCaMPcorr = trialAverage(deltaGCaMPcorr, numTrials, numFrame);
        rStart = baseFrames;
        rEnd = 2*baseFrames;        
 
        arespmask = maskSFDI(t).aff_mask;
        aGCaMPcorrresp = arespmask.*squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
        aGCaMPcorrresp(aGCaMPcorrresp==0) = NaN;
        aHbOresp = arespmask.*squeeze(mean(HbOavg(:,:,rStart:rEnd),3));
        aHbOresp(aHbOresp==0) = NaN;
        aHbRresp = arespmask.*squeeze(mean(HbRavg(:,:,rStart:rEnd),3));
        aHbRresp(aHbRresp==0) = NaN;
        aHbTresp = arespmask.*squeeze(mean(HbTavg(:,:,rStart:rEnd),3));
        aHbTresp(aHbTresp==0) = NaN;

        urespmask = maskSFDI(t).unaff_mask;
        uGCaMPcorrresp = urespmask.*squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
        uGCaMPcorrresp(uGCaMPcorrresp==0) = NaN;
        uHbOresp = urespmask.*squeeze(mean(HbOavg(:,:,rStart:rEnd),3));
        uHbOresp(uHbOresp==0) = NaN;
        uHbRresp = urespmask.*squeeze(mean(HbRavg(:,:,rStart:rEnd),3));
        uHbRresp(uHbRresp==0) = NaN;
        uHbTresp = urespmask.*squeeze(mean(HbTavg(:,:,rStart:rEnd),3));
        uHbTresp(uHbTresp==0) = NaN;
        
        amagHbT(m,t) = mean(aHbTresp(:),'omitnan');
        amagHbO(m,t) = mean(aHbOresp(:),'omitnan');
        amagHbR(m,t) = mean(aHbRresp(:),'omitnan');
        amagG(m,t) = mean(aGCaMPcorrresp(:),'omitnan');
        
        umagHbT(m,t) = mean(uHbTresp(:),'omitnan');
        umagHbO(m,t) = mean(uHbOresp(:),'omitnan');
        umagHbR(m,t) = mean(uHbRresp(:),'omitnan');
        umagG(m,t) = mean(uGCaMPcorrresp(:),'omitnan');
    end
end

magG = 100.*amagG; % convert percent change
magHbO = 1e6.*amagHbO; % convert to micro molar
magHbR = 1e6.*amagHbR;
magHbT = 1e6.*amagHbT;

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
a = 0.01;
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
set(b, 'FaceAlpha', 1)

% For ipsilateral affected 
sigline2([0.69 0.85],[],1.9)
sigline2([0.69 1],[],2.2)
sigline2([0.69 1.15],[],2.5)
sigline2([0.69 1.3],[],2.8)

sigline2([1.69 1.85],[],1.6)
sigline2([1.69 2],[],1.9)
sigline2([1.69 2.15],[],2.2)
sigline2([1.69 2.3],[],2.5)

sigline2([2.69 2.85],[],3.2)
sigline2([2.69 3],[],3.5)
sigline2([2.69 3.15],[],3.8)
sigline2([2.69 3.3],[],4.1)

sigline2([3.69 4],[],-2.4)
sigline2([3.69 4.15],[],-2.7)
sigline2([3.69 4.3],[],-3)
    
group = {'GCaMP','HbT','HbO','HbR'};
ylabel('\DeltaF/F (%)   \DeltaHb (\muM)')
ylim([-4 6])
set(gca, 'XTickLabel', group, 'FontSize', 24);
saveas(gcf,'figures/figure2d_top.png')

magG = 100.*umagG; % convert percent change
magHbO = 1e6.*umagHbO; % convert to micro molar
magHbR = 1e6.*umagHbR;
magHbT = 1e6.*umagHbT;

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
a = 0.01;
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
set(b, 'FaceAlpha', 1)
group = {'GCaMP','HbT','HbO','HbR'};
ylabel('\DeltaF/F (%)   \DeltaHb (\muM)')
ylim([-4 6])
set(gca, 'XTickLabel', group, 'FontSize', 24);
saveas(gcf,'figures/figure2d_bottom.png')

%% Supplementary figure 4b
for t = 1:5
    mmagG(t) = mean(amagG(:,t));
    mmagHbO(t) = mean(amagHbO(:,t));
    mmagHbR(t) = mean(amagHbR(:,t));
    mmagHbT(t) = mean(amagHbT(:,t));
end
mmag = [mmagG; mmagHbT; mmagHbO; mmagHbR];
mmag = 100.*mmag./abs(mmag(:,1));

figure
b = bar(1:4,mmag);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
group = {'GCaMP','HbT','HbO','HbR'};
set(gca, 'XTickLabel', group, 'FontSize', 24);
ylabel({'\DeltaF/F (%)   \DeltaHb (\muM)';'Normalized to pre-stroke'})
set(b, 'FaceAlpha', 1)
ylim([-120 120])
saveas(gcf,'figures/suppfig4b.png')

%% Supplementary figure 5

mousename = 'SS80';
[timepoints, mask] = animals(mousename);
t=1;
load([timepoints{t},'/','act_IOSI_ipsi.mat'])
load([timepoints{t},'/','act_GCaMPcorr_ipsi.mat'])
% Trial average hemodynamics and GCaMP (corrected and uncorrected)
HbOavg = trialAverage(HbO, numTrials, numFrame);
HbRavg = trialAverage(HbR, numTrials, numFrame);
HbTavg = trialAverage(HbT, numTrials, numFrame);
GCaMPcorr = trialAverage(deltaGCaMPcorr, numTrials, numFrame);
count = 0;
fh = figure();
fh.WindowState = 'maximized';

for f = 8:10:120
    count = count+1;
    subplot(8,15,count)
    imagesc(100.*GCaMPcorr(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-5 5])

    subplot(8,15,count+15)
    imagesc(1e6.*HbTavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-4 4])

    subplot(8,15,count+30)
    imagesc(1e6.*HbOavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-7 7])

    subplot(8,15,count+45)
    imagesc(1e6.*HbRavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-3 3])
end
saveas(gcf,'figures/suppfig5_1.png')

t=2;
load([timepoints{t},'/','act_IOSI_ipsi.mat'])
load([timepoints{t},'/','act_GCaMPcorr_ipsi.mat'])
% Trial average hemodynamics and GCaMP (corrected and uncorrected)
HbOavg = trialAverage(HbO, numTrials, numFrame);
HbRavg = trialAverage(HbR, numTrials, numFrame);
HbTavg = trialAverage(HbT, numTrials, numFrame);
GCaMPcorr = trialAverage(deltaGCaMPcorr, numTrials, numFrame);
count = 0;
fh = figure();
fh.WindowState = 'maximized';

for f = 8:10:120
    count = count+1;
    subplot(8,15,count)
    imagesc(100.*GCaMPcorr(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-5 5])

    subplot(8,15,count+15)
    imagesc(1e6.*HbTavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-4 4])

    subplot(8,15,count+30)
    imagesc(1e6.*HbOavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-7 7])

    subplot(8,15,count+45)
    imagesc(1e6.*HbRavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-3 3])
end
saveas(gcf,'figures/suppfig5_2.png')

t=3;
load([timepoints{t},'/','act_IOSI_ipsi.mat'])
load([timepoints{t},'/','act_GCaMPcorr_ipsi.mat'])
% Trial average hemodynamics and GCaMP (corrected and uncorrected)
HbOavg = trialAverage(HbO, numTrials, numFrame);
HbRavg = trialAverage(HbR, numTrials, numFrame);
HbTavg = trialAverage(HbT, numTrials, numFrame);
GCaMPcorr = trialAverage(deltaGCaMPcorr, numTrials, numFrame);
count = 0;
fh = figure();
fh.WindowState = 'maximized';

for f = 8:10:120
    count = count+1;
    subplot(8,15,count)
    imagesc(100.*GCaMPcorr(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-5 5])

    subplot(8,15,count+15)
    imagesc(1e6.*HbTavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-4 4])

    subplot(8,15,count+30)
    imagesc(1e6.*HbOavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-7 7])

    subplot(8,15,count+45)
    imagesc(1e6.*HbRavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-3 3])
end
saveas(gcf,'figures/suppfig5_3.png')

t=5;
load([timepoints{t},'/','act_IOSI_ipsi.mat'])
load([timepoints{t},'/','act_GCaMPcorr_ipsi.mat'])
% Trial average hemodynamics and GCaMP (corrected and uncorrected)
HbOavg = trialAverage(HbO, numTrials, numFrame);
HbRavg = trialAverage(HbR, numTrials, numFrame);
HbTavg = trialAverage(HbT, numTrials, numFrame);
GCaMPcorr = trialAverage(deltaGCaMPcorr, numTrials, numFrame);
count = 0;
fh = figure();
fh.WindowState = 'maximized';

for f = 8:10:120
    count = count+1;
    subplot(8,15,count)
    imagesc(100.*GCaMPcorr(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-5 5])

    subplot(8,15,count+15)
    imagesc(1e6.*HbTavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-4 4])

    subplot(8,15,count+30)
    imagesc(1e6.*HbOavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-7 7])

    subplot(8,15,count+45)
    imagesc(1e6.*HbRavg(:,:,f))
    axis image
    axis off
    colormap jet
    caxis([-3 3])
end
saveas(gcf,'figures/suppfig5_4.png')
