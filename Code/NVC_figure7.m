%% Figure 7: Correlating metrics to behavior
clear
%% Figure 7a
% Behavior over time

behavior = xlsread('GCaMPmaster.xlsx', 'Sheet3');
count = 0;
figure
for m = 1:2:24
    count = count+1;
    pctuse = behavior(m,:)./(behavior(m,:)+behavior(m+1,:));
    asym1(count,:) = (pctuse-pctuse(1))./pctuse(1);
    scatter(1,pctuse(1))
    hold on
end
asym1 = 100.*asym1(:,[2 4 5 6]);
masym = mean(asym1,1);
sasym = std(asym1,1)./sqrt(12);
figure
subplot(1,2,1)
b = bar(1,masym);
b(1).FaceColor = [0.5 0 0];
b(2).FaceColor = [1 0 0];
b(3).FaceColor = [1 0.5 0];
b(4).FaceColor = [1 1 0];
hold on
ngroups = 1;
nbars = 4;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,masym(:,e),sasym(:,e),sasym(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
ylabel({'Percent change in'; 'impaired forelimb use'})
legend('Day 2','Week 1','Week 2','Week 4')
set(gca, 'XTickLabel', '', 'FontSize', 24);

subplot(1,2,2)
scatter(ones(12,1),asym1(:,4),200,'MarkerFaceColor','y','MarkerEdgeColor','k')
xlim([0 2])
ylim([-100 15])
box on
ylabel('Individual recovery at week4')
set(gca, 'XTickLabel', '', 'FontSize', 24);
saveas(gcf,'figures/figure7a.png')


%% Figure 7b
% Overlap between stroke core and pre-stroke forelimb area
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
for m = 1:12
    mask1 = [mouse{m},'brainmaskSFDI.mat'];
    load(mask1)
    stroke_forelimb_overlap(m) = dice(maskSFDI(1).ipsiOutline, double(maskSFDI(3).stroke_mask));
end
figure
b=4;
metric = 100.*stroke_forelimb_overlap;
[p,S] = polyfit(metric,asym1(:,b),1);
f = polyval(p,metric);
h = plot(metric,asym1(:,b),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(metric,f,'-','Color','k','LineWidth',2)
[r,p] = corrcoef(metric,asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2);
pval = p(1,2);
xlabel({'Percent overlap between'; 'forelimb and stroke'})
xlim([0 80])
ylabel({'Percent change in'; 'impaired forelimb use'})
ylim([-110 20])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/figure7b_1.png')

for m = 1:12
    mask1 = [mouse{m},'brainmaskSFDI.mat'];
    load(mask1)
    pixels_stroke(m) = length(maskSFDI(1).aff_prop_mus(maskSFDI(3).stroke_mask));
    area(m) = 52*52*10^-3*10^-3*pixels_stroke(m);
end
figure
b=4;
metric = area;
[p,S] = polyfit(metric,asym1(:,b),1);
f = polyval(p,metric);
h = plot(metric,asym1(:,b),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(metric,f,'-','Color','k','LineWidth',2)
[r,p] = corrcoef(metric,asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2);
pval = p(1,2);
xlabel('Stroke core area (mm^2)')
xlim([0 5])
ylabel({'Percent change in'; 'impaired forelimb use'})
ylim([-110 20])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/figure7b_2.png')

%% Figure 7c
% Behavior and response magnitude
% At week4 and week1

load('respMag.mat')
figure
t=3;
b=4;
metric = respMag(1).impFore(:,t);
[p,S] = polyfit(metric,asym1(:,b),1);
f = polyval(p,metric);
h = plot(metric,asym1(:,b),'o','MarkerFaceColor','g','MarkerEdgeColor','g');
hold on
plot(metric,f,'-','Color','g','LineWidth',2)
[r,p] = corrcoef(metric,asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

metric = respMag(4).impFore(:,t);
[p,S] = polyfit(metric,asym1(:,b),1);
f = polyval(p,metric);
h = plot(metric,asym1(:,b),'o','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(metric,f,'-','Color','r','LineWidth',2)
[r,p] = corrcoef(metric,asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

xlabel({'Response magniture in forelimb at week 1';'GCaMP \DeltaF/F (%) and \DeltaHbT (\muM)'})
xlim([-0.5 3.5])
ylabel({'Percent change in'; 'impaired forelimb use'})
ylim([-110 20])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/figure7c_1.png')


%% Figure 7c
% Behavior and neurovascular coupling correlation
load('GCaMP_Hb_corr_impAff.mat')
b=4;
th = 2;
time=3;
figure
black = rgb('black');
red = rgb('red');
scale = [black;red];
for m = 1:12
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbT(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));    
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbT(m,time,:,th));
    HbTpval(m) = p(1,2);
    HbTrval(m) = r(1,2);
    if HbTpval(m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(1,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(2,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([0 8])
    xlim([-2 8])
     
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbO(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbO(m,time,:,th));
    HbOpval(m) = p(1,2);
    HbOrval(m) = r(1,2);
    
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbR(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbR(m,time,:,th));
    HbRpval(m) = p(1,2);
    HbRrval(m) = r(1,2);
end
xlabel('\DeltaF/F (%)')
ylabel('\DeltaHbT (\muM)')
set(gca,'FontSize',24)
figure
scatter(ones(12,1),HbTrval,200,'MarkerFaceColor','r','MarkerEdgeColor','k')
xlim([0.5 1.5])
ylabel('Correlation coefficient')
ylim([-0.8 0.99])
set(gca, 'FontSize', 24);

asym = asym1(:,b);
figure
metric = HbTrval';
[p,S] = polyfit(metric,asym,1);
f = polyval(p,metric);
h = plot(metric,asym,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(metric,f,'-','Color','k','LineWidth',2)
[r,p] = corrcoef(metric,asym);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)
xlabel({'Evoked responses';'correlation coefficient at week 1'})
xlim([-0.4 1])
ylabel({'Percent change in'; 'impaired forelimb use'})
ylim([-110 20])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/figure7c_2.png')

%% Figure 7d
load('NVC_stats.mat')
figure
t=3;
b=4;
correlation = squeeze(corr_roi(2,t,:));
[p,S] = polyfit(correlation,asym1(:,b),1);
f = polyval(p,correlation);
h = plot(correlation,asym1(:,b),'o','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(correlation,f,'-','Color','k','LineWidth',2)
[r,p] = corrcoef(correlation,asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)
xlabel({'Neurovascular coupling';'correlation coefficient at week 1'})
xlim([0.25 0.75])
ylabel({'Percent change in'; 'impaired forelimb use'})
ylim([-110 20])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/figure7d.png')

%% Figure 7e
% Behavior and oscillations

load('spect_rest_HbT')
start = 73; %freq=0.1-49
stop = 145; %freq=0.3-145
pwrbandG_ipsi = sum(powerG_ipsi(:,:,start:stop),3);
pwrbandH_ipsi = sum(powerHb_ipsi(:,:,start:stop),3);
pwrbandG_contra = sum(powerG_contra(:,:,start:stop),3);
pwrbandH_contra = sum(powerHb_contra(:,:,start:stop),3);

figure
t=3;
b=4;
[p,S] = polyfit(pwrbandG_ipsi(:,t),asym1(:,b),1);
f = polyval(p,pwrbandG_ipsi(:,t));
h = plot(pwrbandG_ipsi(:,t),asym1(:,b),'o','MarkerFaceColor','g','MarkerEdgeColor','g');
hold on
plot(pwrbandG_ipsi(:,t),f,'-','Color','g','LineWidth',2)
[r,p] = corrcoef(pwrbandG_ipsi(:,t),asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)
[p,S] = polyfit(pwrbandG_contra(:,t),asym1(:,b),1);
f = polyval(p,pwrbandG_contra(:,t));
h = plot(pwrbandG_contra(:,t),asym1(:,b),'o','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880]);
hold on
plot(pwrbandG_contra(:,t),f,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
[r,p] = corrcoef(pwrbandG_contra(:,t),asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)
ylabel({'Percent change in'; 'impaired forelimb use'})
ylim([-110 20])
xlim([0.5e-3 5.5e-3])
xlabel('GCaMP Power: 0.15 - 0.3 Hz')
set(gca, 'FontSize', 24);
saveas(gcf,'figures/figure7e_1.png')

figure
[p,S] = polyfit(pwrbandH_ipsi(:,t),asym1(:,b),1);
f = polyval(p,pwrbandH_ipsi(:,t));
h = plot(pwrbandH_ipsi(:,t),asym1(:,b),'o','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(pwrbandH_ipsi(:,t),f,'-','Color','r','LineWidth',2)
[r,p] = corrcoef(pwrbandH_ipsi(:,t),asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)
[p,S] = polyfit(pwrbandH_contra(:,t),asym1(:,b),1);
f = polyval(p,pwrbandH_contra(:,t));
h = plot(pwrbandH_contra(:,t),asym1(:,b),'o','MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerEdgeColor',[0.6350 0.0780 0.1840]);
hold on
plot(pwrbandH_contra(:,t),f,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
[r,p] = corrcoef(pwrbandH_contra(:,t),asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)
ylabel({'Percent change in'; 'impaired forelimb use'})
xlabel('HbT Power at week 1')
ylim([-110 20])
xlim([1e-11 15e-11])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/figure7e_2.png')

%% Behavior and resting state

load('RSFC_analysis_lowHz.mat')
load('bregma_lambda.mat')
mouse = {'../MouseData/SS75/', ...
    '../MouseData/SS76/',...
    '../MouseData/SS77/',...
    '../MouseData/SS78/',...
    '../MouseData/SS79/',...
    '../MouseData/SS80/',...
    '../MouseData/SS81/',...
    '../MouseData/SS82/',...
    '../MouseData/SS83/',...
    '../MouseData/SS84/',...
    '../MouseData/SS85/',...
    '../MouseData/SS93/'};

thresh = 0:0.1:1;
for m = 1:12
    mask = [mouse{m},'brainmaskSFDI.mat'];
    load(mask)
    for t = 1:5        
        mask1 = maskSFDI(t).aff_mask;
        for th = 1:length(thresh)
            tempG = mask1.*LRC_GCaMP{m,t};
            tempG(tempG<=thresh(th)) = NaN;
            areaG(m,t,th) = sum(sum(~isnan(tempG)))./sum(mask1(:));

            tempH = mask1.*LRC_HbO{m,t};
            tempH(tempH<=thresh(th)) = NaN;
            areaH(m,t,th) = sum(sum(~isnan(tempH)))./sum(mask1(:));
            
            tempG(~isnan(tempG)) = 1;
            tempG(isnan(tempG)) = 0;           
            tempH(~isnan(tempH)) = 1;
            tempH(isnan(tempH)) = 0;            
            if isempty(dice(tempG,tempH)) == 1
                overlap_LRC(m,t,th) = 0;
            else
                overlap_LRC(m,t,th) = dice(tempG,tempH);
            end
        end
    end
end

areaG(isnan(areaG)) = 0;
areaH(isnan(areaH)) = 0;
th = 5;
areaG = squeeze(areaG(:,:,th));
areaH = squeeze(areaH(:,:,th));

figure
t=2;
b=4;
metric = areaG(:,t);
[p,S] = polyfit(metric,asym1(:,b),1);
f = polyval(p,metric);
h = plot(metric,asym1(:,b),'o','MarkerFaceColor','g','MarkerEdgeColor','g');
hold on
plot(metric,f,'-','Color','g','LineWidth',2)
[r,p] = corrcoef(metric,asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

metric = areaH(:,t);
[p,S] = polyfit(metric,asym1(:,b),1);
f = polyval(p,metric);
h = plot(metric,asym1(:,b),'o','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(metric,f,'-','Color','r','LineWidth',2)
[r,p] = corrcoef(metric,asym1(:,b));
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

ylabel({'Percent change in'; 'impaired forelimb use'})
xlabel('Proportional area at th = 0.4')
ylim([-110 20])
xlim([0 0.3])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/figure7f.png')

%% Supplementary figure 11
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
for m = 1:12
    mask1 = [mouse{m},'brainmaskSFDI.mat'];
    load(mask1)
    pixels_stroke(m) = length(maskSFDI(1).aff_prop_mus(maskSFDI(3).stroke_mask));
    area(m) = 52*52*10^-3*10^-3*pixels_stroke(m);
end

%% 11a
% Behavior and response magnitude
% At week4 and week1
load('respMag.mat')
figure
t=3;
b=4;
metric = respMag(1).impFore(:,t);
[p,S] = polyfit(area,metric,1);
f = polyval(p,area);
h = plot(area,metric,'o','MarkerFaceColor','g','MarkerEdgeColor','g');
hold on
plot(area,f,'-','Color','g','LineWidth',2)
[r,p] = corrcoef(area,metric);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

metric = respMag(4).impFore(:,t);
[p,S] = polyfit(area,metric,1);
f = polyval(p,area);
h = plot(area,metric,'o','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(area,f,'-','Color','r','LineWidth',2)
[r,p] = corrcoef(area,metric);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

ylabel({'Response magniture in forelimb at week 1';'GCaMP \DeltaF/F (%) and \DeltaHbT (\muM)'})
ylim([-0.5 3.5])
xlabel('Stroke core area (mm^2)')
xlim([0 5])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/suppfig11a.png')

%% 11b

load('GCaMP_Hb_corr_impAff.mat')
b=4;
th = 2;
time=3;
figure
black = rgb('black');
red = rgb('red');
scale = [black;red];
for m = 1:12
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbT(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));    
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbT(m,time,:,th));
    HbTpval(m) = p(1,2);
    HbTrval(m) = r(1,2);
    if HbTpval(m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(1,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(2,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([0 8])
    xlim([-2 8])
     
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbO(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbO(m,time,:,th));
    HbOpval(m) = p(1,2);
    HbOrval(m) = r(1,2);
    
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbR(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbR(m,time,:,th));
    HbRpval(m) = p(1,2);
    HbRrval(m) = r(1,2);
end
xlabel('\DeltaF/F (%)')
ylabel('\DeltaHbT (\muM)')
set(gca,'FontSize',24)
figure
scatter(ones(12,1),HbTrval,200,'MarkerFaceColor','r','MarkerEdgeColor','k')
xlim([0.5 1.5])
ylabel('Correlation coefficient')
ylim([-0.8 0.99])
set(gca, 'FontSize', 24);

metric = HbTrval;
[p,S] = polyfit(area,metric,1);
f = polyval(p,area);
h = plot(area,metric,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(area,f,'-','Color','k','LineWidth',2)
[r,p] = corrcoef(area,metric);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

ylabel({'Evoked responses';'correlation coefficient at week 1'})
ylim([-0.5 1])
xlabel('Stroke core area (mm^2)')
xlim([0 5])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/suppfig11b.png')

%% 11c
load('NVC_stats.mat')
figure
t=3;
b=4;
correlation = squeeze(corr_roi(1,t,:));
metric = correlation;
[p,S] = polyfit(area,metric,1);
f = polyval(p,area);
h = plot(area,metric,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(area,f,'-','Color','k','LineWidth',2)
[r,p] = corrcoef(area,metric);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)
ylabel({'Neurovascular coupling';'correlation coefficient at week 1'})
ylim([0.25 0.75])
xlabel('Stroke core area (mm^2)')
xlim([0 5])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/suppfig11c.png')

%% 11d
load('NVC_stats.mat')
figure
t=5;
b=4;
correlation = squeeze(corr_roi(1,t,:));
metric = correlation;
[p,S] = polyfit(area,metric,1);
f = polyval(p,area);
h = plot(area,metric,'o','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(area,f,'-','Color','k','LineWidth',2)
[r,p] = corrcoef(area,metric);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)
ylabel({'Neurovascular coupling';'correlation coefficient at week 1'})
ylim([0.2 0.8])
xlabel('Stroke core area (mm^2)')
xlim([0 5])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/suppfig11d.png')

%% 11e

load('spect_rest_HbT')
start = 73; %freq=0.1-49
stop = 145; %freq=0.3-145
pwrbandG_ipsi = sum(powerG_ipsi(:,:,start:stop),3);
pwrbandH_ipsi = sum(powerHb_ipsi(:,:,start:stop),3);
pwrbandG_contra = sum(powerG_contra(:,:,start:stop),3);
pwrbandH_contra = sum(powerHb_contra(:,:,start:stop),3);

figure
t = 3;
metric = pwrbandG_ipsi(:,t);
[p,S] = polyfit(area,metric,1);
f = polyval(p,area);
h = plot(area,metric,'o','MarkerFaceColor','g','MarkerEdgeColor','g');
hold on
plot(area,f,'-','Color','g','LineWidth',2)
[r,p] = corrcoef(area,metric);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

metric = pwrbandG_contra(:,t);
[p,S] = polyfit(area,metric,1);
f = polyval(p,area);
h = plot(area,metric,'o','MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880]);
hold on
plot(area,f,'-','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
[r,p] = corrcoef(area,metric);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

ylim([0.5e-3 5.5e-3])
ylabel('GCaMP Power: 0.15 - 0.3 Hz')
xlabel('Stroke core area (mm^2)')
xlim([0 5])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/suppfig11e.png')

%% 11f
figure
t = 3;
metric = pwrbandH_ipsi(:,t);
[p,S] = polyfit(area,metric,1);
f = polyval(p,area);
h = plot(area,metric,'o','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(area,f,'-','Color','r','LineWidth',2)
[r,p] = corrcoef(area,metric);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

metric = pwrbandH_contra(:,t);
[p,S] = polyfit(area,metric,1);
f = polyval(p,area);
h = plot(area,metric,'o','MarkerFaceColor',[0.6350 0.0780 0.1840],'MarkerEdgeColor',[0.6350 0.0780 0.1840]);
hold on
plot(area,f,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',2)
[r,p] = corrcoef(area,metric);
z = 0.5*(log(1+r) - log(1-r));
cc = r(1,2)
pval = p(1,2)

ylim([1e-11 15e-11])
ylabel('HbT Power: 0.15 - 0.3 Hz')
xlabel('Stroke core area (mm^2)')
xlim([0 5])
set(gca, 'FontSize', 24);
saveas(gcf,'figures/suppfig11f.png')
