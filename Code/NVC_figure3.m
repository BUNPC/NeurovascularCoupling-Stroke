%% Stroke NVC figure 3
% Smrithi Sunil
% BOAS Lab

clear
%% figure 3a

mousename = 'SS75';
[timepoints, mask] = animals(mousename);
load(mask)

tp=1;
respmask = maskSFDI(1).ipsiOutline;       
load([timepoints{tp},'/','act_IOSI_ipsi.mat'])
load([timepoints{tp},'/','act_GCaMPcorr_ipsi.mat'])
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

Greshape_new = reshape(Greshape,[size(Greshape,1)*size(Greshape,2)*size(Greshape,3) size(Greshape,4)]);
Greshape_new(any(isnan(Greshape_new), 2), :) = [];
HbOreshape_new = reshape(HbOreshape,[size(HbOreshape,1)*size(HbOreshape,2)*size(HbOreshape,3) size(HbOreshape,4)]);
HbOreshape_new(any(isnan(HbOreshape_new), 2), :) = [];
HbRreshape_new = reshape(HbRreshape,[size(HbOreshape,1)*size(HbOreshape,2)*size(HbOreshape,3) size(HbOreshape,4)]);
HbRreshape_new(any(isnan(HbRreshape_new), 2), :) = [];

xlimG = [-5 25];
ylimG = [-5 10];
xlimH = [-5 25];
ylimH = [-8 12];
time = -5+0.2:0.2:25;
figure
fig = gcf;
options.x_axis = time;
options.handle = figure(fig);
options.error = 'std';
options.color_area = [0.1 0.7 0.1]; 
options.color_line = [0.1 0.7 0.1]; 
options.alpha = 0.3;
options.line_width = 3;
plot_areaerrorbar(100.*Greshape_new,options)
hold on
xlim(xlimG)
ylim(ylimG)
plot([0 5],[0 0],'k-','LineWidth',5)
ylabel('\DeltaF/F (%)')
xlabel('Time (sec)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure3a_1.png')

figure
fig = gcf;
options.handle = figure(fig);
options.color_area = [0.8500 0.3250 0.0980]; 
options.color_line = [0.8500 0.3250 0.0980]; 
plot_areaerrorbar(1e6.*HbOreshape_new,options)
hold on
options.color_area = [0 0.4470 0.7410]; 
options.color_line = [0 0.4470 0.7410]; 
plot_areaerrorbar(1e6.*HbRreshape_new,options)
hold on
xlim(xlimH)
ylim(ylimH)
plot([0 5],[0 0],'k-','LineWidth',5)
ylabel('\DeltaHb (\muM)')
xlabel('Time (sec)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure3a_2.png')

tp=3;
load([timepoints{tp},'/','act_IOSI_ipsi.mat'])
load([timepoints{tp},'/','act_GCaMPcorr_ipsi.mat'])
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

Greshape_new = reshape(Greshape,[size(Greshape,1)*size(Greshape,2)*size(Greshape,3) size(Greshape,4)]);
Greshape_new(any(isnan(Greshape_new), 2), :) = [];
HbOreshape_new = reshape(HbOreshape,[size(HbOreshape,1)*size(HbOreshape,2)*size(HbOreshape,3) size(HbOreshape,4)]);
HbOreshape_new(any(isnan(HbOreshape_new), 2), :) = [];
HbRreshape_new = reshape(HbRreshape,[size(HbOreshape,1)*size(HbOreshape,2)*size(HbOreshape,3) size(HbOreshape,4)]);
HbRreshape_new(any(isnan(HbRreshape_new), 2), :) = [];

time = -5+0.2:0.2:25;
figure
fig = gcf;
options.x_axis = time;
options.handle = figure(fig);
options.color_area = [0.1 0.7 0.1]; 
options.color_line = [0.1 0.7 0.1]; 
options.alpha = 0.3;
options.line_width = 3;
plot_areaerrorbar(100.*Greshape_new,options)
hold on
xlim(xlimG)
ylim(ylimG)
plot([0 5],[0 0],'k-','LineWidth',5)
ylabel('\DeltaF/F (%)')
xlabel('Time (sec)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure3a_3.png')

figure
fig = gcf;
options.handle = figure(fig);
options.color_area = [0.8500 0.3250 0.0980]; 
options.color_line = [0.8500 0.3250 0.0980]; 
plot_areaerrorbar(1e6.*HbOreshape_new,options)
hold on
options.color_area = [0 0.4470 0.7410]; 
options.color_line = [0 0.4470 0.7410]; 
plot_areaerrorbar(1e6.*HbRreshape_new,options)
hold on
xlim(xlimH)
ylim(ylimH)
plot([0 5],[0 0],'k-','LineWidth',5)
ylabel('\DeltaHb (\muM)')
xlabel('Time (sec)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure3a_4.png')

tp=5;
load([timepoints{tp},'/','act_IOSI_ipsi.mat'])
load([timepoints{tp},'/','act_GCaMPcorr_ipsi.mat'])
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

Greshape_new = reshape(Greshape,[size(Greshape,1)*size(Greshape,2)*size(Greshape,3) size(Greshape,4)]);
Greshape_new(any(isnan(Greshape_new), 2), :) = [];
HbOreshape_new = reshape(HbOreshape,[size(HbOreshape,1)*size(HbOreshape,2)*size(HbOreshape,3) size(HbOreshape,4)]);
HbOreshape_new(any(isnan(HbOreshape_new), 2), :) = [];
HbRreshape_new = reshape(HbRreshape,[size(HbOreshape,1)*size(HbOreshape,2)*size(HbOreshape,3) size(HbOreshape,4)]);
HbRreshape_new(any(isnan(HbRreshape_new), 2), :) = [];

time = -5+0.2:0.2:25;
figure
fig = gcf;
options.x_axis = time;
options.handle = figure(fig);
options.color_area = [0.1 0.7 0.1]; 
options.color_line = [0.1 0.7 0.1]; 
options.alpha = 0.3;
options.line_width = 3;
plot_areaerrorbar(100.*Greshape_new,options)
hold on
xlim(xlimG)
ylim(ylimG)
plot([0 5],[0 0],'k-','LineWidth',5)
ylabel('\DeltaF/F (%)')
xlabel('Time (sec)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure3a_5.png')

figure
fig = gcf;
options.handle = figure(fig);
options.color_area = [0.8500 0.3250 0.0980]; 
options.color_line = [0.8500 0.3250 0.0980]; 
plot_areaerrorbar(1e6.*HbOreshape_new,options)
hold on
options.color_area = [0 0.4470 0.7410]; 
options.color_line = [0 0.4470 0.7410]; 
plot_areaerrorbar(1e6.*HbRreshape_new,options)
hold on
xlim(xlimH)
ylim(ylimH)
plot([0 5],[0 0],'k-','LineWidth',5)
ylabel('\DeltaHb (\muM)')
xlabel('Time (sec)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure3a_6.png')


%% Figure 3c
% Similarity in response area

load('GCaMP_Hb_corr_impAff.mat')
figure
C = imfuse(squeeze(thGCaMP{1,1,2,2}),squeeze(thHbO{1,1,2,2}),'ColorChannels',[2 1 0]);
imshow(C)
saveas(gcf,'figures/figure3c_left.png')

for m = 1:12
    for t = 1:5
        magHbO(m,t) = squeeze(mean(simO(m,t,:,2)));
        magHbR(m,t) = squeeze(mean(simR(m,t,:,2)));
        magHbT(m,t) = squeeze(mean(simT(m,t,:,2)));
    end
end

for t = 1:5
    mmagHbO(t) = mean(magHbO(:,t));
    smagHbO(t) = std(magHbO(:,t));
    mmagHbR(t) = mean(magHbR(:,t));
    smagHbR(t) = std(magHbR(:,t));
    mmagHbT(t) = mean(magHbT(:,t));
    smagHbT(t) = std(magHbT(:,t));
end
mmag = [mmagHbT; mmagHbO; mmagHbR];
smag = [smagHbT; smagHbO; smagHbR];
a = 0.05;
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
b = bar(1:3,mmag);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = 3;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mmag(:,e),smag(:,e),smag(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'HbT','HbO','HbR'};
ylabel({'Dice similarity coefficient';'of response area with GCaMP'})
set(b, 'FaceAlpha', 1)
ylim([0 1])

sigline1([0.69 0.85],[],0.66)
sigline1([0.69 1],[],0.7)
sigline1([0.69 1.15],[],0.74)
sigline2([0.69 1.3],[],0.79)
sigline1([1.69 1.85],[],0.73)
sigline1([1.69 2],[],0.77)
sigline2([1.69 2.15],[],0.81)
sigline2([1.69 2.3],[],0.86)
legend('Pre-stroke','Day2','Week1','Week2','Week4')
set(gca, 'XTickLabel', group, 'FontSize', 24);
saveas(gcf,'figures/figure3c_right.png')

%% figure 3d

color{1} = [0.8500 0.3250 0.0980];
color{2} = [0 0.4470 0.7410];
m = 1;
th = 2;

figure
t=1;
[p,S] = polyfit(mGCaMP(m,t,:,th),mHbO(m,t,:,th),1);
f = polyval(p,mGCaMP(m,t,:,th));
h = plot(squeeze(mGCaMP(m,t,:,th)),squeeze(mHbO(m,t,:,th)),'o','MarkerFaceColor',color{1},'MarkerEdgeColor',color{1});
hold on
plot(squeeze(mGCaMP(m,t,:,th)),squeeze(f),'-','Color',color{1},'LineWidth',2)
[r,p] = corrcoef(mGCaMP(m,t,:,th),mHbO(m,t,:,th));
z = 0.5*(log(1+r) - log(1-r));
[p,S] = polyfit(mGCaMP(m,t,:,th),mHbR(m,t,:,th),1);
f = polyval(p,mGCaMP(m,t,:,th));
h = plot(squeeze(mGCaMP(m,t,:,th)),squeeze(mHbR(m,t,:,th)),'o','MarkerFaceColor',color{2},'MarkerEdgeColor',color{2});
hold on
plot(squeeze(mGCaMP(m,t,:,th)),squeeze(f),'-','Color',color{2},'LineWidth',2)
[r,p] = corrcoef(mGCaMP(m,t,:,th),mHbR(m,t,:,th));
z = 0.5*(log(1+r) - log(1-r));  
xlim([0 8])
ylim([-6 10])
ylabel('\DeltaHb (\muM)')
xlabel('\DeltaF/F (%)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure3d_top.png')

figure
t=3;
[p,S] = polyfit(mGCaMP(m,t,:,th),mHbO(m,t,:,th),1);
f = polyval(p,mGCaMP(m,t,:,th));
h = plot(squeeze(mGCaMP(m,t,:,th)),squeeze(mHbO(m,t,:,th)),'o','MarkerFaceColor',color{1},'MarkerEdgeColor',color{1});
hold on
plot(squeeze(mGCaMP(m,t,:,th)),squeeze(f),'-','Color',color{1},'LineWidth',2)
[r,p] = corrcoef(mGCaMP(m,t,:,th),mHbO(m,t,:,th));
z = 0.5*(log(1+r) - log(1-r));
[p,S] = polyfit(mGCaMP(m,t,:,th),mHbR(m,t,:,th),1);
f = polyval(p,mGCaMP(m,t,:,th));
h = plot(squeeze(mGCaMP(m,t,:,th)),squeeze(mHbR(m,t,:,th)),'o','MarkerFaceColor',color{2},'MarkerEdgeColor',color{2});
hold on
plot(squeeze(mGCaMP(m,t,:,th)),squeeze(f),'-','Color',color{2},'LineWidth',2)
[r,p] = corrcoef(mGCaMP(m,t,:,th),mHbR(m,t,:,th));
z = 0.5*(log(1+r) - log(1-r));  
xlim([0 8])
ylim([-6 10])
ylabel('\DeltaHb (\muM)')
xlabel('\DeltaF/F (%)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure3d_middle.png')

figure
t=5;
[p,S] = polyfit(mGCaMP(m,t,:,th),mHbO(m,t,:,th),1);
f = polyval(p,mGCaMP(m,t,:,th));
h = plot(squeeze(mGCaMP(m,t,:,th)),squeeze(mHbO(m,t,:,th)),'o','MarkerFaceColor',color{1},'MarkerEdgeColor',color{1});
hold on
plot(squeeze(mGCaMP(m,t,:,th)),squeeze(f),'-','Color',color{1},'LineWidth',2)
[r,p] = corrcoef(mGCaMP(m,t,:,th),mHbO(m,t,:,th));
z = 0.5*(log(1+r) - log(1-r));
[p,S] = polyfit(mGCaMP(m,t,:,th),mHbR(m,t,:,th),1);
f = polyval(p,mGCaMP(m,t,:,th));
h = plot(squeeze(mGCaMP(m,t,:,th)),squeeze(mHbR(m,t,:,th)),'o','MarkerFaceColor',color{2},'MarkerEdgeColor',color{2});
hold on
plot(squeeze(mGCaMP(m,t,:,th)),squeeze(f),'-','Color',color{2},'LineWidth',2)
[r,p] = corrcoef(mGCaMP(m,t,:,th),mHbR(m,t,:,th));
z = 0.5*(log(1+r) - log(1-r));  
xlim([0 8])
ylim([-6 10])
ylabel('\DeltaHb (\muM)')
xlabel('\DeltaF/F (%)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure3d_bottom.png')

%% Plotting all animals individually

black = rgb('black');
red = rgb('red');
scale = [black;red];

th = 2;
fh = figure();
fh.WindowState = 'maximized';
for time = 1:5
for m = 1:12
    subplot(3,5,time)
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbT(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));    
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbT(m,time,:,th));
    HbTpval(time,m) = p(1,2);
    if HbTpval(time,m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(1,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(2,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([0 8])
    xlim([-2 8])
    set(gca,'FontSize',20)
    
    subplot(3,5,time+5)
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbO(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbO(m,time,:,th));
    HbOpval(time,m) = p(1,2);
    if HbOpval(time,m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(1,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(2,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([0 12])
    xlim([-2 8])
    set(gca,'FontSize',20)
    
    subplot(3,5,time+10)
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbR(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbR(m,time,:,th));
    HbRpval(time,m) = p(1,2);
    if HbRpval(time,m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(1,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(2,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([-8 0])
    xlim([-2 8])
    set(gca,'FontSize',20)
end
end
saveas(gcf,'figures/figure3e.png')

%% Plotting all animals individually by color
% Supplementary figure 6

color = jet(12);

th = 2;
fh = figure();
fh.WindowState = 'maximized';
for time = 1:5
for m = 1:12
    subplot(3,5,time)
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbT(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));    
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbT(m,time,:,th));
    HbTpval(time,m) = p(1,2);
    if HbTpval(time,m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Marker','none','Color',color(m,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Marker','s','MarkerFaceColor',color(m,:),'Color',color(m,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([0 8])
    xlim([-2 8])
    set(gca,'FontSize',20)
    
    subplot(3,5,time+5)
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbO(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbO(m,time,:,th));
    HbOpval(time,m) = p(1,2);
    if HbOpval(time,m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Marker','none','Color',color(m,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Marker','s','MarkerFaceColor',color(m,:),'Color',color(m,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([0 12])
    xlim([-2 8])
    set(gca,'FontSize',20)
    
    subplot(3,5,time+10)
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbR(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbR(m,time,:,th));
    HbRpval(time,m) = p(1,2);
    if HbRpval(time,m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Marker','none','Color',color(m,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Marker','s','MarkerFaceColor',color(m,:),'Color',color(m,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([-8 0])
    xlim([-2 8])
    set(gca,'FontSize',20)
end
end
saveas(gcf,'figures/suppfig6.png')

%% Supplementary figure 7
% Unaffected hemisphere
load('GCaMP_Hb_corr_unimpUnaff.mat')

for m = 1:12
    for t = 1:5
        magHbO(m,t) = squeeze(mean(simO(m,t,:,2)));
        magHbR(m,t) = squeeze(mean(simR(m,t,:,2)));
        magHbT(m,t) = squeeze(mean(simT(m,t,:,2)));
    end
end

for t = 1:5
    mmagHbO(t) = mean(magHbO(:,t));
    smagHbO(t) = std(magHbO(:,t));
    mmagHbR(t) = mean(magHbR(:,t));
    smagHbR(t) = std(magHbR(:,t));
    mmagHbT(t) = mean(magHbT(:,t));
    smagHbT(t) = std(magHbT(:,t));
end
mmag = [mmagHbT; mmagHbO; mmagHbR];
smag = [smagHbT; smagHbO; smagHbR];
a = 0.05;
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
b = bar(1:3,mmag);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = 3;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mmag(:,e),smag(:,e),smag(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'HbT','HbO','HbR'};
ylabel({'Dice similarity coefficient';'of response area with GCaMP'})
set(b, 'FaceAlpha', 1)
ylim([0 1])

legend('Pre-stroke','Day2','Week1','Week2','Week4')
set(gca, 'XTickLabel', group, 'FontSize', 24);
saveas(gcf,'figures/suppfig7a.png')

%% Plotting all animals individually

black = rgb('black');
red = rgb('red');
scale = [black;red];

th = 2;
fh = figure();
fh.WindowState = 'maximized';
for time = 1:5
for m = 1:12
    subplot(3,5,time)
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbT(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));    
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbT(m,time,:,th));
    HbTpval(time,m) = p(1,2);
    if HbTpval(time,m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(1,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(2,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([0 8])
    xlim([-2 8])
    set(gca,'FontSize',20)
    
    subplot(3,5,time+5)
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbO(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbO(m,time,:,th));
    HbOpval(time,m) = p(1,2);
    if HbOpval(time,m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(1,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(2,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([0 12])
    xlim([-2 8])
    set(gca,'FontSize',20)
    
    subplot(3,5,time+10)
    [p,S] = polyfit(mGCaMP(m,time,:,th),mHbR(m,time,:,th),1);
    f = polyval(p,mGCaMP(m,time,:,th));
    [r,p] = corrcoef(mGCaMP(m,time,:,th),mHbR(m,time,:,th));
    HbRpval(time,m) = p(1,2);
    if HbRpval(time,m) < 0.05                               % significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(1,:),'LineWidth',2)
    else                                                     % not significant
        plot(squeeze(mGCaMP(m,time,:,th)),squeeze(f),'Color',scale(2,:),'LineWidth',2)
    end
    hold on
    z = 0.5*(log(1+r) - log(1-r));
    ylim([-8 0])
    xlim([-2 8])
    set(gca,'FontSize',20)
end
end
saveas(gcf,'figures/suppfig7b.png')
