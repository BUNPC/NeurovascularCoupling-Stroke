%% Figure 5
% Power spectrum of non-infarct and contralesion
clear
%% Figure 5a
% Time course of signal, pre-stroke and day 2

% ipsilesional
fs = 5;
lpf = 0.4; % 2; %0.08;
hpf = 0.009; %0.02; %0.009;
wn(2) = lpf/(fs/2);
wn(1) = hpf/(fs/2);
[fb,fa] = butter(2,wn,'bandpass');
animal_number = 'SS76';
[timepoints, mask] = animals(animal_number);
load(mask)
dataDir = '../MouseData/SS76/RestingState/Baseline';
t=1;
load([dataDir,'/','1_v2_rsfc_I.mat'])
[nY,nX,nLam,nT] = size(Iunfilt);
brain_mask = maskSFDI(t).aff_mask + maskSFDI(t).unaff_mask;
brainIdx = find(brain_mask > 0);
yAll = reshape(Iunfilt(:,:,1,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
Chb1 = zeros(nY*nX,nT);
Chb1(brainIdx,:) = ynew';
Chb1 = reshape(Chb1, [nY, nX, nT]);
HbT = Chb1;
rawsigH = reshape(HbT,[nY*nX nT]);
filt_rawsigH = filtfilt(fb,fa,rawsigH')';
filt_rawsigH = reshape(filt_rawsigH,[nY nX nT]);
    
load([dataDir,'/','gcampcorr_v2_rsfc_I.mat'])
[nY,nX,nLam,nT] = size(Iunfilt);
yAll = reshape(Iunfilt(:,:,1,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
Chb1 = zeros(nY*nX,nT);
Chb1(brainIdx,:) = ynew';
Chb1 = reshape(Chb1, [nY, nX, nT]);
deltaGCaMPcorr = Chb1;
rawsigG = reshape(deltaGCaMPcorr,[nY*nX nT]);
filt_rawsigG = filtfilt(fb,fa,rawsigG')';
filt_rawsigG = reshape(filt_rawsigG,[nY nX nT]);
pos = [98 64 25 24];
HbOroi = filt_rawsigH(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
GCaMProi = filt_rawsigG(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
mHbOroi_base = 1e6.*squeeze(mean(mean(HbOroi,1),2));
mGCaMProi_base = 100.*squeeze(mean(mean(GCaMProi,1),2));

dataDir = '../MouseData/SS76/RestingState/Day2';
t=2;
load([dataDir,'/','1_v2_rsfc_I.mat'])
[nY,nX,nLam,nT] = size(Iunfilt);
brain_mask = maskSFDI(t).aff_mask + maskSFDI(t).unaff_mask;
brainIdx = find(brain_mask > 0);
yAll = reshape(Iunfilt(:,:,1,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
Chb1 = zeros(nY*nX,nT);
Chb1(brainIdx,:) = ynew';
Chb1 = reshape(Chb1, [nY, nX, nT]);
HbT = Chb1;
rawsigH = reshape(HbT,[nY*nX nT]);
filt_rawsigH = filtfilt(fb,fa,rawsigH')';
filt_rawsigH = reshape(filt_rawsigH,[nY nX nT]);
    
load([dataDir,'/','gcampcorr_v2_rsfc_I.mat'])
[nY,nX,nLam,nT] = size(Iunfilt);
yAll = reshape(Iunfilt(:,:,1,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
Chb1 = zeros(nY*nX,nT);
Chb1(brainIdx,:) = ynew';
Chb1 = reshape(Chb1, [nY, nX, nT]);
deltaGCaMPcorr = Chb1;
rawsigG = reshape(deltaGCaMPcorr,[nY*nX nT]);
filt_rawsigG = filtfilt(fb,fa,rawsigG')';
filt_rawsigG = reshape(filt_rawsigG,[nY nX nT]);
HbOroi = filt_rawsigH(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
GCaMProi = filt_rawsigG(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
mHbOroi_day2 = 1e6.*squeeze(mean(mean(HbOroi,1),2));
mGCaMProi_day2 = 100.*squeeze(mean(mean(GCaMProi,1),2));

Fs = 5;
T = 1/Fs;
L = size(HbT,3);
timefull = (1:(L))*T;
figure
plot(timefull,mGCaMProi_base,'Color',rgb('dark green'),'LineWidth',2)
hold on
plot(timefull,mHbOroi_base,'Color',rgb('dark red'),'LineWidth',2)
plot(timefull,mGCaMProi_day2,'Color',rgb('green'),'LineWidth',2)
plot(timefull,mHbOroi_day2,'Color',rgb('light red'),'LineWidth',2)
ylabel('\DeltaF/F (%)   \DeltaHbO (\muM)')
xlim([0 300])
ylim([-4 4])
legend('Pre-stroke: GCaMP','Pre-stroke: HbO','Day2: GCaMP','Day2: HbO')
set(gca,'FontSize',24,'XTick','')
saveas(gcf,'figures/figure5a_top.png')

% contralesional
fs = 5;
lpf = 0.4; % 2; %0.08;
hpf = 0.009; %0.02; %0.009;
wn(2) = lpf/(fs/2);
wn(1) = hpf/(fs/2);
[fb,fa] = butter(2,wn,'bandpass');
animal_number = 'SS76';
[timepoints, mask] = animals(animal_number);
load(mask)
dataDir = '../MouseData/SS76/RestingState/Baseline';
t=1;
load([dataDir,'/','1_v2_rsfc_I.mat'])
[nY,nX,nLam,nT] = size(Iunfilt);
brain_mask = maskSFDI(t).aff_mask + maskSFDI(t).unaff_mask;
brainIdx = find(brain_mask > 0);
yAll = reshape(Iunfilt(:,:,1,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
Chb1 = zeros(nY*nX,nT);
Chb1(brainIdx,:) = ynew';
Chb1 = reshape(Chb1, [nY, nX, nT]);
HbT = Chb1;
rawsigH = reshape(HbT,[nY*nX nT]);
filt_rawsigH = filtfilt(fb,fa,rawsigH')';
filt_rawsigH = reshape(filt_rawsigH,[nY nX nT]);
    
load([dataDir,'/','gcampcorr_v2_rsfc_I.mat'])
[nY,nX,nLam,nT] = size(Iunfilt);
yAll = reshape(Iunfilt(:,:,1,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
Chb1 = zeros(nY*nX,nT);
Chb1(brainIdx,:) = ynew';
Chb1 = reshape(Chb1, [nY, nX, nT]);
deltaGCaMPcorr = Chb1;
rawsigG = reshape(deltaGCaMPcorr,[nY*nX nT]);
filt_rawsigG = filtfilt(fb,fa,rawsigG')';
filt_rawsigG = reshape(filt_rawsigG,[nY nX nT]);
pos = [3 64 20 31];
HbOroi = filt_rawsigH(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
GCaMProi = filt_rawsigG(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
mHbOroi_base = 1e6.*squeeze(mean(mean(HbOroi,1),2));
mGCaMProi_base = 100.*squeeze(mean(mean(GCaMProi,1),2));

dataDir = '../MouseData/SS76/RestingState/Day2';
t=2;
load([dataDir,'/','1_v2_rsfc_I.mat'])
[nY,nX,nLam,nT] = size(Iunfilt);
brain_mask = maskSFDI(t).aff_mask + maskSFDI(t).unaff_mask;
brainIdx = find(brain_mask > 0);
yAll = reshape(Iunfilt(:,:,1,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
Chb1 = zeros(nY*nX,nT);
Chb1(brainIdx,:) = ynew';
Chb1 = reshape(Chb1, [nY, nX, nT]);
HbT = Chb1;
rawsigH = reshape(HbT,[nY*nX nT]);
filt_rawsigH = filtfilt(fb,fa,rawsigH')';
filt_rawsigH = reshape(filt_rawsigH,[nY nX nT]);
    
load([dataDir,'/','gcampcorr_v2_rsfc_I.mat'])
[nY,nX,nLam,nT] = size(Iunfilt);
yAll = reshape(Iunfilt(:,:,1,:),[nY*nX,nT]);
y1 = yAll(brainIdx,:);
yMean = mean(y1,1);
a = pinv(yMean*yMean')*yMean*y1';
ynew = y1'-yMean'*a;
Chb1 = zeros(nY*nX,nT);
Chb1(brainIdx,:) = ynew';
Chb1 = reshape(Chb1, [nY, nX, nT]);
deltaGCaMPcorr = Chb1;
rawsigG = reshape(deltaGCaMPcorr,[nY*nX nT]);
filt_rawsigG = filtfilt(fb,fa,rawsigG')';
filt_rawsigG = reshape(filt_rawsigG,[nY nX nT]);
HbOroi = filt_rawsigH(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
GCaMProi = filt_rawsigG(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3),:);
mHbOroi_day2 = 1e6.*squeeze(mean(mean(HbOroi,1),2));
mGCaMProi_day2 = 100.*squeeze(mean(mean(GCaMProi,1),2));

Fs = 5;
T = 1/Fs;
L = size(HbT,3);
timefull = (1:(L))*T;
figure
plot(timefull,mGCaMProi_base,'Color',rgb('dark green'),'LineWidth',2)
hold on
plot(timefull,mHbOroi_base,'Color',rgb('dark red'),'LineWidth',2)
plot(timefull,mGCaMProi_day2,'Color',rgb('green'),'LineWidth',2)
plot(timefull,mHbOroi_day2,'Color',rgb('light red'),'LineWidth',2)
xlabel('Time (sec)')
ylabel('\DeltaF/F (%)   \DeltaHbO (\muM)')
xlim([0 300])
ylim([-4 4])
set(gca,'FontSize',24)
saveas(gcf,'figures/figure5a_bottom.png')

%% Figure 5b
% Combined variance of all animals

load('spect_rest_HbO.mat')
for t = 1:5
    varG_comb = [];
    varH_comb = [];
    for m = 1:12
        temp = varG{m,t};
        temp = reshape(temp,[128*128 1]);
        temp(temp==0)=[];
        temp(isnan(temp))=[];
        varG_comb = cat(1,varG_comb,temp);
        clear temp
        temp = varH{m,t};
        temp = reshape(temp,[128*128 1]);
        temp(temp==0)=[];
        temp(isnan(temp))=[];
        varH_comb = cat(1,varH_comb,temp);
    end
    varG_all{t} = varG_comb;
    varH_all{t} = varH_comb;
end
figure
hold on
histogram(varG_all{1},'FaceColor','k','EdgeColor','k','FaceAlpha',0.5)
histogram(varG_all{2},'FaceColor','b','EdgeColor','b','FaceAlpha',0.5)
histogram(varG_all{5},'FaceColor','c','EdgeColor','c','FaceAlpha',0.5)
xlabel('Variance: GCaMP')
xlim([0 3e-5])
ylim([0 9000])
set(gca,'FontSize',24)
saveas(gcf,'figures/figure5b_left.png')

figure
hold on
histogram(varH_all{1},'FaceColor','k','EdgeColor','k','FaceAlpha',0.5)
histogram(varH_all{2},'FaceColor','b','EdgeColor','b','FaceAlpha',0.5)
histogram(varH_all{5},'FaceColor','c','EdgeColor','c','FaceAlpha',0.5)
legend('Pre-stroke','Day2','Week4')
xlabel('Variance: HbO')
xlim([0 12e-13])
ylim([0 9000])
set(gca,'FontSize',24)
saveas(gcf,'figures/figure5b_right.png')

%% Figure 5c
% Frequency spectrum
load('spect_rest_HbO.mat')
% Remove artifact from GCaMP
for m = 1:12
    for t = 1:5        
        powerG_ipsi(m,t,54) = (powerG_ipsi(m,t,53)+powerG_ipsi(m,t,55))/2;
        powerG_contra(m,t,54) = (powerG_contra(m,t,53)+powerG_contra(m,t,55))/2;
    end
end

color = [0 0 0; 0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
options.alpha = 0.5;
options.line_width = 2; 
options.x_axis = freq;
options.error = 'sem';
for t = 1:5
    figure(1)  
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(powerG_ipsi(:,t,:)),options)
    hold on
    
    figure(2)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(powerHb_ipsi(:,t,:)),options)
    hold on
    
    figure(3)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(powerG_contra(:,t,:)),options)
    hold on
    
    figure(4)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(powerHb_contra(:,t,:)),options)
    hold on
end
figure(1)
xlabel('Frequency (Hz)')
ylabel('Power: GCaMP')
xlim([0.01 0.5])
% ylim([0 0.5])
set(gca,'XScale','log','YScale','log','FontSize',24)
saveas(gcf,'figures/figure5c_1.png')

figure(2)
ylabel('Power: HbO')
xlabel('Frequency (Hz)')
xlim([0.01 0.5])
% ylim([0 0.5])
set(gca,'XScale','log','YScale','log','FontSize',24)
saveas(gcf,'figures/figure5c_2.png')

figure(3)
ylabel('Power: GCaMP')
xlabel('Frequency (Hz)')
xlim([0.01 0.5])
% ylim([0 0.5])
set(gca,'XScale','log','YScale','log','FontSize',24)
saveas(gcf,'figures/figure5c_3.png')

figure(4)
ylabel('Power: HbO')
xlabel('Frequency (Hz)')
xlim([0.01 0.5])
% ylim([0 0.5])
set(gca,'XScale','log','YScale','log','FontSize',24)
saveas(gcf,'figures/figure5c_4.png')

%% Figure 5d
% Power in specific band 0.1-0.3 Hz

load('spect_rest_HbO.mat')
start = 73; 
stop = 121; 
pwrbandG_ipsi = sum(powerG_ipsi(:,:,start:stop),3);
pwrbandH_ipsi = sum(powerHb_ipsi(:,:,start:stop),3);
pwrbandG_contra = sum(powerG_contra(:,:,start:stop),3);
pwrbandH_contra = sum(powerHb_contra(:,:,start:stop),3);

mpwrG_peri = mean(pwrbandG_ipsi,1);
mpwrH_peri = mean(pwrbandH_ipsi,1);
mpwrG_contra = mean(pwrbandG_contra,1);
mpwrH_contra = mean(pwrbandH_contra,1);

mpwr_G = [mpwrG_peri;mpwrG_contra];
mpwr_H = [mpwrH_peri;mpwrH_contra];

spwrG_peri = std(pwrbandG_ipsi,1);
spwrH_peri = std(pwrbandH_ipsi,1);
spwrG_contra = std(pwrbandG_contra,1);
spwrH_contra = std(pwrbandH_contra,1);
spwr_G = [spwrG_peri;spwrG_contra];
spwr_H = [spwrH_peri;spwrH_contra];

a = 0.05;
[hp(1),pG(1)] = ttest2(pwrbandG_ipsi(:,1),pwrbandG_ipsi(:,2),'Alpha',a);
[hp(2),pG(2)] = ttest2(pwrbandG_ipsi(:,1),pwrbandG_ipsi(:,3),'Alpha',a);
[hp(3),pG(3)] = ttest2(pwrbandG_ipsi(:,1),pwrbandG_ipsi(:,4),'Alpha',a);
[hp(4),pG(4)] = ttest2(pwrbandG_ipsi(:,1),pwrbandG_ipsi(:,5),'Alpha',a);
[hc(1),pG(1)] = ttest2(pwrbandG_contra(:,1),pwrbandG_contra(:,2),'Alpha',a);
[hc(2),pG(2)] = ttest2(pwrbandG_contra(:,1),pwrbandG_contra(:,3),'Alpha',a);
[hc(3),pG(3)] = ttest2(pwrbandG_contra(:,1),pwrbandG_contra(:,4),'Alpha',a);
[hc(4),pG(4)] = ttest2(pwrbandG_contra(:,1),pwrbandG_contra(:,5),'Alpha',a);
figure
b = bar(1:2,mpwr_G);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = 2;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mpwr_G(:,e),spwr_G(:,e),spwr_G(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'Ipsilesional','Contralesional'};
set(gca, 'XTickLabel', group, 'FontSize', 24);
ylabel('Power: GCaMP')
set(b, 'FaceAlpha', 1)
ylim([0 4e-3])
sigline1([1.69 1.85],[],3.1e-3)
sigline2([1.69 2],[],3.3e-3)
saveas(gcf,'figures/figure5d_1.png')

a = 0.01;
[hp(1),pG(1)] = ttest2(pwrbandH_ipsi(:,1),pwrbandH_ipsi(:,2),'Alpha',a);
[hp(2),pG(2)] = ttest2(pwrbandH_ipsi(:,1),pwrbandH_ipsi(:,3),'Alpha',a);
[hp(3),pG(3)] = ttest2(pwrbandH_ipsi(:,1),pwrbandH_ipsi(:,4),'Alpha',a);
[hp(4),pG(4)] = ttest2(pwrbandH_ipsi(:,1),pwrbandH_ipsi(:,5),'Alpha',a);
[hc(1),pG(1)] = ttest2(pwrbandH_contra(:,1),pwrbandH_contra(:,2),'Alpha',a);
[hc(2),pG(2)] = ttest2(pwrbandH_contra(:,1),pwrbandH_contra(:,3),'Alpha',a);
[hc(3),pG(3)] = ttest2(pwrbandH_contra(:,1),pwrbandH_contra(:,4),'Alpha',a);
[hc(4),pG(4)] = ttest2(pwrbandH_contra(:,1),pwrbandH_contra(:,5),'Alpha',a);
figure
b = bar(1:2,mpwr_H);
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0 0];
b(3).FaceColor = [1 0 0];
b(4).FaceColor = [1 0.5 0];
b(5).FaceColor = [1 1 0];
hold on
ngroups = 2;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mpwr_H(:,e),spwr_H(:,e),spwr_H(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'Ipsilesional','Contralesional'};
set(gca, 'XTickLabel', group, 'FontSize', 24);
ylabel('Power: HbO')
set(b, 'FaceAlpha', 1)
ylim([0 2.5e-10])
sigline1([0.69 0.85],[],1.9e-10)
sigline1([0.69 1],[],2e-10)
sigline1([1.69 1.85],[],1.6e-10)
sigline1([1.69 2],[],1.7e-10)
saveas(gcf,'figures/figure5d_2.png')

%% Figure 5e
% Spatial power
m=3;
fh = figure();
fh.WindowState = 'maximized';
for t = 1:5
    subplot(2,5,t)
    imagesc(GCaMP_pwr_all{m,t})
    axis image
    colormap bone
    axis off
    caxis([0 0.5e-4])
    colorbar
end
for t = 1:5
    subplot(2,5,t+5)
    imagesc(HbT_pwr_all{m,t})
    axis image
    colormap bone
    axis off
    caxis([0 3e-12])
    colorbar
end
saveas(gcf,'figures/figure5e.png')

%% Figure 5f
% Seperate forepaw power

load('spect_rest_HbO.mat')
start = 2; 
stop = 241; 

% GCaMP
pwrbandG_ipsiforelimb = sum(powerG_ipsiforelimb(:,:,start:stop),3);
pwrbandG_ipsi_nonfore = sum(powerG_ipsi_nonfore(:,:,start:stop),3);
pwrbandG_contraforelimb = sum(powerG_contraforelimb(:,:,start:stop),3);
pwrbandG_contra_nonfore = sum(powerG_contra_nonfore(:,:,start:stop),3);

mpwrG_ipsiforelimb = mean(pwrbandG_ipsiforelimb,1);
mpwrG_ipsi_nonfore = mean(pwrbandG_ipsi_nonfore,1);
mpwrG_contraforelimb = mean(pwrbandG_contraforelimb,1);
mpwrG_contra_nonfore = mean(pwrbandG_contra_nonfore,1);
mpwr_G = [mpwrG_ipsiforelimb;mpwrG_ipsi_nonfore;mpwrG_contraforelimb;mpwrG_contra_nonfore];

spwrG_ipsiforelimb = std(pwrbandG_ipsiforelimb,1);
spwrG_ipsi_nonfore = std(pwrbandG_ipsi_nonfore,1);
spwrG_contraforelimb = std(pwrbandG_contraforelimb,1);
spwrG_contra_nonfore = std(pwrbandG_contra_nonfore,1);
spwr_G = [spwrG_ipsiforelimb;spwrG_ipsi_nonfore;spwrG_contraforelimb;spwrG_contra_nonfore];

a = 0.01;
[hif(1),pG(1)] = ttest2(pwrbandG_ipsiforelimb(:,1),pwrbandG_ipsiforelimb(:,2),'Alpha',a);
[hif(2),pG(2)] = ttest2(pwrbandG_ipsiforelimb(:,1),pwrbandG_ipsiforelimb(:,3),'Alpha',a);
[hif(3),pG(3)] = ttest2(pwrbandG_ipsiforelimb(:,1),pwrbandG_ipsiforelimb(:,4),'Alpha',a);
[hif(4),pG(4)] = ttest2(pwrbandG_ipsiforelimb(:,1),pwrbandG_ipsiforelimb(:,5),'Alpha',a);
[hi(1),pG(1)] = ttest2(pwrbandG_ipsi_nonfore(:,1),pwrbandG_ipsi_nonfore(:,2),'Alpha',a);
[hi(2),pG(2)] = ttest2(pwrbandG_ipsi_nonfore(:,1),pwrbandG_ipsi_nonfore(:,3),'Alpha',a);
[hi(3),pG(3)] = ttest2(pwrbandG_ipsi_nonfore(:,1),pwrbandG_ipsi_nonfore(:,4),'Alpha',a);
[hi(4),pG(4)] = ttest2(pwrbandG_ipsi_nonfore(:,1),pwrbandG_ipsi_nonfore(:,5),'Alpha',a);
[hcf(1),pG(1)] = ttest2(pwrbandG_contraforelimb(:,1),pwrbandG_contraforelimb(:,2),'Alpha',a);
[hcf(2),pG(2)] = ttest2(pwrbandG_contraforelimb(:,1),pwrbandG_contraforelimb(:,3),'Alpha',a);
[hcf(3),pG(3)] = ttest2(pwrbandG_contraforelimb(:,1),pwrbandG_contraforelimb(:,4),'Alpha',a);
[hcf(4),pG(4)] = ttest2(pwrbandG_contraforelimb(:,1),pwrbandG_contraforelimb(:,5),'Alpha',a);
[hc(1),pG(1)] = ttest2(pwrbandG_contra_nonfore(:,1),pwrbandG_contra_nonfore(:,2),'Alpha',a);
[hc(2),pG(2)] = ttest2(pwrbandG_contra_nonfore(:,1),pwrbandG_contra_nonfore(:,3),'Alpha',a);
[hc(3),pG(3)] = ttest2(pwrbandG_contra_nonfore(:,1),pwrbandG_contra_nonfore(:,4),'Alpha',a);
[hc(4),pG(4)] = ttest2(pwrbandG_contra_nonfore(:,1),pwrbandG_contra_nonfore(:,5),'Alpha',a);

fh = figure();
fh.WindowState = 'maximized';
b = bar(1:4,mpwr_G);
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
    er = errorbar(x,mpwr_G(:,e),spwr_G(:,e),spwr_G(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'Ipsi forelimb','Ipsi non-forelimb','Contra forelimb','Contra non-forelimb'};
set(gca, 'XTickLabel', group, 'FontSize', 20);
ylabel('Power: GCaMP')
set(b, 'FaceAlpha', 1)
ylim([0 5e-2])

sigline2([0.69 1.15],[],0.025)
sigline1([0.69 1.3],[],0.027)
sigline1([2.69 2.85],[],0.033)
sigline1([2.69 3],[],0.035)
legend('Pre-stroke','Day2','Week1','Week2','Week4')
saveas(gcf,'figures/figure5f_left.png')


% HbO
pwrbandHb_ipsiforelimb = sum(powerHb_ipsiforelimb(:,:,start:stop),3);
pwrbandHb_ipsi_nonfore = sum(powerHb_ipsi_nonfore(:,:,start:stop),3);
pwrbandHb_contraforelimb = sum(powerHb_contraforelimb(:,:,start:stop),3);
pwrbandHb_contra_nonfore = sum(powerHb_contra_nonfore(:,:,start:stop),3);

mpwrHb_ipsiforelimb = mean(pwrbandHb_ipsiforelimb,1);
mpwrHb_ipsi_nonfore = mean(pwrbandHb_ipsi_nonfore,1);
mpwrHb_contraforelimb = mean(pwrbandHb_contraforelimb,1);
mpwrHb_contra_nonfore = mean(pwrbandHb_contra_nonfore,1);
mpwr_Hb = [mpwrHb_ipsiforelimb;mpwrHb_ipsi_nonfore;mpwrHb_contraforelimb;mpwrHb_contra_nonfore];

spwrHb_ipsiforelimb = std(pwrbandHb_ipsiforelimb,1);
spwrHb_ipsi_nonfore = std(pwrbandHb_ipsi_nonfore,1);
spwrHb_contraforelimb = std(pwrbandHb_contraforelimb,1);
spwrHb_contra_nonfore = std(pwrbandHb_contra_nonfore,1);
spwr_Hb = [spwrHb_ipsiforelimb;spwrHb_ipsi_nonfore;spwrHb_contraforelimb;spwrHb_contra_nonfore];

a = 0.01;
[hif(1),pH(1)] = ttest2(pwrbandHb_ipsiforelimb(:,1),pwrbandHb_ipsiforelimb(:,2),'Alpha',a);
[hif(2),pH(2)] = ttest2(pwrbandHb_ipsiforelimb(:,1),pwrbandHb_ipsiforelimb(:,3),'Alpha',a);
[hif(3),pH(3)] = ttest2(pwrbandHb_ipsiforelimb(:,1),pwrbandHb_ipsiforelimb(:,4),'Alpha',a);
[hif(4),pH(4)] = ttest2(pwrbandHb_ipsiforelimb(:,1),pwrbandHb_ipsiforelimb(:,5),'Alpha',a);
[hi(1),pH(1)] = ttest2(pwrbandHb_ipsi_nonfore(:,1),pwrbandHb_ipsi_nonfore(:,2),'Alpha',a);
[hi(2),pH(2)] = ttest2(pwrbandHb_ipsi_nonfore(:,1),pwrbandHb_ipsi_nonfore(:,3),'Alpha',a);
[hi(3),pH(3)] = ttest2(pwrbandHb_ipsi_nonfore(:,1),pwrbandHb_ipsi_nonfore(:,4),'Alpha',a);
[hi(4),pH(4)] = ttest2(pwrbandHb_ipsi_nonfore(:,1),pwrbandHb_ipsi_nonfore(:,5),'Alpha',a);
[hcf(1),pH(1)] = ttest2(pwrbandHb_contraforelimb(:,1),pwrbandHb_contraforelimb(:,2),'Alpha',a);
[hcf(2),pH(2)] = ttest2(pwrbandHb_contraforelimb(:,1),pwrbandHb_contraforelimb(:,3),'Alpha',a);
[hcf(3),pH(3)] = ttest2(pwrbandHb_contraforelimb(:,1),pwrbandHb_contraforelimb(:,4),'Alpha',a);
[hcf(4),pH(4)] = ttest2(pwrbandHb_contraforelimb(:,1),pwrbandHb_contraforelimb(:,5),'Alpha',a);
[hc(1),pH(1)] = ttest2(pwrbandHb_contra_nonfore(:,1),pwrbandHb_contra_nonfore(:,2),'Alpha',a);
[hc(2),pH(2)] = ttest2(pwrbandHb_contra_nonfore(:,1),pwrbandHb_contra_nonfore(:,3),'Alpha',a);
[hc(3),pH(3)] = ttest2(pwrbandHb_contra_nonfore(:,1),pwrbandHb_contra_nonfore(:,4),'Alpha',a);
[hc(4),pH(4)] = ttest2(pwrbandHb_contra_nonfore(:,1),pwrbandHb_contra_nonfore(:,5),'Alpha',a);

fh = figure();
fh.WindowState = 'maximized';
b = bar(1:4,mpwr_Hb);
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
    er = errorbar(x,mpwr_Hb(:,e),spwr_Hb(:,e),spwr_Hb(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'Ipsi forelimb','Ipsi non-forelimb','Contra forelimb','Contra non-forelimb'};
set(gca, 'XTickLabel', group, 'FontSize', 20);
ylabel('Power: HbO')
set(b, 'FaceAlpha', 1)
ylim([0 1.5e-9])
sigline2([0.69 0.85],[],1.21e-9)
sigline2([0.69 1],[],1.28e-9)
sigline1([1.69 1.85],[],0.8e-9)
sigline1([1.69 2],[],0.85e-9)
sigline1([2.69 2.85],[],0.85e-9)
sigline1([2.69 3],[],0.9e-9)
sigline2([2.69 3.15],[],0.95e-9)
sigline1([3.69 4],[],0.8e-9)
sigline1([3.69 4.15],[],0.85e-9)
saveas(gcf,'figures/figure5f_right.png')


