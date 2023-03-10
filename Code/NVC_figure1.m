%% Stroke NVC figure 1
% Smrithi Sunil
% BOAS Lab

%%
clear

%% figure 1b
clear
data(1) = load('../MouseData/SS75/brainmaskSFDI.mat');

figure
for t = 1:5
    mus = data(1).maskSFDI(t).aff_prop_mus + data(1).maskSFDI(t).unaff_prop_mus;
    subplot(1,5,t)
    imagesc(mus(18:128,80:128))
    colormap jet
    axis image
    axis off
    caxis([0 25])
end
saveas(gcf,'figures/figure1b_top.png')

figure
for t = 1:5
    mus = data(1).maskSFDI(t).aff_prop_mus + data(1).maskSFDI(t).unaff_prop_mus;
    if t == 1 || t == 2 || t == 3
        newMask = imdilate(data(1).maskSFDI(3).stroke_mask, true(30));
        peri = data(1).maskSFDI(t).aff_mask.*abs(newMask - data(1).maskSFDI(3).stroke_mask);
        stroke = data(1).maskSFDI(3).stroke_mask;
        subplot(1,5,t)
        E = mus(18:128,80:128);
        imshow(E) 
        caxis([0 25])
        blue = cat(3, zeros(size(E)), zeros(size(E)), ones(size(E))); 
        hold on
        h = imshow(blue); 
        set(h, 'AlphaData', 0.5*stroke(18:128,80:128)) 
        cyan = cat(3, zeros(size(E)), ones(size(E)), ones(size(E))); 
        h = imshow(cyan); 
        set(h, 'AlphaData', 0.5*peri(18:128,80:128)) 
        hold off
    else
        newMask = imdilate(data(1).maskSFDI(t).stroke_mask, true(30));
        peri = data(1).maskSFDI(t).aff_mask.*abs(newMask - data(1).maskSFDI(t).stroke_mask);
        stroke = data(1).maskSFDI(t).stroke_mask;
        subplot(1,5,t)
        E = mus(18:128,80:128);
        imshow(E) 
        caxis([0 25])
        blue = cat(3, zeros(size(E)), zeros(size(E)), ones(size(E))); 
        hold on
        h = imshow(blue); 
        set(h, 'AlphaData', 0.5*stroke(18:128,80:128)) 
        cyan = cat(3, zeros(size(E)), ones(size(E)), ones(size(E))); 
        h = imshow(cyan); 
        set(h, 'AlphaData', 0.5*peri(18:128,80:128)) 
        hold off
        
    end
end
saveas(gcf,'figures/figure1b_bottom.png')

%% figure 1c

clear
gcamp_uncorr = load('../MouseData/SS82/FunctionalActivation/Baseline/act_GCaMP_ipsi.mat');
gcamp_corr = load('../MouseData/SS82/FunctionalActivation/Baseline/act_GCaMPcorr_ipsi.mat');
hb = load('../MouseData/SS82/FunctionalActivation/Baseline/act_IOSI_ipsi.mat');

wavelength = 470; %input('What wavelength was was imaged first? ');
frameRate = 15; %input('What was the frame rate of the camera? ');
totWave = 3; % total number of wavelengths
numTrials = 20; % total number of trials
trialTime = 30; % length of each trial
baseTime = 5; % length of baseline in seconds
totNumFrames = numTrials*trialTime*frameRate;
numFrame = totNumFrames/totWave;
totalTrialFrames = numFrame/numTrials;
baseFrames = baseTime*frameRate/totWave;
HbOavg = trialAverage(hb.HbO, numTrials, numFrame);
HbRavg = trialAverage(hb.HbR, numTrials, numFrame);
HbTavg = trialAverage(hb.HbT, numTrials, numFrame);
GCaMPuncorr = trialAverage(gcamp_uncorr.deltaGCaMP, numTrials, numFrame);
GCaMPcorr = trialAverage(gcamp_corr.deltaGCaMPcorr, numTrials, numFrame);
time = -(baseTime)+(trialTime/totalTrialFrames):(trialTime/totalTrialFrames):trialTime-baseTime;

rStart = baseFrames;
rEnd = 2*baseFrames;
GCaMPcorrresp = squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
HbOresp = squeeze(mean(HbOavg(:,:,rStart:rEnd),3));
HbRresp = squeeze(mean(HbRavg(:,:,rStart:rEnd),3));
HbTresp = squeeze(mean(HbTavg(:,:,rStart:rEnd),3));
figure
subplot(1,4,1)
imagesc(HbOresp)
colormap jet
axis image
axis off
colorbar
caxis([-6e-6 6e-6])
subplot(1,4,2)
imagesc(HbRresp)
colormap jet
axis image
axis off
colorbar
caxis([-3e-6 3e-6])
subplot(1,4,3)
imagesc(HbTresp)
colormap jet
axis image
axis off
colorbar
caxis([-3e-6 3e-6])
subplot(1,4,4)
imagesc(GCaMPcorrresp)
colormap jet
axis image
axis off
colorbar
caxis([-0.05 0.05])
saveas(gcf,'figures/figure1c_top.png')

roi_GCaMPuncorr = squeeze(mean(mean(GCaMPuncorr(65:80,90:110,:))));
roi_GCaMPcorr = squeeze(mean(mean(GCaMPcorr(65:80,90:110,:))));
roi_HbO = 1e6.*squeeze(mean(mean(HbOavg(65:80,90:110,:))));
roi_HbR = 1e6.*squeeze(mean(mean(HbRavg(65:80,90:110,:))));
roi_HbT = 1e6.*squeeze(mean(mean(HbTavg(65:80,90:110,:))));
figure
hold on
f1 = patch([0 5 5 0],[-8 -8 8 8],[0.8 0.8 0.8],'lineStyle','none');
f2 = plot(time,roi_HbO,'r','Linewidth',2);
f3 = plot(time,roi_HbR,'b','Linewidth',2);
f4 = plot(time,roi_HbT,'g','Linewidth',2);
xlim([-5 25])
ylim([-8 8])
ylabel('GCaMP \DeltaF/F (%)   \DeltaHb (\muM)')
xlabel('Time (sec)')
yyaxis right
f5 = plot(time,100.*roi_GCaMPuncorr,'k--','Linewidth',2);
f6 = plot(time,100.*roi_GCaMPcorr,'k-','Linewidth',2);
xlim([-5 25])
ylim([-8 8])
yticks([])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
legend([f2 f3 f4 f5 f6],'HbO','HbR','HbT','GCaMP(raw)','GCaMP(attn corr)')
set(gca,'FontSize',24)
saveas(gcf,'figures/figure1c_bottom.png')


%% Supplementary figure 1a

data(1) = load('../MouseData/SS75/brainmaskSFDI.mat');
data(2) = load('../MouseData/SS78/brainmaskSFDI.mat');
data(3) = load('../MouseData/SS81/brainmaskSFDI.mat');
data(4) = load('../MouseData/SS82/brainmaskSFDI.mat');
data(5) = load('../MouseData/SS83/brainmaskSFDI.mat');
data(6) = load('../MouseData/SS85/brainmaskSFDI.mat');
% Plot mus at week 1 with stroke boundary
figure
for m = 1:6
    baseline_mus = median(data(m).maskSFDI(1).aff_prop_mus(data(m).maskSFDI(1).aff_mask));
    change_mus = data(m).maskSFDI(3).aff_prop_mus./baseline_mus;
    subplot(1,6,m)
    if m==2 || m==5 || m==6
        img = flip(data(m).maskSFDI(3).aff_prop_mus,2);
        imagesc(img(:,75:128))
        B = bwboundaries(flip(data(m).maskSFDI(3).stroke_mask,2));
        hold on
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2)-75, boundary(:,1), 'k', 'LineWidth', 2)
        end
    else
        imagesc(data(m).maskSFDI(3).aff_prop_mus(:,75:128))
        B = bwboundaries(data(m).maskSFDI(3).stroke_mask);
        hold on
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2)-75, boundary(:,1), 'k', 'LineWidth', 2)
        end
    end
    axis off
    axis image
    colormap jet
    caxis([0 25])
end
saveas(gcf,'figures/suppfig1a.png')

% Supplementary figure 1b
data(1) = load('../MouseData/SS75/brainmaskSFDI.mat');
data(2) = load('../MouseData/SS76/brainmaskSFDI.mat');
data(3) = load('../MouseData/SS77/brainmaskSFDI.mat');
data(4) = load('../MouseData/SS78/brainmaskSFDI.mat');
data(5) = load('../MouseData/SS79/brainmaskSFDI.mat');
data(6) = load('../MouseData/SS80/brainmaskSFDI.mat');
data(7) = load('../MouseData/SS81/brainmaskSFDI.mat');
data(8) = load('../MouseData/SS82/brainmaskSFDI.mat');
data(9) = load('../MouseData/SS83/brainmaskSFDI.mat');
data(10) = load('../MouseData/SS84/brainmaskSFDI.mat');
data(11) = load('../MouseData/SS85/brainmaskSFDI.mat');
data(12) = load('../MouseData/SS93/brainmaskSFDI.mat');

for m = 1:12
    pixels_stroke(m) = length(data(m).maskSFDI(1).aff_prop_mus(data(m).maskSFDI(3).stroke_mask));
    area(m) = 52*52*10^-3*10^-3*pixels_stroke(m);
end

figure
scatter(ones(12,1),area,200,'MarkerFaceColor','r','MarkerEdgeColor','k')
xlim([0 2])
box on
ylabel('Stroke surface area (mm^2)')
set(gca, 'XTickLabel', '', 'FontSize', 24);
saveas(gcf,'figures/suppfig1b.png')

%% Supplementary figure 3

clear
data(1) = load('../MouseData/SS75/brainmaskSFDI.mat');
data(2) = load('../MouseData/SS78/brainmaskSFDI.mat');
data(3) = load('../MouseData/SS81/brainmaskSFDI.mat');
data(4) = load('../MouseData/SS82/brainmaskSFDI.mat');
data(5) = load('../MouseData/SS83/brainmaskSFDI.mat');
data(6) = load('../MouseData/SS85/brainmaskSFDI.mat');
% Plot mu_a at pre stroke
figure
for m = 1:6
    baseline_mua = median(data(m).maskSFDI(1).aff_prop_mua(data(m).maskSFDI(1).aff_mask));
    subplot(1,6,m)
    if m==2 || m==5 || m==6
        img = flip(data(m).maskSFDI(1).aff_prop_mua,2);
        imagesc(img(:,75:128))
    else
        imagesc(data(m).maskSFDI(1).aff_prop_mua(:,75:128))
    end
    axis off
    axis image
    colormap jet
    caxis([0 1])
end
saveas(gcf,'figures/suppfig3a.png')

% Plot mu_s at pre stroke
figure
for m = 1:6
    subplot(1,6,m)
    if m==2 || m==5 || m==6
        img = flip(data(m).maskSFDI(1).aff_prop_mus,2);
        imagesc(img(:,75:128))
    else
        imagesc(data(m).maskSFDI(1).aff_prop_mus(:,75:128))
    end
    axis off
    axis image
    colormap jet
    caxis([0 25])
end
saveas(gcf,'figures/suppfig3b.png')

% Change in mu_a at week 1
figure
for m = 1:6
    baseline_mua = median(data(m).maskSFDI(1).aff_prop_mua(data(m).maskSFDI(1).aff_mask));
    change_mua = data(m).maskSFDI(3).aff_prop_mua./baseline_mua;
    subplot(1,6,m)
    if m==2 || m==5 || m==6
        img = flip(change_mua,2);
        imagesc(img(:,75:128))
    else
        imagesc(change_mua(:,75:128))
    end
    axis off
    axis image
    colormap jet
    caxis([0 2])
end
saveas(gcf,'figures/suppfig3c.png')

% Change in mu_s at week 1
figure
for m = 1:6
    baseline_mus = median(data(m).maskSFDI(1).aff_prop_mus(data(m).maskSFDI(1).aff_mask));
    change_mus = data(m).maskSFDI(3).aff_prop_mus./baseline_mus;
    subplot(1,6,m)
    if m==2 || m==5 || m==6
        img = flip(change_mus,2);
        imagesc(img(:,75:128))
    else
        imagesc(change_mus(:,75:128))
    end
    axis off
    axis image
    colormap jet
    caxis([0 2.5])
end
saveas(gcf,'figures/suppfig3d.png')

