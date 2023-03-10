%% Figure 6
clear
load('RSFC_analysis_highHz.mat')
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

%% Spatial maps
m=3;
fh = figure();
fh.WindowState = 'maximized';
count = 0;

for t = [1 2 5]
    count = count+1;
    subplot(4,3,count)
    imagesc(squeeze(GCaMP_forelimb_conn{m,t}))
    hold on
    plot(stats_ipsi(m).Centroid(1),stats_ipsi(m).Centroid(2),'ko','MarkerFaceColor','k','MarkerSize',10)
    axis image
    colormap jet
    caxis([-1 1])
    axis off
    
    subplot(4,3,count+3)
    imagesc(squeeze(LRC_GCaMP{m,t}))
    axis image
    colormap jet
    caxis([-1 1])
    axis off
    
    subplot(4,3,count+6)
    imagesc(squeeze(BOC_GCaMP{m,t}))
    axis image
    colormap jet
    caxis([-1 1])
    axis off
    
    subplot(4,3,count+9)
    imagesc(squeeze(GCaMP_forelimb_conn_contra{m,t}))
    hold on
    plot(stats_contra(m).Centroid(1),stats_contra(m).Centroid(2),'ko','MarkerFaceColor','k','MarkerSize',10)
    axis image
    colormap jet
    caxis([-1 1])
    axis off
end
saveas(gcf,'figures/figure6_highHz_right.png')

%% Plotting for forelimb connectivity 
thresh = 0:0.1:1;
for m = 1:12
    mask = [mouse{m},'brainmaskSFDI.mat'];
    load(mask)
    for t = 1:5        
        unaff_mask = maskSFDI(t).unaff_mask;
        for th = 1:length(thresh)
            tempG = unaff_mask.*GCaMP_forelimb_conn{m,t};
            tempG(tempG<=thresh(th)) = NaN;
            areaG(m,t,th) = sum(sum(~isnan(tempG)))./sum(unaff_mask(:));

            tempH = unaff_mask.*HbO_forelimb_conn{m,t};
            tempH(tempH<=thresh(th)) = NaN;
            areaH(m,t,th) = sum(sum(~isnan(tempH)))./sum(unaff_mask(:));
            
            tempG(~isnan(tempG)) = 1;
            tempG(isnan(tempG)) = 0;           
            tempH(~isnan(tempH)) = 1;
            tempH(isnan(tempH)) = 0;    
            if isempty(dice(tempG,tempH)) == 1
                overlap_foreconn(m,t,th) = 0;
            else
                overlap_foreconn(m,t,th) = dice(tempG,tempH);
            end
        end
    end
end

areaG(isnan(areaG)) = 0;
areaH(isnan(areaH)) = 0;
fh = figure();
fh.WindowState = 'maximized';
color = [0 0 0; 0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
options.alpha = 0.5;
options.line_width = 2; 
options.x_axis = thresh;
options.error = 'sem';
for t = 1:5
    subplot(1,3,1)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(areaG(:,t,:)),options)
    hold on
    
    subplot(1,3,2)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(areaH(:,t,:)),options)
    hold on
    
    subplot(1,3,3)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(overlap_foreconn(:,t,:)),options)
    hold on
end
subplot(1,3,1)
xlim([0 1])
ylim([0 0.5])
ylabel('Area (normalized)')
xlabel('Threshold')
set(gca, 'FontSize', 24);
subplot(1,3,2)
xlim([0 1])
ylim([0 0.5])
ylabel('Area (normalized)')
xlabel('Threshold')
set(gca, 'FontSize', 24);
subplot(1,3,3)
ylabel('Dice similarity')
xlim([0 1])
ylim([0 1])
ylabel('Area (normalized)')
xlabel('Threshold')
set(gca, 'FontSize', 24);
saveas(gcf,'figures/figure6_highHz_1.png')

%% Supplementary figure 10a
% At 0.4 threshold
th = 5;
areaG = squeeze(areaG(:,:,th));
areaH = squeeze(areaH(:,:,th));
sim = squeeze(overlap_foreconn(:,:,th));
a = 0.05;
[hG(1),pG(1)] = ttest2(areaG(:,1),areaG(:,2),'Alpha',a);
[hG(2),pG(2)] = ttest2(areaG(:,1),areaG(:,3),'Alpha',a);
[hG(3),pG(3)] = ttest2(areaG(:,1),areaG(:,4),'Alpha',a);
[hG(4),pG(4)] = ttest2(areaG(:,1),areaG(:,5),'Alpha',a);
[hH(1),pG(1)] = ttest2(areaH(:,1),areaH(:,2),'Alpha',a);
[hH(2),pG(2)] = ttest2(areaH(:,1),areaH(:,3),'Alpha',a);
[hH(3),pG(3)] = ttest2(areaH(:,1),areaH(:,4),'Alpha',a);
[hH(4),pG(4)] = ttest2(areaH(:,1),areaH(:,5),'Alpha',a);
[hS(1),pG(1)] = ttest2(sim(:,1),sim(:,2),'Alpha',a);
[hS(2),pG(2)] = ttest2(sim(:,1),sim(:,3),'Alpha',a);
[hS(3),pG(3)] = ttest2(sim(:,1),sim(:,4),'Alpha',a);
[hS(4),pG(4)] = ttest2(sim(:,1),sim(:,5),'Alpha',a);

mareaG = mean(areaG,1);
sareaG = std(areaG,1)./sqrt(12);
mareaH = mean(areaH,1);
sareaH = std(areaH,1)./sqrt(12);
msim = mean(sim,1);
ssim = std(sim,1)./sqrt(12);
marea = [mareaG;mareaH;msim];
sarea = [sareaG;sareaH;ssim];
figure
b = bar(1:3,marea);
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
    er = errorbar(x,marea(:,e),sarea(:,e),sarea(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'GCaMP (Area)','HbO (Area)','Dice Similarity'};
set(gca, 'XTickLabel', group, 'FontSize', 24);
set(b, 'FaceAlpha', 1)
ylim([0 0.5])
sigline1([1.69 1.85],[],0.24)
sigline1([2.69 2.85],[],0.61)
sigline2([2.69 3],[],0.65)
saveas(gcf,'figures/suppfig10a_high.png')

%% Interhemispheric connectivity

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
fh = figure();
fh.WindowState = 'maximized';
color = [0 0 0; 0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
options.alpha = 0.5;
options.line_width = 2; 
options.x_axis = thresh;
options.error = 'sem';
for t = 1:5
    subplot(1,3,1)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(areaG(:,t,:)),options)
    hold on
    
    subplot(1,3,2)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(areaH(:,t,:)),options)
    hold on
    
    subplot(1,3,3)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(overlap_LRC(:,t,:)),options)
    hold on
end
subplot(1,3,1)
ylabel('Area (normalized)')
xlabel('Threshold')
xlim([0 1])
ylim([0 1])
set(gca, 'FontSize', 20);
subplot(1,3,2)
xlim([0 1])
ylim([0 1])
ylabel('Area (normalized)')
xlabel('Threshold')
set(gca, 'FontSize', 20);
subplot(1,3,3)
ylabel('Dice similarity')
xlim([0 1])
ylim([0 1])
xlabel('Threshold')
set(gca, 'FontSize', 20);
saveas(gcf,'figures/figure6_highHz_2.png')

%% Supplementary figure 10b
% At 0.4 threshold
th = 5;
areaG = squeeze(areaG(:,:,th));
areaH = squeeze(areaH(:,:,th));
sim = squeeze(overlap_LRC(:,:,th));
a = 0.01;
[hG(1),pG(1)] = ttest2(areaG(:,1),areaG(:,2),'Alpha',a);
[hG(2),pG(2)] = ttest2(areaG(:,1),areaG(:,3),'Alpha',a);
[hG(3),pG(3)] = ttest2(areaG(:,1),areaG(:,4),'Alpha',a);
[hG(4),pG(4)] = ttest2(areaG(:,1),areaG(:,5),'Alpha',a);
[hH(1),pG(1)] = ttest2(areaH(:,1),areaH(:,2),'Alpha',a);
[hH(2),pG(2)] = ttest2(areaH(:,1),areaH(:,3),'Alpha',a);
[hH(3),pG(3)] = ttest2(areaH(:,1),areaH(:,4),'Alpha',a);
[hH(4),pG(4)] = ttest2(areaH(:,1),areaH(:,5),'Alpha',a);
[hS(1),pG(1)] = ttest2(sim(:,1),sim(:,2),'Alpha',a);
[hS(2),pG(2)] = ttest2(sim(:,1),sim(:,3),'Alpha',a);
[hS(3),pG(3)] = ttest2(sim(:,1),sim(:,4),'Alpha',a);
[hS(4),pG(4)] = ttest2(sim(:,1),sim(:,5),'Alpha',a);

mareaG = mean(areaG,1);
sareaG = std(areaG,1)./sqrt(12);
mareaH = mean(areaH,1);
sareaH = std(areaH,1)./sqrt(12);
msim = mean(sim,1);
ssim = std(sim,1)./sqrt(12);
marea = [mareaG;mareaH;msim];
sarea = [sareaG;sareaH;ssim];
figure(8)
b = bar(1:3,marea);
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
    er = errorbar(x,marea(:,e),sarea(:,e),sarea(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'GCaMP (Area)','HbO (Area)','Dice Similarity'};
set(gca, 'XTickLabel', group, 'FontSize', 24);
set(b, 'FaceAlpha', 1)
ylim([0 1])
sigline1([0.69 0.85],[],0.44)
sigline1([1.69 1.85],[],0.59)
sigline1([1.69 2],[],0.63)
sigline1([2.69 2.85],[],0.81)
sigline2([2.69 3],[],0.85)
saveas(gcf,'figures/suppfig10b_high.png')

%% Global connectivity
% Area vs correlation coefficient

thresh = 0:0.1:1;
for m = 1:12
    mask = [mouse{m},'brainmaskSFDI.mat'];
    load(mask)
    for t = 1:5        
        mask1 = maskSFDI(t).aff_mask + maskSFDI(t).unaff_mask;
        for th = 1:length(thresh)
            tempG = mask1.*BOC_GCaMP{m,t};
            tempG(tempG<=thresh(th)) = NaN;
            areaG(m,t,th) = sum(sum(~isnan(tempG)))./sum(mask1(:));

            tempH = mask1.*BOC_HbO{m,t};
            tempH(tempH<=thresh(th)) = NaN;
            areaH(m,t,th) = sum(sum(~isnan(tempH)))./sum(mask1(:));
            
            tempG(~isnan(tempG)) = 1;
            tempG(isnan(tempG)) = 0;           
            tempH(~isnan(tempH)) = 1;
            tempH(isnan(tempH)) = 0;            
            if isempty(dice(tempG,tempH)) == 1
                overlap_BOC(m,t,th) = 0;
            else
                overlap_BOC(m,t,th) = dice(tempG,tempH);
            end
        end
    end
end

areaG(isnan(areaG)) = 0;
areaH(isnan(areaH)) = 0;
fh = figure();
fh.WindowState = 'maximized';
color = [0 0 0; 0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
options.alpha = 0.5;
options.line_width = 2; 
options.x_axis = thresh;
options.error = 'sem';
for t = 1:5
    subplot(1,3,1)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(areaG(:,t,:)),options)
    hold on
    
    subplot(1,3,2)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(areaH(:,t,:)),options)
    hold on
    
    subplot(1,3,3)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(overlap_BOC(:,t,:)),options)
    hold on
end
subplot(1,3,1)
ylabel('Area (normalized)')
xlabel('Threshold')
xlim([0 1])
ylim([0 1])
set(gca, 'FontSize', 20);
subplot(1,3,2)
xlim([0 1])
ylim([0 1])
ylabel('Area (normalized)')
xlabel('Threshold')
set(gca, 'FontSize', 20);
subplot(1,3,3)
ylabel('Dice similarity')
xlabel('Threshold')
xlim([0 1])
ylim([0 1])
set(gca, 'FontSize', 20);
saveas(gcf,'figures/figure6_highHz_3.png')

%% Supplementary figure 10c
% At 0.4 threshold
th = 5;
areaG = squeeze(areaG(:,:,th));
areaH = squeeze(areaH(:,:,th));
sim = squeeze(overlap_BOC(:,:,th));
a = 0.05;
[hG(1),pG(1)] = ttest2(areaG(:,1),areaG(:,2),'Alpha',a);
[hG(2),pG(2)] = ttest2(areaG(:,1),areaG(:,3),'Alpha',a);
[hG(3),pG(3)] = ttest2(areaG(:,1),areaG(:,4),'Alpha',a);
[hG(4),pG(4)] = ttest2(areaG(:,1),areaG(:,5),'Alpha',a);
[hH(1),pG(1)] = ttest2(areaH(:,1),areaH(:,2),'Alpha',a);
[hH(2),pG(2)] = ttest2(areaH(:,1),areaH(:,3),'Alpha',a);
[hH(3),pG(3)] = ttest2(areaH(:,1),areaH(:,4),'Alpha',a);
[hH(4),pG(4)] = ttest2(areaH(:,1),areaH(:,5),'Alpha',a);
[hS(1),pG(1)] = ttest2(sim(:,1),sim(:,2),'Alpha',a);
[hS(2),pG(2)] = ttest2(sim(:,1),sim(:,3),'Alpha',a);
[hS(3),pG(3)] = ttest2(sim(:,1),sim(:,4),'Alpha',a);
[hS(4),pG(4)] = ttest2(sim(:,1),sim(:,5),'Alpha',a);

mareaG = mean(areaG,1);
sareaG = std(areaG,1)./sqrt(12);
mareaH = mean(areaH,1);
sareaH = std(areaH,1)./sqrt(12);
msim = mean(sim,1);
ssim = std(sim,1)./sqrt(12);
marea = [mareaG;mareaH;msim];
sarea = [sareaG;sareaH;ssim];
figure(12)
b = bar(1:3,marea);
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
    er = errorbar(x,marea(:,e),sarea(:,e),sarea(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'GCaMP (Area)','HbO (Area)','Dice Similarity'};
set(gca, 'XTickLabel', group, 'FontSize', 24);
set(b, 'FaceAlpha', 1)
ylim([0 1])
sigline2([0.69 0.85],[],0.63)
sigline2([0.69 1],[],0.67)
sigline2([1.69 2],[],0.63)
saveas(gcf,'figures/suppfig10c_high.png')

%% Plotting for contralateral forelimb connectivity 

thresh = 0:0.1:1;
for m = 1:12
    mask = [mouse{m},'brainmaskSFDI.mat'];
    load(mask)
    for t = 1:5        
        unaff_mask = maskSFDI(t).unaff_mask;
        for th = 1:length(thresh)
            tempG = unaff_mask.*GCaMP_forelimb_conn_contra{m,t};
            tempG(tempG<=thresh(th)) = NaN;
            areaG(m,t,th) = sum(sum(~isnan(tempG)))./sum(unaff_mask(:));

            tempH = unaff_mask.*HbO_forelimb_conn_contra{m,t};
            tempH(tempH<=thresh(th)) = NaN;
            areaH(m,t,th) = sum(sum(~isnan(tempH)))./sum(unaff_mask(:));
            
            tempG(~isnan(tempG)) = 1;
            tempG(isnan(tempG)) = 0;           
            tempH(~isnan(tempH)) = 1;
            tempH(isnan(tempH)) = 0;    
            if isempty(dice(tempG,tempH)) == 1
                overlap_foreconn(m,t,th) = 0;
            else
                overlap_foreconn(m,t,th) = dice(tempG,tempH);
            end
        end
    end
end

areaG(isnan(areaG)) = 0;
areaH(isnan(areaH)) = 0;
fh = figure();
fh.WindowState = 'maximized';
color = [0 0 0; 0.5 0 0; 1 0 0; 1 0.5 0; 1 1 0];
options.alpha = 0.5;
options.line_width = 2; 
options.x_axis = thresh;
options.error = 'sem';
for t = 1:5
    subplot(1,3,1)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(areaG(:,t,:)),options)
    hold on
    
    subplot(1,3,2)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(areaH(:,t,:)),options)
    hold on
    
    subplot(1,3,3)
    fig = gcf;
    options.handle = figure(fig);
    options.color_area = color(t,:); 
    options.color_line = color(t,:); 
    plot_areaerrorbar(squeeze(overlap_foreconn(:,t,:)),options)
    hold on
end
subplot(1,3,1)
xlabel('Threshold')
ylabel('Area (normalized)')
xlim([0 1])
ylim([0 1])
set(gca, 'FontSize', 20);
subplot(1,3,2)
ylabel('Area (normalized)')
xlabel('Threshold')
xlim([0 1])
ylim([0 1])
set(gca, 'FontSize', 20);
subplot(1,3,3)
xlabel('Threshold')
ylabel('Dice similarity')
xlim([0 1])
ylim([0 1])
set(gca, 'FontSize', 20);
saveas(gcf,'figures/figure6_highHz_4.png')

%% Supplementary figure 10d
% At 0.4 threshold
th = 5;
areaG = squeeze(areaG(:,:,th));
areaH = squeeze(areaH(:,:,th));
sim = squeeze(overlap_foreconn(:,:,th));
a = 0.01;
[hG(1),pG(1)] = ttest2(areaG(:,1),areaG(:,2),'Alpha',a);
[hG(2),pG(2)] = ttest2(areaG(:,1),areaG(:,3),'Alpha',a);
[hG(3),pG(3)] = ttest2(areaG(:,1),areaG(:,4),'Alpha',a);
[hG(4),pG(4)] = ttest2(areaG(:,1),areaG(:,5),'Alpha',a);
[hH(1),pG(1)] = ttest2(areaH(:,1),areaH(:,2),'Alpha',a);
[hH(2),pG(2)] = ttest2(areaH(:,1),areaH(:,3),'Alpha',a);
[hH(3),pG(3)] = ttest2(areaH(:,1),areaH(:,4),'Alpha',a);
[hH(4),pG(4)] = ttest2(areaH(:,1),areaH(:,5),'Alpha',a);
[hS(1),pG(1)] = ttest2(sim(:,1),sim(:,2),'Alpha',a);
[hS(2),pG(2)] = ttest2(sim(:,1),sim(:,3),'Alpha',a);
[hS(3),pG(3)] = ttest2(sim(:,1),sim(:,4),'Alpha',a);
[hS(4),pG(4)] = ttest2(sim(:,1),sim(:,5),'Alpha',a);

mareaG = mean(areaG,1);
sareaG = std(areaG,1)./sqrt(12);
mareaH = mean(areaH,1);
sareaH = std(areaH,1)./sqrt(12);
msim = mean(sim,1);
ssim = std(sim,1)./sqrt(12);
marea = [mareaG;mareaH;msim];
sarea = [sareaG;sareaH;ssim];
figure
b = bar(1:3,marea);
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
    er = errorbar(x,marea(:,e),sarea(:,e),sarea(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
group = {'GCaMP (Area)','HbO (Area)','Dice Similarity'};
set(gca, 'XTickLabel', group, 'FontSize', 16);
set(b, 'FaceAlpha', 1)
ylim([0 1])
sigline1([0.69 0.85],[],0.59)
sigline2([1.69 1.85],[],0.45)
saveas(gcf,'figures/suppfig10d_high.png')
