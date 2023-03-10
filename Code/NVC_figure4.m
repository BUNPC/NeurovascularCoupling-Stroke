%% Figure 4
clear

%% Figure 4a,b,c

% Pre-stroke
% IRF within forelimb area
animal_number = 'SS77';
t=1;
[timepoints, mask] = animals(animal_number);
load(mask)
load([timepoints{t},'/act_HRF_ipsi.mat'])
HRF_roi = 1e5.*HRF.*maskSFDI(1).ipsiOutline;
HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
temp = HRF_roi;
temp(isnan(temp)) = 0;
temp = mean(temp,2);
idx = find(temp~=0);
HRF_newroi = HRF_roi(idx,:);
figure
time = -4+0.2:0.2:10;
fig = gcf;
options.x_axis = time;
options.handle = fig;
options.error = 'std';
options.color_area = 'k'; 
options.color_line = 'k'; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(HRF_newroi(:,9:end),options)
hold on
plot([0 0],[-3 3],'k--')
ylabel('Amplitude (arbitrary units)')
ylim([-3 3])
xlim([-3 8])
set(gca,'FontSize',24,'XTick','')
saveas(gcf,'figures/figure4a_top.png')

% Time course of GCaMP, measured HbT, estimated HbT
animal_number = 'SS77';
t=1;
[timepoints, mask] = animals(animal_number);
load(mask)
load([timepoints{t},'/','act_IOSI_ipsi.mat'])
load([timepoints{t},'/','act_GCaMPcorr_ipsi.mat'])
load([timepoints{t},'/','act_HRF_ipsi.mat'])
HRF_roi = 1e5.*HRF.*maskSFDI(1).ipsiOutline;
HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
temp = HRF_roi;
temp(isnan(temp)) = 0;
temp = mean(temp,2);
idx = find(temp~=0);
HRF_roi = HRF_roi(idx,:);
gcamp_roi = deltaGCaMPcorr.*maskSFDI(1).ipsiOutline;
gcamp_roi = reshape(gcamp_roi,[size(gcamp_roi,1)*size(gcamp_roi,2) size(gcamp_roi,3)]);
gcamp_roi = 100.*gcamp_roi(idx,:);
tHbT_roi = HbT.*maskSFDI(1).ipsiOutline;
tHbT_roi = reshape(tHbT_roi,[size(tHbT_roi,1)*size(tHbT_roi,2) size(tHbT_roi,3)]);
tHbT_roi = 1e6.*tHbT_roi(idx,:);
pHbT_roi = pHbT.*maskSFDI(1).ipsiOutline;
pHbT_roi = reshape(pHbT_roi,[size(pHbT_roi,1)*size(pHbT_roi,2) size(pHbT_roi,3)]);
pHbT_roi = 1e6.*pHbT_roi(idx,:);
timefull = 0.2:0.2:600;
start = 26;
stop = size(HbT,3);
figure
plot(timefull, squeeze(mean(gcamp_roi)),'Color',rgb('green'),'LineWidth',2)
hold on
plot(timefull, squeeze(mean(tHbT_roi)),'Color',rgb('red'),'LineWidth',2)
plot(timefull(start:stop), squeeze(mean(pHbT_roi)),'Color',rgb('black'),'LineWidth',2)
ylabel('\DeltaF/F (%)   \DeltaHb (\muM)')
set(gca,'FontSize',24,'XTick','')
x=60;
xlim([x x+120])
ylim([-5 10])
saveas(gcf,'figures/figure4b_top.png')

% Correlation coefficient
figure
imagesc(correlation)
axis image
axis off
colormap turbo
colorbar
caxis([0 1])
set(gca,'FontSize',24)
saveas(gcf,'figures/figure4c_top.png')

%% Figure 4a,b,c
% Day2

% IRF within stroke area
animal_number = 'SS78';
[timepoints, mask] = animals(animal_number);
load(mask)
t=2;
load([timepoints{t},'/act_HRF_ipsi.mat'])
stroke = maskSFDI(3).stroke_mask;
HRF_roi = 1e5.*HRF.*stroke;
HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
temp = HRF_roi;
temp(isnan(temp)) = 0;
temp = mean(temp,2);
idx = find(temp~=0);
HRF_newroi = HRF_roi(idx,:);
figure
time = -4+0.2:0.2:10;
fig = gcf;
options.x_axis = time;
options.handle = fig;
options.error = 'std';
options.color_area = 'k'; 
options.color_line = 'k'; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(HRF_newroi(:,9:end),options)
hold on
plot([0 0],[-3 3],'k--')
ylabel('Amplitude (arbitrary units)')
ylim([-3 3])
xlim([-3 8])
set(gca,'FontSize',24)
saveas(gcf,'figures/figure4a_bottom.png')

% Time course of GCaMP, measured HbT, estimated HbT
load([timepoints{t},'/','act_IOSI_ipsi.mat'])
load([timepoints{t},'/','act_GCaMPcorr_ipsi.mat'])
HRF_roi = 1e5.*HRF.*maskSFDI(1).ipsiOutline;
HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
temp = HRF_roi;
temp(isnan(temp)) = 0;
temp = mean(temp,2);
idx = find(temp~=0);
HRF_roi = HRF_roi(idx,:);
gcamp_roi = deltaGCaMPcorr.*maskSFDI(1).ipsiOutline;
gcamp_roi = reshape(gcamp_roi,[size(gcamp_roi,1)*size(gcamp_roi,2) size(gcamp_roi,3)]);
gcamp_roi = 100.*gcamp_roi(idx,:);
tHbT_roi = HbT.*maskSFDI(1).ipsiOutline;
tHbT_roi = reshape(tHbT_roi,[size(tHbT_roi,1)*size(tHbT_roi,2) size(tHbT_roi,3)]);
tHbT_roi = 1e6.*tHbT_roi(idx,:);
pHbT_roi = pHbT.*maskSFDI(1).ipsiOutline;
pHbT_roi = reshape(pHbT_roi,[size(pHbT_roi,1)*size(pHbT_roi,2) size(pHbT_roi,3)]);
pHbT_roi = 1e6.*pHbT_roi(idx,:);
timefull = 0.2:0.2:600;
start = 26;
stop = size(HbT,3);
figure
plot(timefull, squeeze(mean(gcamp_roi)),'Color',rgb('green'),'LineWidth',2)
hold on
plot(timefull, squeeze(mean(tHbT_roi)),'Color',rgb('red'),'LineWidth',2)
plot(timefull(start:stop), squeeze(mean(pHbT_roi)),'Color',rgb('black'),'LineWidth',2)
ylabel('\DeltaF/F (%)   \DeltaHb (\muM)')
legend('Mesured GCaMP','Measured HbT','Predicted HbT')
xlabel('Time (sec)')
set(gca,'FontSize',24)
x=30;
xlim([x x+120])
ylim([-5 10])
saveas(gcf,'figures/figure4b_bottom.png')

% Correlation coefficient
figure
imagesc(correlation)
axis image
axis off
colormap turbo
colorbar
caxis([0 1])
set(gca,'FontSize',24)
saveas(gcf,'figures/figure4c_bottom.png')

%% Figure 4e
% One example animal in time

animal_number = 'SS78';
[timepoints, mask] = animals(animal_number);
time = -4+0.2:0.2:10;
load(mask)
color = [rgb('red');rgb('blue');rgb('green')];
options.alpha = 0.3;
options.line_width = 2;
options.x_axis = time;
options.error = 'std';
fh = figure();
fh.WindowState = 'maximized';
for t = 1:5
    for h = 1:3
        if h == 1
            load([timepoints{t},'/act_HRFHbO_ipsi.mat'])
        elseif h == 2
            load([timepoints{t},'/act_HRFHbR_ipsi.mat'])
        elseif h == 3
            load([timepoints{t},'/act_HRF_ipsi.mat'])
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
    contra = maskSFDI(1).contraOutline;

    HRF_roi = 1e5.*HRF.*stroke;
    HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
    temp = HRF_roi;
    temp(isnan(temp)) = 0;
    temp = mean(temp,2);
    idx = find(temp~=0);
    HRF_newroi = HRF_roi(idx,:);
    subplot(3,5,t)
    hold on
    fig = gcf;
    options.handle = fig;
    options.color_area = color(h,:); 
    options.color_line = color(h,:); 
    plot_areaerrorbar(HRF_newroi(:,9:end),options)
    hold on
    plot([0 0],[-4 6],'k--')
    axis square
    box on
    ylim([-4 6])
    xlim([-2 5])
    xlabel('Time (sec)')
    set(gca,'FontSize',16)

    HRF_roi = 1e5.*HRF.*peri;
    HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
    temp = HRF_roi;
    temp(isnan(temp)) = 0;
    temp = mean(temp,2);
    idx = find(temp~=0);
    HRF_newroi = HRF_roi(idx,:);
    subplot(3,5,t+5)
    hold on
    fig = gcf;
    options.handle = fig;
    options.color_area = color(h,:); 
    options.color_line = color(h,:); 
    plot_areaerrorbar(HRF_newroi(:,9:end),options)
    hold on
    plot([0 0],[-4 6],'k--')
    axis square
    box on
    ylim([-4 6])
    xlim([-2 5])
    xlabel('Time (sec)')
    set(gca,'FontSize',16)

    HRF_roi = 1e5.*HRF.*contra;
    HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);
    temp = HRF_roi;
    temp(isnan(temp)) = 0;
    temp = mean(temp,2);
    idx = find(temp~=0);
    HRF_newroi = HRF_roi(idx,:);
    subplot(3,5,t+10)
    hold on
    fig = gcf;
    options.handle = fig;
    options.color_area = color(h,:); 
    options.color_line = color(h,:); 
    plot_areaerrorbar(HRF_newroi(:,9:end),options)
    hold on
    plot([0 0],[-4 6],'k--')
    axis square
    box on
    ylim([-4 6])
    xlim([-2 5])
    xlabel('Time (sec)')
    set(gca,'FontSize',16)
    end
end
saveas(gcf,'figures/figure4e.png')

%% Figure 4f
% All animals in time

color = [rgb('red');rgb('blue');rgb('green');rgb('light red');rgb('light blue');rgb('light green')];
time = -4+0.2:0.2:10;
fh = figure();
fh.WindowState = 'maximized';
for roi = 1:3
    for t = 1:5
        if roi == 1
            subplot(3,5,t)
        elseif roi ==2
            subplot(3,5,t+5)
        else
            subplot(3,5,t+10)
        end
        hold on
        for h = 1:3
            if h == 1
                load('NVC_HbO_stats.mat')
            elseif h == 2
                load('NVC_HbR_stats.mat')
            elseif h == 3
                load('NVC_HbT_stats.mat')
            end
            for m = 1:12
                HRFm(m,t,:) = 1e5.*mean(HRF_HbT{roi,t,m},1);        
            end        
            mHRF(t,:) = mean(HRFm(:,t,:),1);   
            for m = 1:12
                plot(time,squeeze(HRFm(m,t,9:end)),'Color',color(3+h,:),'LineWidth',1)
            end
            plot(time,mHRF(t,9:end),'Color',color(h,:),'LineWidth',3)
            hold on
            plot([0 0],[-4 6],'k--')
            box on
            axis square
            ylim([-4 6])
            xlim([-2 5])
            xlabel('Time (sec)')
            set(gca,'FontSize',16)
        end
    end
end
saveas(gcf,'figures/figure4f.png')

%% Figure 4g
% Spatial maps of example animal

m = 7;
load('NVC_HbT_stats.mat')
fh = figure();
fh.WindowState = 'maximized';
count = 0;
for t = 1:5
    if t==1 || t==2 || t==3 || t==5
    count = count+1;
    subplot(3,4,count)
    imshow(Corr_baseirf{m,t})
    axis image
    colormap turbo
    axis off
    caxis([-0.5 1])
    end
end

load('NVC_HbO_stats.mat')
for t = 1:5
    if t==1 || t==2 || t==3 || t==5
    count = count+1;
    subplot(3,4,count)
    imshow(Corr_baseirf{m,t})
    axis image
    colormap turbo
    axis off
    caxis([0 1])
    end
end

load('NVC_HbR_stats.mat')
for t = 1:5
    if t==1 || t==2 || t==3 || t==5
    count = count+1;
    subplot(3,4,count)
    imshow(Corr_baseirf{m,t})
    axis image
    colormap turbo
    axis off
    caxis([0 1])
    end
end
saveas(gcf,'figures/figure4g.png')

%% NVC stats

load('NVC_HbT_stats.mat')
grp = 3;
mcorr_roi = mean(corr_baseirf_roi,3);
stdcorr_roi = std(corr_baseirf_roi,[],3); %./sqrt(size(corr_roi,3));
a = 0.05;
[h_fore(1),p(1)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,2,:),'Alpha',a);
[h_fore(2),p(2)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,3,:),'Alpha',a);
[h_fore(3),p(3)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,4,:),'Alpha',a);
[h_fore(4),p(4)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,5,:),'Alpha',a);
[h_peri(1),p(1)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,2,:),'Alpha',a);
[h_peri(2),p(2)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,3,:),'Alpha',a);
[h_peri(3),p(3)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,4,:),'Alpha',a);
[h_peri(4),p(4)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,5,:),'Alpha',a);
[h_contra(1),p(1)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,2,:),'Alpha',a);
[h_contra(2),p(2)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,3,:),'Alpha',a);
[h_contra(3),p(3)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,4,:),'Alpha',a);
[h_contra(4),p(4)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,5,:),'Alpha',a);
figure
roi = {'Stroke','Peri-infarct','Contralesional'};
x=1:grp;
b = bar(x,mcorr_roi(1:grp,:));
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
set(gca, 'XTick', '');
ylim([-0.5 2])
set(gca,'FontSize',24)
sigline1([0.69 0.85],[],0.83)
sigline1([0.69 1],[],0.89)
sigline1([0.69 1.15],[],0.95)
sigline1([0.69 1.3],[],1.01)
sigline1([1.69 1.85],[],0.83)
sigline1([1.69 2],[],0.89)
sigline2([1.69 2.15],[],0.95)
legend('Pre-Stroke','Day2','Week1','Week2','Week4')
saveas(gcf,'figures/figure4h_top.png')

load('NVC_HbO_stats.mat')
grp = 3;
mcorr_roi = mean(corr_baseirf_roi,3);
stdcorr_roi = std(corr_baseirf_roi,[],3); %./sqrt(size(corr_roi,3));
a = 0.05;
[h_fore(1),p(1)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,2,:),'Alpha',a);
[h_fore(2),p(2)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,3,:),'Alpha',a);
[h_fore(3),p(3)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,4,:),'Alpha',a);
[h_fore(4),p(4)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,5,:),'Alpha',a);
[h_peri(1),p(1)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,2,:),'Alpha',a);
[h_peri(2),p(2)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,3,:),'Alpha',a);
[h_peri(3),p(3)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,4,:),'Alpha',a);
[h_peri(4),p(4)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,5,:),'Alpha',a);
[h_contra(1),p(1)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,2,:),'Alpha',a);
[h_contra(2),p(2)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,3,:),'Alpha',a);
[h_contra(3),p(3)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,4,:),'Alpha',a);
[h_contra(4),p(4)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,5,:),'Alpha',a);
figure
roi = {'Stroke','Peri-infarct','Contralesional'};
x=1:grp;
b = bar(x,mcorr_roi(1:grp,:));
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
set(gca, 'XTick', '');
ylim([0 1.8])
set(gca,'FontSize',24)
sigline1([0.69 0.85],[],0.9)
sigline1([0.69 1],[],0.96)
sigline2([0.69 1.15],[],1.02)
sigline2([0.69 1.3],[],1.08)
sigline1([1.69 1.85],[],0.85)
sigline1([1.69 2],[],0.91)
sigline2([2.69 2.85],[],0.9)
saveas(gcf,'figures/figure4h_middle.png')

load('NVC_HbR_stats.mat')
grp = 3;
mcorr_roi = mean(corr_baseirf_roi,3);
stdcorr_roi = std(corr_baseirf_roi,[],3); %./sqrt(size(corr_roi,3));
a = 0.05;
[h_fore(1),p(1)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,2,:),'Alpha',a);
[h_fore(2),p(2)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,3,:),'Alpha',a);
[h_fore(3),p(3)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,4,:),'Alpha',a);
[h_fore(4),p(4)] = ttest2(corr_baseirf_roi(1,1,:),corr_baseirf_roi(1,5,:),'Alpha',a);
[h_peri(1),p(1)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,2,:),'Alpha',a);
[h_peri(2),p(2)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,3,:),'Alpha',a);
[h_peri(3),p(3)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,4,:),'Alpha',a);
[h_peri(4),p(4)] = ttest2(corr_baseirf_roi(2,1,:),corr_baseirf_roi(2,5,:),'Alpha',a);
[h_contra(1),p(1)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,2,:),'Alpha',a);
[h_contra(2),p(2)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,3,:),'Alpha',a);
[h_contra(3),p(3)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,4,:),'Alpha',a);
[h_contra(4),p(4)] = ttest2(corr_baseirf_roi(3,1,:),corr_baseirf_roi(3,5,:),'Alpha',a);
figure
roi = {'Stroke','Peri-infarct','Contralesional'};
x=1:grp;
b = bar(x,mcorr_roi(1:grp,:));
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
set(gca, 'XTick', '');
ylim([0 1.8])
set(gca,'FontSize',24)
sigline1([0.69 0.85],[],0.9)
sigline1([0.69 1],[],0.96)
sigline1([1.69 1.85],[],0.85)
sigline2([1.69 2],[],0.91)
sigline2([2.69 2.85],[],0.9)
saveas(gcf,'figures/figure4h_bottom.png')

%% Supplementary figure 9a
% NVC stats with post-stroke HRF
% Spatial maps of example animal

m = 7;
load('NVC_HbT_stats_post.mat')
fh = figure();
fh.WindowState = 'maximized';
count = 0;
for t = 1:5
    if t==1 || t==2 || t==3 || t==5
    count = count+1;
    subplot(3,4,count)
    imshow(Corr{m,t})
    axis image
    colormap turbo
    axis off
    caxis([0 1])
    end
end

load('NVC_HbO_stats_post.mat')
for t = 1:5
    if t==1 || t==2 || t==3 || t==5
    count = count+1;
    subplot(3,4,count)
    imshow(Corr{m,t})
    axis image
    colormap turbo
    axis off
    caxis([0 1])
    end
end

load('NVC_HbR_stats_post.mat')
for t = 1:5
    if t==1 || t==2 || t==3 || t==5
    count = count+1;
    subplot(3,4,count)
    imshow(Corr{m,t})
    axis image
    colormap turbo
    axis off
    caxis([0 1])
    end
end
saveas(gcf,'figures/suppfig9a.png')

%% Supplementary figure 9b
load('NVC_HbT_stats_post.mat')
grp = 3;
mcorr_roi = mean(corr_roi,3);
stdcorr_roi = std(corr_roi,[],3); %./sqrt(size(corr_roi,3));
a = 0.01;
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
figure
roi = {'Stroke','Peri-infarct','Contralesional'};
x=1:grp;
b = bar(x,mcorr_roi(1:grp,:));
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
set(gca, 'XTick', '');
ylim([0 1.3])
set(gca,'FontSize',24)
sigline1([0.69 0.85],[],0.9)
sigline1([0.69 1],[],0.96)
sigline1([0.69 1.15],[],1.02)
sigline1([0.69 1.3],[],1.08)
sigline1([1.69 1.85],[],0.85)
sigline1([1.69 2],[],0.91)
legend('Pre-Stroke','Day2','Week1','Week2','Week4')
saveas(gcf,'figures/suppfig9b_top.png')

load('NVC_HbO_stats_post.mat')
grp = 3;
mcorr_roi = mean(corr_roi,3);
stdcorr_roi = std(corr_roi,[],3); %./sqrt(size(corr_roi,3));
a = 0.01;
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
figure
roi = {'Stroke','Peri-infarct','Contralesional'};
x=1:grp;
b = bar(x,mcorr_roi(1:grp,:));
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
set(gca, 'XTick', '');
ylim([0 1.3])
set(gca,'FontSize',24)
sigline1([0.69 0.85],[],0.9)
sigline1([0.69 1],[],0.96)
sigline2([0.69 1.15],[],1.02)
sigline2([0.69 1.3],[],1.08)
sigline1([1.69 1.85],[],0.85)
sigline2([1.69 2],[],0.91)
saveas(gcf,'figures/suppfig9b_middle.png')

load('NVC_HbR_stats_post.mat')
grp = 3;
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
figure
roi = {'Stroke','Peri-infarct','Contralesional'};
x=1:grp;
b = bar(x,mcorr_roi(1:grp,:));
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
set(gca, 'XTick', '');
ylim([0 1.3])
set(gca,'FontSize',24)
sigline1([0.69 0.85],[],0.9)
sigline1([0.69 1],[],0.96)
sigline2([0.69 1.15],[],1.02)
sigline2([0.69 1.3],[],1.08)
sigline1([1.69 1.85],[],0.83)
sigline2([1.69 2],[],0.89)
saveas(gcf,'figures/suppfig9b_bottom.png')
