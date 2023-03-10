%% Estimating the hemodyanmic response function (HRF)
% HRF estimated using Least-Square Deconvolution

%%
clear

pathname{1} = '..MouseData/SS93/RestingState/Baseline';
pathname{2} = '..MouseData/SS93/RestingState/Day2';
pathname{3} = '..MouseData/SS93/RestingState/Week1';
pathname{4} = '..MouseData/SS93/RestingState/Week2';
pathname{5} = '..MouseData/SS93/RestingState/Week4';

for p = 1:length(pathname)
clearvars -except pathname p

dataDir = pathname{p}; %uigetdir('Please select the Data folder');
load([dataDir,'/','act_IOSI_ipsi.mat'])
load([dataDir,'/','act_GCaMPcorr_ipsi.mat'])
sub = dir(dataDir);
dfolders = sub([sub(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));


% Calculate and save IRF for each pixel
% Only run once

Fs = 5;
T = 1/Fs;
L = size(HbT,3);
timefull = (1:(L))*T;

% figure
% imagesc(mean(deltaGCaMPcorr,3))
% axis image
% colormap jet
% caxis([-0.015 0.015])
% axis off
% title('Corrected GCaMP')

deltaGCaMPcorrnew = deltaGCaMPcorr;
HbTnew = HbT;
start = 26;
stop = size(HbT,3);
lambda = 0.1;
frames = 50+start-1;
time = -5+0.2:0.2:10;
HRF = zeros(size(HbT,1),size(HbT,2),frames+3);
pHbT = zeros(size(HbT,1),size(HbT,2),stop-start+1);
peakTime = zeros(size(HbT,1),size(HbT,2));
amplitudePos = zeros(size(HbT,1),size(HbT,2));
amplitudeNeg = zeros(size(HbT,1),size(HbT,2));
fwhm = zeros(size(HbT,1),size(HbT,2));
for y1 = 1:size(HbTnew,1)
    waitbar(y1/size(HbTnew,1))
    for x1 = 1:size(HbTnew,2) 
        if isnan(deltaGCaMPcorrnew(y1,x1,1))
            HRF(y1,x1,:) = zeros(frames+3,1);
            peakTime(y1,x1) = NaN;
            amplitudePos(y1,x1) = NaN;
            amplitudeNeg(y1,x1) = NaN;
            fwhm(y1,x1) = NaN;
            correlation(y1,x1) = NaN;
            SSE(y1,x1) = NaN;
            RMSE(y1,x1) = NaN;  
        else
            tGCaMP = squeeze(deltaGCaMPcorrnew(y1,x1,start:stop));
            tHbT = squeeze(HbTnew(y1,x1,start:stop));
            app = [ones(size(timefull)); timefull; timefull.^2];
            app = app(:,start:stop)';
            I = eye(frames+3);
            temp = tGCaMP;
            X = convmtx(temp,length(temp));
            X = X(start:stop,1:frames);
            X = [app, X];
            Y = tHbT;
            h = (((X'*X) + (lambda.*I))^-1) * (X'*Y);
            HRF(y1,x1,:) = h; 
            
            idx = find(HRF(y1,x1,29:53) == max(HRF(y1,x1,29:53)));
            if length(idx)>1
                idx = idx(1);
            end
            peakTime(y1,x1) = time(25+idx);
            amplitudePos(y1,x1) = max(HRF(y1,x1,29:53));
            amplitudeNeg(y1,x1) = min(HRF(y1,x1,34:53));
            halfMax = amplitudePos(y1,x1)/2;      
            index1 = find(HRF(y1,x1,29:53) >= halfMax, 1, 'first');
            index2 = find(HRF(y1,x1,29:53) >= halfMax, 1, 'last');
            if isempty(index1) || isempty(index2)
                fwhm(y1,x1) = NaN;
            else
                fwhm(y1,x1) = time(index2) - time(index1);  
            end

            ptHbT = conv(squeeze(HRF(y1,x1,4:end)), tGCaMP);
            ptHbT = ptHbT(start:stop) +  app*squeeze(HRF(y1,x1,1:3));
            correlation(y1,x1) = corr(ptHbT,tHbT);
            SSE(y1,x1) = sum((ptHbT-tHbT).^2);
            RMSE(y1,x1) = sqrt(SSE(y1,x1)/length(ptHbT));
            pHbT(y1,x1,:) = ptHbT;
        end
    end
end
HRF(HRF==0) = NaN;

% figure
% histogram(correlation)
% xlabel('Correlation coefficient')
% ylabel('Pixels')
% xlim([0 1])
figure
imagesc(correlation)
axis off
axis image
colormap jet
caxis([0 1])
colorbar

save([dataDir,'/','act_HRF.mat'],'HRF','peakTime','amplitudePos','amplitudeNeg','fwhm','pHbT','correlation','SSE','RMSE')

end

%% Calculation based on IRF
% 
dataDir = uigetdir('Please select the Data folder');
load([dataDir,'/','act_IOSI_ipsi.mat'])
load([dataDir,'/','act_GCaMPcorr_ipsi.mat'])
sub = dir(dataDir);
dfolders = sub([sub(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
load([dataDir,'/','act_HRF.mat'])
for y = 1:size(HRF,1)
    for x = 1:size(HRF,2)
        sse(y,x) = norm(squeeze(HRF(y,x,10:25)));
    end
end

fh = figure(101);
fh.WindowState = 'maximized';
subplot(3,3,3)
imagesc(sse)
axis image
colormap jet
axis off
colorbar
caxis([0 3*10^-5])
title('Norm from -3 to 0 sec')

timefull = 0.2:0.2:600;
start = 26;
stop = 3000;

subplot(3,3,1)
imagesc(mean(deltaGCaMPcorr,3))
axis image
colormap jet
caxis([-0.015 0.015])
axis off
title('Mean corrected GCaMP')
colorbar

subplot(3,3,2)
imagesc(correlation)
axis image
colormap jet
caxis([0 1])
axis off
title('Correlation coefficient: GCaMP to HbT')
hold on
colorbar
color = ['k','b','m'];
for t = 1:3
    figure(101)
    subplot(3,3,2)
    h = imrect;
    fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
    setPositionConstraintFcn(h,fcn);
    pos = wait(h);
    delete(h);
    subplot(3,3,2)
    plot([pos(1),pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)],[pos(2),pos(2)+pos(4),pos(2)+pos(4),pos(2),pos(2)],color(t),'LineWidth',3);
    x(1) = int32(pos(1));
    x(2) = x(1)+int32(pos(3));
    y(1) = int32(pos(2));
    y(2) = y(1)+int32(pos(4));
    HRF_roi = HRF(y(1):y(2),x(1):x(2),:);
    HRF_roi = reshape(HRF_roi,[size(HRF_roi,1)*size(HRF_roi,2) size(HRF_roi,3)]);

    gcamp_roi = deltaGCaMPcorr(y(1):y(2),x(1):x(2),:);
    gcamp_roi = reshape(gcamp_roi,[size(gcamp_roi,1)*size(gcamp_roi,2) size(gcamp_roi,3)]);
    
    tHbT_roi = HbT(y(1):y(2),x(1):x(2),:);
    tHbT_roi = reshape(tHbT_roi,[size(tHbT_roi,1)*size(tHbT_roi,2) size(tHbT_roi,3)]);
    
    pHbT_roi = pHbT(y(1):y(2),x(1):x(2),:);
    pHbT_roi = reshape(pHbT_roi,[size(pHbT_roi,1)*size(pHbT_roi,2) size(pHbT_roi,3)]);
    
    corr_roi = correlation(y(1):y(2),x(1):x(2));
    peakTime_roi = peakTime(y(1):y(2),x(1):x(2));
    amplitude_roi = amplitudePos(y(1):y(2),x(1):x(2));
    fwhm_roi = fwhm(y(1):y(2),x(1):x(2));
    time = -4+0.2:0.2:10;
    
    figure(101)
    subplot(3,3,t+3)
    fig = gcf;
    options.x_axis = time;
    options.handle = figure(101);
    options.error = 'std';
    options.color_area = color(t); 
    options.color_line = color(t); 
    options.alpha = 0.5;
    options.line_width = 2;
    plot_areaerrorbar(HRF_roi(:,9:end),options)
    xlabel('Time (sec)')
    title('IRF')
    ylim([-1.5*10^-5 3*10^-5])
    hold on
    plot([0 0],[-1.5*10^-5 3*10^-5],'k--')
    
    mcorr_roi = mean(mean(corr_roi));
    mpeakTime = mean(mean(peakTime_roi));
    mamplitude = mean(mean(amplitude_roi));
    mfwhm = mean(mean(fwhm_roi));
%     str = {['Time to peak = ',num2str(mpeakTime),'s'],['Amplitude = ',num2str(mamplitude)],['Full-width half-max = ',num2str(mfwhm),'s']};
%     dim = [0.52 0.6 0.3 0.3];
%     annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    subplot(3,3,t+6)
    plot(timefull, squeeze(mean(tHbT_roi)))
    hold on
    plot(timefull(start:stop), squeeze(mean(pHbT_roi)))
    ylim([-5*1e-6 10*1e-6])
    ylabel('Change in HbT (\muM)')
    yyaxis right
    plot(timefull, squeeze(mean(gcamp_roi)))
    legend('Mesured HbT','Predicted HbT','Measured GCaMP')
    ylim([-0.05 0.15])
    xlabel('Time (sec)')
    ylabel('GCaMP (\DeltaF/F)')

end

[filepath1,name,ext] = fileparts(dataDir);
name1 = regexprep(name,'_',' ');
[filepath,name,ext] = fileparts(filepath1);
name2 = regexprep(name,'_',' ');
str = [name1 ' ' name2];
dim = [.01 .7 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',16);
saveas(gcf,[filepath1,'/','IRF ',name1,'.tif'])

%% HRF properties 

figure
subplot(1,4,1)
imagesc(correlation)
axis image
colormap jet
caxis([0 1])
axis off
colorbar
title('Correlation coefficient of fit')
subplot(1,4,2)
imagesc(amplitudePos)
axis image
colormap jet
% caxis([0 1])
axis off
colorbar
title('Amplitude')
subplot(1,4,3)
imagesc(peakTime)
axis image
colormap jet
caxis([0 3])
axis off
colorbar
title('Time to peak (sec)')
subplot(1,4,4)
imagesc(fwhm)
axis image
colormap jet
caxis([0 4])
axis off
colorbar
title('Width at half-max (sec)')


%% HRF through Least-Square Deconvolution 

timefull = 0.2:0.2:600;
figure(100)
for t = 1:2
subplot(4,3,t)
imagesc(mean(HbT,3))
axis image
colormap jet
caxis([-2e-6 2e-6])
axis off
title('HbT')
[x1,y1] = ginput(1);
x1 = round(x1);
y1 = round(y1);
hold on
plot(x1,y1,'ko','MarkerFaceColor','k');
hold off

start = 26;
time = -5+0.2:0.2:12;
stop = 3000;
tGCaMP = squeeze(deltaGCaMPcorr(y1,x1,start:stop));
% tGCaMP = movmean(tGCaMP,3);
tHbT = squeeze(HbT(y1,x1,start:stop));

lambda = 0.1;
frames = 50+start-1;
app = [ones(size(timefull)); timefull; timefull.^2; timefull.^3];
app = app(:,start:stop)';
I = eye(frames+4);
temp = tGCaMP;
X = convmtx(temp,length(temp));
X = X(start:stop,1:frames);
X = [app, X];
Y = tHbT;
h = (((X'*X) + (lambda.*I))^-1) * (X'*Y);
HRF = h; 

ptHbT = conv(HRF(5:end), tGCaMP);
ptHbT = ptHbT(start:stop) +  app*HRF(1:4);
% ptHbT = app*HRF(1:4);
correlation = corr(ptHbT,tHbT);
SSE = sum((ptHbT-tHbT).^2);
RMSE = sqrt(SSE/length(ptHbT));

figure(100)
if t == 1
    subplot(4,3,4)
    plot(time(1:frames),HRF(5:end))
    xlabel('Time')
    title('HRF')
    xlim([-5 10])
    subplot(4,3,[7 8 9])
    plot(timefull(1:length(tHbT)),tHbT)
    hold on
    plot(timefull(1:length(ptHbT)),ptHbT)
    xlabel('Time')
    ylabel('Change in HbT')
    title(['All 20 trials (r=',num2str((correlation),'%4.2f'),')',' (RMSE=',num2str(RMSE),')'])
    ylim([-6e-6 8e-6])
    stim = 1:30:600;
    x1 = [stim+5; stim+5; stim; stim];
    y1 = [-6e-6*ones(length(stim),1)'; 8e-6*ones(length(stim),1)'; 8e-6*ones(length(stim),1)'; -6e-6*ones(length(stim),1)'];
    patch(x1,y1,'black','EdgeAlpha',0.1,'FaceAlpha',0.1)
    yyaxis right
    plot(timefull(1:length(tGCaMP)),tGCaMP)    
    ylim([-0.1 0.1])
    legend('Measured HbT','Deconvolution model','','Measured GCaMP')

%     subplot(4,3,10)
%     plot(timefull(1:150),tHbT(1:150))
%     hold on
%     plot(timefull(1:150),ptHbT(1:150))
%     xlabel('Time')
%     title('Trial 1')
%     ylabel('Change in HbT')
%     ylim([-6e-6 8e-6])
%     stim = 5;
%     x1 = [stim+5; stim+5; stim; stim];
%     y1 = [-6e-6*ones(length(stim),1)'; 8e-6*ones(length(stim),1)'; 8e-6*ones(length(stim),1)'; -6e-6*ones(length(stim),1)'];
%     patch(x1,y1,'black','EdgeAlpha',0.1,'FaceAlpha',0.1)
elseif t == 2
    subplot(4,3,5)
    plot(time(1:frames),HRF(5:end))
    xlim([-5 10])
    xlabel('Time')
    title('HRF')
    subplot(4,3,[10 11 12])
    plot(timefull(1:length(tHbT)),tHbT)
    hold on
    plot(timefull(1:length(ptHbT)),ptHbT)
    xlabel('Time')
    ylabel('Change in HbT')
    title(['All 20 trials (r=',num2str((correlation),'%4.2f'),')',' (RMSE=',num2str(RMSE),')'])
    ylim([-6e-6 8e-6])
    stim = 1:30:600;
    x1 = [stim+5; stim+5; stim; stim];
    y1 = [-6e-6*ones(length(stim),1)'; 8e-6*ones(length(stim),1)'; 8e-6*ones(length(stim),1)'; -6e-6*ones(length(stim),1)'];
    patch(x1,y1,'black','EdgeAlpha',0.1,'FaceAlpha',0.1)
    yyaxis right
    plot(timefull(1:length(tGCaMP)),tGCaMP)    
    ylim([-0.1 0.1])
    legend('Measured HbT','Deconvolution model','Measured GCaMP')
%     subplot(4,3,11)
%     plot(timefull(1:150),tHbT(1:150))
%     hold on
%     plot(timefull(1:150),ptHbT(1:150))
%     xlabel('Time')
%     ylabel('Change in HbT')
%     ylim([-1.5e-6 4e-6])
%     stim = 5;
%     x1 = [stim+5; stim+5; stim; stim];
%     y1 = [-1.5*1e-6*ones(length(stim),1)'; 4*1e-6*ones(length(stim),1)'; 4*1e-6*ones(length(stim),1)'; -1.5*1e-6*ones(length(stim),1)'];
%     patch(x1,y1,'black','EdgeAlpha',0.1,'FaceAlpha',0.1)
%     title('Trial 1')
elseif t == 3
    subplot(4,3,6)
    plot(timefull(1:length(HRF)),HRF)
    xlabel('Time')
    title('HRF')
    subplot(4,3,9)
    plot(timefull(1:length(tHbT)),tHbT)
    hold on
    plot(timefull(1:length(ptHbT)),ptHbT)
    xlabel('Time')
    ylabel('Change in HbT')
    title(['All 10 trials (r=',num2str((correlation),'%4.2f'),')'])
    ylim([-6e-6 8e-6])
    stim = 5:30:600;
    x1 = [stim+5; stim+5; stim; stim];
    y1 = [-6e-6*ones(length(stim),1)'; 8e-6*ones(length(stim),1)'; 8e-6*ones(length(stim),1)'; -6e-6*ones(length(stim),1)'];
    patch(x1,y1,'black','EdgeAlpha',0.1,'FaceAlpha',0.1)
    subplot(4,3,12)
    plot(timefull(1:600),tHbT(1:600))
    hold on
    plot(timefull(1:600),ptHbT(1:600))
    xlabel('Time')
    ylabel('Change in HbT')
    ylim([-1.5e-6 4e-6])
    stim = 5;
    x1 = [stim+5; stim+5; stim; stim];
    y1 = [-1.5*1e-6*ones(length(stim),1)'; 4*1e-6*ones(length(stim),1)'; 4*1e-6*ones(length(stim),1)'; -1.5*1e-6*ones(length(stim),1)'];
    patch(x1,y1,'black','EdgeAlpha',0.1,'FaceAlpha',0.1)
    title('Trial 1')
end
end

