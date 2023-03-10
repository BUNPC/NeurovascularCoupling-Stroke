%% Spectrogram analyzer for resting state

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
fs = 5;
lpf = 0.4; % 2; %0.08;
hpf = 0.12; %0.02; %0.009;
wn(2) = lpf/(fs/2);
wn(1) = hpf/(fs/2);
[fb,fa] = butter(2,wn,'bandpass');
for m = 1:length(mouse)
m
baseline = [mouse{m},'RestingState/Baseline'];
day2 = [mouse{m},'RestingState/Day2'];
week1 = [mouse{m},'RestingState/Week1'];
week2 = [mouse{m},'RestingState/Week2'];
week4 = [mouse{m},'RestingState/Week4'];
mask = [mouse{m},'brainmaskSFDI.mat'];
pathname = {baseline, day2, week1, week2, week4};
mask1 = [mouse{m},'brainmaskSFDI.mat'];
load(mask1)
for t = 1:length(pathname)
    t
    time = -5+0.2:0.2:10;
    dataDir = pathname{t};
    
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
    HbTGS(m,t,:) = yMean;
    HbT = Chb1;
    rawsigH = reshape(HbT,[nY*nX nT]);
    filt_rawsigH = filtfilt(fb,fa,rawsigH')';
    filt_rawsigH = reshape(filt_rawsigH,[nY nX nT]);
    varH{m,t} = var(filt_rawsigH,[],3);
     
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
    deltaGCaMPcorrGS(m,t,:) = yMean;
    deltaGCaMPcorr = Chb1;
    rawsigG = reshape(deltaGCaMPcorr,[nY*nX nT]);
    filt_rawsigG = filtfilt(fb,fa,rawsigG')';
    filt_rawsigG = reshape(filt_rawsigG,[nY nX nT]);
    varG{m,t} = var(filt_rawsigG,[],3);
    
    HbT_p = reshape(HbT,[nY*nX,nT]);
    pixpow = bandpower(HbT_p',fs,[0.002 0.5]);
    pixpow = reshape(pixpow,[nY nX]);
    HbT_pwr_all{m,t} = pixpow;
    clear pixpow
    G_p = reshape(deltaGCaMPcorr,[nY*nX,nT]);
    pixpow = bandpower(G_p',fs,[0.002 0.5]);
    pixpow = reshape(pixpow,[nY nX]);
    GCaMP_pwr_all{m,t} = pixpow;
    clear pixpow
    HbT_p = reshape(HbT,[nY*nX,nT]);
    pixpow = bandpower(HbT_p',fs,[0.009 0.08]);
    pixpow = reshape(pixpow,[nY nX]);
    HbT_pwr_low{m,t} = pixpow;
    clear pixpow
    G_p = reshape(deltaGCaMPcorr,[nY*nX,nT]);
    pixpow = bandpower(G_p',fs,[0.009 0.08]);
    pixpow = reshape(pixpow,[nY nX]);
    GCaMP_pwr_low{m,t} = pixpow;
    clear pixpow
    HbT_p = reshape(HbT,[nY*nX,nT]);
    pixpow = bandpower(HbT_p',fs,[0.12 0.4]);
    pixpow = reshape(pixpow,[nY nX]);
    HbT_pwr_high{m,t} = pixpow;
    clear pixpow
    G_p = reshape(deltaGCaMPcorr,[nY*nX,nT]);
    pixpow = bandpower(G_p',fs,[0.12 0.4]);
    pixpow = reshape(pixpow,[nY nX]);
    GCaMP_pwr_high{m,t} = pixpow;
    clear pixpow
    
    if t == 1 || t == 2
        newMask = imdilate(maskSFDI(3).stroke_mask, true(30));
        peri = abs(newMask - maskSFDI(3).stroke_mask);       
        stroke = maskSFDI(3).stroke_mask;
    else        
        newMask = imdilate(maskSFDI(t).stroke_mask, true(30));
        peri = abs(newMask - maskSFDI(t).stroke_mask);
        stroke = maskSFDI(t).stroke_mask;
    end
    ipsi = maskSFDI(t).aff_mask;
    contra = maskSFDI(t).unaff_mask;
    contraforelimb = maskSFDI(1).contraOutline;
    ipsiforelimb = maskSFDI(1).ipsiOutline;
    ipsi_nonfore = abs(maskSFDI(t).aff_mask - maskSFDI(1).ipsiOutline);
    contra_nonfore = abs(maskSFDI(t).unaff_mask - maskSFDI(1).contraOutline);
    
    HbT_core = HbT.*stroke;
    HbT_core = reshape(HbT_core,[size(HbT_core,1)*size(HbT_core,2) size(HbT_core,3)]);
    gcamp_core = deltaGCaMPcorr.*stroke;
    gcamp_core = reshape(gcamp_core,[size(gcamp_core,1)*size(gcamp_core,2) size(gcamp_core,3)]);
    mHbT_roi = mean(HbT_core,2);
    idx = find(mHbT_roi~=0);
    HbT_core = HbT_core(idx,:);
    gcamp_core = gcamp_core(idx,:);
    N = size(HbT_core,2);
    freq = 0:fs/N:fs/2;
    spectHb_core = zeros(length(idx),N);
    spectG_core = zeros(length(idx),N);
    for f = 1:length(idx)
        spectHb_core(f,:) = fft(HbT_core(f,:));
        spectG_core(f,:) = fft(gcamp_core(f,:));        
    end
    PSDH = (1/(fs*N)).*(abs(spectHb_core(:,1:N/2+1)).^2);
    PSDH(:,2:end-1) = 2.*PSDH(:,2:end-1);
    powerHb_core(m,t,:) = mean(PSDH,1);
    PSDG = (1/(fs*N)).*(abs(spectG_core(:,1:N/2+1)).^2);
    PSDG(:,2:end-1) = 2.*PSDG(:,2:end-1);
    powerG_core(m,t,:) = mean(PSDG,1);
    clear mHbT_roi PSDG PSDH
    
    HbT_peri = HbT.*peri;
    HbT_peri = reshape(HbT_peri,[size(HbT_peri,1)*size(HbT_peri,2) size(HbT_peri,3)]);
    gcamp_peri = deltaGCaMPcorr.*peri;
    gcamp_peri = reshape(gcamp_peri,[size(gcamp_peri,1)*size(gcamp_peri,2) size(gcamp_peri,3)]);
    mHbT_roi = mean(HbT_peri,2);
    idx = find(mHbT_roi~=0);
    HbT_peri = HbT_peri(idx,:);
    gcamp_peri = gcamp_peri(idx,:);
    N = size(HbT_peri,2);
    spectHb_peri = zeros(length(idx),N);
    spectG_peri = zeros(length(idx),N);
    for f = 1:length(idx)
        spectHb_peri(f,:) = fft(HbT_peri(f,:));
        spectG_peri(f,:) = fft(gcamp_peri(f,:));
        pixpow(f,:) = bandpower(HbT_peri(f,:),fs,[0 fs/2]);
    end    
    PSDH = (1/(fs*N)).*(abs(spectHb_peri(:,1:N/2+1)).^2);
    PSDH(:,2:end-1) = 2.*PSDH(:,2:end-1);
    powerHb_peri(m,t,:) = mean(PSDH,1);
    PSDG = (1/(fs*N)).*(abs(spectG_peri(:,1:N/2+1)).^2);
    PSDG(:,2:end-1) = 2.*PSDG(:,2:end-1);
    powerG_peri(m,t,:) = mean(PSDG,1);
    clear mHbT_roi PSDG PSDH
    
    HbT_ipsi = HbT.*ipsi;
    HbT_ipsi = reshape(HbT_ipsi,[size(HbT_ipsi,1)*size(HbT_ipsi,2) size(HbT_ipsi,3)]);
    gcamp_ipsi = deltaGCaMPcorr.*ipsi;
    gcamp_ipsi = reshape(gcamp_ipsi,[size(gcamp_ipsi,1)*size(gcamp_ipsi,2) size(gcamp_ipsi,3)]);
    mHbT_roi = mean(HbT_ipsi,2);
    idx = find(mHbT_roi~=0);
    HbT_ipsi = HbT_ipsi(idx,:);
    gcamp_ipsi = gcamp_ipsi(idx,:);
    N = size(HbT_ipsi,2);
    spectHb_ipsi = zeros(length(idx),N);
    spectG_ipsi = zeros(length(idx),N);
    for f = 1:length(idx)
        spectHb_ipsi(f,:) = fft(HbT_ipsi(f,:));
        spectG_ipsi(f,:) = fft(gcamp_ipsi(f,:));
        pixpow(f,:) = bandpower(HbT_ipsi(f,:),fs,[0 fs/2]);
    end    
    PSDH = (1/(fs*N)).*(abs(spectHb_ipsi(:,1:N/2+1)).^2);
    PSDH(:,2:end-1) = 2.*PSDH(:,2:end-1);
    powerHb_ipsi(m,t,:) = mean(PSDH,1);
    PSDG = (1/(fs*N)).*(abs(spectG_ipsi(:,1:N/2+1)).^2);
    PSDG(:,2:end-1) = 2.*PSDG(:,2:end-1);
    powerG_ipsi(m,t,:) = mean(PSDG,1);
    clear mHbT_roi PSDG PSDH
    
    HbT_contra = HbT.*contra;
    HbT_contra = reshape(HbT_contra,[size(HbT_contra,1)*size(HbT_contra,2) size(HbT_contra,3)]);
    gcamp_contra = deltaGCaMPcorr.*contra;
    gcamp_contra = reshape(gcamp_contra,[size(gcamp_contra,1)*size(gcamp_contra,2) size(gcamp_contra,3)]);
    mHbT_roi = mean(HbT_contra,2);
    idx = find(mHbT_roi~=0);
    HbT_contra = HbT_contra(idx,:);
    gcamp_contra = gcamp_contra(idx,:);
    N = size(HbT_contra,2);
    spectHb_contra = zeros(length(idx),N);
    spectG_contra = zeros(length(idx),N);
    for f = 1:length(idx)
        spectHb_contra(f,:) = fft(HbT_contra(f,:));
        spectG_contra(f,:) = fft(gcamp_contra(f,:));
    end
    PSDH = (1/(fs*N)).*(abs(spectHb_contra(:,1:N/2+1)).^2);
    PSDH(:,2:end-1) = 2.*PSDH(:,2:end-1);
    powerHb_contra(m,t,:) = mean(PSDH,1);
    PSDG = (1/(fs*N)).*(abs(spectG_contra(:,1:N/2+1)).^2);
    PSDG(:,2:end-1) = 2.*PSDG(:,2:end-1);
    powerG_contra(m,t,:) = mean(PSDG,1);
    clear mHbT_roi PSDG PSDH
    
    HbT_contraforelimb = HbT.*contraforelimb;
    HbT_contraforelimb = reshape(HbT_contraforelimb,[size(HbT_contraforelimb,1)*size(HbT_contraforelimb,2) size(HbT_contraforelimb,3)]);
    gcamp_contraforelimb = deltaGCaMPcorr.*contraforelimb;
    gcamp_contraforelimb = reshape(gcamp_contraforelimb,[size(gcamp_contraforelimb,1)*size(gcamp_contraforelimb,2) size(gcamp_contraforelimb,3)]);
    mHbT_roi = mean(HbT_contraforelimb,2);
    idx = find(mHbT_roi~=0);
    HbT_contraforelimb = HbT_contraforelimb(idx,:);
    gcamp_contraforelimb = gcamp_contraforelimb(idx,:);
    N = size(HbT_contraforelimb,2);
    spectHb_contraforelimb = zeros(length(idx),N);
    spectG_contraforelimb = zeros(length(idx),N);
    for f = 1:length(idx)
        spectHb_contraforelimb(f,:) = fft(HbT_contraforelimb(f,:));
        spectG_contraforelimb(f,:) = fft(gcamp_contraforelimb(f,:));
    end
    PSDH = (1/(fs*N)).*(abs(spectHb_contraforelimb(:,1:N/2+1)).^2);
    PSDH(:,2:end-1) = 2.*PSDH(:,2:end-1);
    powerHb_contraforelimb(m,t,:) = mean(PSDH,1);
    PSDG = (1/(fs*N)).*(abs(spectG_contraforelimb(:,1:N/2+1)).^2);
    PSDG(:,2:end-1) = 2.*PSDG(:,2:end-1);
    powerG_contraforelimb(m,t,:) = mean(PSDG,1);
    clear mHbT_roi PSDG PSDH
    
    HbT_ipsiforelimb = HbT.*ipsiforelimb;
    HbT_ipsiforelimb = reshape(HbT_ipsiforelimb,[size(HbT_ipsiforelimb,1)*size(HbT_ipsiforelimb,2) size(HbT_ipsiforelimb,3)]);
    gcamp_ipsiforelimb = deltaGCaMPcorr.*ipsiforelimb;
    gcamp_ipsiforelimb = reshape(gcamp_ipsiforelimb,[size(gcamp_ipsiforelimb,1)*size(gcamp_ipsiforelimb,2) size(gcamp_ipsiforelimb,3)]);
    mHbT_roi = mean(HbT_ipsiforelimb,2);
    idx = find(mHbT_roi~=0);
    HbT_ipsiforelimb = HbT_ipsiforelimb(idx,:);
    gcamp_ipsiforelimb = gcamp_ipsiforelimb(idx,:);
    N = size(HbT_ipsiforelimb,2);
    spectHb_ipsiforelimb = zeros(length(idx),N);
    spectG_ipsiforelimb = zeros(length(idx),N);
    for f = 1:length(idx)
        spectHb_ipsiforelimb(f,:) = fft(HbT_ipsiforelimb(f,:));
        spectG_ipsiforelimb(f,:) = fft(gcamp_ipsiforelimb(f,:));
    end
    PSDH = (1/(fs*N)).*(abs(spectHb_ipsiforelimb(:,1:N/2+1)).^2);
    PSDH(:,2:end-1) = 2.*PSDH(:,2:end-1);
    powerHb_ipsiforelimb(m,t,:) = mean(PSDH,1);
    PSDG = (1/(fs*N)).*(abs(spectG_ipsiforelimb(:,1:N/2+1)).^2);
    PSDG(:,2:end-1) = 2.*PSDG(:,2:end-1);
    powerG_ipsiforelimb(m,t,:) = mean(PSDG,1);
    clear mHbT_roi PSDG PSDH
    
    HbT_ipsi_nonfore = HbT.*ipsi_nonfore;
    HbT_ipsi_nonfore = reshape(HbT_ipsi_nonfore,[size(HbT_ipsi_nonfore,1)*size(HbT_ipsi_nonfore,2) size(HbT_ipsi_nonfore,3)]);
    gcamp_ipsi_nonfore = deltaGCaMPcorr.*ipsi_nonfore;
    gcamp_ipsi_nonfore = reshape(gcamp_ipsi_nonfore,[size(gcamp_ipsi_nonfore,1)*size(gcamp_ipsi_nonfore,2) size(gcamp_ipsi_nonfore,3)]);
    mHbT_roi = mean(HbT_ipsi_nonfore,2);
    idx = find(mHbT_roi~=0);
    HbT_ipsi_nonfore = HbT_ipsi_nonfore(idx,:);
    gcamp_ipsi_nonfore = gcamp_ipsi_nonfore(idx,:);
    N = size(HbT_ipsi_nonfore,2);
    spectHb_ipsi_nonfore = zeros(length(idx),N);
    spectG_ipsi_nonfore = zeros(length(idx),N);
    for f = 1:length(idx)
        spectHb_ipsi_nonfore(f,:) = fft(HbT_ipsi_nonfore(f,:));
        spectG_ipsi_nonfore(f,:) = fft(gcamp_ipsi_nonfore(f,:));
    end
    PSDH = (1/(fs*N)).*(abs(spectHb_ipsi_nonfore(:,1:N/2+1)).^2);
    PSDH(:,2:end-1) = 2.*PSDH(:,2:end-1);
    powerHb_ipsi_nonfore(m,t,:) = mean(PSDH,1);
    PSDG = (1/(fs*N)).*(abs(spectG_ipsi_nonfore(:,1:N/2+1)).^2);
    PSDG(:,2:end-1) = 2.*PSDG(:,2:end-1);
    powerG_ipsi_nonfore(m,t,:) = mean(PSDG,1);
    clear mHbT_roi PSDG PSDH
    
    HbT_contra_nonfore = HbT.*contra_nonfore;
    HbT_contra_nonfore = reshape(HbT_contra_nonfore,[size(HbT_contra_nonfore,1)*size(HbT_contra_nonfore,2) size(HbT_contra_nonfore,3)]);
    gcamp_contra_nonfore = deltaGCaMPcorr.*contra_nonfore;
    gcamp_contra_nonfore = reshape(gcamp_contra_nonfore,[size(gcamp_contra_nonfore,1)*size(gcamp_contra_nonfore,2) size(gcamp_contra_nonfore,3)]);
    mHbT_roi = mean(HbT_contra_nonfore,2);
    idx = find(mHbT_roi~=0);
    HbT_contra_nonfore = HbT_contra_nonfore(idx,:);
    gcamp_contra_nonfore = gcamp_contra_nonfore(idx,:);
    N = size(HbT_contra_nonfore,2);
    spectHb_contra_nonfore = zeros(length(idx),N);
    spectG_contra_nonfore = zeros(length(idx),N);
    for f = 1:length(idx)
        spectHb_contra_nonfore(f,:) = fft(HbT_contra_nonfore(f,:));
        spectG_contra_nonfore(f,:) = fft(gcamp_contra_nonfore(f,:));
    end
    PSDH = (1/(fs*N)).*(abs(spectHb_contra_nonfore(:,1:N/2+1)).^2);
    PSDH(:,2:end-1) = 2.*PSDH(:,2:end-1);
    powerHb_contra_nonfore(m,t,:) = mean(PSDH,1);
    PSDG = (1/(fs*N)).*(abs(spectG_contra_nonfore(:,1:N/2+1)).^2);
    PSDG(:,2:end-1) = 2.*PSDG(:,2:end-1);
    powerG_contra_nonfore(m,t,:) = mean(PSDG,1);
    clear mHbT_roi PSDG PSDH
end
end

%%
save('spect_rest_HbO.mat','mouse','HbTGS','deltaGCaMPcorrGS','HbT_pwr_all','GCaMP_pwr_all',...
    'HbT_pwr_low','GCaMP_pwr_low','HbT_pwr_high','GCaMP_pwr_high','varG','varH',...
    'powerHb_core','powerG_core','powerHb_peri','powerG_peri','powerHb_ipsi','powerG_ipsi',...
    'powerHb_contra','powerG_contra','powerHb_contraforelimb','powerG_contraforelimb',...
    'powerHb_ipsiforelimb','powerG_ipsiforelimb','powerHb_ipsi_nonfore','powerG_ipsi_nonfore',...
    'powerHb_contra_nonfore','powerG_contra_nonfore','freq','fs','N','-v7.3')

%% Spatial power

for m = 1:12
    mask1 = [mouse{m},'brainmaskSFDI.mat'];
    load(mask1)
    for t = 1:5
        forelimb_ipsi = maskSFDI(1).ipsiOutline;
        powerG_aff(m,t) = sum(sum(GCaMP_pwr_high{m,t}.*forelimb_ipsi))./sum(sum(forelimb_ipsi));
        powerH_aff(m,t) = sum(sum(HbT_pwr_high{m,t}.*forelimb_ipsi))./sum(sum(forelimb_ipsi));
        forelimb_contra = maskSFDI(1).contraOutline;
        powerG_unaff(m,t) = sum(sum(GCaMP_pwr_high{m,t}.*forelimb_contra))./sum(sum(forelimb_contra));
        powerH_unaff(m,t) = sum(sum(HbT_pwr_high{m,t}.*forelimb_contra))./sum(sum(forelimb_contra));
    end
end

mpwrG_peri = mean(powerG_aff,1);
mpwrH_peri = mean(powerH_aff,1);
mpwrG_contra = mean(powerG_unaff,1);
mpwrH_contra = mean(powerH_unaff,1);
mpwr_G = [mpwrG_peri;mpwrG_contra];
mpwr_H = [mpwrH_peri;mpwrH_contra];

spwrG_peri = std(powerG_aff,1);
spwrH_peri = std(powerH_aff,1);
spwrG_contra = std(powerG_unaff,1);
spwrH_contra = std(powerH_unaff,1);
spwr_G = [spwrG_peri;spwrG_contra];
spwr_H = [spwrH_peri;spwrH_contra];

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
set(gca, 'XTickLabel', group, 'FontSize', 16);
ylabel('Power: GCaMP')
set(b, 'FaceAlpha', 1)



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
set(gca, 'XTickLabel', group, 'FontSize', 16);
ylabel('Power: HbO')
set(b, 'FaceAlpha', 1)
legend('Pre-stroke','Day2','Week1','Week2','Week4')

%%

figure
f = freq;
for m = 1:length(mouse)
    subplot(3,4,m)
    hold on
    for t = 1:3
        if t == 1
        semilogx(f,squeeze((powerHb_peri(m,t,:))),'k-')
        ylabel('Power: HbT')
        xlim([0.02 1])
        elseif t == 2
        semilogx(f,squeeze((powerHb_peri(m,t,:))),'r-')
        ylabel('Power: HbT')
        xlim([0.02 1])
        elseif t == 3
        semilogx(f,squeeze((powerHb_peri(m,t,:))),'g-')
        ylabel('Power: HbT')
        xlim([0.02 1])
        end
    end
end



%% Global signal

fs = 5;
HbTGS = HbTGS(:,:,1:N/2);
mHbTGSBaseline = nanmean(squeeze(HbTGS(:,1,:)));
mHbTGSDay2 = nanmean(squeeze(HbTGS(:,2,:)));
mHbTGSWeek1 = nanmean(squeeze(HbTGS(:,3,:)));
mHbTGSWeek2 = nanmean(squeeze(HbTGS(:,4,:)));
mHbTGSWeek4 = nanmean(squeeze(HbTGS(:,5,:)));
[h2,p2] = ttest2(HbTGS(:,1,:),HbTGS(:,2,:),'Alpha',0.05);
[h3,p3] = ttest2(HbTGS(:,1,:),HbTGS(:,3,:),'Alpha',0.05);
[h4,p4] = ttest2(HbTGS(:,1,:),HbTGS(:,3,:),'Alpha',0.05);

figure
semilogx(f,mHbTGSBaseline,'k-','LineWidth',2)
hold on
semilogx(f,mHbTGSDay2,'r-','LineWidth',2)
semilogx(f,mHbTGSWeek1,'g-','LineWidth',2)
semilogx(f,mHbTGSWeek2,'b-','LineWidth',2)
semilogx(f,mHbTGSWeek4,'m-','LineWidth',2)
% semilogx(f(find(h2 == 1)),0.9e-6,'r*')
% semilogx(f(find(h3 == 1)),0.8e-6,'g*')
legend('Pre-stroke','Day2','Week1')
ylabel('Power: HbT Global Signal')
xlabel('Frequency (Hz)')
% ylim([0 7e-8])
xlim([0.02 1])
set(gca,'XScale','log','FontSize',16)

deltaGCaMPcorrGS = deltaGCaMPcorrGS(:,:,1:N/2);
mdeltaGCaMPcorrGSBaseline = nanmean(squeeze(deltaGCaMPcorrGS(:,1,:)));
mdeltaGCaMPcorrGSDay2 = nanmean(squeeze(deltaGCaMPcorrGS(:,2,:)));
mdeltaGCaMPcorrGSWeek1 = nanmean(squeeze(deltaGCaMPcorrGS(:,3,:)));
mdeltaGCaMPcorrGSWeek2 = nanmean(squeeze(deltaGCaMPcorrGS(:,4,:)));
mdeltaGCaMPcorrGSWeek4 = nanmean(squeeze(deltaGCaMPcorrGS(:,5,:)));
[h2,p2] = ttest2(deltaGCaMPcorrGS(:,1,:),deltaGCaMPcorrGS(:,2,:),'Alpha',0.05);
[h3,p3] = ttest2(deltaGCaMPcorrGS(:,1,:),deltaGCaMPcorrGS(:,3,:),'Alpha',0.05);

figure
semilogx(f,mdeltaGCaMPcorrGSBaseline,'k-','LineWidth',2)
hold on
semilogx(f,mdeltaGCaMPcorrGSDay2,'r-','LineWidth',2)
semilogx(f,mdeltaGCaMPcorrGSWeek1,'g-','LineWidth',2)
semilogx(f,mdeltaGCaMPcorrGSWeek2,'b-','LineWidth',2)
semilogx(f,mdeltaGCaMPcorrGSWeek4,'m-','LineWidth',2)
% semilogx(f(find(h2 == 1)),0.015,'r*')
% semilogx(f(find(h3 == 1)),0.013,'g*')
legend('Pre-stroke','Day2','Week1')
ylabel('Power: GCaMP Global Signal')
xlabel('Frequency (Hz)')
% ylim([0 7e-8])
xlim([0.02 1])
set(gca,'XScale','log','FontSize',16)

%%
powerHb_peri = powerHb_peri(:,:,1:N/2);
mPowerHbBaseline_peri = nanmean(squeeze(powerHb_peri(:,1,:)));
mPowerHbDay2_peri = nanmean(squeeze(powerHb_peri(:,2,:)));
mPowerHbWeek1_peri = nanmean(squeeze(powerHb_peri(:,3,:)));
mPowerHbWeek2_peri = nanmean(squeeze(powerHb_peri(:,4,:)));
mPowerHbWeek4_peri = nanmean(squeeze(powerHb_peri(:,5,:)));
a = 0.01;
[h2,p2] = ttest2(powerHb_peri(:,1,:),powerHb_peri(:,2,:),'Alpha',a);
[h3,p3] = ttest2(powerHb_peri(:,1,:),powerHb_peri(:,3,:),'Alpha',a);
[h4,p4] = ttest2(powerHb_peri(:,1,:),powerHb_peri(:,4,:),'Alpha',a);
[h5,p5] = ttest2(powerHb_peri(:,1,:),powerHb_peri(:,5,:),'Alpha',a);

figure
semilogx(f,mPowerHbBaseline_peri,'Color',[0 0 0],'LineWidth',2)
hold on
semilogx(f,mPowerHbDay2_peri,'Color',[0.5 0 0],'LineWidth',2)
semilogx(f,mPowerHbWeek1_peri,'Color',[1 0 0],'LineWidth',2)
semilogx(f,mPowerHbWeek2_peri,'Color',[1 0.5 0],'LineWidth',2)
semilogx(f,mPowerHbWeek4_peri,'Color',[1 1 0],'LineWidth',2)
semilogx(f(find(h2 == 1)),3.8e-11,'r*')
semilogx(f(find(h3 == 1)),3.6e-11,'g*')
semilogx(f(find(h4 == 1)),3.4e-11,'b*')
semilogx(f(find(h5 == 1)),3.2e-11,'m*')
legend('Pre-stroke','Day2','Week1','Week2','Week4')
ylabel('Power: HbT')
xlabel('Frequency (Hz)')
% ylim([0 7e-8])
xlim([0.02 1])
set(gca,'XScale','log','FontSize',16)

powerHb_contra = powerHb_contra(:,:,1:N/2);
mPowerHbBaseline_contra = nanmean(squeeze(powerHb_contra(:,1,:)));
mPowerHbDay2_contra = nanmean(squeeze(powerHb_contra(:,2,:)));
mPowerHbWeek1_contra = nanmean(squeeze(powerHb_contra(:,3,:)));
mPowerHbWeek2_contra = nanmean(squeeze(powerHb_contra(:,4,:)));
mPowerHbWeek4_contra = nanmean(squeeze(powerHb_contra(:,5,:)));
[h2,p2] = ttest2(powerHb_contra(:,1,:),powerHb_contra(:,2,:),'Alpha',a);
[h3,p3] = ttest2(powerHb_contra(:,1,:),powerHb_contra(:,3,:),'Alpha',a);
[h4,p4] = ttest2(powerHb_contra(:,1,:),powerHb_contra(:,4,:),'Alpha',a);
[h5,p5] = ttest2(powerHb_contra(:,1,:),powerHb_contra(:,5,:),'Alpha',a);

figure
semilogx(f,mPowerHbBaseline_contra,'k-','LineWidth',2)
hold on
semilogx(f,mPowerHbDay2_contra,'r-','LineWidth',2)
semilogx(f,mPowerHbWeek1_contra,'g-','LineWidth',2)
semilogx(f,mPowerHbWeek2_contra,'b-','LineWidth',2)
semilogx(f,mPowerHbWeek4_contra,'m-','LineWidth',2)
semilogx(f(find(h2 == 1)),3.8e-11,'r*')
semilogx(f(find(h3 == 1)),3.6e-11,'g*')
semilogx(f(find(h4 == 1)),3.4e-11,'b*')
semilogx(f(find(h5 == 1)),3.2e-11,'m*')
legend('Pre-stroke','Day2','Week1','Week2','Week4')
ylabel('Power: HbT')
xlabel('Frequency (Hz)')
% ylim([0 7e-8])
xlim([0.02 1])
set(gca,'XScale','log','FontSize',16)

% sPowerHbBaseline = std(squeeze(powerHb_peri(:,1,:)))/sqrt(size(powerHb_peri,1));
% sPowerHbDay2 = std(squeeze(powerHb_peri(:,2,:)))/sqrt(size(powerHb_peri,1));
% sPowerHbWeek1 = std(squeeze(powerHb_peri(:,3,:)))/sqrt(size(powerHb_peri,1));
% sPowerHbWeek2 = std(squeeze(powerHb_peri(:,4,:)))/sqrt(size(powerHb_peri,1));
figure
fig = gcf;
options.x_axis = f;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = 'k'; 
options.color_line = 'k'; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(powerHb_peri(:,1,:)),options)
hold on
options.color_area = 'r'; 
options.color_line = 'r'; 
hold on
plot_areaerrorbar(squeeze(powerHb_peri(:,2,:)),options)
options.color_area = 'g'; 
options.color_line = 'g'; 
hold on
plot_areaerrorbar(squeeze(powerHb_peri(:,3,:)),options)
options.color_area = 'b'; 
options.color_line = 'b'; 
hold on
plot_areaerrorbar(squeeze(powerHb_peri(:,4,:)),options)
options.color_area = 'm'; 
options.color_line = 'm'; 
hold on
plot_areaerrorbar(squeeze(powerHb_peri(:,5,:)),options)
hold on
plot(f(find(h2 == 1)),3.8e-11,'r*')
plot(f(find(h3 == 1)),3.6e-11,'g*')
plot(f(find(h4 == 1)),3.4e-11,'b*')
plot(f(find(h5 == 1)),3.2e-11,'m*')
% ylim([0 3e-8])
xlim([0.02 1])
ylabel('Power: HbT')
xlabel('Frequency (Hz)')
set(gca,'FontSize',16)

%%

fs = 5;
powerG_peri = powerG_peri(:,:,1:N/2);
mPowerGBaseline_peri = nanmean(squeeze(powerG_peri(:,1,:)));
mPowerGDay2_peri = nanmean(squeeze(powerG_peri(:,2,:)));
mPowerGWeek1_peri = nanmean(squeeze(powerG_peri(:,3,:)));
mPowerGWeek2_peri = nanmean(squeeze(powerG_peri(:,4,:)));
mPowerGWeek4_peri = nanmean(squeeze(powerG_peri(:,5,:)));
[h2,p2] = ttest2(powerG_peri(:,1,:),powerG_peri(:,2,:),'Alpha',0.05);
[h3,p3] = ttest2(powerG_peri(:,1,:),powerG_peri(:,3,:),'Alpha',0.05);

figure
semilogx(f,mPowerGBaseline_peri,'k-','LineWidth',2)
hold on
semilogx(f,mPowerGDay2_peri,'r-','LineWidth',2)
semilogx(f,mPowerGWeek1_peri,'g-','LineWidth',2)
semilogx(f,mPowerGWeek2_peri,'b-','LineWidth',2)
semilogx(f,mPowerGWeek4_peri,'m-','LineWidth',2)
% semilogx(f(find(h2 == 1)),3.9,'r*')
% semilogx(f(find(h3 == 1)),3.8,'g*')
legend('Pre-stroke','Day2','Week1')
ylabel('Power: GCaMP')
xlabel('Frequency (Hz)')
% ylim([0 4])
xlim([0.02 1])
set(gca,'XScale','log','FontSize',16)

powerG_contra = powerG_contra(:,:,1:N/2);
mPowerGBaseline_contra = nanmean(squeeze(powerG_contra(:,1,:)));
mPowerGDay2_contra = nanmean(squeeze(powerG_contra(:,2,:)));
mPowerGWeek1_contra = nanmean(squeeze(powerG_contra(:,3,:)));
mPowerGWeek2_contra = nanmean(squeeze(powerG_contra(:,4,:)));
mPowerGWeek4_contra = nanmean(squeeze(powerG_contra(:,5,:)));
[h2,p2] = ttest2(powerG_contra(:,1,:),powerG_contra(:,2,:),'Alpha',0.05);
[h3,p3] = ttest2(powerG_contra(:,1,:),powerG_contra(:,3,:),'Alpha',0.05);

figure
semilogx(f,mPowerGBaseline_contra,'k-','LineWidth',2)
hold on
semilogx(f,mPowerGDay2_contra,'r-','LineWidth',2)
semilogx(f,mPowerGWeek1_contra,'g-','LineWidth',2)
semilogx(f,mPowerGWeek2_contra,'b-','LineWidth',2)
semilogx(f,mPowerGWeek4_contra,'m-','LineWidth',2)
% semilogx(f(find(h2 == 1)),3.9,'r*')
% semilogx(f(find(h3 == 1)),3.8,'g*')
legend('Pre-stroke','Day2','Week1')
ylabel('Power: GCaMP')
xlabel('Frequency (Hz)')
% ylim([0 4])
xlim([0.02 1])
set(gca,'XScale','log','FontSize',16)
% 
figure
fig = gcf;
options.x_axis = f;
options.handle = figure(fig);
options.error = 'sem';
options.color_area = 'k'; 
options.color_line = 'k'; 
options.alpha = 0.5;
options.line_width = 2;
plot_areaerrorbar(squeeze(powerG_peri(:,1,:)),options)
hold on
options.color_area = 'r'; 
options.color_line = 'r'; 
hold on
plot_areaerrorbar(squeeze(powerG_peri(:,2,:)),options)
options.color_area = 'g'; 
options.color_line = 'g'; 
hold on
plot_areaerrorbar(squeeze(powerG_peri(:,3,:)),options)
options.color_area = 'b'; 
options.color_line = 'b'; 
hold on
plot_areaerrorbar(squeeze(powerG_peri(:,4,:)),options)
options.color_area = 'm'; 
options.color_line = 'm'; 
hold on
plot_areaerrorbar(squeeze(powerG_peri(:,5,:)),options)
xlim([0.02 1])
% ylim([0 4])
ylabel('Power: GCaMP')
xlabel('Frequency (Hz)')
set(gca,'FontSize',16)

%% Peak oscillation

powerHb_peri = powerHb_peri(:,:,1:N/2);
lowHz = find(f == 0.15);
highHz = find(f == 0.4);
for m = 1:length(mouse)
    for t = 1:3
        meanPowerHb_peri(m,t) = mean(powerHb_peri(m,t,lowHz:highHz));
    end
end
figure
bar(meanPowerHb_peri)
normMeanPowerHb_peri = meanPowerHb_peri./meanPowerHb_peri(:,1);
figure
subplot(1,2,1)
bar(1:3,mean(meanPowerHb_peri,1))
hold on
er = errorbar(1:3,mean(meanPowerHb_peri),std(meanPowerHb_peri)./sqrt(10),std(meanPowerHb_peri)./sqrt(10));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
subplot(1,2,2)
bar(1:3,mean(normMeanPowerHb_peri,1))
hold on
er = errorbar(1:3,mean(normMeanPowerHb_peri),std(normMeanPowerHb_peri)./sqrt(10),std(normMeanPowerHb_peri)./sqrt(10));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

powerG_peri = powerG_peri(:,:,1:N/2);
lowHz = find(f == 0.15);
highHz = find(f == 0.4);
for m = 1:length(mouse)
    for t = 1:3
        meanPowerG_peri(m,t) = mean(powerG_peri(m,t,lowHz:highHz));
    end
end
figure
bar(meanPowerG_peri)
normMeanPowerG_peri = meanPowerG_peri./meanPowerG_peri(:,1);
figure
subplot(1,2,1)
bar(1:3,mean(meanPowerG_peri,1))
hold on
er = errorbar(1:3,mean(meanPowerG_peri),std(meanPowerG_peri)./sqrt(10),std(meanPowerG_peri)./sqrt(10));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
subplot(1,2,2)
bar(1:3,mean(normMeanPowerG_peri,1))
hold on
er = errorbar(1:3,mean(normMeanPowerG_peri),std(normMeanPowerG_peri)./sqrt(10),std(normMeanPowerG_peri)./sqrt(10));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';



%%

powerHb_contra = powerHb_contra(:,:,1:N/2);
lowHz = find(f == 0.15);
highHz = find(f == 0.4);
for m = 1:length(mouse)
    for t = 1:3
        meanPowerHb_contra(m,t) = mean(powerHb_contra(m,t,lowHz:highHz));
    end
end
figure
bar(meanPowerHb_contra)
normMeanPowerHb_contra = meanPowerHb_contra./meanPowerHb_contra(:,1);
figure
subplot(1,2,1)
bar(1:3,mean(meanPowerHb_contra,1))
hold on
er = errorbar(1:3,mean(meanPowerHb_contra),std(meanPowerHb_contra)./sqrt(10),std(meanPowerHb_contra)./sqrt(10));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
subplot(1,2,2)
bar(1:3,mean(normMeanPowerHb_contra,1))
hold on
er = errorbar(1:3,mean(normMeanPowerHb_contra),std(normMeanPowerHb_contra)./sqrt(10),std(normMeanPowerHb_contra)./sqrt(10));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

powerG_contra = powerG_contra(:,:,1:N/2);
lowHz = find(f == 0.15);
highHz = find(f == 0.4);
for m = 1:length(mouse)
    for t = 1:3
        meanPowerG_contra(m,t) = mean(powerG_contra(m,t,lowHz:highHz));
    end
end
figure
bar(meanPowerG_contra)
normMeanPowerG_contra = meanPowerG_contra./meanPowerG_contra(:,1);
figure
subplot(1,2,1)
bar(1:3,mean(meanPowerG_contra,1))
hold on
er = errorbar(1:3,mean(meanPowerG_contra),std(meanPowerG_contra)./sqrt(10),std(meanPowerG_contra)./sqrt(10));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
subplot(1,2,2)
bar(1:3,mean(normMeanPowerG_contra,1))
hold on
er = errorbar(1:3,mean(normMeanPowerG_contra),std(normMeanPowerG_contra)./sqrt(10),std(normMeanPowerG_contra)./sqrt(10));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';


% figure
% hold on
% for m = 1:10
%     for t = 1:2
%         if t == 1
%         f = (0:n/2-1)*(fs/n);     % frequency range
%         semilogx(f,squeeze((powerHb(m,t,1:n/2))),'k-')
%         ylabel('Power: HbT')
%         ylim([0 10e-8])
%         xlim([0.02 1])
%         elseif t == 2
%         f = (0:n/2-1)*(fs/n);     % frequency range
%         semilogx(f,squeeze((powerHb(m,t,1:n/2))),'r-')
%         ylabel('Power: HbT')
%         ylim([0 10e-8])
%         xlim([0.02 1])
%         end
%     end
% end
% xlabel('Frequency (Hz)')
% legend('Baseline','Day2')
% set(gca,'XScale','log','FontSize',16)

%%
bSS75 = -88;
bSS76 = 5;
bSS77 = -30;
bSS78 = -10;
bSS79 = -23;
bSS80 = -59;
bSS81 = -6;
bSS82 = -95;
bSS83 = -13;
bSS84 = -4;
bSS85 = -21;
bSS93 = -10;
behavior = [bSS75 bSS76 bSS77 bSS78 bSS79 bSS80 bSS81 bSS82 bSS83 bSS84 bSS85 bSS93];

figure
subplot(3,2,1)
plot(meanPowerHb_peri(:,1),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerHb_peri(:,1),behavior);
title(['Baseline, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
subplot(3,2,2)
plot(meanPowerHb_peri(:,2),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerHb_peri(:,2),behavior);
title(['Day2, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
subplot(3,2,3)
plot(meanPowerHb_peri(:,3),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerHb_peri(:,3),behavior);
title(['Week1, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
% subplot(3,2,4)
% plot(meanPowerHb_peri(:,4),behavior, 'k.','MarkerSize',20)
% [r, p] = corrcoef(meanPowerHb_peri(:,4),behavior);
% title(['Week2, r = ', num2str(r(2)),', p = ', num2str(p(2))])
% ylim([-100 20])
% % xlim([0 1])
% subplot(3,2,5)
% plot(meanPowerHb_peri(:,5),behavior, 'k.','MarkerSize',20)
% [r, p] = corrcoef(meanPowerHb_peri(:,5),behavior);
% title(['Week4, r = ', num2str(r(2)),', p = ', num2str(p(2))])
% ylim([-100 20])
% % xlim([0 1])

figure
subplot(3,2,1)
plot(meanPowerG_peri(:,1),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerG_peri(:,1),behavior);
title(['Baseline, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
subplot(3,2,2)
plot(meanPowerG_peri(:,2),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerG_peri(:,2),behavior);
title(['Day2, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
subplot(3,2,3)
plot(meanPowerG_peri(:,3),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerG_peri(:,3),behavior);
title(['Week1, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
% subplot(3,2,4)
% plot(meanPowerG_peri(:,4),behavior, 'k.','MarkerSize',20)
% [r, p] = corrcoef(meanPowerG_peri(:,4),behavior);
% title(['Week2, r = ', num2str(r(2)),', p = ', num2str(p(2))])
% ylim([-100 20])
% % xlim([0 1])
% subplot(3,2,5)
% plot(meanPowerG_peri(:,5),behavior, 'k.','MarkerSize',20)
% [r, p] = corrcoef(meanPowerG_peri(:,5),behavior);
% title(['Week4, r = ', num2str(r(2)),', p = ', num2str(p(2))])
% ylim([-100 20])
% % xlim([0 1])


figure
subplot(3,2,1)
plot(meanPowerHb_contra(:,1),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerHb_contra(:,1),behavior);
title(['Baseline, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
subplot(3,2,2)
plot(meanPowerHb_contra(:,2),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerHb_contra(:,2),behavior);
title(['Day2, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
subplot(3,2,3)
plot(meanPowerHb_contra(:,3),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerHb_contra(:,3),behavior);
title(['Week1, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
% subplot(3,2,4)
% plot(meanPowerHb_contra(:,4),behavior, 'k.','MarkerSize',20)
% [r, p] = corrcoef(meanPowerHb_contra(:,4),behavior);
% title(['Week2, r = ', num2str(r(2)),', p = ', num2str(p(2))])
% ylim([-100 20])
% % xlim([0 1])
% subplot(3,2,5)
% plot(meanPowerHb_contra(:,5),behavior, 'k.','MarkerSize',20)
% [r, p] = corrcoef(meanPowerHb_contra(:,5),behavior);
% title(['Week4,r = ', num2str(r(2)),', p = ', num2str(p(2))])
% ylim([-100 20])
% % xlim([0 1])

figure
subplot(3,2,1)
plot(meanPowerG_contra(:,1),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerG_contra(:,1),behavior);
title(['Baseline, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
subplot(3,2,2)
plot(meanPowerG_contra(:,2),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerG_contra(:,2),behavior);
title(['Day2, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
subplot(3,2,3)
plot(meanPowerG_contra(:,3),behavior, 'k.','MarkerSize',20)
[r, p] = corrcoef(meanPowerG_contra(:,3),behavior);
title(['Week1, r = ', num2str(r(2)),', p = ', num2str(p(2))])
ylim([-100 20])
% xlim([0 1])
% subplot(3,2,4)
% plot(meanPowerG_contra(:,4),behavior, 'k.','MarkerSize',20)
% [r, p] = corrcoef(meanPowerG_contra(:,4),behavior);
% title(['Week2, r = ', num2str(r(2)),', p = ', num2str(p(2))])
% ylim([-100 20])
% % xlim([0 1])
% subplot(3,2,5)
% plot(meanPowerG_contra(:,5),behavior, 'k.','MarkerSize',20)
% [r, p] = corrcoef(meanPowerG_contra(:,5),behavior);
% title(['Week4, r = ', num2str(r(2)),', p = ', num2str(p(2))])
% ylim([-100 20])
% % xlim([0 1])

