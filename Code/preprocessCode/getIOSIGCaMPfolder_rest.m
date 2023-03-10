%% Simultaneous hemodynamic and calcium imaging
% Pre-process code for hemodynamics and GCaMP
% Hemodynamics wavelengths: 530nm and 625nm
% GCaMP excitation: 470nm, emission: 525nm
%
% getIOSIGCaMPFolder - Reads raw tif files recorded by HCImage (Hamamatsu
% software) from the selected folder and calculates optical density changes
% and oxy, deoxy, and total hemoglobin changes.
% 
% Other m-files required: pathlengths.m, GetExtinctions.m
%
% Data is saved as a mat file within the selected raw images folder. 
%
% Smrithi Sunil
% email: ssunil@bu.edu
% BOAS Lab, Boston University

%% Load data and crop
clear

pathname{1} = '..MouseData/SS93/RestingState/Baseline';
pathname{2} = '..MouseData/SS93/RestingState/Day2';
pathname{3} = '..MouseData/SS93/RestingState/Week1';
pathname{4} = '..MouseData/SS93/RestingState/Week2';
pathname{5} = '..MouseData/SS93/RestingState/Week4';

for p = 1:length(pathname)
clearvars -except pathname p

% Read the Data (images)
wavelength = 470; %input('What wavelength was was imaged first? ');
freq = 15;
if wavelength == 470
    lambda = [470 530 625];
elseif wavelength == 530
    lambda = [530 625 470];
elseif wavelength == 625
    lambda = [625 470 530];
else
    error('Enter a valid wavelength')
end

dataDir = pathname{p}; %uigetdir('Please select the Data folder');
dataFiles = dir([dataDir '/*.tif']);
I = double(imread([dataDir '/' dataFiles(2).name]));
I = imresize(I,0.25);
Ly = size(I,1);
Lx = size(I,2);


numFrames = 9000;
data_stack = zeros(Ly,Lx,numFrames);
for i = 1:numFrames
    waitbar(i/numFrames)
    temp_tiff = double(imread([dataDir '/' dataFiles(i).name]));
    temp_tiff = imresize(temp_tiff,0.25);
    data_stack(:,:,i) = temp_tiff;
end

s = regexp(dataDir,'/','split');
load([s{1},'/',s{2},'/','brainmaskSFDI.mat'])
mask = maskSFDI(p).aff_mask + maskSFDI(p).unaff_mask;

if lambda(1) == 470
    led_470 = data_stack(:,:,1:3:end);
    led_530 = data_stack(:,:,2:3:end);
    led_625 = data_stack(:,:,3:3:end);
elseif lambda(1) == 530
    led_530 = data_stack(:,:,1:3:end);
    led_625 = data_stack(:,:,2:3:end);
    led_470 = data_stack(:,:,3:3:end);
elseif lambda(1) == 625
    led_625 = data_stack(:,:,1:3:end);
    led_470 = data_stack(:,:,2:3:end);
    led_530 = data_stack(:,:,3:3:end);
end
led_470 = led_470(:,:,301:2700);
led_530 = led_530(:,:,301:2700);
led_625 = led_625(:,:,301:2700);

dataIOSI(1,:,:,:) = temporalDetrend(led_530);
dataIOSI(2,:,:,:) = temporalDetrend(led_625);
dataGCaMP(1,:,:,:) = temporalDetrend(led_470);

dataIOSI(1,:,:,:) = squeeze(dataIOSI(1,:,:,:)) .* mask;
dataIOSI(2,:,:,:) = squeeze(dataIOSI(2,:,:,:)) .* mask;
dataGCaMP(1,:,:,:) = squeeze(dataGCaMP(1,:,:,:)) .* mask;

save([dataDir,'/','rest_LED_v2.mat'],'dataIOSI','dataGCaMP','numFrames','dataDir','freq','lambda','-v7.3');


%% IOSI pre-processing
numFrame = size(dataIOSI,4);
lambdaIOSI = [530 625];
numWave = size(dataIOSI,1);
Lx = size(dataIOSI,2);
Ly = size(dataIOSI,3);

H = fspecial('gaussian',3,1.5); %For spacial filter of 5x5 pixels with 1.3 pixel standard deviation
opticalDensity = zeros([numWave numFrame Ly Lx]);
PL  = pathlengths(lambdaIOSI);

for i = 1:1:numWave %iterate through each wavelength
    waitbar(single(i)/single(numWave))
    intensity0 = mean(squeeze(dataIOSI(i,:,:,:)),3);
    for f = 1:numFrame
        opticalDensity_unfilt = -log(double(squeeze(dataIOSI(i,:,:,f)))./intensity0)./PL(i);
        foo = conv2(opticalDensity_unfilt, H, 'same'); % spatially smooth optical density
        opticalDensity(i,f,:,:) = reshape(foo, [1,1,Ly,Lx]);
    end
end

% Calculate changes in oxy and deoxy hemoglobin 
e = GetExtinctions(lambdaIOSI);
opticalDensity = reshape(opticalDensity, [numWave,numFrame,Ly*Lx]);
Hb = zeros(numFrame,Ly*Lx,3);
for f = 1:numFrame
    waitbar(single(f)/single(numFrame))
    for i = 1:Ly*Lx
        Hb(f,i,1:2) = (e\[squeeze(opticalDensity(1,f,i)); squeeze(opticalDensity(2,f,i))])';
    end
end
Hb(:,:,3) = Hb(:,:,1)+Hb(:,:,2);
Hb = reshape(Hb, [numFrame,Ly,Lx,3]);
Hb = permute(Hb, [2,3,4,1]);
Hb(isnan(Hb)) = 0;
opticalDensity = reshape(opticalDensity, [numWave,numFrame,Ly,Lx]);

numFrames = numFrame*3;
exp_times = zeros(numFrames,2);
exp_times(:,1) = 0:1/freq:(numFrames-1)/freq;
for u = 1:length(lambda)
     idx = u:3:numFrames;
     exp_times(idx,2) = lambda(u);
end
tRS = exp_times(:,1);
tRS = tRS(1:3:end);
nT = length(tRS);

fs = 5;
lpf = 0.08; % 2; %0.08;
hpf = 0.009; %0.02; %0.009;

wn(2) = lpf/(fs/2);
wn(1) = hpf/(fs/2);


for u = 1:3
    u
    y = reshape(Hb(:,:,u,:),[Lx*Ly nT]);
    [fb,fa] = butter(2,wn,'bandpass');
    yfilt = filtfilt(fb,fa,y')';
    yunfilt = reshape(y,[Lx*Ly nT]);
    I=reshape(yfilt,[Lx Ly 1 nT]);
    Iunfilt=reshape(yunfilt,[Lx Ly 1 nT]);
    img = double(imread([dataDir '/' dataFiles(2).name]));
    I0 = img(1:2:end,1:2:end);
    I0 = imresize(I0,0.5,'bilinear');
    I0_unsampled = I0;
    save([dataDir,'/',num2str(u) '_v2_rsfc_I.mat'],'I','Iunfilt','tRS','I0', 'I0_unsampled', '-v7.3');
end


%% GCaMP pre-processing (correction method 1)

GCaMPfull = squeeze(dataGCaMP(1,:,:,:));
numFrame = size(dataGCaMP,4);
H = fspecial('gaussian',3 ,1.5); %For spacial filter of 5x5 pixels with 1.3 pixel standard deviation
intensity0 = mean(GCaMPfull,3);
deltaGCaMP = zeros(Ly,Lx,numFrame);
for f = 1:numFrame
    temp = log(double(GCaMPfull(:,:,f))./intensity0);
    deltaGCaMP(:,:,f) = conv2(temp, H, 'same'); % spatially smooth optical density
end
deltaGCaMP(isnan(deltaGCaMP))=0;
numFrames = numFrame*3;
exp_times = zeros(numFrames,2);
exp_times(:,1) = 0:1/freq:(numFrames-1)/freq;
for u = 1:length(lambda)
     idx = u:3:numFrames;
     exp_times(idx,2) = lambda(u);
end
tRS = exp_times(:,1);
tRS = tRS(1:3:end);
nT = length(tRS);
y = reshape(deltaGCaMP(:,:,:),[Lx*Ly nT]);
fs = 5;
lpf = 0.08; % 2; %0.08;
hpf = 0.009; %0.02; %0.009;
wn(2) = lpf/(fs/2);
wn(1) = hpf/(fs/2);
[fb,fa] = butter(2,wn,'bandpass');
yfilt = filtfilt(fb,fa,y')';
yunfilt = reshape(y,[Lx*Ly nT]);
I=reshape(yfilt,[Lx Ly 1 nT]);
Iunfilt=reshape(yunfilt,[Lx Ly 1 nT]);
img = double(imread([dataDir '/' dataFiles(2).name]));
I0 = img(1:2:end,1:2:end);
I0 = imresize(I0,0.5,'bilinear');
I0_unsampled = I0;
save([dataDir,'/','gcamp_v2_rsfc_I.mat'],'I','Iunfilt','tRS','I0', 'I0_unsampled', '-v7.3');

clearvars -except dataDir pathname p mask
HbO = load([dataDir '/1_v2_rsfc_I.mat']);
HbR = load([dataDir '/2_v2_rsfc_I.mat']);
GCaMP = load([dataDir '/gcamp_v2_rsfc_I.mat']);
lambda = [473; 530];
e = GetExtinctions(lambda);
Xex = 0.56/10; % Change to centimeters
Xem = 0.57/10;
HbO.I = permute(HbO.I,[1 2 4 3]);
HbO.Iunfilt = permute(HbO.Iunfilt,[1 2 4 3]);
HbR.I = permute(HbR.I,[1 2 4 3]);
HbR.Iunfilt = permute(HbR.Iunfilt,[1 2 4 3]);
GCaMP.I = permute(GCaMP.I,[1 2 4 3]);
GCaMP.Iunfilt = permute(GCaMP.Iunfilt,[1 2 4 3]);
mua_ex = (e(1,1).*HbO.I) + (e(1,2).*HbR.I);
mua_em = (e(2,1).*HbO.I) + (e(2,2).*HbR.I);
temp = (mua_ex.*Xex) + (mua_em.*Xem) + GCaMP.I;
GCaMPcorr.I = exp(temp)-1;
mua_ex = (e(1,1).*HbO.Iunfilt) + (e(1,2).*HbR.Iunfilt);
mua_em = (e(2,1).*HbO.Iunfilt) + (e(2,2).*HbR.Iunfilt);
temp = (mua_ex.*Xex) + (mua_em.*Xem) + GCaMP.Iunfilt;
GCaMPcorr.Iunfilt = exp(temp)-1;
I = GCaMPcorr.I;
I = permute(I, [1 2 4 3]);
Iunfilt = GCaMPcorr.Iunfilt;
Iunfilt = permute(Iunfilt, [1 2 4 3]);
tRS = GCaMP.tRS;
I0 = GCaMP.I0;
I0_unsampled = GCaMP.I0_unsampled;
save([dataDir,'/','gcampcorr_v2_rsfc_I.mat'],'I','Iunfilt','tRS','I0', 'I0_unsampled', '-v7.3');

end
