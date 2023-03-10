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

pathname{1} = '..MouseData/SS93/RestingState/Baseline';
pathname{2} = '..MouseData/SS93/RestingState/Day2';
pathname{3} = '..MouseData/SS93/RestingState/Week1';
pathname{4} = '..MouseData/SS93/RestingState/Week2';
pathname{5} = '..MouseData/SS93/RestingState/Week4';

for p = 1:length(pathname)
clearvars -except pathname p
    p
% Read the Data (images)
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

data_stack = zeros(Ly,Lx,totNumFrames);
for i = 1:totNumFrames
    waitbar(i/totNumFrames)
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

dataIOSI(1,:,:,:) = temporalDetrend(led_530);
dataIOSI(2,:,:,:) = temporalDetrend(led_625);
dataGCaMP(1,:,:,:) = temporalDetrend(led_470);

dataIOSI(1,:,:,:) = squeeze(dataIOSI(1,:,:,:)) .* mask;
dataIOSI(2,:,:,:) = squeeze(dataIOSI(2,:,:,:)) .* mask;
dataGCaMP(1,:,:,:) = squeeze(dataGCaMP(1,:,:,:)) .* mask;

save([dataDir,'/','act_LED.mat'],'baseFrames','totalTrialFrames','numTrials','numFrame','dataIOSI','dataGCaMP','-v7.3');

%% IOSI pre-processing

lambdaIOS = [530 625];
numWave = size(dataIOSI,1);

H = fspecial('gaussian',3,1.5); %For spacial filter of 5x5 pixels with 1.3 pixel standard deviation
opticalDensity = zeros([numWave numFrame Ly Lx]);
PL  = pathlengths(lambdaIOS);
for i = 1:1:numWave %iterate through each wavelength
    waitbar(single(i)/single(numWave))
    intensity0 = zeros(Ly,Lx,numTrials);
    for t = 1:numTrials
        intensity0(:,:,t) = mean(squeeze(dataIOSI(i,:,:,1+(totalTrialFrames*(t-1)):baseFrames-5+(totalTrialFrames*(t-1)))),3);        
    end
    intensity0 = mean(intensity0,3);
    for f = 1:numFrame
        opticalDensity_unfilt = -log(double(squeeze(dataIOSI(i,:,:,f)))./intensity0)./PL(i);
        foo = conv2(opticalDensity_unfilt, H, 'same'); % spatially smooth optical density
        opticalDensity(i,f,:,:) = reshape(foo, [1,1,Ly,Lx]);
    end
end

% Calculate changes in oxy and deoxy hemoglobin 
e = GetExtinctions(lambdaIOS);
opticalDensity = reshape(opticalDensity, [numWave,numFrame,Ly*Lx]);
Hb = zeros(numFrame,Ly*Lx,2);
for f = 1:numFrame
    waitbar(single(f)/single(numFrame))
    for i = 1:Ly*Lx
        Hb(f,i,:) = (e\[squeeze(opticalDensity(1,f,i)); squeeze(opticalDensity(2,f,i))])';
    end
end

Hb = reshape(Hb, [numFrame,Ly,Lx,2]);
opticalDensity = reshape(opticalDensity, [numWave,numFrame,Ly,Lx]);
HbO = Hb(:,:,:,1);
HbR = Hb(:,:,:,2);
HbT = HbO + HbR;

opticalDensity = permute(opticalDensity, [1,3,4,2]);
HbO = permute(HbO, [2,3,1]);
HbR = permute(HbR, [2,3,1]);
HbT = permute(HbT, [2,3,1]);

save([dataDir,'/','act_IOSI_ipsi.mat'],'baseFrames','totalTrialFrames','numTrials','numFrame','HbO','HbR','HbT','-v7.3');

%% GCaMP pre-processing

dataGCaMP = squeeze(dataGCaMP(1,:,:,:));
GCaMPfull = dataGCaMP; %temporalDetrendAdam(dataGCaMP);
H = fspecial('gaussian',3,1.5); %For spacial filter of 5x5 pixels with 1.3 pixel standard deviation
intensity0 = zeros(Ly,Lx,numTrials);
for t = 1:numTrials
    intensity0(:,:,t) = mean(squeeze(GCaMPfull(:,:,1+(totalTrialFrames*(t-1)):baseFrames-5+(totalTrialFrames*(t-1)))),3);        
end
intensity0 = mean(intensity0,3);
deltaGCaMP = zeros(Ly,Lx,numFrame);
for f = 1:numFrame
    deltaGCaMP_unfilt = log(double(GCaMPfull(:,:,f))./intensity0);
    deltaGCaMP(:,:,f) = conv2(deltaGCaMP_unfilt, H, 'same'); % spatially smooth optical density
end

save([dataDir,'/','act_GCaMP.mat'],'baseFrames','totalTrialFrames','numTrials','numFrame','deltaGCaMP','-v7.3');


%% GCaMP pre-processing (correction method 1)

clearvars -except dataDir pathname p maskSFDI
load([dataDir,'/act_IOSI_ipsi.mat'])
load([dataDir,'/act_GCaMP.mat'])
load([dataDir,'/act_LED.mat'])
s = regexp(dataDir,'/','split');
load([s{1},'/',s{2},'/',s{3},'/',s{4},'/','PL_grid.mat'])

% GCaMP correction for hemodynamic cross-talk
lambda = [473; 530];
e = GetExtinctions(lambda);
mua_ex = (e(1,1).*HbO) + (e(1,2).*HbR);
mua_em = (e(2,1).*HbO) + (e(2,2).*HbR);
temp = ((mua_ex.*XexScl_grid) + (mua_em.*Xem_grid)) + deltaGCaMP;
deltaGCaMPcorr = exp(temp)-1;
save([dataDir,'/','act_GCaMPcorr_ipsi.mat'],'baseFrames','totalTrialFrames','numTrials','numFrame',...
    'deltaGCaMPcorr','Xex_grid','Xem_grid','-v7.3');

end