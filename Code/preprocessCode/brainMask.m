%% Save mask of cranial window for each mouse

pathname{1} = '..MouseData/SS93/RestingState/Baseline';
pathname{2} = '..MouseData/SS93/RestingState/Day2';
pathname{3} = '..MouseData/SS93/RestingState/Week1';
pathname{4} = '..MouseData/SS93/RestingState/Week2';
pathname{5} = '..MouseData/SS93/RestingState/Week4';

for p = 1:length(pathname)
    fh = figure;
    fh.WindowState = 'maximized';
    dataDir = pathname{p}; 
    load(dataDir)
    subplot(1,2,1)
    imagesc(prop_mua)
    axis image
    colormap jet
    axis off
    caxis([0 1])
    subplot(1,2,2)
    imagesc(prop_mus./10)
    axis image
    colormap jet
    axis off
    caxis([0 3])
    for r = 1:2
        [M,Xi,Yi] = roipoly;
        hold on
        plot(Xi,Yi,'color','k')
        if r == 1
            maskSFDI(p).aff_mask = M;
            maskSFDI(p).aff_prop_mua = prop_mua.*M;
            maskSFDI(p).aff_prop_mus = prop_mus.*M;
        elseif r == 2
            maskSFDI(p).unaff_mask = M;
            maskSFDI(p).unaff_prop_mua = prop_mua.*M;
            maskSFDI(p).unaff_prop_mus = prop_mus.*M;
        end
    end
end

s = regexp(dataDir,'\','split');
save([s{1},'\',s{2},'\','brainmaskSFDI','.mat'],'maskSFDI','-v7.3');
close all

%% Save stroke outline in brainmaskSFDI
clear
mouse = {SS84};
load([mouse{1},'brainmaskSFDI.mat'])
maskSFDI(1).stroke_mask = logical(zeros(128));
maskSFDI(2).stroke_mask = logical(zeros(128));

for t = 3:5
    figure
    imagesc(maskSFDI(t).aff_prop_mus)
    colormap jet
    axis image
    axis off
    caxis([0 25])
    [M,Xi,Yi] = roipoly;
    hold on
    plot(Xi,Yi,'color','k')
    maskSFDI(t).stroke_mask = M;
end
save([mouse{1},'brainmaskSFDI.mat'],'maskSFDI','-v7.3')

%% Save activation outlines in brain mask
% Ipsi and contra forelimb activations

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
baseFrames = 25;
rStart = baseFrames;
rEnd = 2*baseFrames;
for m = 1:length(mouse)
    m
    [timepoints, img, mask] = animals(mousename{m});
    load(mask)
    load([timepoints{1},'\act_trialAvg1_v2_2.mat'])
    GCaMPcorrresp = squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
    cutGCaMPcorrresp = maskSFDI(1).aff_mask.*GCaMPcorrresp;
    peakG = prctile(cutGCaMPcorrresp(:),99);
    cutGCaMPcorrresp(cutGCaMPcorrresp<=0.75*peakG) = 0;
    cutGCaMPcorrresp(cutGCaMPcorrresp>0.75*peakG) = 1;
    cutGCaMPcorrresp(isnan(cutGCaMPcorrresp)) = 0;
    maskSFDI(1).ipsiOutline = cutGCaMPcorrresp;

    [timepoints, img, mask] = animals_contra(mousename{m});
    load([timepoints{1},'\act_trialAvg1_v2_2.mat'])
    GCaMPcorrresp = squeeze(mean(GCaMPcorr(:,:,rStart:rEnd),3));
    cutGCaMPcorrresp = maskSFDI(1).unaff_mask.*GCaMPcorrresp;
    peakG = prctile(cutGCaMPcorrresp(:),99);
    cutGCaMPcorrresp(cutGCaMPcorrresp<=0.75*peakG) = 0;
    cutGCaMPcorrresp(cutGCaMPcorrresp>0.75*peakG) = 1;
    cutGCaMPcorrresp(isnan(cutGCaMPcorrresp)) = 0;
    maskSFDI(1).contraOutline = cutGCaMPcorrresp;
    
    save([mouse{m},'brainmaskSFDI.mat'],'maskSFDI','-v7.3')
end


%% Draw and save stroke, affected, and unaffected hemisphere outlines 

pathname{1} = '..MouseData/SS93/RestingState/Baseline';
pathname{2} = '..MouseData/SS93/RestingState/Day2';
pathname{3} = '..MouseData/SS93/RestingState/Week1';
pathname{4} = '..MouseData/SS93/RestingState/Week2';
pathname{5} = '..MouseData/SS93/RestingState/Week4';

% baseline

figure
dataDir = pathname{1}; 
dataFiles = dir([dataDir '/*.tif']);
I = double(imread([dataDir '/' dataFiles(2).name]));
I = imresize(I,0.25);
imagesc(I)
axis image
colormap jet
axis off
for r = 1:2
    [M,Xi,Yi] = roipoly;
    hold on
    plot(Xi,Yi,'color','k')
    if r == 1
        mask.affected = M;
    elseif r == 2
        mask.unaffected = M;
    end
    mask.stroke = logical(zeros(size(I)));
end
save([dataDir,'/','mask.mat'],'mask','-v7.3');

% day2
figure
dataDir = pathname{2}; 
dataFiles = dir([dataDir '/*.tif']);
I = double(imread([dataDir '/' dataFiles(4).name]));
I = imresize(I,0.25);
imagesc(I)
axis image
colormap jet
axis off
for r = 1:3
    [M,Xi,Yi] = roipoly;
    hold on
    plot(Xi,Yi,'color','k')
    if r == 1
        mask.stroke = M;
    elseif r == 2
        mask.affected = M;
    elseif r == 3
        mask.unaffected = M;
    end
end
save([dataDir,'/','mask.mat'],'mask','-v7.3');

% week1-week4
for p = 3:5
    figure
    dataDir = pathname{p}; 
    dataFiles = dir([dataDir '/*.tif']);
    I = double(imread([dataDir '/' dataFiles(2).name]));
    I = imresize(I,0.25);
    imagesc(I)
    axis image
    colormap jet
    axis off
    for r = 1:3
        [M,Xi,Yi] = roipoly;
        hold on
        plot(Xi,Yi,'color','k')
        if r == 1
            mask.stroke = M;
        elseif r == 2
            mask.affected = M;
        elseif r == 3
            mask.unaffected = M;
        end
    end
    save([dataDir,'/','mask.mat'],'mask','-v7.3');
end





