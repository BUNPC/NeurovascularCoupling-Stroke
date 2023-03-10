%% Register to mouse brain atlas

%% Registration to atlas for baseline
clear
% close all
dataDir = uigetdir('Please select the Data folder');
sub = dir(dataDir);
folders = sub([sub(:).isdir]);
folderIOSI = regexp({folders(:).name},'IOSI','match');
dfolders = folders(find(~cellfun('isempty',folderIOSI)));
baseFrames = 25;
rStart = baseFrames;
rEnd = 2*baseFrames;
prompt1 = sprintf('Which is the ipsi-lesional hemisphere?');
button1 = questdlg(prompt1,'Select region','Right','Left','No');
if strcmpi(button1,'Right')
    atlas_coordinates = [249.7575  420.0912; 249.7639  235.8948; 310 225];
else
    atlas_coordinates = [249.7639  235.8948;  249.7575  420.0912; 310 432];
end
atlas = imread('miceatlas.png');

fixed = double(atlas(:,:,1));
figure(1)
colormap jet
subplot(2,2,1) 
imagesc(fixed)
axis image
caxis([prctile(fixed(:),5), prctile(fixed(:),95)]); 
for i = 1:3
filename = 'act_trialavg.mat';
pathname = [dataDir,'\',dfolders(i).name];
[~,~,ext] = fileparts(filename);
data(i) = load([pathname,'\', filename]);
I0 = data(i).GCaMPcorr;
figure(1)
subplot(2,2,1)
hold on
plot(atlas_coordinates(i,2),atlas_coordinates(i,1),'ko','MarkerFaceColor','k','MarkerSize',15)

moving = double(mean(I0(:,:,rStart:rEnd),3));     
subplot(2,2,i+1)
imagesc(moving)
caxis([prctile(moving(:),2), prctile(moving(:),98)]); 
axis image
colormap jet
[x2(i), y2(i)] = ginput(1);
hold on
plot(x2(i), y2(i),'ko','MarkerFaceColor','k','MarkerSize',15)

end
saveas(gcf,[dataDir,'\','atlas_ref','.tif'])
[D,Z,TRANSFORM] = procrustes([y2' x2'], atlas_coordinates);

figure
for i = 1:3
filename = 'act_trialavg.mat';
pathname = [dataDir,'\',dfolders(i).name];
[~,~,ext] = fileparts(filename);
data(i) = load([pathname,'\', filename]);
I0 = data(i).GCaMPcorr;
moving = double(mean(I0(:,:,rStart:rEnd),3));     

movingImg = zeros(size(fixed,1),size(fixed,2));
for u = 1:size(fixed,1)
    for v = 1:size(fixed,2)
        idx = round(TRANSFORM.b * [u v] * TRANSFORM.T + TRANSFORM.c(1,:));
        if idx(1) < 1 || idx(1) > size(moving,1) || idx(2) < 1 || idx(2) > size(moving,2)
            movingImg(u,v) = 0;
        else
            movingImg(u,v) = moving(idx(1),idx(2));
        end
    end
end
movingImg(movingImg==0) = NaN;

regImg = movingImg;
regImg(regImg<=(nanstd(nanstd(movingImg))*7)) = NaN;
subplot(1,3,i)
imshowpair(fixed,regImg)
end

save([dataDir,'\','tform_atlas.mat'],'D','Z','TRANSFORM','fixed','-v7.3');

%% Register for time points after baseline

clear
% close all
dataDir = uigetdir('Please select the Data folder');
sub = dir(dataDir);
folders = sub([sub(:).isdir]);
folderIOSI = regexp({folders(:).name},'IOSI','match');
dfolders = folders(find(~cellfun('isempty',folderIOSI)));
prompt1 = sprintf('Which is the ipsi-lesional hemisphere?');
button1 = questdlg(prompt1,'Select region','Right','Left','No');
if strcmpi(button1,'Right')
    atlas_coordinates = [249.7575  420.0912; 249.7639  235.8948; 310 225];
else
    atlas_coordinates = [249.7639  235.8948;  249.7575  420.0912; 310 432];
end
atlas = imread('miceatlas.png');

fixed = double(atlas(:,:,1));
figure(1)
colormap jet
subplot(1,2,1) 
imagesc(fixed)
axis image
caxis([prctile(fixed(:),5), prctile(fixed(:),95)]); 

for i = 1:3
pathname = [dataDir,'\',dfolders(i).name];
dataFiles = dir([pathname '\*.tif']);
data(i).img = double(imread([pathname '\' dataFiles(2).name]));
data(i).img = imresize(data(i).img,0.25);
I0 = data(i).img;

figure(1)
subplot(1,2,1)
hold on
plot(atlas_coordinates(i,2),atlas_coordinates(i,1),'ko','MarkerFaceColor','k','MarkerSize',15)

moving = data(i).img;
subplot(1,2,2)
imagesc(moving)
caxis([prctile(moving(:),2), prctile(moving(:),98)]); 
axis image
colormap jet
[x2(i), y2(i)] = ginput(1);
end

[D,Z,TRANSFORM] = procrustes([y2' x2'], atlas_coordinates);

figure
for i = 1:3
pathname = [dataDir,'\',dfolders(i).name];
dataFiles = dir([pathname '\*.tif']);
data(i).img = double(imread([pathname '\' dataFiles(2).name]));
data(i).img = imresize(data(i).img,0.25);
I0 = data(i).img;   

movingImg = zeros(size(fixed,1),size(fixed,2));
for u = 1:size(fixed,1)
    for v = 1:size(fixed,2)
        idx = round(TRANSFORM.b * [u v] * TRANSFORM.T + TRANSFORM.c(1,:));
        if idx(1) < 1 || idx(1) > size(moving,1) || idx(2) < 1 || idx(2) > size(moving,2)
            movingImg(u,v) = 0;
        else
            movingImg(u,v) = moving(idx(1),idx(2));
        end
    end
end
movingImg(movingImg==0) = NaN;

regImg = movingImg;
regImg(regImg<=(nanstd(nanstd(movingImg))*2)) = NaN;
subplot(1,3,i)
imshowpair(fixed,regImg)
end

save([dataDir,'\','tform_atlas.mat'],'D','Z','TRANSFORM','fixed','-v7.3');

%% Transform rest data

clear
[filename,pathname] = uigetfile({'*.mat;*.tiff;*.tif'},'Please select transformation');
load([pathname filename]);

for f = 1:2
[filename,pathname] = uigetfile({'*.mat;*.tiff;*.tif'},'Please select resting state data');
load([pathname filename]);

moving = LRC;     
movingImg = zeros(size(fixed,1),size(fixed,2));
for u = 1:size(fixed,1)
    for v = 1:size(fixed,2)
        idx = round(TRANSFORM.b * [u v] * TRANSFORM.T + TRANSFORM.c(1,:));
        if idx(1) < 1 || idx(1) > size(moving,1) || idx(2) < 1 || idx(2) > size(moving,2)
            movingImg(u,v) = 0;
        else
            movingImg(u,v) = moving(idx(1),idx(2));
        end
    end
end
movingImg(movingImg==0) = NaN;
figure
imshowpair(fixed,movingImg)
regLRC = movingImg;
str_1 = filename(1:strfind(filename,'.')-1);
save([pathname,'\',str_1, '_atlas.mat'],'regLRC','D','Z','TRANSFORM','fixed','-v7.3');
end

%% Select forelimb area and plot correlation 

figure
imagesc(fixed)
axis image
colormap jet
axis off
[mask_forepaw,Xi,Yi] = roipoly;
hold on;
plot(Xi,Yi,'color','k');
hold off;

%%

baseline = 'E:\SS77\RestingState\Baseline\RestingState_IOSI\1_rsfc_I_LRC_atlas.mat';
day2 = 'E:\SS77\RestingState\Day2\RestingState_IOSI\1_rsfc_I_LRC_atlas.mat';
week1 = 'E:\SS77\RestingState\Week1\RestingState_IOSI\1_rsfc_I_LRC_atlas.mat';
week2 = 'E:\SS77\RestingState\Week2\RestingState_IOSI\1_rsfc_I_LRC_atlas.mat';
week4 = 'E:\SS77\RestingState\Week4\RestingState_IOSI\1_rsfc_I_LRC_atlas.mat';
load mask_forepaw
load mask_hindpaw
load mask_whisker
load mask_motor
timepoint = {baseline; day2; week1; week2; week4};
count=1;
figure
for t = 1:5
    load(timepoint{t})
    
    forepaw_r = mask_forepaw.*regLRC;
    forepaw_r(forepaw_r==0) = NaN;
    subplot(5,4,count)
    imagesc(forepaw_r)
    axis image
    colormap jet
    caxis([-1 1])
    axis off
    mforepaw(t) = nanmean(nanmean(forepaw_r));
    stdforepaw(t) = nanstd(nanstd(forepaw_r));
    count = count+1;
    
    hindpaw_r = mask_hindpaw.*regLRC;
    hindpaw_r(hindpaw_r==0) = NaN;
    subplot(5,4,count)
    imagesc(hindpaw_r)
    axis image
    colormap jet
    caxis([-1 1])
    axis off
    mhindpaw(t) = nanmean(nanmean(hindpaw_r));
    stdhindpaw(t) = nanstd(nanstd(hindpaw_r));
    count = count+1;
    
    whisker_r = mask_whisker.*regLRC;
    whisker_r(whisker_r==0) = NaN;
    subplot(5,4,count)
    imagesc(whisker_r)
    axis image
    colormap jet
    caxis([-1 1])
    axis off
    mwhisker(t) = nanmean(nanmean(whisker_r));
    stdwhisker(t) = nanstd(nanstd(whisker_r));
    count = count+1;

    motor_r = mask_motor.*regLRC;
    motor_r(motor_r==0) = NaN;
    subplot(5,4,count)
    imagesc(motor_r)
    axis image
    colormap jet
    caxis([-1 1])
    axis off
    mmotor(t) = nanmean(nanmean(motor_r));
    stdmotor(t) = nanstd(nanstd(motor_r));
    count = count+1;
end

atlas_region = {'Forelimb','Hindlimb','Whisker','Motor'};
mean_region1 = [mforepaw; mhindpaw; mwhisker; mmotor];
error = [stdforepaw; stdhindpaw; stdwhisker; stdmotor];
figure
x=1:4;
bar(x,mean_region1)
hold on
ngroups = 4;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mean_region1(:,e),error(:,e),error(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
hold off
ylim([-0.3 1])
set(gca, 'XTickLabel', atlas_region,'FontSize',16);
legend('Pre-Stroke','Day2','Week1','Week2','Week4')

%
baseline = 'E:\SS77\RestingState\Baseline\RestingState_IOSI\gcampcorr2_rsfc_I_LRC_atlas.mat';
day2 = 'E:\SS77\RestingState\Day2\RestingState_IOSI\gcampcorr2_rsfc_I_LRC_atlas.mat';
week1 = 'E:\SS77\RestingState\Week1\RestingState_IOSI\gcampcorr2_rsfc_I_LRC_atlas.mat';
week2 = 'E:\SS77\RestingState\Week2\RestingState_IOSI\gcampcorr2_rsfc_I_LRC_atlas.mat';
week4 = 'E:\SS77\RestingState\Week4\RestingState_IOSI\gcampcorr2_rsfc_I_LRC_atlas.mat';
load mask_forepaw
load mask_hindpaw
load mask_whisker
load mask_motor
timepoint = {baseline; day2; week1; week2; week4};

for t = 1:5
    load(timepoint{t})
    
    forepaw_r = mask_forepaw.*regLRC;
    forepaw_r(forepaw_r==0) = NaN;
    mforepaw(t) = nanmean(nanmean(forepaw_r));
    stdforepaw(t) = nanstd(nanstd(forepaw_r));

    hindpaw_r = mask_hindpaw.*regLRC;
    hindpaw_r(hindpaw_r==0) = NaN;
    mhindpaw(t) = nanmean(nanmean(hindpaw_r));
    stdhindpaw(t) = nanstd(nanstd(hindpaw_r));

    whisker_r = mask_whisker.*regLRC;
    whisker_r(whisker_r==0) = NaN;
    mwhisker(t) = nanmean(nanmean(whisker_r));
    stdwhisker(t) = nanstd(nanstd(whisker_r));

    motor_r = mask_motor.*regLRC;
    motor_r(motor_r==0) = NaN;
    mmotor(t) = nanmean(nanmean(motor_r));
    stdmotor(t) = nanstd(nanstd(motor_r));
end

atlas_region = {'Forelimb','Hindlimb','Whisker','Motor'};
mean_region2 = [mforepaw; mhindpaw; mwhisker; mmotor];
error = [stdforepaw; stdhindpaw; stdwhisker; stdmotor];
figure
x=1:4;
bar(x,mean_region2)
hold on
ngroups = 4;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    er = errorbar(x,mean_region2(:,e),error(:,e),error(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
hold off
ylim([-0.3 1])
set(gca, 'XTickLabel', atlas_region,'FontSize',16);
legend('Pre-Stroke','Day2','Week1','Week2','Week4')


figure
x=1:4;
bar(x,mean_region1./mean_region1(:,1))
hold on
ngroups = 4;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    z = mean_region1./mean_region1(:,1);
    er = errorbar(x,z(:,e),error(:,e),error(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
set(gca, 'XTickLabel', atlas_region,'FontSize',16);
legend('Pre-Stroke','Day2','Week1','Week2','Week4')
ylim([-1 1.5])


figure
x=1:4;
bar(x,mean_region2./mean_region2(:,1))
hold on
ngroups = 4;
nbars = 5;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for e = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*e-1) * groupwidth / (2*nbars);
    z = mean_region2./mean_region2(:,1);
    er = errorbar(x,z(:,e),error(:,e),error(:,e));    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
end
set(gca, 'XTickLabel', atlas_region,'FontSize',16);
legend('Pre-Stroke','Day2','Week1','Week2','Week4')
ylim([-1 1.5])

