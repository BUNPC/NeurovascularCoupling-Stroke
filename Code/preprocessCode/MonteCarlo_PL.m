%% Monte Carlo simulation for pathlength 


pathname{1} = '..MouseData/SS93/RestingState/Baseline';
pathname{2} = '..MouseData/SS93/RestingState/Day2';
pathname{3} = '..MouseData/SS93/RestingState/Week1';
pathname{4} = '..MouseData/SS93/RestingState/Week2';
pathname{5} = '..MouseData/SS93/RestingState/Week4';

p=3;
dataDir = pathname{p}; %uigetdir('Please select the Data folder');
s = regexp(dataDir,'/','split');
load([s{1},'/',s{2},'/','brainmaskSFDI.mat'])
mua = maskSFDI(p).aff_prop_mua + maskSFDI(p).unaff_prop_mua;
mus = maskSFDI(p).aff_prop_mus + maskSFDI(p).unaff_prop_mus;
figure
imagesc(mus)
axis image
colormap jet
caxis([0 25])
hold on
h = imrect;
pos = wait(h);
pos = round(pos);
delete(h);

for p = 1:length(pathname)
clearvars -except pathname p pos

dataDir = pathname{p}; %uigetdir('Please select the Data folder');
s = regexp(dataDir,'/','split');
load([s{1},'/',s{2},'/','brainmaskSFDI.mat'])

mua = maskSFDI(p).aff_prop_mua + maskSFDI(p).unaff_prop_mua;
mus = maskSFDI(p).aff_prop_mus + maskSFDI(p).unaff_prop_mus;

% optical property values from SFDI baseline of all mice
mean_bmua = 0.58;
std_bmua = 0.08;
mean_bmus = 10.401;
std_bmus = 1.6034;
mus_thresh = mean_bmus + (2*std_bmus);
mask = maskSFDI(p).aff_mask + maskSFDI(p).unaff_mask;

% GCaMP correction for hemodynamic cross-talk
lambda = [473; 530];
e = GetExtinctions(lambda);
Xex = 0.54/10; % Change to centimeters
Xem = 0.84/10;

mus_crop = mus(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3));
mua_crop = mua(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3));

if p == 1
    Xex_grid = Xex.*mask.*ones(size(mus,1),size(mus,2));
    Xem_grid = Xem.*mask.*ones(size(mus,1),size(mus,2));
    idx = 0;
    idy = 0;
    cfg = 0;
else

    Xem_grid_crop = ones(size(mus_crop,1),size(mus_crop,2));
    [idx,idy] = find(mus_crop>mus_thresh & mus_crop<100);
    
    for i = 1:length(idx)
        i
        waitbar(single(i)/single(length(idx)))
        cfg.nphoton=1e6;
        cfg.vol=uint8(ones(60,60,60));
        n = 1.33;
        g = 0.9;
        cfg.prop = [0 0 1 1; mua_crop(idx(i),idy(i)) mus_crop(idx(i),idy(i)) g n];
        cfg.tstart = 0;
        cfg.tend = 5e-9;
        cfg.tstep = 1e-10;
        cfg.srcpos = [30 30 1];
        cfg.srcdir = [0 0 1];
        cfg.detpos = [30 30 1 1];
        cfg.srctype = 'disk';
        cfg.srcparam1 = [60 0 0 0];
        cfg.gpuid = 1;
        cfg.autopilot = 1;
        cfg.isreflect = 0;
        cfg.debuglevel = 'P';
        [fluxs_disk,detp_disk] = mcxlab(cfg);
        prop = cfg.prop;
        avgpath_disk = mcxmeanpath(detp_disk,prop);
        Xem_grid_crop(idx(i),idy(i)) = avgpath_disk/2/10; % Change to cm and only take one direction
    end

Xex_grid = Xex.*mask.*ones(size(mus,1),size(mus,2));
Xem_grid = mask.*ones(size(mus,1),size(mus,2));
Xem_grid(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3)) = Xem_grid_crop;
Xem_grid(Xem_grid==1) = Xem;
% Xem_grid(mus<=mus_thresh) = Xem;
% Xem_grid(mus == 100) = Xem;

end
    
    
Xem_grid = mask.*Xem_grid;
totPL = Xex_grid + Xem_grid;
[needScalex,needScaley] = find(totPL < 0.13 & totPL > 0.01);
XexScl_grid = Xex_grid;
for i = 1:length(needScalex)
    XexScl_grid(needScalex(i),needScaley(i)) = 0.138 - Xem_grid(needScalex(i),needScaley(i));
end
save([dataDir,'/','PL_grid.mat'],'Xem_grid','Xex_grid','XexScl_grid','idx','idy','cfg','-v7.3');
end

