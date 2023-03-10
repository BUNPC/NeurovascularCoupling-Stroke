
clear
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

figure
for m = 1:12
    imagesc(data(m).maskSFDI(3).aff_prop_mus)
    axis image
    colormap jet
    caxis([0 25])
    hold on
    h = imrect;
    pos = wait(h);
    pos = round(pos);
    delete(h);
    for t = 1:4
    data(m).maskSFDI(t).aff_prop_mus = data(m).maskSFDI(t).aff_prop_mus(pos(2):1:pos(2)+pos(4),pos(1):1:pos(1)+pos(3));
    end
end

% Get baseline statistics
for m = 1:length(data)    
    temp = data(m).maskSFDI(1).aff_prop_mua; % + data(m).maskSFDI(1).unaff_prop_mua;
    temp(temp == 100) = 0;
    mua(m).base = temp;
    mua(m).mbase = mean(nonzeros(mua(m).base(:)));
    mua(m).stdbase = std(nonzeros(mua(m).base(:)));
    
    temp = data(m).maskSFDI(1).aff_prop_mus; % + data(m).maskSFDI(1).unaff_prop_mus;
    temp(temp == 100) = 0;
    mus(m).base = temp;
    mus(m).mbase = mean(nonzeros(mus(m).base(:)));
    mus(m).stdbase = std(nonzeros(mus(m).base(:)));
    data(m).mouse = m;
    
    temp = data(m).maskSFDI(3).aff_prop_mua; % + data(m).maskSFDI(1).unaff_prop_mua;
    temp(temp == 100) = 0;
    mua(m).week1 = temp;
    mua(m).mweek1 = mean(nonzeros(mua(m).week1(:)));
    mua(m).stdweek1 = std(nonzeros(mua(m).week1(:)));
    
    temp = data(m).maskSFDI(3).aff_prop_mus; % + data(m).maskSFDI(1).unaff_prop_mus;
    temp(temp == 100) = 0;
    mus(m).week1 = temp;
    mus(m).mweek1 = mean(nonzeros(mus(m).week1(:)));
    mus(m).stdweek1 = std(nonzeros(mus(m).week1(:)));
    data(m).mouse = m;
end

figure
count = 1;
subplot(2,2,1)
for m = 1:length(data)    
    bar([count;count+1],[mua(m).mbase;mua(m).mweek1])
    hold on
    er = errorbar([count;count+1],[mua(m).mbase;mua(m).mweek1],[mua(m).stdbase;mua(m).stdweek1],[mua(m).stdbase;mua(m).stdweek1]);
    er.Color = [0 0 0];   
    er.LineStyle = 'none';  
    count = count+2;
end
ylim([0 1])
ylabel('Absorption coefficient (mm^-^1)')
title('Individual mice')
set(gca,'FontSize',16,'XTickLabel',{''})
subplot(2,2,2)
mbase = mean([mua.mbase]);
stdbase = std([mua.mbase]);
mweek1 = mean([mua.mweek1]);
stdweek1 = std([mua.mweek1]);
bar([mbase;mweek1])
hold on
er = errorbar(1:2,[mbase;mweek1],[stdbase;stdweek1],[stdbase;stdweek1]);
er.Color = [0 0 0]; 
er.LineStyle = 'none';  
ylim([0 1])
title('Average of all mice')
set(gca,'FontSize',16,'XTickLabel',{''})  

subplot(2,2,3)
for m = 1:length(data)    
    bar([count;count+1],[mus(m).mbase;mus(m).mweek1])
    hold on
    er = errorbar([count;count+1],[mus(m).mbase;mus(m).mweek1],[mus(m).stdbase;mus(m).stdweek1],[mus(m).stdbase;mus(m).stdweek1]);
    er.Color = [0 0 0];   
    er.LineStyle = 'none';  
    count = count+2;
end
ylim([0 25])
ylabel('Scattering coefficient (mm^-^1)')
title('Individual mice')
set(gca,'FontSize',16,'XTickLabel',{''})
subplot(2,2,4)
mbase = mean([mus.mbase]);
stdbase = std([mus.mbase]);
mweek1 = mean([mus.mweek1]);
stdweek1 = std([mus.mweek1]);
bar([mbase;mweek1])
hold on
er = errorbar(1:2,[mbase;mweek1],[stdbase;stdweek1],[stdbase;stdweek1]);
er.Color = [0 0 0]; 
er.LineStyle = 'none';  
ylim([0 20])
title('Average of all mice')
set(gca,'FontSize',16,'XTickLabel',{''})  

%% GFP mouse
clear
data(1) = load('F:\SS86\brainmaskSFDI.mat');

for m = 1:length(data)  
    for b = 1:3
    temp = data(m).maskSFDI(b).aff_prop_mua + data(m).maskSFDI(b).unaff_prop_mua;
    temp(temp == 100) = 0;
    mua(b).base = temp;
    mua(b).mbase = mean(nonzeros(mua(b).base(:)));
    mua(b).stdbase = std(nonzeros(mua(b).base(:)));
    
    temp = data(m).maskSFDI(b).aff_prop_mus + data(m).maskSFDI(b).unaff_prop_mus;
    temp(temp == 100) = 0;
    mus(b).base = temp;
    mus(b).mbase = mean(nonzeros(mus(b).base(:)));
    mus(b).stdbase = std(nonzeros(mus(b).base(:)));
    end
end
figure
subplot(2,2,1)
for m = 1:3   
    bar(m,mua(m).mbase)
    hold on
    er = errorbar(m,mua(m).mbase,mua(m).stdbase,mua(m).stdbase);
    er.Color = [0 0 0];                            
end
ylim([0 1])
ylabel('Absorption coefficient')
title('Individual mice')
set(gca,'FontSize',16,'XTickLabel',{''})
subplot(2,2,2)
mbase = mean([mua.mbase]);
stdbase = std([mua.mbase]);
bar(mbase)
hold on
er = errorbar(1,mbase,stdbase,stdbase);
er.Color = [0 0 0]; 
ylim([0 1])
title('Average of all mice')
set(gca,'FontSize',16,'XTickLabel',{''})  

subplot(2,2,3)
for m = 1:3   
    bar(m,mus(m).mbase)
    hold on
    er = errorbar(m,mus(m).mbase,mus(m).stdbase,mus(m).stdbase);
    er.Color = [0 0 0];                            
end
ylim([0 20])
ylabel('Scattering coefficient')
title('Individual mice')
set(gca,'FontSize',16,'XTickLabel',{''})
subplot(2,2,4)
mbase = mean([mus.mbase]);
stdbase = std([mus.mbase]);
bar(mbase)
hold on
er = errorbar(1,mbase,stdbase,stdbase);
er.Color = [0 0 0]; 
ylim([0 20])
title('Average of all mice')
set(gca,'FontSize',16,'XTickLabel',{''})  

% affmua = cat(1, affmua{:});
% affmus = cat(1, affmus{:});
% unaffmua = cat(1, unaffmua{:});
% unaffmus = cat(1, unaffmus{:});

% mua.base = affmua{1}+unaffmua{1};
% mua.mbase = mean(nonzeros(mua.base(:)));
% mua.stdbase = std(nonzeros(mua.base(:)));
% mua.SD1 = mua.mbase + (1*mua.stdbase);
% mua.SD2 = mua.mbase + (2*mua.stdbase);
% mua.SD3 = mua.mbase + (3*mua.stdbase);
% 
% mua.day2 = affmua{2}+unaffmua{2};
% mua.day2(mua.day2 > mua.SD2) = 100;
% mua.week1 = affmua{3}+unaffmua{3};
% mua.week1(mua.week1 > mua.SD3) = 100;
% 
% mus.base = affmus{1}+unaffmus{1};
% mus.mbase = mean(nonzeros(mus.base(:)));
% mus.stdbase = std(nonzeros(mus.base(:)));
% mus.SD3 = mus.mbase + (3*mus.stdbase);
% mus.SD6 = mus.mbase + (6*mus.stdbase);
% mus.SD9 = mus.mbase + (9*mus.stdbase);
% 
% mus.day2 = affmus{2}+unaffmus{2};
% mus.day2(mus.day2 > mus.SD2) = 100;
% mus.week1 = affmus{3}+unaffmus{3};
% mus.week1(mus.week1 > mus.SD9) = 100;


%%
figure
subplot(1,4,1)
imagesc(unaffmua)
axis image
caxis([0 1])
colorbar 
colormap jet
subplot(1,4,2)
imagesc(affmua)
axis image
caxis([0 1])
colorbar 
colormap jet
subplot(1,4,3)
imagesc(unaffmus)
axis image
caxis([0 30])
colorbar 
colormap jet
subplot(1,4,4)
imagesc(affmus)
axis image
caxis([0 30])
colorbar 
colormap jet


