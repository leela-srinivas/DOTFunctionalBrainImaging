% combined_DOT_full_pipeline.m
% —————————————————————————————————————————————
% Computes A on‐the‐fly, selects NN==2 & WL==2 channels, extracts & averages
% ΔOD epochs, rejects noisy channels, applies optional filtering,
% performs Tikhonov‐SVD reconstruction, and saves a mid‐slice video.

clear; clc; close all;

%% --- Step 1: Load data & define voxel grid ---
load('NeuroDOT_Data_Sample_CCW1.mat','data','info');

% Voxel grid (mm)
xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];
mmX = 2; mmY = 2; mmZ = 2;
[Ygrid, Xgrid, Zgrid] = meshgrid( ...
    yBnds(1):mmY:yBnds(2), ...
    xBnds(1):mmX:xBnds(2), ...
    zBnds(1):mmZ:zBnds(2) ...
);
[nx, ny, nz] = size(Xgrid);
voxCrd = [Xgrid(:), Ygrid(:), Zgrid(:)];
numVoxels = size(voxCrd,1);

% Optical properties
mua    = 0.0192;
musp   = 0.6726;
D      = 1/(3*(mua + musp));
mu_eff = sqrt(mua / D);

%% --- Step 2: Read measured source–detector pairs & build Afull ---
pairsTable = info.pairs(:,1:2);
varNames   = pairsTable.Properties.VariableNames;
numPairs   = height(pairsTable);
pairs      = zeros(numPairs,2);
for k = 1:2
    col = pairsTable.(varNames{k});
    if ~isnumeric(col)
        pairs(:,k) = cellfun(@str2double, col);
    else
        pairs(:,k) = col;
    end
end

Afull = zeros(numPairs, numVoxels);
for p = 1:numPairs
    sIdx   = pairs(p,1);
    dIdx   = pairs(p,2);
    srcPos = info.optodes.spos3(sIdx,:);
    detPos = info.optodes.dpos3(dIdx,:);
    
    r_src = sqrt(sum((voxCrd - srcPos).^2, 2));
    r_det = sqrt(sum((voxCrd - detPos).^2, 2));
    r_sd  = norm(srcPos - detPos);
    
    G_src = exp(-mu_eff*r_src)./(4*pi*D*r_src);
    G_det = exp(-mu_eff*r_det)./(4*pi*D*r_det);
    G_sd  = exp(-mu_eff*r_sd) /(4*pi*D*r_sd);
    
    Afull(p,:) = (1/D) * (G_src .* G_det / G_sd)';
end

%% --- Step 3: Channel selection (NN==2 & WL==2) ---
sel      = info.pairs.NN==2 & info.pairs.WL==2;
A_sel    = Afull(sel, :);
data_sel = data(sel, :);
nCh_sel  = sum(sel);

%% --- Step 4: Epoch extraction & ΔOD per run (robust to end‐of‐data) ---
syn    = info.paradigm.synchpts;    % sync points
winLen = 4000;                       % samples per epoch
T0     = size(data_sel, 2);

% Only keep those runs whose full window fits inside the data length T0
runStarts = syn(2:end);
validRuns = runStarts(runStarts + winLen - 1 <= T0);
nRuns     = numel(validRuns);

runsDod = zeros(nCh_sel, winLen, nRuns);
for kRun = 1:nRuns
    idx = validRuns(kRun) : validRuns(kRun) + winLen - 1;
    Y   = data_sel(:, idx);
    % compute ΔOD relative to this run’s own mean:
    runsDod(:,:,kRun) = -log(bsxfun(@rdivide, Y, mean(Y,2)));
end

% Average ΔOD across valid runs
avgDod = mean(runsDod, 3);

if nRuns < numel(syn)-1
    warning('Dropped %d incomplete run(s) at end because window exceeded data length.', ...
            (numel(syn)-1) - nRuns);
end

%% --- Step 5: Std‐dev based channel rejection ---
v    = std(avgDod, [], 2);
keep = v <= 0.05;
avgDod = avgDod(keep, :);
A      = A_sel(keep, :);

%% --- Step 6: Optional signal processing (detrend + bandpass) ---
fs     = info.system.framerate;
[b,a]  = butter(4, [0.01 0.2] / (fs/2));
for ch = 1:size(avgDod,1)
    avgDod(ch,:) = filtfilt(b,a,avgDod(ch,:));
end

%% --- Step 7: SVD + Tikhonov inversion ---
lambda = 1e-12;
[U,S,V] = svd(A, 'econ');
s       = diag(S);
invF    = s ./ (s.^2 + lambda^2);

T    = size(avgDod, 2);
Xrec = zeros(size(V,1), T);
for t = 1:T
    Uy        = U' * avgDod(:,t);
    Xrec(:,t) = V * (invF .* Uy);
end

% Reshape into 4D volume [nx×ny×nz×T]
vol4D = reshape(Xrec, [nx, ny, nz, T]);

%% --- Step 8: Mid‐slice & median filter over time ---
zmid     = round(nz/2);
volSlice = squeeze(vol4D(:,:,zmid,:));
for i = 1:nx
    for j = 1:ny
        volSlice(i,j,:) = medfilt1(squeeze(volSlice(i,j,:)), 5);
    end
end

cmin = min(volSlice(:));
cmax = max(volSlice(:));

%% --- Step 9: Animate & save video ---
vObj = VideoWriter('Combined_DOT_Recon.avi');
vObj.FrameRate = fs;
open(vObj);

fig  = figure('Name','Combined DOT Reconstruction','NumberTitle','off');
hImg = imagesc(volSlice(:,:,1), [cmin cmax]);
axis image off; colormap(jet); colorbar;
title('t = 0.00 s');

tVec = (0:T-1) / fs;
for t = 1:T
    set(hImg, 'CData', volSlice(:,:,t));
    title(sprintf('t = %.2f s', tVec(t)));
    frame = getframe(gcf);
    writeVideo(vObj, frame);
end

close(vObj);
fprintf('✅ Video saved: Combined_DOT_Recon.avi\n');
