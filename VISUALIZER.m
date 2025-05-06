% combined_DOT_with_Lcurve_and_Zchoice.m

clear; clc; close all;

%% --- Step 1: Load data & define voxel grid ---
load('NeuroDOT_Data_Sample_CCW1.mat','data','info');

% Grid
xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];
mmX = 2; mmY = 2; mmZ = 2;
[Ygrid, Xgrid, Zgrid] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2));
[nx, ny, nz] = size(Xgrid);
voxCrd = [Xgrid(:), Ygrid(:), Zgrid(:)];
numVoxels = size(voxCrd,1);

% Optical properties
mua = 0.0192; musp = 0.6726;
D = 1/(3*(mua + musp)); mu_eff = sqrt(mua / D);

%% --- Step 2: Build A matrix from measured pairs ---
pairsTable = info.pairs(:,1:2);
varNames = pairsTable.Properties.VariableNames;
numPairs = height(pairsTable);
pairs = zeros(numPairs,2);
for k = 1:2
    col = pairsTable.(varNames{k});
    if ~isnumeric(col), pairs(:,k) = cellfun(@str2double, col);
    else, pairs(:,k) = col; end
end

Afull = zeros(numPairs, numVoxels);
for p = 1:numPairs
    sIdx = pairs(p,1); dIdx = pairs(p,2);
    src = info.optodes.spos3(sIdx,:); det = info.optodes.dpos3(dIdx,:);
    r_src = sqrt(sum((voxCrd - src).^2, 2));
    r_det = sqrt(sum((voxCrd - det).^2, 2));
    r_sd = norm(src - det);
    G_src = exp(-mu_eff*r_src)./(4*pi*D*r_src);
    G_det = exp(-mu_eff*r_det)./(4*pi*D*r_det);
    G_sd = exp(-mu_eff*r_sd)/(4*pi*D*r_sd);
    Afull(p,:) = (1/D)*(G_src .* G_det / G_sd)';
end

%% --- Step 3: Channel selection (NN==2 & WL==2) ---
sel = info.pairs.NN==2 & info.pairs.WL==2;
A_sel = Afull(sel,:);
data_sel = data(sel,:);
nCh_sel = sum(sel);

%% --- Step 4: Epoch extraction & ΔOD (robust windowing) ---
syn = info.paradigm.synchpts;
winLen = 388;
T0 = size(data_sel,2);
fs = info.system.framerate;
runStarts = syn(2:end);
validRuns = runStarts(runStarts + winLen - 1 <= T0);
nRuns = numel(validRuns);

runsDod = zeros(nCh_sel, winLen, nRuns);
for k = 1:nRuns
    idx = validRuns(k):validRuns(k)+winLen-1;
    Y = data_sel(:, idx);
    runsDod(:,:,k) = -log(bsxfun(@rdivide, Y, mean(Y,2)));
end

avgDod = mean(runsDod, 3);
if nRuns < numel(syn)-1
    warning('Dropped %d incomplete runs.', (numel(syn)-1) - nRuns);
end

%% --- Step 5: Std-dev channel rejection ---
v = std(avgDod, [], 2);
keep = v <= 0.05;
avgDod = avgDod(keep,:);
A = A_sel(keep,:);

%% --- Step 6: Band-pass filter (no detrend) ---
[b,a] = butter(2, [0.01 0.2]/(fs/2));
for ch = 1:size(avgDod,1)
    avgDod(ch,:) = filtfilt(b, a, avgDod(ch,:));
end

%% --- Step 7: L-curve to choose optimal lambda ---
disp('Performing L-curve analysis...');
[U,S,V] = svd(A, 'econ');
s = diag(S);
T = size(avgDod, 2);
Uy = U' * avgDod;

lambdas = logspace(-15, -10, 100);
residuals = zeros(size(lambdas));
solutions = zeros(size(lambdas));
for i = 1:length(lambdas)
    lam = lambdas(i);
    filt = s ./ (s.^2 + lam^2);
    X = V * (filt .* Uy);
    residuals(i) = norm(A * X - avgDod, 'fro');
    solutions(i) = norm(X, 'fro');
end

loglog(residuals, solutions, '-o'); grid on;
xlabel('Residual norm ||Ax - y||'); ylabel('Solution norm ||x||');
title('L-Curve'); drawnow;

[~, elbowIdx] = min(abs(log10(residuals ./ max(residuals)) + log10(solutions ./ max(solutions))));
lambda_opt = lambdas(elbowIdx);
fprintf('✅ Optimal λ from L-curve: %.2e\n', lambda_opt);

%% --- Step 8: Reconstruct
invF = s ./ (s.^2 + lambda_opt^2);
Xrec = zeros(size(V,1), T);
for t = 1:T
    Xrec(:,t) = V * (invF .* (U' * avgDod(:,t)));
end
vol4D = reshape(Xrec, [nx, ny, nz, T]);

%% --- Step 9: Choose z-slice to visualize ---
fprintf('Available z slices: z = [%s]\n', num2str(unique(voxCrd(:,3))'));
z_values = unique(voxCrd(:,3));
z_idx = input(sprintf('Which z-depth to display? Choose from %s: ', mat2str(z_values')));
z_slice_idx = find(abs(z_values - z_idx) < 1e-3, 1);
if isempty(z_slice_idx)
    error('Invalid z value selected.');
end
z_vox_idx = round(z_slice_idx);

volSlice = squeeze(vol4D(:,:,z_vox_idx,:));
for i = 1:nx
    for j = 1:ny
        volSlice(i,j,:) = medfilt1(squeeze(volSlice(i,j,:)), 5);
    end
end

cmin = min(volSlice(:));
cmax = max(volSlice(:));

%% --- Step 10: Animate and save video with accurate timing ---
vObj = VideoWriter(sprintf('DOT_Recon_z%d.avi', z_idx));
vObj.FrameRate = fs;
open(vObj);

fig = figure('Name','DOT Video (Z-slice)','NumberTitle','off');
hImg = imagesc(volSlice(:,:,1), [cmin cmax]);
axis image off; colormap(jet); colorbar;
title('t = 0.00 s');

t0 = validRuns(1) / fs;  % use sync timing for accuracy
for t = 1:T
    set(hImg, 'CData', volSlice(:,:,t));
    title(sprintf('t = %.2f s', t0 + (t-1)/fs));
    frame = getframe(gcf);
    writeVideo(vObj, frame);
end

close(vObj);
fprintf('🎥 Video saved: DOT_Recon_z%d.avi\n', z_idx);
