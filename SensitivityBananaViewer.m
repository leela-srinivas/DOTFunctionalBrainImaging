%% Forward Problem Light Propagation Model
% Created March 31st, 2025
% Last Updated: April 9th, 2025

close all; clear
debug = false;
load('NeuroDOT_Data_Sample_CCW1.mat');

%% Global Variable Declarations: properties of the brain we're trying to image
global mua musp nu D xBnds yBnds zBnds mmX mmY mmZ X Y Z voxCrd;

mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
nu = 1.4; 
D = 1/(3*(mua+musp));

% Define bounds on medium
xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];  
mmX = 2; mmY = 2; mmZ = 2; 
voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

[Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 

%% SD Pair Generation

srcs = info.optodes.spos2;
srcs = cat(2, srcs, zeros(24,1));
dets = info.optodes.dpos2;
dets = cat(2, dets, zeros(28,1));
pairs = info.pairs; pairs = table2array(pairs(1:672,1:2));

figure();
hold on
scatter(srcs(:, 1), srcs(:, 2));
scatter(dets(:, 1), dets(:, 2));
hold off

if debug; disp(srcs); end

%% SD Pair Generation, Sparse DOT Array Model

raw_srcs = info.optodes.spos2;
srcs = raw_srcs(1:2:end, :);
srcs = cat(2, srcs, zeros(length(srcs), 1));
numSrcs = length(srcs);

raw_dets = info.optodes.dpos2;
dets_count = 0;
for i=1:length(raw_dets)
    if mod(i, 8) ~= 2 && mod(i, 8) ~= 4 && mod(i, 8) ~= 5 && mod(i, 8) ~= 7
        
        dets_count = dets_count + 1;
        dets(dets_count) = raw_dets(i);
    end
end

dets = cat(2, dets, zeros(length(dets),1));
numDets = length(dets);

raw_pairs = info.pairs; 
pairs = table2array(raw_pairs(1:2:numDets * numSrcs,1:2));
pairs(pairs(:, 1) > 12, :) = [];
%pairs(mod(pairs(:,2),2) == 0, :) = [];

figure();
hold on
scatter(srcs(:, 1), srcs(:, 2));
scatter(dets(:, 1), dets(:, 2));
hold off
legend("Sources","Detectors");

%clear raw_dets raw_srcs;

%% Sensitivity Matrix Generation

for i = 1:2
    tmpSrc2Voxels = greensSrc(srcs(i,:));
end
for i = 1:2
    tmpVoxels2Dets = greensDet(dets(i,:));
end

numVoxCrds = length(tmpSrc2Voxels);
numSDPairs = length(pairs);
sensitivityMatrix = zeros(numSDPairs, numVoxCrds);
lightFalloffVals = zeros(numSDPairs, 2);

for i = 1:size(pairs,1)
    j = pairs(i,1);
    k = pairs(i,2);
    tmpSrc2Voxels = greensSrc(srcs(j,:));
    tmpVoxels2Dets = greensDet(dets(k,:));
    [tmpSrc2Det, dist] = greensSrc2Det(srcs(j,:),dets(k,:));
    lightFalloffVals(i, 1) = dist;
    lightFalloffVals(i, 2) = tmpSrc2Det;


    sensitivityMatrix(i, :) = 1/D * tmpSrc2Voxels .* tmpVoxels2Dets.' / tmpSrc2Det;
end

%% visualizing A matrix

% sensitivity banana is all x's for one y
% y is all the source-detector pairs (measurements). x is all the voxel crds.

figure();
hold on
colormap summer
for i= 1:numSDPairs
    reshaped_A = reshape(sensitivityMatrix(i, :), [71 31 5]);
    % sgtitle("Sensitivity Matrix: Source " + pairs(i, 1) + ", Detector " + pairs(i, 2));
    for j = 1:2
        imagesc((reshaped_A(:, :, j)));
        title("Sensitivity Matrix: Source " + pairs(i, 1) + ", Detector " + pairs(i, 2) + ", Z = " + j);
        colorbar();
        pause(0.2);
    end
    pause(0.1);
end

hold off

%% Convolving HRF -- not used

load('hrf_DOT3.mat');
%conv(hrf, nu);

%% Simulating Light Falloff

figure();
scatter(lightFalloffVals(1:numSDPairs/32, 1), log(lightFalloffVals(1:numSDPairs/32, 2))); % only plotting first 32nd, it's all the same anyway
ylabel("log of SD-pair intensity");
xlabel("Distance between SD-pair (nm probably)");

% labeling for funsies
for i = 1:numSDPairs/32
    labelstr = "(" + pairs(i, 1) + ", " + pairs(i, 2) + ")";
    
    text(lightFalloffVals(i, 1), log(lightFalloffVals(i, 2)), labelstr);
end


%% Using Light Falloff to determine which SD-Pairs we use

% take out anything above 40 nm
dist_threshold = 40;

newLightFallOffVals = lightFalloffVals(lightFalloffVals(:,1)<=dist_threshold,:);
boolLightFalloffs = lightFalloffVals(:,1)<=dist_threshold;
j = 1;
maxSDPairs = numSDPairs;
for i = 1:maxSDPairs
    if boolLightFalloffs(i) == 0
        sensitivityMatrix(j, :) = [];

        j = j - 1;
        maxSDPairs = maxSDPairs - 1;
    end
    j = j + 1;
end

%% visualize after taking out extraneous SD pairs
visualizeAMatrix(sensitivityMatrix);


function reshaped_A = visualizeAMatrix(A)
    
    numSDPairs = size(A, 1);

    figure();
    hold on
    colormap summer
    for i= 1:numSDPairs
        reshaped_A = reshape(A(i, :), [71 31 5]);
        % sgtitle("Sensitivity Matrix: Source " + pairs(i, 1) + ", Detector " + pairs(i, 2));
        for j = 1:5
            imagesc((reshaped_A(:, :, j)));
            title("Sensitivity Matrix");
            colorbar();
            pause(0.2);
        end
        pause(0.1);
    end
end


% Code from InfiniteGreensFunctionSlab.m
function GsAnalytic = greensSrc(pos)

    mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
    musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
    nu = 1.4; 
    D = 1/(3*(mua+musp));
    mu_eff = sqrt(mua/D); 

    % Define bounds on medium
    xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];  
    mmX = 2; mmY = 2; mmZ = 2; 
    
    [Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 
    voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector
    nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels

    r = pdist2(pos,voxCrd); % distance from source to each voxel
    
    GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, each voxel
    
    tmp = reshape(GsAnalytic,nX,nY,nZ); % reshape as 3D volume

end

function GsAnalytic = greensDet(pos)

    mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
    musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
    nu = 1.4; 
    D = 1/(3*(mua+musp));
    mu_eff = sqrt(mua/D); 

    % Define bounds on medium
    xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];  
    mmX = 2; mmY = 2; mmZ = 2; 
    
    [Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 
    voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector
    nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels
    
    %srcPos = [0 0 0]; % source position, in 3D coordinates 
    
    r = pdist2(voxCrd, pos); % distance from detector to each voxel
    
    GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, each voxel
    
    tmp = reshape(GsAnalytic,nX,nY,nZ); % reshape as 3D volume

end

function [GsAnalytic, r] = greensSrc2Det(srcPos, detPos)

    mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
    musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
    nu = 1.4; 
    D = 1/(3*(mua+musp));
    mu_eff = sqrt(mua/D); 
    
    % Define bounds on medium
    xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];  
    mmX = 2; mmY = 2; mmZ = 2; 
    
    [Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 
    voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector
    nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels
    
    %srcPos = [0 0 0]; % source position, in 3D coordinates 
    
    r = pdist2(srcPos, detPos); % distance from detector to each voxel
    
    GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, each voxel
end
