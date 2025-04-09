%% Forward Problem Light Propagation Model
% Created March 31st, 2025
% Last Updated: April 7th, 2025

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

NPerWall = 3; % number of sources/detectors (will probably need separate code generating these pairs later on
srcs = info.optodes.spos2;
srcs = cat(2, srcs, zeros(24,1));
dets = info.optodes.dpos2;
dets = cat(2, dets, zeros(28,1));
pairs = info.pairs; pairs = table2array(pairs(1:672,1:2));

if debug; disp(srcs); end


%srcsMatrix(:, 3) = linspace(zBnds(1), zBnds(2), N);


% srcs = [0 30 0; -45 0 0];
% dets = [0 -30 0; 45 0 0];


%% Sensitivity Matrix Generation

N = 4 * NPerWall;
numMeasurements = N * N - 1;
%tmpSrc2Det = zeros(1, numMeasurements);
%tmpSrc2Voxels = zeros(numMeasurements, 24);
%tmpVoxels2Dets = zeros(24, numMeasurements);

for i = 1:2
    tmpSrc2Voxels = greensSrc(srcs(i,:));
end
for i = 1:2
    tmpVoxels2Dets = greensDet(dets(i,:));
end

for i = 1:size(pairs,1)
    j = pairs(i,1);
    k = pairs(i,2);
    tmpSrc2Voxels = greensSrc(srcs(j,:));
    tmpVoxels2Dets = greensDet(dets(k,:));
    tmpSrc2Det = greensSrc2Det(srcs(j,:),dets(k,:));

    sensitivityMatrix(i, :) = 1/D * tmpSrc2Voxels .* tmpVoxels2Dets.' / tmpSrc2Det;
end

% for i = 1:N
%     for j = 1:N - 1
%         tmpSrc2Voxels(((i - 1) * N + j), :) = greensSrc(srcs(i, :));
%         tmpVoxels2Dets(:, ((i - 1) * N + j)) = greensDet(dets(i, :));
%     end
%     tmpSrc2Det(i) = greensSrc2Det(srcs(i, :), dets(i, :));
% end
% tmpSrc2Voxels = [greensSrc(srcs(1, :)); greensSrc(srcs(2, :))];
% tmpVoxels2Det = [greensDet(dets(1, :)) greensDet(dets(2, :))];

%N = length(tmpSrc2Voxels(:, 1));

% tmpSrc2Det = zeros(1, N);
% tmpSrc2Det = [greensSrc2Det(srcs(1, :), dets(1, :)), greensSrc2Det(srcs(2, :), dets(2, :))];
% % 

%sensitivityMatrix = zeros(numMeasurements, length(tmpSrc2Voxels));

% for k=1:(24*28)
%     %tmpSrc2Dets(i) = greensSrc2Det(srcPos1, detPos1);
%     disp(k)
%     sensitivityMatrix(k, :) = 1/D * tmpSrc2Voxels(:, k) .* tmpVoxels2Dets(k, :).' / tmpSrc2Det(k);
% end

%% visualizing A matrix

% sensitivity banana is all x's for one y
% y is all the source-detector pairs (measurements). x is all the voxel
% crds.

numSDPairs = 672;
numVoxCrds = 11005;

figure();
hold on
colormap parula
for i= 1:numSDPairs
    reshaped_A = reshape(sensitivityMatrix(i, :), [71 31 5]);

    for j = 1:5
        imagesc((reshaped_A(:, :, j)));
        title("Source " + pairs(i, 1) + ", Detector " + pairs(i, 2) + ", Z = " + j);
        colorbar();
        pause(0.2);
    end
    pause(0.1);
end

hold off



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
    
    %srcPos = [0 0 0]; % source position, in 3D coordinates 
    
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

function GsAnalytic = greensSrc2Det(srcPos, detPos)

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
