%% Forward Problem Light Propagation Model
% Created March 31st, 2025
% Last Updated: April 2nd, 2025

close all; clear

%% Global Variable Declarations: properties of the brain we're trying to image
global mua musp nu D xBnds yBnds zBnds mmX mmY mmZ X Y Z voxCrd;

mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
nu = 1.4; 
D = 1/(3*(mua+musp));

% Define bounds on medium
xBnds = [-30 30]; yBnds = [-45 45]; zBnds = [1 30];  
mmX = 2; mmY = 2; mmZ = 2; 
voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

[Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 

%% Visualize in 3D



%% Sensitivity Matrix Generation
N = 2; % number of SD pairs (will probably need separate code generating these pairs later on

srcs = [0 30 0; -45 0 0];
dets = [0 -30 0; 45 0 0];

tmpSrc2Voxels = [greensSrc(srcs(1, :)); greensSrc(srcs(2, :))];
tmpVoxels2Det = [greensDet(dets(1, :)) greensDet(dets(2, :))];

tmpSrc2Det = zeros(1, N);
tmpSrc2Det = [greensSrc2Det(srcs(1, :), dets(1, :)), greensSrc2Det(srcs(2, :), dets(2, :))];
% 
sensitivityMatrix = zeros(N, length(tmpSrc2Voxels));

for i=1:N
    %tmpSrc2Dets(i) = greensSrc2Det(srcPos1, detPos1);
    sensitivityMatrix(i, :) = 1/D * tmpSrc2Voxels(i, :) .* tmpVoxels2Det(:, i).' / tmpSrc2Det(i);
end

% Code from InfiniteGreensFunctionSlab.m
function GsAnalytic = greensSrc(pos)

    mua = .0192; % flags.op.mua_gray=[0.0180,0.0192];
    musp = 0.6726; % flags.op.musp_gray=[0.8359,0.6726];
    nu = 1.4; 
    D = 1/(3*(mua+musp));
    mu_eff = sqrt(mua/D); 

    % Define bounds on medium
    xBnds = [-30 30]; yBnds = [-45 45]; zBnds = [1 30];  
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
    xBnds = [-30 30]; yBnds = [-45 45]; zBnds = [1 30];  
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
    xBnds = [-30 30]; yBnds = [-45 45]; zBnds = [1 30];  
    mmX = 2; mmY = 2; mmZ = 2;

    [Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); % generate coordinates for slab 
    nX = size(X,1); nY = size(X,2); nZ = size(X,3); % volume is nX x nY x nZ voxels
    
    %srcPos = [0 0 0]; % source position, in 3D coordinates 
    
    r = pdist2(srcPos, detPos); % distance from detector to each voxel
    
    GsAnalytic = 1./(4*pi*D*r).*exp(-mu_eff*r); % Green's function, each voxel
end
