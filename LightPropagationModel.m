%% Forward Problem Light Propagation Model
% Created March 31st, 2025
% Last Updated: April 1st, 2025
% Note: run after Infinite Greens Function Slab.m

close all; clear

%% Visualize in 3D

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

%% idk yet
N = 1; % number of SD pairs
srcPos = [0 10 0]; % example, not committed to this
detPos = [0 -10 0]; % example, not committed to this

tmpSrc2Voxels = greensSrc(srcPos);
tmpVoxels2Det = greensDet(detPos);
tmpSrc2Det = greensSrc2Det(srcPos, detPos);

matrix = 1/D * tmpSrc2Voxels .* tmpVoxels2Det;

% Code from InfiniteGreensFunctionSlab.m
function tmp = greensSrc(pos)

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

function tmp = greensDet(pos)

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
