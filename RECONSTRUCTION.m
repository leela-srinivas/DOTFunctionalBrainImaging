%% ReconstructImageVolume.m
% This script loads your DOT/fNIRS dataset, builds the sensitivity matrix,
% performs Tikhonov-regularized image reconstruction for one time frame,
% and opens the interactive Volume Viewer to display the 3D volume.

clear; close all; clc;

%% Step 1: Load Data
% Ensure the file 'NeuroDOT_Data_Sample_CCW1.mat' is in your MATLAB working directory.
load('NeuroDOT_Data_Sample_CCW1.mat'); 

%% Step 2: Build the Sensitivity Matrix (A)
[A, gridSize] = buildSensitivityMatrix(info);
% 'A' is the sensitivity matrix; 'gridSize' gives the dimensions for reshaping the reconstruction.

%% Step 3: Prepare the Measurement Data
% Compute the differential log-ratio: y = -log(I / (mean(I over time)))
mean_data = mean(data, 2);
y = -log(data ./ repmat(mean_data, 1, size(data,2)));

% Select a single time frame for reconstruction (modify the index as needed)
timeIndex = 50;
y_frame = y(:, timeIndex);

%% Step 4: Regularized Image Reconstruction
lambda = 0.1;   % Adjust this regularization parameter if needed
x_recon = myImageReconstruction(y_frame, A, lambda, gridSize);

%% Step 5: Visualize the Reconstructed Image in 3D Using volumeViewer
% Use volumeViewer interactively to view the 3D volume (absolute value).
volumeViewer(abs(x_recon));
title('Reconstructed Image (3D Interactive View)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = myImageReconstruction(y, A, lambda, gridSize)
    % Performs image reconstruction via Tikhonov regularization.
    % Solves x_hat = (A'*A + lambda*I) \ (A'*y)
    A_reg = (A' * A + lambda * eye(size(A,2))) \ (A');
    x_hat = A_reg * y;
    x = reshape(x_hat, gridSize);
end

function [A, gridSize] = buildSensitivityMatrix(info)
    % Constructs the sensitivity matrix A for DOT imaging.
    % Defines a voxel grid over the imaging medium and uses optode positions from info.optodes.
    
    % Define voxel grid boundaries and resolution (adjust as required)
    xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];
    mmX = 2; mmY = 2; mmZ = 2;
    [Ygrid, Xgrid, Zgrid] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2));
    gridSize = size(Xgrid);
    
    % Get voxel coordinates (each row = [x y z])
    voxCrd = [Xgrid(:), Ygrid(:), Zgrid(:)];
    
    % Retrieve source and detector 3D positions from info.optodes
    if isfield(info.optodes, 'spos3')
        srcs = info.optodes.spos3;
    else
        error('Source positions not found: info.optodes.spos3 is missing.');
    end
    if isfield(info.optodes, 'dpos3')
        dets = info.optodes.dpos3;
    else
        error('Detector positions not found: info.optodes.dpos3 is missing.');
    end

    % Extract source-detector pairs from info.pairs. Convert non-numeric types if necessary.
    if isfield(info, 'pairs')
        T = info.pairs(:, 1:2); 
        varNames = T.Properties.VariableNames;
        numPairs = height(T);
        pairs = zeros(numPairs, 2);
        for k = 1:2
            col = T.(varNames{k});
            if ~isnumeric(col)
                pairs(:, k) = cellfun(@(c) str2double(c), col);
            else
                pairs(:, k) = col;
            end
        end
    else
        % If not provided, assume all combinations.
        [iSrc, iDet] = ndgrid(1:size(srcs,1), 1:size(dets,1));
        pairs = [iSrc(:), iDet(:)];
    end

    % Initialize the sensitivity matrix A (numPairs x numVoxels)
    numPairs = size(pairs, 1);
    numVoxels = size(voxCrd, 1);
    A = zeros(numPairs, numVoxels);
    
    % Set optical properties (adjust if needed)
    mua = 0.0192; 
    musp = 0.6726;  
    D = 1/(3*(mua + musp));
    mu_eff = sqrt(mua / D);
    
    % Loop over each pair to compute sensitivity at each voxel.
    for p = 1:numPairs
        sIdx = pairs(p, 1);
        dIdx = pairs(p, 2);
        srcPos = srcs(sIdx, :);
        detPos = dets(dIdx, :);
        
        % Compute distances from source and detector to each voxel.
        r_src = sqrt(sum((voxCrd - srcPos).^2, 2));
        r_det = sqrt(sum((voxCrd - detPos).^2, 2));
        
        % Calculate Green's functions from source and to detector.
        G_src = 1./(4 * pi * D * r_src) .* exp(-mu_eff * r_src);
        G_det = 1./(4 * pi * D * r_det) .* exp(-mu_eff * r_det);
        
        % Normalize by direct source-detector coupling.
        r_sd = norm(srcPos - detPos);
        G_sd = 1/(4 * pi * D * r_sd) * exp(-mu_eff * r_sd);
        
        % Compute sensitivity: A(p,:) = (1/D) * (G_src .* G_det / G_sd)
        A(p, :) = (1/D) * (G_src .* G_det / G_sd)';
    end
end
