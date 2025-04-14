%% ReconstructImageVolume_LogColor_NoiseControlled.m
% This script reconstructs a 3D image from DOT/fNIRS data with noise added
% according to Week 11 instructions: Gaussian noise added to log-ratio measurements.
% It uses Tikhonov regularization and displays log-scaled 3D images with volumeViewer.

clear; close all; clc;

%% Step 1: Load Data
load('NeuroDOT_Data_Sample_CCW1.mat'); 

%% Step 2: Build the Sensitivity Matrix (A)
[A, gridSize] = buildSensitivityMatrix(info);

%% Step 3: Prepare Measurement Data (Log-Ratio) and Add Gaussian Noise
% Compute the clean log-ratio data
mean_data = mean(data, 2);
y_clean = -log(data ./ mean_data);  % clean log-ratio signal

% Choose time point for reconstruction
timeIndex = 50;
y_clean_t = y_clean(:, timeIndex);

% Set target NSR = 1/SNR
target_nsr = 0.05;  % try 0.01 to 0.1 (5% noise)
signal_norm = norm(y_clean_t);
noise_std = target_nsr * signal_norm / sqrt(numel(y_clean_t));

% Add Gaussian noise to the measurement vector
y_noisy_t = y_clean_t + noise_std * randn(size(y_clean_t));

%% Step 4: Tikhonov Regularized Image Reconstruction
lambda = 1;  % Regularization parameter (adjustable)
x_recon = myImageReconstruction(y_noisy_t, A, lambda, gridSize);

%% Step 5: Visualize 3D Volume (Log Scale)
log_x = log10(abs(x_recon) + eps);  % Log-scale the absolute image
volumeViewer(log_x);
title('Log-scaled Reconstructed Image with NSR = 0.05');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = myImageReconstruction(y, A, lambda, gridSize)
    A_reg = (A' * A + lambda * eye(size(A,2))) \ (A');
    x_hat = A_reg * y;
    x = reshape(x_hat, gridSize);
end

function [A, gridSize] = buildSensitivityMatrix(info)
    xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];
    mmX = 2; mmY = 2; mmZ = 2;
    [Ygrid, Xgrid, Zgrid] = meshgrid(yBnds(1):mmY:yBnds(2), ...
                                     xBnds(1):mmX:xBnds(2), ...
                                     zBnds(1):mmZ:zBnds(2));
    gridSize = size(Xgrid);
    voxCrd = [Xgrid(:), Ygrid(:), Zgrid(:)];
    
    srcs = info.optodes.spos3;
    dets = info.optodes.dpos3;

    % Extract source-detector pairs
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

    % Optical properties
    mua = 0.0192; musp = 0.6726;
    D = 1 / (3 * (mua + musp));
    mu_eff = sqrt(mua / D);

    % Build A matrix
    numVoxels = size(voxCrd, 1);
    A = zeros(numPairs, numVoxels);
    for p = 1:numPairs
        sIdx = pairs(p, 1);
        dIdx = pairs(p, 2);
        srcPos = srcs(sIdx, :);
        detPos = dets(dIdx, :);
        r_src = sqrt(sum((voxCrd - srcPos).^2, 2));
        r_det = sqrt(sum((voxCrd - detPos).^2, 2));
        G_src = 1./(4*pi*D*r_src) .* exp(-mu_eff*r_src);
        G_det = 1./(4*pi*D*r_det) .* exp(-mu_eff*r_det);
        r_sd = norm(srcPos - detPos);
        G_sd = 1/(4*pi*D*r_sd) * exp(-mu_eff*r_sd);
        A(p, :) = (1/D) * (G_src .* G_det / G_sd)';
    end
end
