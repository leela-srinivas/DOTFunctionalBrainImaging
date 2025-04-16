function [x] = imageReconstruction(y, A, lambda, noise_factor)

    signal_norm = norm(y);
    noise_std = noise_factor * signal_norm / sqrt(numel(y));
    y_with_noise = y + noise_std * randn(size(y));

    disp("Inverting sensitivity matrix");
    A_reg = inv(A'*A + lambda*eye(size(A,2)))*A';

    disp("Getting x_hat");
    x_hat = A_reg*y_with_noise;

    % Define bounds on medium
    xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];  
    mmX = 2; mmY = 2; mmZ = 2; 
    [Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); 
    voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

    x = reshape(x_hat, size(X));
    sliceViewer(x,'Colormap',hot(256))
end

