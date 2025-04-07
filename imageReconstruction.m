function [x] = imageReconstruction(y, A, lambda)
    A_reg = inv(A'*A + lambda*eye(size(A,2)))*A';
    x_hat = A_reg*y;

    % Define bounds on medium
    xBnds = [-70 70]; yBnds = [-30 30]; zBnds = [1 10];  
    mmX = 2; mmY = 2; mmZ = 2; 
    [Y X Z] = meshgrid(yBnds(1):mmY:yBnds(2), xBnds(1):mmX:xBnds(2), zBnds(1):mmZ:zBnds(2)); 
    voxCrd = double([X(:) Y(:) Z(:)]); % coordinates (mm) for each voxel, reshaped as 1D vector

    x = reshape(x_hat, size(X));
end