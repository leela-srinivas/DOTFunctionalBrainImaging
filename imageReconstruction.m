function [x] = imageReconstruction(y, A, lambda)
    A_reg = inv(A'*A + lambda*eye(size(A,2)))*A';
    x = A_reg*y;
end