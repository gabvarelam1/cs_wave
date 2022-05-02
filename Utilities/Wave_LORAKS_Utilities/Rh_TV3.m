function [ y ] = Rh_TV3( x, nx, ny, nz )
%R_TV Summary of this function goes here
%   Detailed explanation goes here

    x = x(:);
    x1 = reshape(x(1:nx*ny*nz), [nx ny nz]);
    x2 = reshape(x(nx*ny*nz+1 : 2*nx*ny*nz), [nx ny nz]);
    x3 = reshape(x(2*nx*ny*nz+1 : 3*nx*ny*nz), [nx ny nz]);
    
    
    v1 = x1 - circshift(x1, [1 0 0]);
    v2 = x2 - circshift(x2, [0 1 0]);
    v3 = x3 - circshift(x3, [0 0 1]);
    
    y = v1(:)+v2(:)+v3(:);


    
end

