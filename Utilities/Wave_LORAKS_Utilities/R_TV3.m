function [ y ] = R_TV3( x, nx, ny, nz )
%R_TV Summary of this function goes here
%   Detailed explanation goes here


    x = reshape(x, [nx ny nz]);
    v1 = x - circshift(x, [-1 0 0]);
    v2 = x - circshift(x, [0 -1 0]);
    v3 = x - circshift(x, [0 0 -1]);
    
    y = [ vect(v1); vect(v2); vect(v3)];
%     y = [vect(x-circshift(x,[-1 0 0])); vect(x-circshift(x,[0 -1 0])); vect(x-circshift(x,[0 0 -1]))];
end

