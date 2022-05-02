function [ data ] = ft3( data )
%2DFT Summary of this function goes here
%   Detailed explanation goes here


data = fftshift(fft(fft2(ifftshift(data)), [], 3));

% datasize = size(data);
% n = 1;
% if datasize>3
%     n = prod(datasize(4:end));
% end
% data = ifftshift(data);   
% for k = 1:n
%     data(:,:,:,k) = fftn(data(:,:,:,k));
% end
% 
% data = fftshift(data);

end

