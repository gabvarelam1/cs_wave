function [ data ] = ift3( data )
%2DFT Summary of this function goes here
%   Detailed explanation goes here

%need to implement faster 3d fft
data = fftshift(ifft(ifft2(ifftshift(data)), [], 3));
% datasize = size(data);
% n = 1;
% if datasize>3
%     n = prod(datasize(4:end));
% end
% data = ifftshift(data);   
% for k = 1:n
%     data(:,:,:,k) = ifftn(data(:,:,:,k));
% end
% 
% data = fftshift(data);

end

    