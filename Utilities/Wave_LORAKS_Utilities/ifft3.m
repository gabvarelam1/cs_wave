function [out] = ifft3(data)
%FFT3 Summary of this function goes here
%   Detailed explanation goes here

  
out = ifft(ifft2(data), [], 3) ;
end

