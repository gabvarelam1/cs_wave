function [out] = fft3(data)
%FFT3 Summary of this function goes here
%   Detailed explanation goes here
out = fft(fft2(data), [], 3);
end

