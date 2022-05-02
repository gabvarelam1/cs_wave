function res = fft3c(x)
%%% This function does the FFT3 fourier transform.
%%% The idea is to be explicit and avoid non-fourier-related dimensions
%%% from the data.

%%% by gvm

S = size(x);
fctr = S(1)*S(2)*S(3);

x = reshape(x,S(1),S(2),S(3),prod(S(4:end)));

res = zeros(size(x));
for n=1:size(x,4)
	res(:,:,:,n) = 1/sqrt(fctr)*fftshift(fftn(ifftshift(x(:,:,:,n))));
end

res = reshape(res,S);



end