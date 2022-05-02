function res = ifft3c(x)
%%% This function does the FFT3 fourier transform.
%%% The idea is to be explicit and avoid non-fourier-related dimensions
%%% from the data. This is specially thought to number of coils usually
%%% found in the 4th dimension.

%%% by gvm

S = size(x);
fctr = S(1)*S(2)*S(3);

x = reshape(x,S(1),S(2),S(3),prod(S(4:end)));

res = zeros(size(x));
for n=1:size(x,4)
	res(:,:,:,n) = (fctr/sqrt(fctr))*fftshift(ifftn(ifftshift(x(:,:,:,n))));
end

res = reshape(res,S);



end