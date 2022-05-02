function Nic = filtfilt_3D(ncc, opt, N1, N2, N3, Nc, R)
% Fast computation of zero phase filtering (for alg=4)
fltlen = size(ncc,2)/Nc;    % filter length
numflt = size(ncc,1);       % number of filters

% LORAKS kernel is circular.
% Following indices account for circular elements in a square patch
% [in1,in2, in3] = meshgrid(-R:R,-R:R, -R:R);
% idx = find(in1.^2+in2.^2<=R^2);
% in1 = in1(idx)';
% in2 = in2(idx)';

[in1,in2,in3] = meshgrid(-R:R,-R:R,-R:R);
idx = find(in1.^2+in2.^2+in3.^2<=R^2);
in1 = in1(idx)';
in2 = in2(idx)';
in3 = in3(idx)';


ind = sub2ind([2*R+1, 2*R+1, 2*R+1],R+1+in1,R+1+in2,R+1+in3);

filtfilt = zeros([(2*R+1)*(2*R+1)*(2*R+1),Nc,numflt],'like',ncc);
filtfilt(ind,:,:) = reshape(permute(ncc,[2,1]),[fltlen,Nc,numflt]);
filtfilt = reshape(filtfilt,(2*R+1),(2*R+1),(2*R+1),Nc,numflt);

cfilt = conj(filtfilt);

if opt == 'S'       % for S matrix
    ffilt = conj(filtfilt);
else                % for C matrix
    ffilt = flip(flip(flip(filtfilt,1),2),3);
end

ccfilt = fft(fft2(cfilt,4*R+1, 4*R+1),4*R+1,3);
fffilt = fft(fft2(ffilt,4*R+1, 4*R+1),4*R+1,3);

clear('cfilt','ffilt','filtfilt');

patch = ifft3(sum(reshape(ccfilt,4*R+1,4*R+1,4*R+1,1,Nc,numflt) ...
    .* reshape(fffilt,4*R+1,4*R+1,4*R+1,Nc,1,numflt),6));

clear('ccfilt','fffilt');


if opt == 'S'       % for S matrix
%     Nic = fft3(circshift(padarray(patch, [N1-1-2*R N2-1-2*R N3-1-2*R],'post'),[-4*R-rem(N1,2) -4*R-rem(N2,2) -4*R-rem(N3,2)]));
    Nic = circshift(padarray(patch, [N1-1-2*R N2-1-2*R N3-1-2*R],'post'),[-4*R-rem(N1,2) -4*R-rem(N2,2) -4*R-rem(N3,2)]);
    clear('patch');
    for k = 1:Nc
        Nic(:,:,:,:,k) = fft3(Nic(:,:,:,:,k));
    end

    
else                % for C matrix
    Nic = circshift(padarray(patch, [N1-1-2*R N2-1-2*R N3-1-2*R], 'post'),[-2*R -2*R -2*R]);
    clear('patch');
    for k = 1:Nc
        Nic(:,:,:,:,k) = fft3(Nic(:,:,:,:,k));
    end
%     Nic = fft3(circshift(padarray(patch, [N1-1-2*R N2-1-2*R N3-1-2*R], 'post'),[-2*R -2*R -2*R]));
end
end
