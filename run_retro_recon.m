clear all; close all; clc;
addpath(genpath('./Utilities'));
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
% Data Parameters %

data_number = 1;
usp_case = 4;  % up to 4
cswlevels = 1; % automatic
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%%
%%% Data
cd('Data')
load(strcat('data',num2str(data_number)));
brain_mask = single(niftiread(strcat('data',num2str(data_number),'_brain_mask.nii.gz')));
cd ..
%%
% reference
ideal = cartesian_img_reference;
[sx,sy,sz] = size(ideal);
maxx = max(abs(ideal(:)))/4; % for visualization
% Add noise
% Add some noise
SNR = 15; % or []

if ~isempty(SNR)
    
    % The original ksp is without additional noise and it can be easily
    % replaced if you want to add more noise in here.
    
    
    ksp_cartesian = fft3c(ideal);
    mu = mean(abs(ksp_cartesian(:)));
    index = 1;
    limit_readout = round(.8 * WaveParams.crop_size/4);
    background_block = zeros(size(ideal));
    background_block(1:limit_readout,1:sy,1:sz) = 1; 
                                         % by looking a data1 number 135
                                          % but take advatange of the
                                          % oversampling for this purpose
    std_noise_range = .1 : 0.01 : 3; % std_noise amplifies the mean of ksp for significant noise level
 
    figure;
    tol = .2; % SNR 10 vs SNR 9.8 -> the same thing
    for std_noise = std_noise_range
        
        noisyKSP = ksp_cartesian + std_noise*mu*( randn(size(ksp_cartesian)) + 1i*randn(size(ksp_cartesian)) );
        noisyIMG = ifft3c(noisyKSP);
        snr(index) = mean(abs(noisyIMG(logical(brain_mask)))) / std(abs(noisyIMG(logical(background_block))));
        
        if abs(snr(index) - SNR) < tol
            break;
        end
        
        plot(snr);
        drawnow;
        index = index+1;
        
    end
    
    [~,pos] = min(abs(SNR-snr));
    std_noise = std_noise_range(pos);
    noisyKSP = ksp_cartesian + std_noise*mu*( randn(size(ksp_cartesian)) + 1i*randn(size(ksp_cartesian)) );
    noisyIMG = ifft3c(noisyKSP);
    snr(index) = mean(abs(noisyIMG(logical(brain_mask)))) / std(abs(noisyIMG(logical(background_block))));

    % k space data
    
    img_test = noisyIMG;
    
else
    img_test = ideal;
end
%%
oversampling_factor = size(WaveParams.PSF,1)/sx;


%% here we create the complete wave-based k-space
sc = size(C,4);
S = ones([sx*oversampling_factor,sy,sz,sc]);
S2d = squeeze(S(1,:,:,1));

figure, imshow(S2d,[])
% 
crop_size = WaveParams.crop_size;
number_channels = WaveParams.number_channels;
pad_size = WaveParams.pad_size;

fft1Dc = @(Ksp, Dim) fftshift( fft( ifftshift(Ksp , Dim) ,[],Dim) ,Dim) / sqrt(size(Ksp,Dim));
ifft1Dc = @(Im, Dim) fftshift( ifft( ifftshift(Im , Dim) ,[],Dim) ,Dim) * sqrt(size(Im, Dim));

iF1d = @(ksp_,dim_) ifft1Dc(ksp_,dim_);
F1d = @(ksp_,dim_) fft1Dc(ksp_,dim_);

% Forward
Cs = @(sensitivity_coils , im_) padarray(repmat(im_,[1,1,1,number_channels]) .* sensitivity_coils, [pad_size,0,0,0]);
PFx = @(Cs_ , wave_point_spread_function) bsxfun(@times, F1d(Cs_,1) , wave_point_spread_function);
Fyz = @(im_ , S_) S_ .* ( F1d( F1d(im_, 2) ,3 ) ) ;


%  Inverse
iFyz = @(ksp_, S_) iF1d( iF1d(ksp_ .* S_, 3) , 2 );
iPFx = @(Im_yz_kspace_x , wave_point_spread_function) ...
    iF1d( bsxfun(@times,Im_yz_kspace_x , wave_point_spread_function) , 1 );
Ct = @(sensitivity_coils , im_) sum( ...
conj(sensitivity_coils) .* im_(1+end/2-crop_size/2:end/2+crop_size/2,:,:,:), 4);

% k space data
y = Fyz(PFx(Cs(C,img_test),WaveParams.PSF),S);

% image from wave
x = Ct( C , iPFx( iFyz( y , S ) , conj(WaveParams.PSF) ) );
%%% NOTE: it should be almost identical to the Cartesian image "ideal"

[mid] = ceil((size(ideal)+1)/2);
    
%%% magnitude
figure, imshow([abs(squeeze(x(mid(1),:,:))),...
    abs(squeeze(ideal(mid(1),:,:)))],[],'init',200);
figure, imshow([abs(squeeze(x(:,mid(2),:))),...
    abs(squeeze(ideal(:,mid(2),:)))],[],'init',200);
figure, imshow([abs(squeeze(x(:,:,mid(3)))),...
    abs(squeeze(ideal(:,:,mid(3))))],[],'init',200); 

%%% phase
figure, imshow([angle(squeeze(x(mid(1),:,:))),...
    angle(squeeze(ideal(mid(1),:,:)))],[],'init',200);
figure, imshow([angle(squeeze(x(:,mid(2),:))),...
    angle(squeeze(ideal(:,mid(2),:)))],[],'init',200);
figure, imshow([angle(squeeze(x(:,:,mid(3)))),...
    angle(squeeze(ideal(:,:,mid(3))))],[],'init',200);
%%
%% Now, we create our undersampled data
%%%

switch usp_case
    case 1
        folder_name = 'uniform_4x3';
        method_usf = 'gr'; % {'gr','cs'} 
        usf_x = 4;
        usf_y = 3;
        calibration_lines = 0;
        elliptical_flag = 0;

        %%%
        %%% Design Undersampling Pattern
        % USPdesign(Method,Ry,Rz,Ncalib, FlagElliptical , Size_ksp)
        [S,S2d,usf_real,dir2save] = USPdesign(method_usf,usf_x,usf_y,...
                                    calibration_lines,elliptical_flag,...
                                    size(ksp));
    case 2
         % uniform Poisson
         folder_name = 'uPoisson_12';
         load mask_poissonVDOFF_Rt12
         figure, imshow(mask,[])
         S = repmat(permute(mask,[3,1,2]),[size(WaveParams.PSF,1),1,1,number_channels]);
     
    case 3
        
         % uniform Poisson
         folder_name = 'vdPoisson_12';
         load Poisson_12x_256x192
         figure, imshow(mask,[])
         S = repmat(permute(mask,[3,1,2]),[size(WaveParams.PSF,1),1,1,number_channels]);
    case 4
        folder_name = 'uniform_3x3';
        method_usf = 'gr'; % {'gr','cs'} 
        usf_x = 3;
        usf_y = 3;
        calibration_lines = 0;
        elliptical_flag = 0;

        %%%
        %%% Design Undersampling Pattern
        % USPdesign(Method,Ry,Rz,Ncalib, FlagElliptical , Size_ksp)
        [S,S2d,usf_real,dir2save] = USPdesign(method_usf,usf_x,usf_y,...
                                    calibration_lines,elliptical_flag,...
                                    size(ksp));
end
%%
mask = squeeze(S(1,:,:,1));
y_us = y.*S;

x_sub = Ct( C , iPFx( iFyz( y_us , S ) , conj(WaveParams.PSF) ) );

%%% magnitude

figure, imshow([abs(squeeze(x_sub(mid(1),:,:))),...
    abs(squeeze(ideal(mid(1),:,:)))],[],'init',200);
figure, imshow([abs(squeeze(x_sub(:,mid(2),:))),...
    abs(squeeze(ideal(:,mid(2),:)))],[],'init',200);
figure, imshow([abs(squeeze(x_sub(:,:,mid(3)))),...
    abs(squeeze(ideal(:,:,mid(3))))],[],'init',200); 

%%% phase
figure, imshow([angle(squeeze(x_sub(mid(1),:,:))),...
    angle(squeeze(ideal(mid(1),:,:)))],[],'init',200);
figure, imshow([angle(squeeze(x_sub(:,mid(2),:))),...
    angle(squeeze(ideal(:,mid(2),:)))],[],'init',200);
figure, imshow([angle(squeeze(x_sub(:,:,mid(3)))),...
    angle(squeeze(ideal(:,:,mid(3))))],[],'init',200);

%%%  Reconstruction Algorithms
%%% 
%%
%% Run CS-Wave
%%
%% Contribution from Gabriel Varela-Mattatall
%% Please cite the following if you refer to this section
%% Varela-Mattatall et al (MRM, 2021) 
%% "Automatic determination of the regularization weighting for wavelet-based
%% compressed sensning MRI reconstructions "

%%% Please consider citing as well:

%%% Ong, (MRM, 2018)
%%% "General Phase Regularized Reconstruction Using Phase Cycling"

wname = 'db4';
%%% NEED TO ANALZY THE INPUT
%%% WHAT HAPPENS TO K-SPACE and all stuff....
[recon, recons_all, optKit] = WaveCS_mFISTA3(y_us,S,C,x_sub,WaveParams,...
                                             wname,'Lambda','multiple','Levels',1);

csWave = struct;
csWave.ReconFinal = recon;
csWave.ReconAll = recons_all;
csWave.Regularization_weighting = optKit.Lambda;
csWave.Reconstruction = recons_all(:,:,:,1);
csWave.ReconTime = sum(optKit.Iter_time(1:50));
csWave.OptKit = optKit;
%% this is to check results and its going to be replaced at the end of the script

NRMSE = @(recon) norm(recon(:)-ideal(:)) / norm(ideal(:));
nrmse_vec_image = [NRMSE(csWave.Reconstruction),NRMSE(csWave.ReconFinal)]

vec = @(a) a(:);
brain_mask = logical(brain_mask);
NRMSE = @(recon) norm(recon(brain_mask)-ideal(brain_mask)) / norm(ideal(brain_mask));
nrmse_vec_brain = [NRMSE(csWave.Reconstruction),NRMSE(csWave.ReconFinal)]


%% Run Wave-CAIPI and Wave-LORAKS

%% Contribution from Berkin Bilgic and Tae Hyung Kim
%% Please cite the following if you refer to these sections
%% Bilgic et al (MRM, 2015) "Wave-CAIPI for highly accelerated 3D imaging"
%% Kim et al (MRM, 2019) "Wave-LORAKS: 
%% Combining Wave encoding with low-rank matrix modeling for 
%% more highly accelerated 3D imaging"
%%
%%% Please consider citing as well the seminar paper of LORAKS

%%% Haldar, (IEEE on Medical Imaging, 2014)
%%% "Low-Rank Modeling of Local k-Space Neighborhoods (LORAKS) 
%%% for Constrained MRI"

%%
%%%% THINGS YOU WANT TO CHANGE
%%%%%%%%%
%%%%%%%%%

Radius = 4;
rank_S = 2000; 
lambda = 1e-1;

%%%%%%%%%
%%%%%%%%%
%%%%
os_factor = oversampling_factor;      % oversampling factor in readout

[nx,ny,nz,~]    =   size(ideal);

nx_os = nx * os_factor;

psf_wav = WaveParams.PSF;
receive = C;                    clear C;

%--------------------------------------------------------------------------
%% undersampled wave k-space data
%--------------------------------------------------------------------------

% find: img, fftcn, ifftc, nx_os, nx, psf_wav, receive, num_chan, M

%kspace_wave = fftc(fftc(fftc(padarray(img, (nx_os-nx)/2), 1) .* repmat(psf_wav, [1,1,1,size(receive,4)]), 2), 3);
kspace_wave = y_us;

[M(1), M(2), M(3), num_chan] = size(kspace_wave);

M3d = repmat(permute(mask, [3,1,2]), [M(1),1,1]);



%%
OS = os_factor; % oversampling factor in readout


% load data  
kData = kspace_wave; clear('kspace_wave');
psf = single(psf_wav); clear('psf_wav');
coil_sens = receive; clear('receive');
kMask = M3d; clear('M3d')


%% Reference from fully sampled k-space data
[nx, ny, nz, nc] = size(kData);
nx = nx/OS;


%% Define Fourier Transforms
ft_x = @(x) fftshift(fft(ifftshift(x, 1), [], 1), 1) / sqrt(size(x,1));
ft_yz = @(x) fftshift(fft(fft(ifftshift(x),[],2),[],3)) / sqrt(size(x,2) * size(x,3));
ft3 = @(x) fftshift(fft(fft2(ifftshift(x)), [], 3)) / sqrt(size(x,1)*size(x,2)*size(x,3));

ift_x = @(x) fftshift(ifft(ifftshift(x, 1), [], 1), 1) * sqrt(size(x,1));
ift_yz = @(x) fftshift(ifft(ifft(ifftshift(x),[],2),[],3)) * sqrt(size(x,2) * size(x,3));
ift3 = @(x) fftshift(ifft(ifft2(ifftshift(x)), [], 3))  * sqrt(size(x,1)*size(x,2)*size(x,3));



%%
% zeropadding 
ZP = zeros([nx*OS ny nz nc], 'single');
ZP(1+end/2-end/2/OS:end/2+end/2/OS,:,:,:,:) = 1;
Smtx.type='()';
Smtx.subs{:} = find(vect(ZP));
  
fig_num = 0;

imagesc3d2(ideal, [nx,ny,nz]/2, fig_num+1, [180,0,0], [0,3e-7], 0, 'R1')


%% Wave-CAIPI

% Wave-Encoding Forward models
% kMask = repmat(kMask,[1,1,1,num_chan]);
%Ah = @(x)  vect(sum(conj(coil_sens).*reshape(subsref(ift_x(ift_yz(reshape(x,[nx*OS ny nz nc]) .* repmat(kMask,[1,1,1,nc])).*repmat(conj(psf),[1,1,1,nc])),S),[nx ny nz nc]),4)); 
Ah = @(x)  vect(sum(conj(coil_sens).*reshape(subsref(ift_x(ift_yz(reshape(x,[nx*OS ny nz nc]) .* kMask).*conj(psf)),Smtx),[nx ny nz nc]),4)); 
AhA = @(x) vect(sum(conj(coil_sens).*reshape(subsref(ift_x(ift_yz( ft_yz(ft_x(subsasgn(ZP,Smtx,coil_sens.*reshape(x,[nx ny nz 1]))).*psf).*kMask).*conj(psf)),Smtx),[nx ny nz nc]),4));


Ahd = Ah(kData(:)); 
clear('kData');

%% WAVE-CAIPI (CG iterations -- Same as above block using pcg function, CG algorithm is hard-coded)

% The following steps are equivalent to 
% x = pcg(AhA, Ahd)


max_iter = 20;
tol = 1e-4;

x = vect(zeros(nx, ny, nz));


% Conjugate gradient algorithm
tic;
r = Ahd - AhA(x);
p = r;
rsold = r'*r;
%%
for k = 1:max_iter
    z = AhA(p);
    alpha = rsold / (p'*z);
    x = x + alpha * p;
    r = r - alpha *z;
    rsnew = r'*r;
    if sqrt(rsnew) < tol
        break;
    end
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
    
    NRMSE = (norm(x(:) - ideal(:)) / norm(ideal(:)));
    disp(['iter: ' num2str(k) ', err: ' num2str(norm(x(:) - ideal(:)) / norm(ideal(:))) ', residue: ' num2str(sqrt(rsnew))])
    imagesc3d2(reshape(x, [nx ny nz]), [nx,ny,nz]/2, fig_num+20, [180,0,0], [0,3e-5], 0, ['WAVE-CAIPI iter' num2str(k)])
    

end
recontime_wave_caipi = toc;
recon = reshape(x, [nx ny nz]);
wave_caipi = recon;

wCaipi = struct;
wCaipi.Reconstruction = wave_caipi;
wCaipi.ReconTime = recontime_wave_caipi;

imagesc3d2(wave_caipi, [nx,ny,nz]/2, fig_num+2, [180,0,0], [0,3e-4], 0, 'WAVE-CAIPI')
NRMSE = (norm(x(:) - ideal(:)) / norm(ideal(:)));

% save('Wave_CAIPI.mat','wave_caipi', 'NRMSE');


%% Calibration data for Wave-LORAKS


nx_ac = 64;      %ACS region size x
ny_ac = 64;      %ACS region size y
nz_ac = 64;      %ACS region size z

x_ac_idx = (floor(nx/2)+1 - floor(nx_ac/2)) : (floor(nx/2)+1 + floor(nx_ac/2) - ~rem(nx_ac,2));
y_ac_idx = (floor(ny/2)+1 - floor(ny_ac/2)) : (floor(ny/2)+1 + floor(ny_ac/2) - ~rem(ny_ac,2)) ;
z_ac_idx = (floor(nz/2)+1 - floor(nz_ac/2)) : (floor(nz/2)+1 + floor(nz_ac/2) - ~rem(nz_ac,2)) ;

% Radius for LORAKS kernel

[P_Cac, Ph_Cac, ~, ~, P_Sac Ph_Sac ~, ~, ~, sizeCa, ~, sizeSa] = generate_LORAKS_operators_3D(nx_ac, ny_ac, nz_ac, Radius);


% Calibration data from WAVE-CAIPI recon (central region is relatively accurate)
k_calib = ft3(coil_sens.*reshape(wave_caipi,[nx ny nz]));
acData = double(k_calib(x_ac_idx, y_ac_idx, z_ac_idx, :));

ZD = @(x) padarray(reshape(x,[nx ny nz nc]),[2*Radius, 2*Radius, 2*Radius], 'post');
ZD_H = @(x) x(1:nx,1:ny,1:nz,:);
% 
B = @(x) vect(ft3(coil_sens.*reshape(x, [nx ny nz])));           % SENSE encoding
Bh = @(x)  vect(sum(conj(coil_sens).*ift3(reshape(x,[nx ny nz nc])),4));         % Adjoint of SENSE encoding

% lambda = 1e-1;  % adjust

% Regularization 
LhL = @(x) Bh(ZD_H(ifft3(squeeze(sum(Nic.*fft3(ZD(B(x))),4)))));


%% Calibration for Wave-LORAKS (S-version: has phase constraint)
clear('Nic','Nis','phi','LhL')


tic
if rank_S>0
%     Cac = [];
    Sac = zeros(sizeSa(1)*nc, sizeSa(2),'single');
    for k = 1:nc
        Sac(sizeSa(1)*(k-1)+1:sizeSa(1)*k,:) = single(P_Sac(decomplexify(vect(acData(:,:,:,k)))));
    end
    Us = svd_left(Sac); clear('Sac');
    nss = Us(:,rank_S+1:end)';%clear('Us');

    nf = size(nss,1);
    patchSize= size(nss,2) / (2*nc);
    nss = reshape(nss,[nf, patchSize, 2*nc]);
    nss_h = reshape(nss(:,:,1:2:end)+1j*nss(:,:,2:2:end),nf, []);

    Nic = filtfilt_3D(nss_h, 'C', nx, ny, nz, nc, Radius);
    Nis = filtfilt_3D(nss_h, 'S', nx, ny, nz, nc, Radius);

end
calibtime = toc;
disp(calibtime);

%% Wave-LORAKS (S-version: has phase constraint)

% Regularization 
LhL = @(x) 2*Bh( (ZD_H(ifft3(squeeze(sum(Nic.*fft3(ZD(B(x))),4))))) ...
    -(ZD_H(ifft3(squeeze(sum(Nis.*conj(fft3(ZD(B(x)))),4))))));

% Wave-LORAKS model (lambda tuning)
phi = @(x) AhA(x) + lambda*LhL(x);

max_iter = 40;
tol = 1e-4;

% x = vect(zeros(nx, ny, nz));
x = vect(wave_caipi);        % Initialization with WAVE-CAIPI recon


tic;
r = Ahd - phi(x);
p = r;
rsold = r'*r;

for k = 1:max_iter
    z = phi(p);
    alpha = rsold / (p'*z);
    x = x + alpha * p;
    r = r - alpha *z;
    rsnew = r'*r;
    if sqrt(rsnew) < tol
        break;
    end
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
    
    NRMSE = (norm(x(:) - ideal(:)) / norm(ideal(:)));
    disp(['iter: ' num2str(k) ', err: ' num2str(NRMSE) ', residue: ' num2str(sqrt(rsnew))])
    imagesc3d2(reshape(x, [nx ny nz]), [nx,ny,nz]/2, fig_num+20, [180,0,0], [0,3e-4], 0, ['WAVE-LORAKS (S) iter' num2str(k)])
    
        
end
recontime_loraks_s = toc;
recon3 = reshape(x, [nx ny nz]);
recon_LORAKS_S = recon3;

imagesc3d2(recon_LORAKS_S, [nx,ny,nz]/2, fig_num+4, [180,0,0], [0,3e-4], 0, 'WAVE-LORAKS (S)')
NRMSE = (norm(x(:) - ideal(:)) / norm(ideal(:)));

% save('LORAKS_S.mat','LORAKS_S', 'NRMSE','r_S');
%%
wLoraks = struct;
wLoraks.Radius = Radius;
wLoraks.Rank_S_matrix = rank_S;
wLoraks.Regularization_weighting = lambda;
wLoraks.Reconstruction = recon_LORAKS_S;
wLoraks.ReconTime = recontime_loraks_s;

%% 
cd('Reconstructions');

mkdir(folder_name);
cd(folder_name);
date_exe = date;

NRMSE = @(recon) norm(recon(:)-ideal(:)) / norm(ideal(:));
nrmse_vec_order = {'loraks','caipi','cs_15iter','cs_50iter'};
nrmse_vec_image = [NRMSE(wLoraks.Reconstruction),NRMSE(wCaipi.Reconstruction),NRMSE(csWave.Reconstruction),NRMSE(csWave.ReconFinal)]

vec = @(a) a(:);
brain_mask = logical(brain_mask);
NRMSE = @(recon) norm(recon(brain_mask)-ideal(brain_mask)) / norm(ideal(brain_mask));
nrmse_vec_brain = [NRMSE(wLoraks.Reconstruction),NRMSE(wCaipi.Reconstruction),NRMSE(csWave.Reconstruction),NRMSE(csWave.ReconFinal)]
%%
save('Reconstructions','date_exe','wCaipi','wLoraks','csWave','cartesian_img_reference','mask','nrmse_vec_brain','nrmse_vec_image','nrmse_vec_order','-v7.3');
cd ..
cd ..

figure, imshow([angle(squeeze(wLoraks.Reconstruction(:,:,97))),angle(squeeze(wCaipi.Reconstruction(:,:,97))),angle(squeeze(csWave.Reconstruction(:,:,97))),angle(squeeze(csWave.ReconFinal(:,:,97))),angle(squeeze(ideal(:,:,97)))],[ ])
%%
figure, imshow([abs(squeeze(wLoraks.Reconstruction(:,:,97))),abs(squeeze(wCaipi.Reconstruction(:,:,97))),abs(squeeze(csWave.Reconstruction(:,:,97))),abs(squeeze(csWave.ReconFinal(:,:,97))),abs(squeeze(ideal(:,:,97)))],[0,max(abs(ideal(:)))/4])
figure, imshow(mask,[],'init',400);
%%

