%% PI-CS reconstruction example for book chapter in Advanced Neuro MR Techniques and Applications
% 
% It uses a Total Generalized Variation Regularizer described in [1].
% Pseudo-random sampling masks are provided for R=4 and 26 reference
% lines at the center of k-space and R=6 and 16 lines. Coil sensitivities
% were pre-estimated using ESPIRIT [2].
% 
% [1] Knoll F, Bredies K, Pock T, Stollberger, R.
% Second Order Total Generalized Variation (TGV) for MRI.
% Magnetic Resonance in Medicine. 65: 480-491 (2011).
% 
% [2] Uecker M, Lai P, Murphy MJ, Virtue P, Elad M, Pauly JM, Vasanawala SS, Lustig M.
% ESPIRiT - An eigenvalue approach to autocalibrating parallel MRI: Where SENSE meetsGRAPPA.
% Magnetic Resonance in Medicine 71:, 990?1001 (2014).
% 
% June 2020
% Florian Knoll (florian.knoll@nyumc.org)

clear all; close all; clc;

addpath('./utils');
addpath('./tgv');

%% Load Data
datapath = './data/';
load([datapath, 'rawdata17']);
load([datapath, 'espirit17']);
load([datapath, 'mask_random4_768_396_nRef26']); acc = 4; refLines = 26;
% load([datapath, 'mask_random6_768_396_nRef16']); acc = 6; refLines = 16;
nCh = size(rawdata,3);
img = ifft2c(rawdata);
img_sos = sosComb(img);

%% TGV Parameters
reduction = 2^(-8);     % usually there is no need to change this
alpha = 5e-6/reduction; % usually there is no need to change this 
maxit = 1000;            % 500 use 1000 Iterations for optimal image quality
senseIter = 100;        % iterations for sensitivity estimation

%% Generate mask
K  = @(x) 1/sqrt(length(x(:))) .* mask .* (fftshift(fft2(fftshift(x))));
Kh = @(x) sqrt(length(x(:)))   .* fftshift(ifft2(fftshift(mask.*x)));
K_ref  = @(x) 1/sqrt(length(x(:))) .* mask_ref .* (fftshift(fft2(fftshift(x))));
Kh_ref = @(x) sqrt(length(x(:)))   .* fftshift(ifft2(fftshift(mask_ref.*x)));

%% Normalize data
rawdata = rawdata./max(abs(rawdata(:)));

for ii = 1:nCh
    img_full(:,:,ii) = fftshift(ifft2(fftshift(rawdata(:,:,ii))));
end
img_full = sqrt(sum(abs(img_full).^2,3));

%% TGV2 Reconstruction
tic
imgTGV2 = tgv2_l2_2D_pd(sensitivities, rawdata, K, Kh, 2*alpha, alpha, maxit, reduction);
disp(['Time: ', num2str(toc/60), ' min']);

%% Flip and remove OS
imgTGV2 = fliplr(flipud(imgTGV2));
imgTGV2 = imgTGV2(1+nFE/4:3*nFE/4,:);

%% Display Results
imgTGV2 = normalize01(imgTGV2);
figure,imshow(brighten(imgTGV2,0.2),[0,0.75]);

/
