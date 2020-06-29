%% SENSE reconstruction example for book chapter in Advanced Neuro MR Techniques and Applications
% 
% Uses an implementation of Klaas Pruessman's iterative SENSE method [1].
% Regular equidistant sampling masks are provided for R=4 and 26 reference
% lines at the center of k-space and R=6 and 16 lines. Coil sensitivities
% were pre-estimated using ESPIRIT [2].
% 
% [1] Pruessmann KP, Weiger M, Boernert P, Boesiger P.
% Advances in sensitivity encoding with arbitrary k-space trajectories.
% Magnetic Resonance in Medicine 46: 638-651 (2001).
% 
% [2] Uecker M, Lai P, Murphy MJ, Virtue P, Elad M, Pauly JM, Vasanawala SS, Lustig M.
% ESPIRiT - An eigenvalue approach to autocalibrating parallel MRI: Where SENSE meetsGRAPPA.
% Magnetic Resonance in Medicine 71:, 990?1001 (2014).

% June 2020
% Florian Knoll (florian.knoll@nyumc.org)

clear all; close all; clc;

addpath('./utils');
addpath('./sense');

%% Load data
datapath = './data/';
load([datapath, 'rawdata17']);
load([datapath, 'espirit17']);
load([datapath, 'mask_ipat4_768_396_nRef26']); acc = 4; refLines = 26;
% load([datapath, 'mask_ipat6_768_396_nRef16']); acc = 6; refLines = 16;
nCh = size(rawdata,3);
img = ifft2c(rawdata);
img_sos = sosComb(img);

%% CG SENSE reconstuction parameters
alpha = 0;      % Tikhonov regularization: Unused
tol   = 1e-4;   % CG tolerance
maxitCG = 1000; % Number of CG iterations: Obsolete, iterations are defined via tolerance

%% Subsampling
rawdata_subs = rawdata.*repmat(mask,[1,1,nCh]);
R_true = size(rawdata_subs,1)*size(rawdata_subs,2) / size(find(rawdata_subs(:,:,1)),1);

%% CGSENSE Reconstruction
tic
cgsenseRecon = pmri_cgsense(rawdata_subs,mask,sensitivities,zeros(nPE,nFE),alpha,tol,maxitCG,1);
disp(['Elapsed time: ', num2str(toc/60), ' min']);

%% Flip and remove readout oversampling
cgsenseRecon = fliplr(flipud(cgsenseRecon));
cgsenseRecon = cgsenseRecon(1+nFE/4:3*nFE/4,:);

%% Display results
cgsenseRecon = normalize01(cgsenseRecon);
figure,imshow(brighten(cgsenseRecon,0.2),[0,0.75]);


