%% GRAPPA reconstruction example for book chapter in Advanced Neuro MR Techniques and Applications
% 
% GRAPPA reconstruction based on [1]. Regular equidistant sampling masks are provided for
% R=4 and 26 reference lines at the center of k-space and R=6 and 16 lines.
% This didactic example is based on the Grappa implementation that is
% provided as part of Miki Lustig's code package for SPIRiT [2].
% 
% [1] Griswold MA, Jakob PM, Heidemann RM, Nittka M, Jellus V, Wang J,
% Kiefer B, Haase A. Generalized Autocalibrating Partially Parallel
% Acquisitions (GRAPPA). Magnetic Resonancein Medicine 47: 1202-1210 (2002).
% 
% [2] Lustig M, Pauly J. SPIRiT: Iterative Self-Consistent Parallel Imaging
% Reconstruction from Arbitrary k-Space Sampling. Magnetic Resonance in Medicine.
% 64:457-471 (2010). 
% 
% The code is available from Miki's webpage:
% https://people.eecs.berkeley.edu/~mlustig/softwaregrappa/SPIRiT_v0.1.tar.gz
% The following functions from this package need to be added to the Matlab path:
% calibrate.m, corrMatrix.m, getCalibSize.m, GRAPPA.m
% 
% June 2020
% Florian Knoll (florian.knoll@nyumc.org)
 
clear all; close all; clc;

addpath('./utils');
addpath('./grappa/');

%% Load data
datapath = './data/';
load([datapath, 'rawdata17']);
load([datapath, 'mask_ipat4_768_396_nRef26']); acc = 4; refLines = 26;
% load([datapath, 'mask_ipat6_768_396_nRef16']); acc = 6; refLines = 16;
pattern = mask; clear mask;
nCh = size(rawdata,3);
img = ifft2c(rawdata);
img_sos = sosComb(img);

%% GRAPPA reconstruction Parameters
kSize = [5,5]; % GRAPPA kernel size
CalibTyk = 0;  % Tikhonov regularization in the calibration: 0.01

%% Subsampling
rawdata_subs = rawdata.*repmat(pattern,[1,1,nCh]);
R_true = size(rawdata_subs,1)*size(rawdata_subs,2) / size(find(rawdata_subs(:,:,1)),1);
% figure, imshow(pattern);
[CalibSize, dcomp] = getCalibSize(pattern); % get size of calibration area from mask

%% Scale data such that zero-filled density compensated k-space norm is 1 
DATAcomp = rawdata_subs.*repmat(dcomp,[1,1,nCh]);
scale_fctr = norm(DATAcomp(:))/sqrt(nCh)/20;
rawdata_subs = rawdata_subs/scale_fctr;
% figure,kshow(rawdata_subs(:,:,3));

%% GRAPPA reconstruction
disp('-------------------------');
disp('GRAPPA reconstruction');
disp('-------------------------');
tic;
kCalib = crop(rawdata_subs,[CalibSize,nCh]);
res_grappa = GRAPPA(rawdata_subs,kCalib,kSize,CalibTyk);
disp(['Time: ', num2str(toc/60), ' min']); 
im_grappa = ifft2c(res_grappa);
reconGrappaSos = sosComb(im_grappa);

%% Flip and remove readout oversampling
reconGrappaSos = fliplr(flipud(reconGrappaSos));
reconGrappaSos = reconGrappaSos(1+nFE/4:3*nFE/4,:);

%% Display results
reconGrappaSos = normalize01(reconGrappaSos);
figure,imshow(brighten(reconGrappaSos,0.2),[0,0.75]);
