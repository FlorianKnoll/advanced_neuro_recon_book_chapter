%% Partial Fourier POCS brain reconstruction example for book chapter in Advanced Neuro MR Techniques and Applications
% 
% The script simulates a 5/8 half Fourier acquisition and shows reconstructions
% from simple zero filling, brute force conjugate symmetry reconstruction
% as well as an iterative phase constrained reconstruction using the
% projection on convex sets (POCS) algorithm [1].
% 
% [1] Cuppen J, van Est A. Reducing MR imaging time by one-sided reconstruction. Magnetic Resonance Imaging. 5:526-527 (1987).
% 
% June 2020
% Florian Knoll (florian.knoll@nyumc.org)

clear all; close all; clc;

addpath('./utils');
addpath('./pocs');

%% Load data
datapath = './data/';
load([datapath, 'espirit17']); clear sensitivities;
[nFE,nPE]=size(reference);
reference = reference(1+nFE/4:3*nFE/4,:);
[nFE,~]=size(reference);
reference = reference(:,(nPE-nFE)/2+1:end-(nPE-nFE)/2);
[nFE,nPE]=size(reference);

%% generate a synthetic phase
% This is done for didactic reasons. If you want to see the partial Fourier
% reconstruction without this additional corruption, you can just comment this out.
SNR = 10;
[xx,yy] = meshgrid(linspace(-pi/2,pi/2,nPE));
synthphase = double(sin(4*yy.^2))';
NoiseToAdd = (norm(synthphase(:))/(size(synthphase,1)+size(synthphase,2))/2)/SNR; 
disp(['NoiseToAdd: ',num2str(NoiseToAdd)]);
synthphase = synthphase + NoiseToAdd*randn(size(synthphase));
img_temp =  reference .* exp(-1j*synthphase);
rawdata = fft2c(img_temp);

%% PF
PF = 5/8;
PF_crop = nPE*PF;

rawdata(:,PF_crop:end)=0;
figure; kshow(abs(rawdata));
title('Zero-filled partial Fourier data');
drawnow;

I_zfill = ifft2c(rawdata);
rawdata(:,PF_crop:end)=[];

%% Direct partial Fourier reconstruction
rawdata_crop = rawdata(:,1:nPE-PF_crop+1);
rawdata_conj = flipud(fliplr(conj(rawdata_crop)));
figure; kshow(rawdata_conj);
rawdata_conj = cat(2,rawdata,rawdata_conj);
I_conj = ifft2c(rawdata_conj);

%% Phase estimate
P = EstPhase(rawdata,ceil(nPE*PF)-1-nPE/2); % input data and size of extra lines
figure; imshow(fliplr(flipud(P)),[]); 
title('Phase estimate');

%% POCS
pocsIter = 10;
Ipocs = pf_pocs(rawdata,PF,pocsIter);

%% Calculate difference maps
ErrorZF = (abs(reference) - abs(I_zfill))/max(abs(reference(:)));
ErrorC  = (abs(reference) - abs(I_conj))/max(abs(reference(:)));
ErrorP  = (abs(reference) - abs(Ipocs))/max(abs(reference(:)));

scaling_error_plot = 0.5;

figure; imshow(fliplr(flipud(ErrorZF)),[0,scaling_error_plot*max(ErrorZF(:))]); title('Zero-filling error image'); % colorbar; drawnow;
figure; imshow(fliplr(flipud(ErrorC)),[0,scaling_error_plot*max(ErrorZF(:))]); title('Conjugate-symmetry error image'); % colorbar; drawnow;
figure; imshow(fliplr(flipud(ErrorP)),[0,scaling_error_plot*max(ErrorZF(:))]); title('POCS error image'); % colorbar; drawnow;

figure; imshow(fliplr(flipud(ErrorZF)),[0,scaling_error_plot*max(ErrorZF(:))]); title('Zero-filling error image'); h = colorbar; set(h, 'ylim', [0 0.15]); drawnow;
figure; imshow(fliplr(flipud(ErrorC)),[0,scaling_error_plot*max(ErrorZF(:))]); title('Conjugate-symmetry error image'); h = colorbar; set(h, 'ylim', [0 0.15]); drawnow;
figure; imshow(fliplr(flipud(ErrorP)),[0,scaling_error_plot*max(ErrorZF(:))]); title('POCS error image'); h = colorbar; set(h, 'ylim', [0 0.15]); drawnow;

%% Calculate quantitative evaluation
reference = abs(reference);
reference = normalize01(reference);
I_zfill = normalize01(I_zfill);
I_conj = normalize01(I_conj);
Ipocs = normalize01(Ipocs);

[ssimval_zfill,~] = ssim(reference,I_zfill);
[ssimval_conj,~] = ssim(reference,I_conj);
[ssimval_pocs,~] = ssim(reference,Ipocs);
disp(['SSIM ZF: ', num2str(ssimval_zfill)]);
disp(['SSIM CONJ: ', num2str(ssimval_conj)]);
disp(['SSIM POCS: ', num2str(ssimval_pocs)]);

%% Display reconstructions
fsize = 12;
brightening = 0.4;

figure; imshow(fliplr(flipud(brighten(reference,brightening))),[]); 
title('Reference');
drawnow;

figure; imshow(fliplr(flipud(brighten(I_zfill,brightening))),[]);
title('Zero-filled reconstruction');
h = text(15,nFE-25,sprintf('%.3f',ssimval_pocs));
set(h,'FontSize',fsize,'Horizontalalignment','left','VerticalAlignment','middle','Color', [1 1 1],'FontWeight','bold');
drawnow;

figure; imshow(fliplr(flipud(brighten(I_conj,brightening))),[]);
title('Conjugate-symmetry reconstruction');
h = text(15,nFE-25,sprintf('%.3f',ssimval_conj));
set(h,'FontSize',fsize,'Horizontalalignment','left','VerticalAlignment','middle','Color', [1 1 1],'FontWeight','bold');
drawnow;

figure; imshow(fliplr(flipud(brighten(Ipocs,brightening))),[]); title('POCS');
h = text(15,nFE-25,sprintf('%.3f',ssimval_pocs));
set(h,'FontSize',fsize,'Horizontalalignment','left','VerticalAlignment','middle','Color', [1 1 1],'FontWeight','bold');
drawnow;

