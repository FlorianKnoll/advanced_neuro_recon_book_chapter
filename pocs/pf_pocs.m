function [I] = pf_pocs(kdata,PF,Nite)
% POCS reconstruction for partial Fourier (PF) MRI
% kdata: asymmetric k-space data
% N: size of the reconstructed PF dimension % Nite: number of iterations
% I: reconstructed magnitude image
[x1 y1] = size(kdata);
N = ceil(y1/PF)+1;
% N = y1/PF;
a = y1-N/2;
P = EstPhase(kdata,a);
Sk0 = zerofill(kdata,PF);
I0 = ifft2c(Sk0);
It = I0;
for t = 1:Nite
    It = abs(It).*exp(1i*P);
    Skt = fft2c(It);
    Skt(:,1:y1) = kdata;
    It = ifft2c(Skt);
end
I = It;
end