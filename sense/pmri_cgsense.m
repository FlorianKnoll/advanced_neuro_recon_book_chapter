function  us = pmri_cgsense(data,pattern,c,x0,alpha,tol,maxit,matlabCG,display)
%
% us = pmri_cgsense_random(data,pattern,c,x0,alpha,tol,maxit,matlabCG,display)
% reconstruct subsampled PMRI data using CG SENSE [1]
%
% INPUT
% data:    3D array of coil images (in k-space)
% pattern: mask of acquired Fourier coefficients
% c:       3D array of coil sensitivities
% x0:      initial guess for image
%
% OUTPUT
% us:      reconstructed image
%
% Christian Clason and Florian Knoll, June 2014
% 
% [1] Pruessmann, K. P.; Weiger, M.; Boernert, P. and Boesiger, P.
% Advances in sensitivity encoding with arbitrary k-space trajectories.
% Magn Reson Med 46: 638-651 (2001)
%
%
% =========================================================================

%% set up parameters and operators
[nx,ny,nc] = size(data);

% centered FFT, IFFT
cfft  = @(x) fftshift(fft2(fftshift(x)))./sqrt(nx.*ny);
cifft = @(x) fftshift(ifft2(fftshift(x))).*sqrt(nx.*ny);

% sampling operator
F  = @(x) pattern.*cfft(x);
FH = @(x) cifft(pattern.*x);

cscale = 0;
for i = 1:nc
    cscale = cscale+abs((c(:,:,i))).^2;
end
cscale = sqrt(cscale);

%% CG iteration
% precompute complex conjugates
cbar = conj(c);

% right hand side: -K^*residual 
y  = zeros(nx,ny);
for i = 1:nc
    y = y + FH(data(:,:,i)).*cbar(:,:,i);
end

% system matrix: F'^T*F' + alpha I
M  = @(x) applyM(F,FH,c,cbar,x) + alpha*x;

%% CG iterations
if matlabCG
    % Use Matlab CG
    x = pcg(M,y(:),tol,maxit);
else
    % Own CG
    x = 0*y(:); r=y(:); p = r; rr = r'*r; it=0;
    while(rr>tol*norm(r))&&(it<maxit)
        Ap = M(p);
        a = rr/(p'*Ap);
        x = x + a*p;
        rnew = r - a*Ap;
        b = (rnew'*rnew)/rr;
        r=rnew;
        rr = r'*r;
        p = r + b*p;
        it=it+1;
        
        if display
            u_it = reshape(x,nx,ny);
            kspace_it = cfft(u_it);
            u_it = cifft(kspace_it);
            figure(99);
            subplot(1,2,1),imshow(abs(u_it),[]); % colorbar;
            title(['Image CG iteration ' num2str(it)]);
            subplot(1,2,2),kshow(kspace_it);
            title(['k-space iteration ' num2str(it)]);
            drawnow;
        end
    end
    disp(['Own CG: iter ', num2str(it), ' tol: ', num2str(rr/norm(r))]);
end

% update iterates
us  = reshape(x,nx,ny);
us = cscale.*us;

% end main function

%% Derivative evaluation

function y = applyM(F,FH,c,cconj,x)
[nx,ny,nc] = size(c);
dx = reshape(x,nx,ny);

y  = zeros(nx,ny);
for i = 1:nc
    y = y + cconj(:,:,i).*FH(F(c(:,:,i).*dx));
end
y = y(:);
% end function applyM
