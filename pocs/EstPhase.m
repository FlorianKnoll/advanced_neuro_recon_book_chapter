function [ Y ] = EstPhase(X,a)
% X: kdata
% a: # of extra lines.
[m n]=size(X);
Phase = X(:,n-2*a+1:end);
HW = hamming(m,'symmetric')*hamming(2*a,'symmetric')';
PH = Phase.*HW;
ny = (n-a)*2;
Y = zeros(m,ny);
Y(:,ny/2-a+1:ny/2+a) = PH;
Y = angle(ifft2c(Y));
end
