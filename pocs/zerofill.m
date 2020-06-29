function [ zfill ] = zerofill( kdata,PF )
%Zero fill for PF data
[nkx nky]=size(kdata);
nky1 = ceil(nky/PF)+1;
zfill = zeros(nkx, nky1);
zfill (:,1:nky) = kdata;
end

