function [sosImg] = sosComb(arrayCoilImg)
%  SOS combination of coil array data
%========================================================================
%  [SOSIMG] = SOSCOMB(ARRAYCOILIMG)
%========================================================================
%  INPUT: arrayCoilImg: 3D Matrix with array coil data
%         
%  OUTPUT: sosImg: Sum of squares combination

sosImg = sqrt(sum(abs(arrayCoilImg).^2,3));


