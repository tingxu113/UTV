function [DyD1,DyTDyD1] =diff2(D1)
% compute fixed quantities
 [m1,n1] = size(D1);
 DyD1(1:m1-1,:)=diff(D1,1,1);
 DyD1(m1,:)=D1(1,:)-D1(m1,:);
 DyTDyD1(:,1:n1-1)=diff(D1,1,2);
 DyTDyD1(:,n1)=D1(:,1)-D1(:,n1);
%  DyT(:,1:n1-1)=diff(D1,1,2);
%  DyT(:,n1)=D1(:,1)-D1(:,n1);
% otfDx = psf2otf([1,-1],sizeI);
% otfDy = psf2otf([1;-1],sizeI);
% conjoDx = conj(otfDx);
% conjoDy = conj(otfDy);
% Denom1 = abs(otfDx).^2;
% Denom2 = abs(otfDy).^2;
end
