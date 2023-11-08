function [DyT] = getC(D1)
[m1,n1] = size(D1);
% compute fixed quantities
 DyT(:,1:n1-1)=diff(D1,1,2);
 DyT(:,n1)=D1(:,1)-D1(:,n1);
end