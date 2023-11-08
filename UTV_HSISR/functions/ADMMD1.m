function AA=ADMMD1(Y,E0, B,Z, sf,s0,fft_B,A,mu1,mu2,lambda1,tol)
AA=zeros(size(E0));
Dd1=zeros(size(E0,1)-1,size(E0,2));
Dy=diffy(E0);
relerr=[];
for i=1:100
    A_old=AA;
    Vv1=soft(Dy*AA-Dd1,lambda1/2*mu2);
    AA=DIC_CG1( Y,E0,AA, B,Z, sf,s0,fft_B,A,mu1,mu2,Vv1,Dd1,Dy);
    Dd1=Dd1+Vv1-Dy*AA;
    relerr(i)=abs(norm(AA(:)-A_old(:))/(norm(A_old)+eps));
    relerr(1)=1;
    if relerr(i)<tol
        break
    end
end
end
    
    