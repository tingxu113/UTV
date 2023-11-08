function E=DIC_CG1(  Y,E0,E,B,Z, sf,s0,fft_B,A,mu1,mu2,V111,D111,Dy )

% optimize the dictionaries via Conjugate gradient 
 maxIter =40;
%  P1=P'*P;

A1=A*A';
B1=B*B'+mu1*eye(size(B,1));
H=PTx( Z,fft_B,s0,sf )*A'+mu1*E0+Y*B'+mu2*Dy'*(V111+D111);

r0=H-E*B1-PTPx( E,fft_B,s0,sf )*A1-mu2*Dy'*Dy*E;
p0=r0;
 for i=1:maxIter
    pp= (p0*B1+PTPx( p0,fft_B,s0,sf )*A1);
    pp1=p0(:)'*pp(:);

    a=(r0(:)')*r0(:)/pp1;
    E=E+a*p0;
    r1=r0-a*pp;
         
    b1=(r1(:)'*r1(:))/(r0(:)'*r0(:));
    p1=r1+b1*p0;
    p0=p1;
    r0=r1;
 end
 
 
