function HR_HSI= UTV_SR(HSI,MSI,T,BW,BH,downsampling_scale,par,s0,S)
beta=par.beta;
mu1=par.mu1;
mu2=par.mu2;
mu3=par.mu3;
mu4=par.mu4;
lambda1=par.lambda1;
lambda2=par.lambda2;
lambda3=par.lambda3;
lambda4=par.lambda4;

 %%  simulate LR-HSI
 Y_h_bar=hyperConvert2D(HSI);
 HSI1=Unfold(HSI,size(HSI),1);
 HSI2=Unfold(HSI,size(HSI),2);
 HSI3=Unfold(HSI,size(HSI),3);
   %%  simulate HR-MSI
 MSI1=Unfold(MSI,size(MSI),1);
 MSI2=Unfold(MSI,size(MSI),2);
 MSI3=Unfold(MSI,size(MSI),3);
%% inilization D1 D2 D3 C
  [m n]=size(MSI1);
  D=MSI1;
  params.Tdata = 2;            % Number of target sparse vectors
  params.dictsize =par.n1;     
  params.iternum = 100;
  params.DUC =1; 
        
  D1 = trainD(D,MSI1,[],[],params);
       
  params.dictsize =par.n2; 
  D2 = trainD(MSI2,MSI2,[],[],params);
        
  rng(10,'twister')
  D3=vca(Y_h_bar,par.n3);

  D_1=ifft(fft(D1).*repmat(BW,[1 par.n1]));
  D_1=D_1(s0:downsampling_scale:end,:); 
  
  D_2=ifft(fft(D2).*repmat(BH,[1 par.n2]));
  D_2=D_2(s0:downsampling_scale:end,:);
  D_3=T*D3;
  
  D11{1}=D_1;
  D11{2}=D_2;
  D11{3}=D3;
  D22{1}=D1;
  D22{2}=D2;
  D22{3}=D_3;
  C=zeros(size(D1,2),size(D2,2),size(D3,2));
  tol=1e-2;

  C=sparse_tucker2( D11,D22, HSI,MSI, lambda1,C,beta ,mu4,tol);

for i=1:60
%%  update U1
CC=ttm(tensor(C),{D2,D_3},[2 3]); 
CC2=Unfold(double(CC),size(CC),1);
CC3=ttm(tensor(C),{D_2,D3},[2 3]);
CC4=Unfold(double(CC3),size(CC3),1);
D1  = ADMMD1( MSI1, D1, CC2,HSI1,downsampling_scale,s0,BW,CC4,beta,mu1,lambda2,tol);
D_1=ifft(fft(D1).*repmat(BW,[1 par.n1]));
D_1=D_1(s0:downsampling_scale:end,:);
 %%  update U2
CC=ttm(tensor(C),{D1,D_3},[1 3]);
CC2=Unfold(double(CC),size(CC),2);
CC3=ttm(tensor(C),{D_1,D3},[1 3]);
CC4=Unfold(double(CC3),size(CC3),2);
D2= ADMMD1( MSI2, D2, CC2,HSI2,downsampling_scale,s0,BH,CC4,beta,mu2,lambda3,tol);
D_2=ifft(fft(D2).*repmat(BH,[1 par.n2]));
D_2=D_2(s0:downsampling_scale:end,:); 
%%  update U3
CC=ttm(tensor(C),{D_1,D_2},[1 2]);
CC2=Unfold(double(CC),size(CC),3);
CC3=ttm(tensor(C),{D1,D2},[1 2]);
CC4=Unfold(double(CC3),size(CC3),3);
D3 = ADMMD3( HSI3, D3, CC2,MSI3,T,CC4,beta,mu3,lambda4,tol);
D_3=T*D3; 
 %%  update G
 D11{1}=D_1;
 D11{2}=D_2;
 D11{3}=D3;
 D22{1}=D1;
 D22{2}=D2;
 D22{3}=D_3;
 C =  sparse_tucker2( D11,D22, HSI,MSI, lambda1,C,beta,mu4,tol );
    
HR_HSI=ttm(tensor(C),{D1,D2,D3},[1 2 3]);
HR_HSI=double(HR_HSI);

HR_HSI=ttm(tensor(C),{D1,D2,D3},[1 2 3]);
HR_HSI=double(HR_HSI);
end
       


