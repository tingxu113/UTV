% The method for Hyperspectral image super-resolution
% Copyright(c) 2020 T. Xu
% All Rights Reserved.
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for Hyperspectral image super-
% resolution from a pair of low-resolution hyperspectral image and a high-
% resolution multispectral image.
% 
% if you use this code, Please cite the following paper:
%
% T. Xu, T.Z. Huang, L.J. Deng, X.L. Zhao, and J.Huang, Hyperspectral Image
% Superresoluition Using Unidirectional Total Variation With Tucker
% Decomposition, IEEE J. Sel. Topics Appl. Earth Observ. Remote Sens.,doi:
% 10.1109/JSTARS.2020.3012566.

%--------------------------------------------------------------------------

clear
clc

addpath('functions', 'ImprovedDL', 'tensor_toolbox_2.6','sisal')
S=load('paints.mat');
S=S.paints;
%% HR-HSI
S=double(S);
S=S(1:512,1:512,:);
S=S/max(S(:));
[M,N,L] = size(S);
%% the downsampling matrix of the  spectral mode
F=create_F;
%%  simulate LR-HSI
S_bar = hyperConvert2D(S);
downsampling_scale=16;
s0=downsampling_scale/2;
BW=ones(16,1)/16;
BW1=psf2otf(BW,[M 1]);
S_w=ifft(fft(S).*repmat(BW1,1,N,L)); 

BH=ones(16,1)/16;
BH1=psf2otf(BH,[N 1]);
aa=fft(permute(S_w,[2 1 3]));
S_h=(aa.*repmat(BH1,1,M,L));
S_h= permute(ifft(S_h),[2 1 3]);  
 
Y_h=S_h(s0:downsampling_scale:end,s0:downsampling_scale:end,:);% uniform downsamping
Y_h_bar=hyperConvert2D(Y_h);
HSI=hyperConvert3D(Y_h_bar,M/downsampling_scale, N/downsampling_scale );

%% simulate HR-MS
Y = F*S_bar; 
MSI=hyperConvert3D(Y,M,N); 
%% The proposed Super-resolution method
HSI(HSI<0) = 0;
MSI(MSI<0) = 0;
par.n1 = 280;par.n2 = 300;par.n3 = 12;par.beta =1e-2;
par.lambda1 =1e-7;par.lambda2 = 1e-8; par.lambda3 = 1e-4; par.lambda4 = 1e2;   
par.mu1=1e-2;par.mu2=1e-2;par.mu3=1e-2;par.mu4=1e-2;%The penalty parameters mu1,mu2,mu3,and mu4 are fixed.

fprintf('It will take few minu. for a 512x512x3 image...\n')
t=clock;
Z4= UTV_SR(HSI,MSI,F,BW1,BH1,downsampling_scale,par,s0,S);
t4=etime(clock,t);

[psnr4,rmse4, ergas4, sam4, uiqi4,ssim4,DD4,CC4] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z4)), 0, 1.0/downsampling_scale);

 










