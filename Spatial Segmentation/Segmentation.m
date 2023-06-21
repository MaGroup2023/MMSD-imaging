
%%
clear;clc;close all;
R = 263; %行数
C1 = 114; %列数
C2 = 114;
N = 300; %m/z个数
[mz,ionMapsLE,ionMapsHE,TICLE,TICHE] = MS1Convert("20220624-004-MBSI-neg-ms1-slic-im3-tran40.txt",R,C1,C2);
%%
[Smz,Sid] = sort(mz,'ascend');
ionMapsLE = ionMapsLE(:,:,Sid)./TICLE;
ionMapsHE = ionMapsHE(:,:,Sid)./TICHE;
DataMatrixHE = reshape(ionMapsHE,R*(C1+C2),[]);
DataMatrixLE = reshape(ionMapsLE,R*(C1+C2),[]);

%%
% load("Mask.mat");
Mask = ones(R,C1+C2);
MList = reshape(Mask,R*(C1+C2),[]);
DataMatrixHE(MList==0,:)=0;
DataMatrixLE(MList==0,:)=0;
[Coeff,Score,Latent,~,Exp] = pca(DataMatrixLE);
S = Score(:,1:4);
S = reshape(S,R,(C1+C2),[]);
%% Spatial aware segmentation
r =3;
alpha = zeros(2*r+1,2*r+1);
for i = -r:r
    for j = -r:r
        alpha(i+r+1,j+r+1) = 1/((i^2+j^2)/2+1);
    end
end
% alpha = [0.5,0.7,0.5;0.7,1,0.7;0.5,0.7,0.5];
% Mask = ones(R,C1+C2);
K= 4;
[ClassID,CentroidSpectrum,SADistance,Iter,Loss] = SAKMeans(S,K,Mask,r,alpha,45);
%%
save ClassID.mat ClassID -mat
%%
figure(5)
plot(log10(Loss))
set(gca,'TickDir','out');
box off
grid off
%%
figure(6)
% Mask = (ClassID~=0);
Mask = (ClassID==1)|(ClassID==2)|(ClassID==4)|(ClassID==6)|(ClassID==7)|(ClassID==10);
imshow(Mask,[]),colormap gray;
save("Mask.mat","Mask","-mat");

