% this code is based on the paper 
% An improved FOG removal method for the traffic monitoring 
%image by Wenli Zhang 
% Author of the Code : B.Abraham David 

clc
clear all;
close all;
cd IMAGES
[J P]=uigetfile('*.*','select the source file');
I=imread(strcat(P,J));
cd ..
figure,imshow(I);title('Input Image');
I=double(I);
%==== applying Drak channel ===================
[height, width, ~] = size(I);
patchSize = 15;                             %the patch size is set to be 15 x 15
padSize = 7;                                % half the patch size to pad the image  
JD=zeros(height, width);                    % the dark channel
tim = padarray(I,[padSize padSize], Inf); 
for j = 1:height
    for i = 1:width
        patch =tim(j:(j+patchSize-1),i:(i+patchSize-1),:);
        JD(j,i) = min(patch(:));
     end
end
%======================================
% atmospheric Light adjustment
A=atmLight(double(I),JD);
%=======================================
omega = 0.95;           % the amount of haze we're keeping
im3 = zeros(size(I));
for ind = 1:3 
    im3(:,:,ind) = I(:,:,ind)./A(ind);
end
[height, width, ~] = size(im3);
patchSize = 15;                             %the patch size is set to be 15 x 15
padSize = 7;                                % half the patch size to pad the image  
JD=zeros(height, width);                    % the dark channel
timm = padarray(im3,[padSize padSize], Inf); 
for j = 1:height
    for i = 1:width
        patch =timm(j:(j+patchSize-1),i:(i+patchSize-1),:);
        JDD(j,i) = min(patch(:));
     end
end
out=1-omega*JDD;
t0 = 0.1;

J = zeros(size(I));
for ind = 1:3
   J(:,:,ind) = A(ind) + (I(:,:,ind) - A(ind))./max(out,t0); 
end
J = J./(max(max(max(J))));
figure,imshow(J);title('Dehazed Image');
%======================================================
%========= Retinex Algorithm =========================
[Y,X]=meshgrid(1:height,1:width);
c=500;
Fnok = exp(-((X.^2)+(Y.^2))./(c.^2));
K = 1/(sum(sum(Fnok)));
F = K.*Fnok; 
IR = double(I(:,:,1));
IG = double(I(:,:,2));
IB = double(I(:,:,3));
FF  = fftshift(fft2(F'));
IFR = fftshift(fft2(IR)); IFR=FF.*IFR; IFR=real(ifft2(ifftshift(IFR)));
IFG = fftshift(fft2(IG)); IFG=FF.*IFG; IFG=real(ifft2(ifftshift(IFG))); 
IFB = fftshift(fft2(IB)); IFB=FF.*IFB; IFB=real(ifft2(ifftshift(IFB))); 
RR = log10(double(IR))-log10(IFR);
RG = log10(double(IG))-log10(IFG);
RB = log10(double(IB))-log10(IFB);
OUT(:,:,1)=uint8(255*RR/(max(max([RR RG RB]))));
OUT(:,:,2)=uint8(255*RG/(max(max([RR RG RB]))));
OUT(:,:,3)=uint8(255*RB/(max(max([RR RG RB]))));
figure,imshow(OUT);title('Retinex Image');
%=========================================================
[a h v d]=dwt2(OUT(:,:,1),'haar');
R1(:,:,1)=idwt2(a,2*h, 2*v,2*d,'haar');
[a h v d]=dwt2(OUT(:,:,2),'haar');
R1(:,:,2)=idwt2(a,2*h, 2*v,2*d,'haar');
[a h v d]=dwt2(OUT(:,:,3),'haar');
R1(:,:,3)=idwt2(a,2*h, 2*v,2*d,'haar');
figure,imshow(uint8(R1));title('R+WT method');
%====================================

