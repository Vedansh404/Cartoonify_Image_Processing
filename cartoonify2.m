clc;
clear; 
close all;
pkg load image;

#addpath('add the path of the folder where the gs function is stored');
% Read the image 
a=imread('baboon.png');
% a = imresize(a,2);
b = size(a);
b1 = size(a);
%imshow(a)
% Convert to Ycbcr 
a1 = rgb2ycbcr(a);

%Choose the First plane
a = a1(:,:,1);
a = double(a);

%Initialize the parameters
n = 11;                             %Filter Size
n1=ceil(n/2);
vars =150;                           %Spacial Variance
varr = 150;                          %Pixel Value Variance
c=0;
c1=0;
r=1;                                %No. Of Reptitions
msg = 'Almost done';
x = 0;
f = waitbar(x,msg);
%%

function out = gs(p,q,var)
%     d = size(p);
    out = 1/((2*pi)*var^2)*exp(-(p-q)*(p-q)'/(2*var^2));
    return

    

%Bilateral Filter loop
for i1 = 1:r
for i=n1:b(1)-n1
    for j=n1:b(2)-n1
        for k=1:n
            for l=1:n
            c=c+gs(sqrt((-n1+k)^2+(-n1+l)^2),0,vars)*gs(a(i-n1+k,j-n1+l),a(i,j),varr)*a(i-n1+k,j-n1+l);
            c1=c1+gs(sqrt((-n1+k)^2+(-n1+l)^2),0,vars)*gs(a(i-n1+k,j-n1+l),a(i,j),varr);
            end
        end
        
        d(i-n1+1,j-n1+1)=c/c1;
        c=0;
        c1=0;
    end
x = i/(b(1)-n1);
waitbar(x,f)  
end

a=d;
clear d
b = size(a);
end
close(f)
%%

%Merge the Y plane back
d1 = uint8(a);
d1(:,:,2)= a1(r*n1-(r-1):b1(1)-(r*n1),r*n1:b1(2)-r*n1+r-1,2);
d1(:,:,3)= a1(r*n1-(r-1):b1(1)-(r*n1),r*n1:b1(2)-r*n1+r-1,3);

figure;
subplot(1,2,1)
imshow(uint8(ycbcr2rgb(a1)));
title('Original Image')

subplot(1,2,2)
imshow(ycbcr2rgb(d1));
title('Bilateral Filter Output Image')