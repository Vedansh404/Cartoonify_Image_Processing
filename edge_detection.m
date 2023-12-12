clc;close all;clear all;pkg load image;

img = imread('parrot.png');
img = double(img);

function [edge] = edgedetector(input_img)

img = double(input_img);

R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

[Rx,Ry] = gradient(R);
[Gx,Gy] = gradient(G);
[Bx,By] = gradient(B);


Sx = Rx.^2+Gx.^2+Bx.^2;
Sy = Ry.^2+Gy.^2+By.^2;

D = sqrt((Sx+Sy)/2);  
 
edge  = D;

end

edges = edgedetector(img);
edges = edges/max(edges(:));

subplot(2,1,1)
imshow(mat2gray(img));
title('Original')
subplot(2,1,2)
imshow(edges);
title('Edges')