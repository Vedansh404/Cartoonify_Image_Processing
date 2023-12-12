clc;close all;clear all;pkg load image;
tic();
% filename = '../data/flower.png';
% filename = '../data/dolphin.jpeg';
% for ind = 3:12
% filename = strcat(['../data/im_' int2str(ind) '.bmp']);
%filename = '../data/im_1.bmp';

img = imread('parrot.png');
img = double(img);

a = 24;    %Quantization Factor

img_smooth = img;
%%  Apply Median Filter to remove Salt and Pepper noise
img_smooth(:,:,1) = medfilt2(img(:,:,1),[7,7]);
img_smooth(:,:,2) = medfilt2(img(:,:,2),[7,7]);
img_smooth(:,:,3) = medfilt2(img(:,:,3),[7,7]);



function [filteredImage ]= myBilateralFiltering(corruptedImage,sd_sp,sd_int,windowSize)


[x,y] = meshgrid(-windowSize:windowSize,-windowSize:windowSize);
spaceGaussian = exp(-(x.^2+y.^2)/(2*(sd_sp^2)));


[ m, n,d ] = size(corruptedImage);


for i = 1:m
    for j = 1:n
        leftb = max(i-windowSize,1);
        rightb = min(i+windowSize,m);
        topb = max(j-windowSize,1);
        botb = min(j+windowSize,n);
        
        for k = 1:d      
            submatrix = corruptedImage(leftb:rightb,topb:botb,k);
            intensityGaussian = exp(-(corruptedImage(i,j,k)-submatrix).^2/(2*sd_int^2));
            weightMatrix = intensityGaussian.*spaceGaussian((leftb:rightb)-i+windowSize+1,(topb:botb)-j+windowSize+1);
            filteredImage(i,j,k) = sum(sum(weightMatrix.*submatrix))/sum(sum(weightMatrix));
        end
    end
end

end






function [edge] = edgedetector(input_img)

img = double(input_img);

R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

[Rx,Ry] = gradient(R);
[Gx,Gy] = gradient(G);
[Bx,By] = gradient(B);

% Sx = sqrt(Rx.^2+Gx.^2+Bx.^2);
% Sy = sqrt( Ry.^2+Gy.^2+By.^2);
% Sxy = sqrt(Rx.*Ry+Gx.*Gy+Bx.*By);

Sx = Rx.^2+Gx.^2+Bx.^2;
Sy = Ry.^2+Gy.^2+By.^2;
Sxy = Rx.*Ry+Gx.*Gy+Bx.*By;


D = sqrt(abs((Sx+Sy).^2-4*(1)*(Sx.*Sy-Sxy.^2))); % Discriminant of the Characteristic Equation of Image structure matrix
eigVal1 = sqrt((Sx+Sy+D)/2);  % Solutions of Characteristic Equation are the Eigen Values
% eigVal2 = (Sx+Sy-D)/2;  

edge  = eigVal1;

end




%%
for i=1:3
    img_smooth = myBilateralFiltering(img_smooth,10,20,3);
end
filtered = img_smooth;

edges = edgedetector(img);
edges = edges/max(edges(:));
cartoon_img = filtered;


%%
for i = 1:3
    % Quantize the Values
    t = a*floor(filtered(:,:,i)./a);
    t(edges>0.18) = 0;
    cartoon_img(:,:,i) = t;
end
% file_name = strcat(['../Results/' filename(9:length(filename)-4) '_toon.' 'png']);
% imwrite(mat2gray(cartoon_img),file_name)
% end
figure
subplot(2,2,4)
imshow(mat2gray(cartoon_img));
title('Cartoon')
subplot(2,1,1)
imshow(mat2gray(img));
title('Original')
subplot(2,2,3)
imshow(edges);
title('Edges')

toc()