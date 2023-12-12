clc; close all; clear all;
pkg load image;
tic();

filename = 'man.png';
img = imread(filename);
img = double(img);
img_smooth = img;

%% Apply Median Filter to remove Salt and Pepper noise
img_smooth(:,:,1) = medfilt2(img(:,:,1),[7,7]);
img_smooth(:,:,2) = medfilt2(img(:,:,2),[7,7]);
img_smooth(:,:,3) = medfilt2(img(:,:,3),[7,7]);

function [filteredImage ]= BilateralFiltering(corruptedImage,sd_sp,sd_int,windowSize)

#sd_sp      -> Standard deviation for spatial domain
#sd_int     -> Standard deviation for intensity domain
#windowSize -> Size of square mask

[x,y] = meshgrid(-windowSize:windowSize,-windowSize:windowSize);  # to get matrices of x and y to be used below
spaceGaussian = exp(-(x.^2+y.^2)/(2*(sd_sp^2)));                  # formula for gaussian filter

[m,n,d] = size(corruptedImage);

for i = 1:m
    for j = 1:n
        leftb  = max(i-windowSize,1);
        rightb = min(i+windowSize,m);
        topb   = max(j-windowSize,1);
        botb   = min(j+windowSize,n);
        
        for k = 1:d      
            submatrix = corruptedImage(leftb:rightb,topb:botb,k);
            intensityGaussian = exp(-(corruptedImage(i,j,k)-submatrix).^2/(2*sd_int^2));
            weightMatrix = intensityGaussian.*spaceGaussian((leftb:rightb)-i+windowSize+1,(topb:botb)-j+windowSize+1);
            filteredImage(i,j,k) = sum(sum(weightMatrix.*submatrix)) / sum(sum(weightMatrix));
        end
    end
end

end
##

%%

img_smooth = BilateralFiltering(img_smooth,10,20,3);
filtered = img_smooth;

##
function [edge] = edgedetector(input_img)  
input_image = uint8(input_img);
input_image = rgb2gray(input_image);
% Convert the image to double
input_image = double(input_image);
% Pre-allocate the filtered_image matrix with zeros
filtered_image = zeros(size(input_image));
% Sobel Operator Mask 
##Mx = [-1 0 1; -2 0 2; -1 0 1]; 
##My = [-1 -2 -1; 0 0 0; 1 2 1];
% Canny Operator Mask
Mx = [-1 -1 -1; 0 0 0; 1 1 1];
My = [-1 0 1; -1 0 1; -1 0 1];

  for i = 1:size(input_image, 1) - 2 
      for j = 1:size(input_image, 2) - 2 
  
        % Gradient approximations
        Gx = sum(sum(Mx.*input_image(i:i+2, j:j+2)));
        Gy = sum(sum(My.*input_image(i:i+2, j:j+2)));
                 
        % Calculate magnitude of vector
        filtered_image(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
         
      end
  end
filtered_image = uint8(filtered_image); 
edge=filtered_image;

end
##
edges       = edgedetector(img);
edges_bw    = 255-edges;
edges1      = edges/max(edges(:));
cartoon_img = filtered;
%%
a = 24;  % Quantization Factor
for i = 1:3
    % Quantize the Values
    t = a*floor(filtered(:,:,i)./a);
    t(edges1>0.18) = 0;
    cartoon_img(:,:,i) = t;
end
% file_name = strcat(['../Results/' filename(9:length(filename)-4) '_toon.' 'png']);
% imwrite(mat2gray(cartoon_img),file_name)
% end

figure
subplot(2,1,1)
imshow(mat2gray(img)); title('Original');
subplot(2,2,3)
imshow(edges); title('Edges');
subplot(2,2,4)
imshow(mat2gray(cartoon_img)); title('Cartoon');
toc()