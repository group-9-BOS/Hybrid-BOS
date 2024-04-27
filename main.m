%% This code implements cross-correlation algorithm for the coarse displacement field 
% and then refines the field using optical flow algorithm
clear all
close all
%% input
I1     = imread('shockbg.tif');                                             % Load image 1
I1     = I1(:,:,1);  

I2     = imread('shockdis.tif');                                             % Load image 2
I2 = I2(:,:,1);  
% 
% % Convert to 8-bit
% I1 = im2uint8(I1);
% I2 = im2uint8(I2);

% Display the bit depth and class of the images
info1 = imfinfo('shockbg.tif');
info2 = imfinfo('shockdis.tif');

bitDepth1 = info1.BitDepth;
bitDepth2 = info2.BitDepth;

disp(['Bit Depth of Image 1: ', num2str(bitDepth1)]);
disp(['Bit Depth of Image 2: ', num2str(bitDepth2)]);
disp(size(I1));
disp(class(I1));

% Convert color images to grayscale
% if size(I1, 3) > 1
%    I1 = rgb2gray(I1);
% end
% if size(I2, 3) > 1
%     I2 = rgb2gray(I2);
% end

% %storing original images
% I1Orig = imread('I1_vortexpair_dt0p03.tif');                                          
% I1Orig     = I1Orig(:,:,1);  
% 
% I2Orig = imread('I2_vortexpair_dt0p03.tif');                                          
% I2Orig     = I2Orig(:,:,1);  


%% cross-correlation algorithm and parameters

[ux_cc,uy_cc]= cc_main(I1,I2);

%% optical flow algorithm

%Input to the optical flow
Im1=imread('image_with_FD.tif'); %using FD image based on accuracy

Im2=I2;

[ux_of, uy_of]=opticalflow_main(Im1, Im2);

%% final displacement values after implementation of both algorithms

[u_final,v_final]= added(ux_of,uy_of,ux_cc,uy_cc);
%%
esc = 7692.31;                  %esc = how much pixel is in 1 meter in image
h = 0.02;                    %daimeter of nozzel
L = 0.285;                    %distance between background and subject
Mag = 2.7;                   %Magnification 
n_0=1.0002921;                %Refractive Index of sorrounding gas (air)
rho_0=1.225;                  %Density of sorrounding gas (air)
G = 2.2649e-4;                %Gladstone-Dale constant for air
scale= (1/esc); 

[density] = density_est(u_final,v_final); %Density estimation
%% displaying results and plots

% Grid for x and y coordinates in pixels
x = (1:size(u_final, 2)) ;
y = (1:size(u_final, 1)) ;
[X Y]=meshgrid(x,y);

% Contour Plot of Ux
figure(1);
imagesc(x, y, u_final);
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
title('Contour Plot of U');
colorbar;
colormap('parula');
hold off;

% Contour Plot of Uy
figure(2);
imagesc(x, y, v_final);
xlabel('x (pixels)');
ylabel('y (pixels)');
axis image;
set(gca, 'YDir', 'reverse');
title('Contour Plot of V');
colorbar;
colormap('parula');
hold off;

figure11 = figure;
axes1 = axes('Parent',figure11);
mesh_meters(density,1,scale)
shading interp
view(0,-90)
colorbar('peer',axes1);
xlabel('X (m)');
ylabel('Y (m)');
title('Contour plot of Density (kg/m^3)');
%% post-processing code
% [x,y,z]=peaks(100);
% z=.4*z;
% tri=delaunay(x,y);
% vtkwrite
% vtkwrite('density.vtk','structured_points','density',density,'x',x,'y',y)
