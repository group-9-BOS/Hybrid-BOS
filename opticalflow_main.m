%%This code is the implementation of Horn-Schunck based optical flow algorithm for BOS

function [ux_of,uy_of] = opticalflow_main(Im1,Im2)

% For Image Pre-Processing

% For local illumination intensity adjustment
size_average=13; % in pixels

% Gausian filter size for removing random noise in images
size_filter=10; % in pixels

% correction of the global and local intensity change in images
I1=Im1;
I2=Im2;
[m1,n1]=size(I1);
window_shifting=[1;n1;1;m1]; % [x1,x2,y1,y2] defines a rectangular window for global correction
[I1,I2]=illumination_correction(I1,I2,window_shifting,size_average);

%filtering to reduce random noise and downsampling images if displacements are large
%  [I1,I2] = filtering(I1,I2,scale_im,size_filter);



%% optical flow calculation for a refined field 

[ux,uy,ux_horn,uy_horn,error1]=opticalflow(I1,I2);
ux_of=ux;    
uy_of=uy;    

end
