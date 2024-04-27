function [ux,uy,ux_horn,uy_horn,error1]=opticalflow(I1,I2)
%% Horn's solution as an initial approximation of u and v  

%coefficient matrices for finite difference scheme
D1 = [0, 0, 0; 0,-1,-1;0,1,1]/2; %1st order 
F1 = [0, 0, 0; 0,1,1;0,1,1]/4; %2nd order

%applying convolution using coefficient matrices
Ix = imfilter((I1+I2)/2, D1, 'symmetric',  'same'); 
Iy = imfilter((I1+I2)/2, D1', 'symmetric',  'same');
It = imfilter(I2-I1, F1, 'symmetric',  'same');

%iteration parameters
maxiter_horn=500; 
tol_horn = 10^(-12);


[u,v] = hs_estimation(Ix, Iy, It, tol_horn, maxiter_horn);
ux_horn = v;
uy_horn = u;


%% new approximation of u and v using physics developed by liu-shen for optical flow

Dm=0*10^(-3);
h=1; %step size

%using central diffference formula for calculating laplacian 
H = [1, 1, 1; 1,0,1;1,1,1]; %coefficient matrix
gradu = -u.*imfilter(ones(size(u)), H/(h*h), 'same') + imfilter(u, H/(h*h), 'same');
f=Dm*gradu;

%iteration parameters
maxiter=100;
tol = 10^(-8);

dx=1; 
dt=1; % unit time


[u,v,error1] = liu_shen_estimator(I1, I2, f, dx, dt, tol, maxiter, uy_horn, ux_horn);

ux=v;
uy=u;

end




