function [u,v,error] = liu_shen_estimator(I0, I1, f, dx, dt, tol, maxnum, u0, v0)
% Ix,Iy,It: represent the spatial and temporal derivatives of intensity
lambda_phy = 2000; %regularization parameter

%coefficient matrices for derivatives 
Der = [0, -1, 0; 0,0,0; 0,1,0]/2; %%% partial derivative 
Mix = [1, 0, -1; 0,0,0;-1,0,1]/4; %%% mixed partial derivatives
A = [0, 1, 0; 0,0,0;0,1,0]; %%% average
Der2 =  [0, 1, 0; 0,-2,0;0,1,0]; %%% partial derivative
H = [1, 1, 1; 1,0,1;1,1,1]; 
h=1;

IIx = I0.*imfilter(I0, Der/dx, 'replicate',  'same');
IIy = I0.*imfilter(I0, Der'/dx, 'replicate',  'same');
II = I0.*I0;
Ixt = I0.*imfilter((I1-I0)/dt-f, Der/dx, 'replicate',  'same');
Iyt = I0.*imfilter((I1-I0)/dt-f, Der'/dx, 'replicate',  'same');

k=0;
total_error=100000000;
u=u0;
v=v0;

[r,c]=size(I1);


%Matrix manipulation
[r,c]=size(I0);
cmtx = imfilter(ones(size(I0)), H/(h*h), 'same');

A11 = I0.*(imfilter(I0, Der2/(h*h), 'replicate',  'same')-2*I0/(h*h)) - lambda_phy*cmtx; 
A22 = I0.*(imfilter(I0, Der2'/(h*h), 'replicate',  'same')-2*I0/(h*h)) - lambda_phy*cmtx; 
A12 = I0.*imfilter(I0, Mix/(h*h), 'replicate',  'same'); 
    
DetA = A11.*A22-A12.*A12;

B11 = A22./DetA;
B12 = -A12./DetA;
B22 = A11./DetA;

error=[];
while total_error > tol & k < maxnum
    total_error;
    bu = 2*IIx.*imfilter(u, Der/dx, 'replicate',  'same')+ ... 
        IIx.*imfilter(v, Der'/dx, 'replicate',  'same')+ ...
        IIy.*imfilter(v, Der/dx, 'replicate',  'same') + ... 
        II.*imfilter(u, A/(dx*dx), 'replicate',  'same')+ ... 
        II.*imfilter(v, Mix/(dx*dx), 'replicate',  'same') + ... 
        lambda_phy*imfilter(u, H/(dx*dx), 'same')+Ixt;
    
    bv = IIy.*imfilter(u, Der/dx, 'replicate',  'same') + ...
         IIx.*imfilter(u, Der'/dx, 'replicate',  'same') + ...
         2*IIy.*imfilter(v, Der'/dx, 'replicate',  'same')+ ...
         II.*imfilter(u, Mix/(dx*dx), 'replicate',  'same') + ...
         II.*imfilter(v, A'/(dx*dx), 'replicate',  'same')+ ... 
         lambda_phy*imfilter(v, H/(dx*dx), 'same')+Iyt;
     
    unew = -(B11.*bu+B12.*bv);
    vnew = -(B12.*bu+B22.*bv);
    total_error = (norm(unew-u,'fro')+norm(vnew-v,'fro'))/(r*c);
    u = unew;
    v = vnew;
    error=[error; total_error];
    k=k+1 ;
end