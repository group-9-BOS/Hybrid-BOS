function [density] = density_est(u_final,v_final)

 u= u_final;
 v = v_final;
[x,y] = size(u_final);    % size of image
[grid_X,grid_Y] = meshgrid(1:y,1:x); % creating grid points

figure,quiver(u,v); % plotting x and y displacement
%%
% esc = 47167;                  %esc = how much pixel is in 1 meter in image
% h = 0.004;                    %daimeter of nozzel
% L = 0.147;                    %distance between background and subject
% Mag = 0.94;                   %Magnification 
% n_0=1.0002921;                %Refractive Index of sorrounding gas (air)
% rho_0=1.204;                  %Density of sorrounding gas (air)
% G = 2.2649e-4;                %Gladstone-Dale constant for air
% scale= (1/esc); 
 %%
%1/4/2024
esc = 7692.31;                  %esc = how much pixel is in 1 meter in image
h = 0.02;                    %daimeter of nozzel
L = 0.285;                    %distance between background and subject
Mag = 2.7;                   %Magnification 
n_0=1.0002921;                %Refractive Index of sorrounding gas (air)
rho_0=1.225;                  %Density of sorrounding gas (air)
G = 2.2649e-4;                %Gladstone-Dale constant for air
scale= (1/esc); 
%%

grid_X = grid_X*scale;        % changing grid oints to meter value
grid_Y = grid_Y*scale;

% grad_x = gradient(u);    % finding x gradient
% grad_y = gradient(v);    % finding y gradient
[grad] = gradd(u,v);       % calling function gradd to calculate sum of the gradient
% rhs = grad_x + grad_y;
rhs = grad;             
k= (n_0)/(G*Mag*L*h);      %Constant term on Poisson equation
RHS=real((k)*rhs/esc^2);

%%
[uu] = jacobi(x-1,y-1,RHS,scale);   % calling function jacobi to solve for poisson equation
%%
density=uu;
%plotting the results thus obtained
figure11 = figure;
axes1 = axes('Parent',figure11);
mesh_meters(density,1,scale)
shading interp
view(0,90)
colorbar('peer',axes1);
xlabel('X (m)');
ylabel('Y (m)');
title('Contour plot of Density (kg/m^3)');
% % figure,contourf(uu)
% imsave(uu)
% vtkwrite( 'try1.vtk','structured_points','uu',uu,'grid_X',grid_X,'grid_Y',grid_Y);

end
