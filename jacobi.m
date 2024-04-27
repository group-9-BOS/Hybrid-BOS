function [uu] = jacobi(x,y,RHS,scale,rho)
b = RHS;
P = ones(x,y); %initializing value
uu = ones(x,y);
dx =1% x*scale/(x-1);   %grid width 1 pixel as resolution of image is equal to displacement vectors 
dy = 1%y*scale/(y-1);
iter = 0; %initial itteration
maxiter = 3000; % maximum no of itteration
        %boundary condition initialization
           %Neuman boundary condition
           uu(:,1) =  uu(:,2); 
           uu(:,end)=    uu(:,end-1); 
           %Dirichlet boundaey condition
%          uu(1,:) =    1.225;
%           uu(end,:) =  1.225;
%%
i = 2:x-1;
 j = 2:y-1;
         tolerance = 0.0001;
      
         error = inf; %initial residual
         while error > tolerance
             iter = iter+1;
          uu(i,j) = -(dx^2/4)*(-b(i,j) - (P(i+1,j)+P(i-1,j)+P(i,j+1)+P(i,j-1))/dx^2);  % central difference second order formula
         
          %boundary condition inside loop
          %Neuman boundary condition
          uu(:,1) =  uu(:,2); 
          uu(:,end)=    uu(:,end-1); 
           %Dirichlet boundaey condition
%           uu(1,:) =    1.225;
%           uu(end,:) =  1.225;
          
          error = max(max(abs((uu-P)./uu)));
          P = uu;
           disp(['Maximum error is  ',num2str(error)]);
           disp(['iteration  ',num2str(iter)]);
         if iter > maxiter
             break;
         end
         end
end

