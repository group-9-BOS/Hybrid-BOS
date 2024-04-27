function [grad] = gardd(u,v)

%u   % dispacement in x direction
%v  % displacement in y direction

[row,column] = size(u);
% memory allocation for soting du and dv

du = zeros(row-2,column-2);
dv = zeros(row-2,column-2);

for i = 2 : row - 2
    for j = 2 : column - 2
        du(i,j) = (u(i,j+1)-u(i,j-1))/(1); % central difference formula
        dv(i,j) = (v(i+1,j)-v(i-1,j))/(1);
    end
end

grad = (du + dv);
%Rhs(isnan(grad)) = 0 ;


