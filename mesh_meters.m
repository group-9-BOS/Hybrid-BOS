function mesh_meters(data, step,scale)

[rows, cols] = size(data);
divisions = floor(rows/step); %number of intervals/divisions in PIVLab plot

xStart = 1;
dx = 1;
%N = rows;
N = cols;
x = xStart + (0:N-1)*dx; %list of divisions in px units
yStart = 1;
dy = 1;
% M = cols;
M = rows;
y = yStart + (0:M-1)*dy; %list of divisions in px units

X = x*step*scale;%convert px to m
Y = y*step*scale;%convert px to m

% ylim([0 0.105]);
% xlim([0 0.077]);

mesh(X,Y, data);
%surf(X,Y, data)