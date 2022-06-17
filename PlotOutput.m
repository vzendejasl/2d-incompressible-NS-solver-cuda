clear;clc;close all;

% Grid points need to match original computation
Nx = 2048;
Lx = 2*pi;
Ny = Nx;
Ly = Lx;
dx = Lx/Nx;
dy = Ly/Ny;

% Staggered grid coordinates for U-momentum equation
x_umom = linspace(0,Lx-dx,Nx);
y_umom = linspace(dy/2,Ly-dy/2,Ny);

% Staggered grid coordinates for V-momentum equation
x_vmom = linspace(dx/2,Lx-dx/2,Nx);
y_vmom = linspace(0,Ly-dy,Ny);

% Cell-centered X,Y grid coordinates (for curl plots)
[X,Y] = meshgrid(x_vmom,y_umom);
X = X';
Y = Y';

%CPU visualization
curl_num = importdata("CPU_output.txt");

figure(1);
contourf(X',Y',curl_num');
xlabel('x'); ylabel('y');
colorbar(); caxis([-2 2]);
title("CPU Output")

%GPU visualization
curl_num = importdata("GPU_output.txt");

figure(2);
contourf(X',Y',curl_num');
xlabel('x'); ylabel('y');
colorbar(); caxis([-2 2]);
title("GPU Output")