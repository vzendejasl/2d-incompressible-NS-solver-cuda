
%% Appendix F: Test of Divergence-Free Condition for Projection Method

clear
close all

global dx dy Nx Ny Re K1_prime K2_prime

Nx = 32;
Lx = 2*pi;
Ny = Nx;
Ly = Lx;
Re = 1;

dx = Lx/Nx;
dy = Ly/Ny;

% Compute modified wavenumbers
K1 = 2.*pi./Lx.*(-Nx/2:Nx/2-1).*ones(Nx);
K2 = (2.*pi./Ly.*(-Ny/2:Ny/2-1).*ones(Ny))';
K1_prime = 1./dx.*sqrt(2-2.*cos(K1.*dx));
K2_prime = 1./dy.*sqrt(2-2.*cos(K2.*dy));

% Define grid
x = linspace(0,Lx-dx,Nx);
y = linspace(0,Ly-dy,Ny);
[X,Y] = meshgrid(x,y);
X = X';
Y = Y';

% Make randomly generated velocity field divergence-free
U0 = rand(Nx,Ny);
V0 = rand(Nx,Ny);
[U,V] = compute_projection_step(U0,V0);

% Compute divergence of updated flowfield
dudx = zeros(Nx,Ny);
dudx(1:end-1,:) = (U(2:end,:)-U(1:end-1,:))./dx;
dudx(end,:) = (U(1,:)-U(end,:))./dx;

dvdy = zeros(Nx,Ny);
dvdy(:,1:end-1) = (V(:,2:end)-V(:,1:end-1))./dy;
dvdy(:,end) = (V(:,1)-V(:,end))./dy;

div = dudx + dvdy;

disp("Maximum divergence = " + num2str(max(max(abs(div)))));
