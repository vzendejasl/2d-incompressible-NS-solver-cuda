
%% Appendix I: Taylor Green Vortex

clear
close all

global dx dy Nx Ny Re K1_prime K2_prime

Nx = 2048;
Lx = 2*pi;
Ny = Nx;
Ly = Lx;
Re = 5000;
perturb = 0; % 1 for perturbations
analytical = 0; % Compare to analytical solution for A=a=b=1, B=-1 case

dx = Lx/Nx;
dy = Ly/Ny;

% Staggered grid coordinates for U-momentum equation
x_umom = linspace(0,Lx-dx,Nx);
y_umom = linspace(dy/2,Ly-dy/2,Ny);
[X_umom,Y_umom] = meshgrid(x_umom,y_umom);
X_umom = X_umom';
Y_umom = Y_umom';

% Staggered grid coordinates for V-momentum equation
x_vmom = linspace(dx/2,Lx-dx/2,Nx);
y_vmom = linspace(0,Ly-dy,Ny);
[X_vmom,Y_vmom] = meshgrid(x_vmom,y_vmom);
X_vmom = X_vmom';
Y_vmom = Y_vmom';

% Cell-centered X,Y grid coordinates (for curl plots)
[X,Y] = meshgrid(x_vmom,y_umom);
X = X';
Y = Y';

% Modified wavenumbers
K1 = 2.*pi./Lx.*(-Nx/2:Nx/2-1).*ones(Nx);
K2 = (2.*pi./Ly.*(-Ny/2:Ny/2-1).*ones(Ny))';
K1_prime = 1./dx.*sqrt(2-2.*cos(K1.*dx));
K2_prime = 1./dy.*sqrt(2-2.*cos(K2.*dy));

% Initial Condition
A = 1;
a = 1;
B = -1;
b = 1;

% Perturbation term
if perturb == 1
    V_perturb = (-1+2.*rand(Nx,Ny)).*1e-2;
else
    V_perturb = zeros(Nx,Ny);
end
U_num = A.*sin(a.*X_umom).*cos(b.*Y_umom);
V_num = B.*cos(a.*X_vmom).*sin(b.*Y_vmom)+V_perturb;

% Iterate in time
if analytical == 1
    dt = 0.001;
    t_final = 0.1;
else
    dt = 0.02;
    t_plot = [1 25 50 60];
    t_final = t_plot(end);
    t_plot_index = 1;
    figure;
end

t0 = clock;
t = dt;
count = 0;
while t <= 50*dt
    % Numerical solution
    [U_num,V_num] = time_advance_RK3(U_num,V_num,dt);
    curl_num = compute_curl(U_num,V_num);
%     
%     if (analytical == 0) && (t >= t_plot(t_plot_index))
%         subplot(2,2,t_plot_index);
%         contourf(X',Y',curl_num');
%         xlabel('x'); ylabel('y'); title("t = "+num2str(t_plot(t_plot_index)));
%         colorbar(); caxis([-2 2]);
%         t_plot_index = t_plot_index+1;
%     end

%     % For making video:
%     figure(1); 
%     contourf(X',Y',curl_num');
%     xlabel('x'); ylabel('y');
%     title("\omega(x,y), t = "+num2str(t,'%.1f'));
%     colorbar(); caxis([-2 2]);

    % Update time
    t = t + dt;
    count = count +1
end
ms = round(etime(clock,t0) * 1000)

% Make plots
if analytical == 1
    
    U_exact = sin(X_umom).*cos(Y_umom).*exp(-2.*t_final./Re);
    V_exact = -cos(X_vmom).*sin(Y_vmom).*exp(-2.*t_final./Re);
    curl_exact = compute_curl(U_exact,V_exact);
    
    % Contour Plots at t = 0.1
    figure;
    subplot(2,2,1); contourf(X_umom',Y_umom',U_num'); 
    xlabel('x'); ylabel('y'); title('U(x,y)'); colorbar();
    subplot(2,2,2); contourf(X_vmom',Y_vmom',V_num');
    xlabel('x'); ylabel('y'); title('V(x,y)');colorbar();
    subplot(2,2,3); contourf(X',Y',curl_num');
    xlabel('x'); ylabel('y'); title('\omega(x,y)'); colorbar();
    subplot(2,2,4); contourf(X',Y',curl_num'-curl_exact');
    xlabel('x'); ylabel('y'); title('|\omega_{numerical}(x,y) - \omega_{exact}(x,y)|'); colorbar();
    
    % Error plot
    % *NOTE*: These errors were computed by running this main script with
    % varying N values from 64 to 256 and manually computing the max error as:
    % max(max(abs(exact sol. - numerical sol.)))
    figure;
    N = [64 128 192 256];
    U_max_error = [3.2077e-08 8.0285e-09 3.5690e-09 2.0077e-09];
    V_max_error = [3.2077e-08 8.0285e-09 3.5690e-09 2.0077e-09];
    curl_max_error = [6.3898e-08 1.6041e-08 7.1348e-09 4.0144e-09];
    hold on
    plot(N,U_max_error,'-or','LineWidth',2);
    plot(N,V_max_error,'--xb','LineWidth',3);
    plot(N,curl_max_error,'-ok','LineWidth',2);
    hold off
    grid on
    xlabel('N'); ylabel('max_{x,y}(|Exact Sol.- Numerical Sol.|)');
    set(gca,'XScale','log','YScale','log');
    legend('U(x,y)','V(x,y)','\omega(x,y)','Location','Southwest');
    
end
